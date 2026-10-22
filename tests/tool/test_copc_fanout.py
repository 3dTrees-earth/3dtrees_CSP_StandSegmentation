"""One spatial decode, exact point order per instance, and bounded pipe failure."""
from concurrent.futures import ThreadPoolExecutor
from io import BytesIO
import json
import os
from pathlib import Path
import sys
import tempfile
import threading
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'workflow'))
sys.path.insert(0, str(ROOT / 'tests/tool'))
from copc_fanout import InstanceFanout, dispatch
from test_copc_owned_inventory import configuration
import copc_owned_inventory as owned


class FanoutTests(unittest.TestCase):
    def test_single_read_preserves_all_owned_points_and_instance_order(self):
        config = configuration()
        config.update(core_size=100, buffer=20)
        index = owned.build_index(config)
        for layout in index['layouts']:
            for number, window in enumerate(layout['windows']):
                serial = BytesIO()
                reference_audit = owned.stream_window(config, index, layout['slug'], number, serial)
                width = len(layout['columns'])
                id_column = layout['columns'].index(layout['instance'])
                expected = np.frombuffer(serial.getvalue(), dtype='<f8').reshape(-1, width)
                with tempfile.TemporaryDirectory() as directory, ThreadPoolExecutor() as pool:
                    groups, futures = [], []
                    for worker, ids in enumerate(np.array_split(window['ids'], min(3, len(window['ids'])))):
                        path = Path(directory)/f'{worker}.fifo'
                        os.mkfifo(path, 0o600)
                        groups.append(dict(ids=[int(i) for i in ids], pipe=str(path)))
                        futures.append(pool.submit(path.read_bytes))
                    audit_path = Path(directory)/'audit.json'
                    audit = dispatch(config, index, layout['slug'], number, groups, audit_path)
                    self.assertTrue(audit['complete'])
                    self.assertEqual(audit['decoded_points'], reference_audit['decoded_points'])
                    self.assertEqual(audit['observed_counts'], reference_audit['observed_counts'])
                    self.assertTrue(json.loads(audit_path.read_text())['complete'])
                    for group, future in zip(groups, futures):
                        actual = np.frombuffer(future.result(timeout=5), dtype='<f8').reshape(-1, width)
                        self.assertTrue(set(actual[:, id_column]).issubset(group['ids']))
                        for key in group['ids']:
                            np.testing.assert_array_equal(actual[actual[:, id_column] == key],
                                                          expected[expected[:, id_column] == key])
                    self.assertFalse(any(path.suffix == '.bin' for path in Path(directory).iterdir()))

    def test_prefetch_crosses_batches_while_another_worker_is_blocked(self):
        # Worker 1 deliberately consumes nothing until worker 2 gets a later
        # batch. A synchronous per-batch drain deadlocks in this scenario.
        with tempfile.TemporaryDirectory() as directory, ThreadPoolExecutor() as pool:
            paths = [Path(directory)/f'{i}.fifo' for i in range(2)]
            for path in paths:
                os.mkfifo(path)
            later_received = threading.Event()
            def stalled_reader():
                with paths[0].open('rb') as stream:
                    if not later_received.wait(5):
                        raise AssertionError('Later spatial batch could not advance')
                    return stream.read()
            def later_reader():
                with paths[1].open('rb') as stream:
                    result = stream.read(8)
                    later_received.set()
                    return result + stream.read()
            first = pool.submit(stalled_reader)
            second = pool.submit(later_reader)
            payload = np.ones(32768, dtype='<f8').tobytes()
            with InstanceFanout(['id'], 'id',
                                [dict(ids=[i+1], pipe=str(p)) for i,p in enumerate(paths)],
                                prefetch_bytes=len(payload) + 8) as sink:
                sink.write(payload)
                sink.write(np.array([2], dtype='<f8').tobytes())
                sink.flush()
                self.assertLessEqual(sink.peak_queued_bytes, sink.capacity)
                self.assertEqual(sink.bytes_sent, [len(payload), 8])
            self.assertEqual(first.result(timeout=5), payload)
            self.assertEqual(second.result(timeout=5), np.array([2], dtype='<f8').tobytes())

    def test_small_queue_applies_backpressure_without_losing_rows(self):
        with tempfile.TemporaryDirectory() as directory, ThreadPoolExecutor() as pool:
            path = Path(directory)/'bounded.fifo'
            os.mkfifo(path)
            future = pool.submit(path.read_bytes)
            payload = np.ones(32768, dtype='<f8').tobytes()
            with InstanceFanout(['id'], 'id', [dict(ids=[1], pipe=str(path))],
                                prefetch_bytes=1024) as sink:
                sink.write(payload)
                sink.flush()
                self.assertLessEqual(sink.peak_queued_bytes, 1024)
            self.assertEqual(future.result(timeout=5), payload)

    def test_regular_file_cannot_be_used_as_point_transport(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'not-a-pipe'
            path.write_bytes(b'keep')
            with self.assertRaisesRegex(ValueError, 'FIFO'):
                with InstanceFanout(['id'], 'id', [dict(ids=[1], pipe=str(path))]):
                    pass
            self.assertEqual(path.read_bytes(), b'keep')

    def test_missing_consumer_has_bounded_startup_wait(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'unopened.fifo'
            os.mkfifo(path)
            with self.assertRaises(TimeoutError):
                with InstanceFanout(['id'], 'id', [dict(ids=[1], pipe=str(path))], startup_timeout=0.05):
                    pass

    def test_dead_consumer_fails_without_a_success_audit(self):
        config = configuration()
        config.update(core_size=100, buffer=20)
        index = owned.build_index(config)
        layout = index['layouts'][0]
        with tempfile.TemporaryDirectory() as directory, ThreadPoolExecutor() as pool:
            path = Path(directory)/'dead.fifo'
            os.mkfifo(path)
            def fail():
                with path.open('rb') as stream:
                    stream.read(1)
            future = pool.submit(fail)
            audit = Path(directory)/'audit.json'
            with self.assertRaises(BrokenPipeError):
                dispatch(config, index, layout['slug'], 0,
                         [dict(ids=layout['windows'][0]['ids'], pipe=str(path))], audit)
            future.result(timeout=5)
            self.assertFalse(audit.exists())

    def test_duplicate_worker_ids_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'disjoint'):
            InstanceFanout(['id'], 'id', [dict(ids=[1]), dict(ids=[1])])


if __name__ == '__main__':
    unittest.main()
