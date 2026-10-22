"""Validate experimental worker reduction against the existing serial index."""
import copy
from io import BytesIO
import sys
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'workflow'))
from test_copc_owned_inventory import configuration
import copc_owned_inventory as owned
import copc_parallel_index as parallel


def canonical(index):
    result = {key: index[key] for key in ('input_points', 'layouts', 'source_identity', 'available_dimensions', 'decompressed_groups')}
    result['stand_hull'] = sorted(map(tuple, index['stand_hull']))
    result['raw_summary_sha256'] = index['audit']['raw_summary_sha256']
    return result


class ParallelTests(unittest.TestCase):
    def test_serial_and_workers_preserve_raw_summaries_hull_and_windows(self):
        config = configuration()
        expected = canonical(owned.build_index(copy.deepcopy(config)))
        for workers in (1, 2, 4, 8):
            with self.subTest(workers=workers):
                current = dict(copy.deepcopy(config), index_workers=workers, index_queue_mb=32)
                actual = owned.build_index(current)
                self.assertEqual(canonical(actual), expected)
                self.assertLessEqual(actual['audit']['peak_inflight_batches'], workers)
                self.assertLessEqual(actual['audit']['peak_queued_compressed_bytes'], 32*1024**2)
                for layout in actual['layouts']:
                    for number in range(len(layout['windows'])):
                        audit = owned.stream_window(current, actual, layout['slug'], number, BytesIO())
                        self.assertTrue(audit['complete'])

    def test_duplicate_ranges_header_count_and_truncated_input_rejected(self):
        config = dict(configuration(), index_workers=2)
        load = owned.hierarchy_nodes
        for kind, error in [('duplicate', 'Overlapping'), ('missing', 'header')]:
            def broken(*args, **kwargs):
                nodes = list(load(*args, **kwargs))
                return iter(nodes + nodes[:1] if kind=='duplicate' else nodes[:-1])
            with patch.object(owned, 'hierarchy_nodes', broken):
                with self.assertRaisesRegex(ValueError, error):
                    owned.build_index(copy.deepcopy(config))
        with self.assertRaisesRegex(ValueError, 'Truncated'):
            parallel.read_batch(BytesIO(b'abc'), np.array([[0, 4, 1]]), 4)

    def test_large_compressed_node_respects_batch_and_queue_limits(self):
        mib = 1024**2
        nodes = np.array([[0, 40*mib, 500000], [40*mib, 20*mib, 500000],
                          [60*mib, 20*mib, 500000]], dtype=np.uint64)
        batches = list(parallel.batches(nodes, 2000000))
        self.assertEqual([size for _, _, size in batches], [60*mib, 20*mib])
        self.assertEqual(sum(count for _, count, _ in batches), 1500000)
        with self.assertRaisesRegex(ValueError, 'queue budget'):
            list(parallel.batches(nodes, 2000000, 32*mib))

    def test_worker_failure_does_not_return_index(self):
        config = dict(configuration(), index_workers=2)
        # Corrupt each real compressed payload; the spawned worker must raise.
        read = parallel.read_batch
        def corrupt(*args):
            payload, table = read(*args)
            return payload[:1], table
        with patch.object(parallel, 'read_batch', corrupt):
            with self.assertRaises(Exception):
                owned.build_index(config)

    def test_order_independent_integer_merge_and_collinear_hull(self):
        parts = [dict(FM={str(i % 2): dict(id=i % 2, count=1, sum_x_raw=2**62,
                 sum_y_raw=-2**62, mins_raw=[i, i, i], maxs_raw=[i, i, i])}) for i in range(8)]
        a, b = {'FM': {}}, {'FM': {}}
        for p in parts: parallel.merge_summaries(a, copy.deepcopy(p))
        for p in reversed(parts): parallel.merge_summaries(b, copy.deepcopy(p))
        self.assertEqual(a, b)
        self.assertEqual(a['FM']['0']['sum_x_raw'], 2**64)
        hull = parallel.update_hull(np.empty((0,2)), np.arange(10), np.arange(10))
        np.testing.assert_array_equal(hull, [[0,0], [9,9]])
        merged = parallel.update_hull(hull, np.array([-1, 4]), np.array([-1, 4]))
        np.testing.assert_array_equal(merged, [[-1,-1], [9,9]])

    def test_source_mutation_and_metadata_budget_rejected(self):
        config = dict(configuration(), index_workers=1)
        original = owned.identity(config['input'])
        with patch.object(owned, 'identity', side_effect=[original, dict(original, mtime_ns=0)]):
            with self.assertRaisesRegex(ValueError, 'changed during'):
                owned.build_index(copy.deepcopy(config))
        with patch.object(parallel, 'MAX_NODES', 1):
            with self.assertRaisesRegex(ValueError, 'metadata budget'):
                owned.build_index(copy.deepcopy(config))

if __name__ == '__main__':
    unittest.main()
