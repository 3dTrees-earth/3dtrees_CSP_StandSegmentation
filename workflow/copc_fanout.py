"""Bounded, backpressured delivery of one COPC read to disjoint ID workers.

FIFOs contain only kernel buffers: no point payload is written to the filesystem.
The writer opens existing FIFOs only and never creates or truncates regular files.
"""
import errno
import json
import os
from pathlib import Path
from collections import deque
import select
import threading
import stat
import time

import numpy as np


class InstanceFanout:
    def __init__(self, columns, instance, groups, *, startup_timeout=60, prefetch_bytes=256 * 1024**2):
        if prefetch_bytes < 1:
            raise ValueError("Prefetch capacity must be positive")
        self.capacity = int(prefetch_bytes)
        self.queued_bytes = self.peak_queued_bytes = 0
        self.queues = [deque() for _ in groups]
        self.condition = threading.Condition()
        self.stopping = False
        self.error = None
        self.thread = None
        self.width = len(columns)
        self.id_column = columns.index(instance)
        pairs = sorted((int(key), number) for number, group in enumerate(groups) for key in group['ids'])
        if not pairs or len({key for key, _ in pairs}) != len(pairs):
            raise ValueError('Worker IDs must be nonempty and disjoint')
        self.ids = np.asarray([key for key, _ in pairs], dtype=np.float64)
        self.owners = np.asarray([owner for _, owner in pairs])
        self.groups = groups
        self.timeout = startup_timeout
        self.fds = []
        self.bytes_sent = [0] * len(groups)

    def __enter__(self):
        deadline = time.monotonic() + self.timeout
        try:
            for group in self.groups:
                path = group['pipe']
                if not stat.S_ISFIFO(os.lstat(path).st_mode):
                    raise ValueError('Point transport must be a FIFO, never a regular file')
                while True:
                    try:
                        fd = os.open(path, os.O_WRONLY | os.O_NONBLOCK | os.O_NOFOLLOW)
                        break
                    except OSError as error:
                        if error.errno != errno.ENXIO:
                            raise
                        if time.monotonic() >= deadline:
                            raise TimeoutError('Inventory worker did not open its point pipe') from error
                        time.sleep(0.02)
                self.fds.append(fd)
                if not stat.S_ISFIFO(os.fstat(fd).st_mode):
                    raise ValueError('Point transport changed while opening')
        except BaseException:
            self.__exit__(None, None, None)
            raise
        self.thread = threading.Thread(target=self._drain, name="copc-fanout", daemon=True)
        self.thread.start()
        return self

    def __exit__(self, *_):
        with self.condition:
            self.stopping = True
            self.condition.notify_all()
        if self.thread is not None:
            self.thread.join()
        for fd in self.fds:
            os.close(fd)
        self.fds.clear()
        self.queues.clear()

    def _check_error(self):
        if self.error is not None:
            raise self.error

    def _drain(self):
        # Independent writer overlaps FIFO delivery with decoding subsequent
        # spatial batches. Account whole byte buffers until fully released.
        try:
            while True:
                with self.condition:
                    self.condition.wait_for(lambda: self.stopping or self.queued_bytes)
                    if self.stopping:
                        return
                    active = {fd: i for i, fd in enumerate(self.fds) if self.queues[i]}
                _, ready, _ = select.select([], list(active), [], 0.1)
                for fd in ready:
                    i = active[fd]
                    with self.condition:
                        data, offset = self.queues[i][0]
                    try:
                        sent = os.write(fd, memoryview(data)[offset:offset + 65536])
                    except BlockingIOError:
                        continue
                    if not sent:
                        raise BrokenPipeError('Inventory worker stopped consuming points')
                    with self.condition:
                        self.bytes_sent[i] += sent
                        if offset + sent == len(data):
                            self.queues[i].popleft()
                            self.queued_bytes -= len(data)
                            self.condition.notify_all()
                        else:
                            self.queues[i][0] = (data, offset + sent)
        except BaseException as error:
            with self.condition:
                self.error = error
                self.condition.notify_all()

    def write(self, payload):
        rows = np.frombuffer(payload, dtype='<f8').reshape(-1, self.width)
        position = np.searchsorted(self.ids, rows[:, self.id_column])
        if np.any(position >= len(self.ids)) or np.any(self.ids[position] != rows[:, self.id_column]):
            raise ValueError('Point ID has no inventory worker')
        owners = self.owners[position]
        # Selection preserves original order within each instance. Besides the
        # bounded queue, only this decoded batch and one routing copy are live.
        chunk_bytes = min(self.capacity, 1024**2)
        for i in range(len(self.fds)):
            selected = rows[owners == i].tobytes()
            for start in range(0, len(selected), chunk_bytes):
                data = selected[start:start + chunk_bytes]
                with self.condition:
                    self.condition.wait_for(lambda: self.error is not None or
                                            self.queued_bytes + len(data) <= self.capacity)
                    self._check_error()
                    self.queues[i].append((data, 0))
                    self.queued_bytes += len(data)
                    self.peak_queued_bytes = max(self.peak_queued_bytes, self.queued_bytes)
                    self.condition.notify_all()
        return len(payload)

    def flush(self):
        with self.condition:
            self.condition.wait_for(lambda: self.error is not None or not self.queued_bytes)
            self._check_error()


def dispatch(config, index, slug, window, groups, audit_path):
    from copc_owned_inventory import stream_window
    layout = next(item for item in index['layouts'] if item['slug'] == slug)
    expected = set(layout['windows'][window]['ids'])
    if {key for group in groups for key in group['ids']} != expected:
        raise ValueError('Worker assignment does not cover the complete owned window')
    audit_path = Path(audit_path)
    if audit_path.exists():
        raise ValueError('Refusing to reuse an existing window success audit')
    with InstanceFanout(layout['columns'], layout['instance'], groups,
                        prefetch_bytes=int(config.get('prefetch_mb', 256)) * 1024**2) as sink:
        audit = stream_window(config, index, slug, window, sink)
        audit['prefetch_capacity_bytes'] = sink.capacity
        audit['prefetch_peak_bytes'] = sink.peak_queued_bytes
        audit['worker_count'] = len(groups)
        audit['worker_stream_bytes'] = sink.bytes_sent
        # Publish success before EOF so each consumer can verify global success.
        audit_path.write_text(json.dumps(audit, allow_nan=False) + '\n')
    return audit
