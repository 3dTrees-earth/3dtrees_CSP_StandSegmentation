"""Experimental bounded COPC indexing; one ordered reader, summary workers.

Only compressed bytes cross process boundaries. Decoded records stay in their
worker, which returns integer summaries and hull vertices. The normal index
path remains the default until comparative benchmarks justify changing it.
"""
from array import array
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import multiprocessing
import os
import sys
import time

import laspy
import lazrs
import numpy as np
from stream_inventory_partitions import update_hull, decompression_groups

MAX_NODES = 2_000_000  # compact table plus sort/copy workspace: at most 128 MB
MAX_HULL_VERTICES = 1_000_000
MAX_COMPRESSED_BYTES = 64 * 1024**2
_WORKER = None


def ordered_nodes(reader, audit):
    from copc_owned_inventory import hierarchy_nodes
    data = array('Q')
    source_bytes = os.fstat(reader.source.fileno()).st_size
    for node in hierarchy_nodes(reader, None, max_pages=100_000, max_page_bytes=8 * 1024**2):
        if len(data) // 3 >= MAX_NODES:
            raise ValueError('COPC hierarchy exceeds indexing metadata budget')
        if node.offset < 0 or not 0 < node.byte_size <= MAX_COMPRESSED_BYTES:
            raise ValueError('Invalid or oversized COPC data range')
        if node.offset + node.byte_size > source_bytes:
            raise ValueError('COPC data range extends beyond source')
        data.extend((node.offset, node.byte_size, node.point_count))
    nodes = np.frombuffer(data, dtype=np.uint64).reshape(-1, 3)
    nodes = nodes[np.argsort(nodes[:, 0], kind='stable')]
    if len(nodes) and np.any(nodes[1:, 0] < nodes[:-1, 0] + nodes[:-1, 1]):
        raise ValueError('Overlapping or duplicate COPC nodes')
    total = sum(map(int, nodes[:, 2]))
    if total != reader.header.point_count:
        raise ValueError('COPC hierarchy count disagrees with the LAS header')
    audit.update(hierarchy_points=total, node_count=len(nodes),
                 metadata_table_bytes=nodes.nbytes, max_hierarchy_nodes=MAX_NODES,
                 flat_metadata_workspace_budget_bytes=64 * MAX_NODES,
                 max_hierarchy_pages=100_000, max_retained_hierarchy_page_payload_bytes=8 * 1024**2)
    return nodes


def batches(nodes, point_limit, compressed_limit=MAX_COMPRESSED_BYTES):
    start = count = size = 0
    for i, (_, byte_size, points) in enumerate(nodes):
        points, byte_size = int(points), int(byte_size)
        if byte_size > compressed_limit:
            raise ValueError('COPC node exceeds compressed queue budget')
        if points > point_limit:
            raise ValueError(f'COPC node has {points} points, above batch limit {point_limit}')
        if i > start and (count + points > point_limit or size + byte_size > compressed_limit or i-start >= 4096):
            yield nodes[start:i], count, size
            start, count, size = i, 0, 0
        count += points
        size += byte_size
    if len(nodes):
        yield nodes[start:], count, size


def read_batch(source, nodes, size):
    """Coalesce adjacent chunks and reject short reads before dispatch."""
    payload = bytearray(size)
    cursor = 0
    groups = []
    for offset, length, _ in nodes:
        offset, length = int(offset), int(length)
        if groups and groups[-1][0] + groups[-1][1] == offset:
            groups[-1][1] += length
        else:
            groups.append([offset, length])
    for offset, length in groups:
        source.seek(offset)
        view = memoryview(payload)[cursor:cursor+length]
        received = 0
        while received < length:
            n = source.readinto(view[received:])
            if not n:
                raise ValueError('Truncated COPC compressed data')
            received += n
        cursor += length
    table = [(int(n), int(b)) for _, b, n in nodes]
    return payload, table


def initialize_worker(header, vlr, selection, layouts, non_tree):
    global _WORKER
    _WORKER = header, vlr, selection.to_lazrs(), layouts, non_tree


def summarize_compressed(payload, table, expected):
    from copc_owned_inventory import summarize_batch, MAX_IDS
    header, vlr, selection, layouts, non_tree = _WORKER
    if len(payload) != sum(size for _, size in table) or expected != sum(n for n, _ in table):
        raise ValueError("Incomplete compressed COPC batch")
    decoded = np.zeros(expected * header.point_format.size, dtype=np.uint8)
    lazrs.decompress_points_with_chunk_table(payload, vlr, decoded, table, selection)
    record = laspy.PackedPointRecord.from_buffer(decoded, header.point_format)
    points = laspy.ScaleAwarePointRecord(record.array, record.point_format, header.scales, header.offsets)
    if len(points) != expected:
        raise ValueError('COPC node decompression returned an incomplete batch')
    states = {}
    for slug, dimension in layouts:
        states[slug] = {}
        summarize_batch(states[slug], points, dimension, non_tree)
        if sum(map(len, states.values())) > MAX_IDS:
            raise ValueError('Too many distinct instance IDs across layouts')
    hull = update_hull(np.empty((0, 2)), np.asarray(points.x), np.asarray(points.y))
    if len(hull) > MAX_HULL_VERTICES:
        raise ValueError('Stand hull exceeds indexing metadata budget')
    return states, hull, expected


def merge_summaries(target, partial):
    from copc_owned_inventory import MAX_IDS
    for slug, instances in partial.items():
        for key, incoming in instances.items():
            if key not in target[slug]:
                target[slug][key] = incoming
            else:
                state = target[slug][key]
                for field in ('count', 'sum_x_raw', 'sum_y_raw'):
                    state[field] += incoming[field]  # Python integers cannot overflow
                state['mins_raw'] = list(map(min, zip(state['mins_raw'], incoming['mins_raw'])))
                state['maxs_raw'] = list(map(max, zip(state['maxs_raw'], incoming['maxs_raw'])))
        if sum(map(len, target.values())) > MAX_IDS:
            raise ValueError('Too many distinct instance IDs across layouts')


def index_summaries(reader, config, dimension_names):
    workers = int(config['index_workers'])
    budget = int(config.get('index_queue_mb', 256)) * 1024**2
    if not 1 <= workers <= 8 or budget < 32 * 1024**2 or budget > 1024**3:
        raise ValueError('Index workers must be 1..8 and queue budget 32..1024 MiB')
    if not 1 <= config['batch_points'] <= 10_000_000:
        raise ValueError('COPC batch_points must be between 1 and 10000000')
    # Spawned processes initialize native pools from these values. No nested
    # decoder/BLAS parallelism: at most workers + one reader/merger CPU.
    for name in ('RAYON_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[name] = '1'
    audit = dict(decoded_points=0, max_batch_points=0, index_workers=workers,
                 queue_budget_bytes=budget, peak_queued_compressed_bytes=0,
                 peak_inflight_batches=0, read_seconds=0., merge_seconds=0.,
                 compressed_bytes=0, native_threads_per_worker=1)
    nodes = ordered_nodes(reader, audit)
    limit = min(config['batch_points'], max(1, 256 * 1024**2 // reader.header.point_format.size))
    layouts = [(x['slug'], dimension_names[x['instance']]) for x in config['layouts']]
    initialization = (reader.header, reader.laszip_vlr.record_data,
                      decompression_groups(reader.header, {"X", "Y", "Z"} | {name for _, name in layouts})[0],
                      layouts, np.atleast_1d(config['non_tree_ids']))
    states = {slug: {} for slug, _ in layouts}
    hull = np.empty((0, 2))
    def merge(result):
        nonlocal hull
        started = time.monotonic()
        partial, vertices, count = result
        merge_summaries(states, partial)
        hull = update_hull(hull, vertices[:, 0], vertices[:, 1])
        if len(hull) > MAX_HULL_VERTICES:
            raise ValueError('Stand hull exceeds indexing metadata budget')
        audit['decoded_points'] += count
        audit['max_batch_points'] = max(audit['max_batch_points'], count)
        audit['merge_seconds'] += time.monotonic() - started
        print(f"copc_index={audit['decoded_points']}/{reader.header.point_count}", file=sys.stderr, flush=True)
    def read(nodes, size):
        started = time.monotonic()
        result = read_batch(reader.source, nodes, size)
        audit['read_seconds'] += time.monotonic() - started
        audit['compressed_bytes'] += size
        return result
    if workers == 1:
        initialize_worker(*initialization)
        for batch, count, size in batches(nodes, limit, min(MAX_COMPRESSED_BYTES, budget)):
            payload, table = read(batch, size)
            audit['peak_queued_compressed_bytes'] = max(audit['peak_queued_compressed_bytes'], size)
            audit['peak_inflight_batches'] = 1
            merge(summarize_compressed(payload, table, count))
    else:
        pending, queued = {}, 0
        with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context('spawn'),
                                 initializer=initialize_worker, initargs=initialization) as pool:
            def drain():
                nonlocal queued
                done, _ = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    queued -= pending.pop(future)
                    merge(future.result())
            for batch, count, size in batches(nodes, limit, min(MAX_COMPRESSED_BYTES, budget)):
                while pending and (len(pending) >= workers or queued + size > budget):
                    drain()
                payload, table = read(batch, size)
                future = pool.submit(summarize_compressed, payload, table, count)
                pending[future] = size
                queued += size
                audit['peak_queued_compressed_bytes'] = max(audit['peak_queued_compressed_bytes'], queued)
                audit['peak_inflight_batches'] = max(audit['peak_inflight_batches'], len(pending))
                del payload
            while pending:
                drain()
    return states, hull, audit
