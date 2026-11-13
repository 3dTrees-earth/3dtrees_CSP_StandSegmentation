#!/usr/bin/env python3
"""Index complete instances, then stream owned COPC windows without point files.

The first pass stores counts, global centroids and extents, never point payloads.
Each instance belongs to exactly one core. A read window contains its 20 m
halo AND the full extents of all owned instances. COPC nodes are decompressed
in bounded batches; only owned IDs are sent to R as little-endian float64 rows.

Laspy 2.6.1's public query() materializes an entire query. The private decoder
call below is deliberately confined to iter_batches and covered by fixture
tests. Every hierarchy level is read; no LOD/resolution reduction is allowed.
"""
from __future__ import annotations

import argparse
from array import array
from collections import namedtuple
import json
import math
import os
from pathlib import Path
import resource
import sys
import time

import laspy
from laspy.copc import Bounds, CopcReader, HierarchyPage
import numpy as np

from stream_inventory_partitions import decompression_groups, update_hull, requested_dimensions



def identity(path):
    stat = Path(path).stat()
    return {"bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}


Node = namedtuple("Node", "offset byte_size point_count")


def hierarchy_nodes(reader, bounds):
    """Visit flat hierarchy pages without constructing the whole octree.

    Each hierarchy entry already includes its voxel key. Page references own
    their subtree; the referenced page replaces that reference. Only the DFS
    stack of pages and compact visited page offsets are retained.
    """
    root_bounds = Bounds(reader.copc_info.center - reader.copc_info.halfsize,
                         reader.copc_info.center + reader.copc_info.halfsize)
    pages = [iter(reader.root_page.entries.values())]
    seen_pages = {reader.copc_info.hierarchy_root_offset}
    while pages:
        entry = next(pages[-1], None)
        if entry is None:
            pages.pop()
            continue
        if bounds is not None and not entry.key.bounds(root_bounds).overlaps(bounds):
            continue
        if entry.point_count == -1:
            if entry.offset in seen_pages:
                raise ValueError("Duplicate or cyclic COPC hierarchy page")
            if entry.byte_size <= 0 or entry.byte_size % 32:
                raise ValueError("Invalid COPC hierarchy page")
            seen_pages.add(entry.offset)
            reader.source.seek(entry.offset)
            payload = reader.source.read(entry.byte_size)
            if len(payload) != entry.byte_size:
                raise ValueError("Truncated COPC hierarchy page")
            page = HierarchyPage.from_bytes(payload)
            pages.append(iter(page.entries.values()))
        elif entry.point_count > 0:
            yield Node(entry.offset, entry.byte_size, entry.point_count)
        elif entry.point_count < 0:
            raise ValueError("Invalid COPC point count")


def iter_batches(reader, bounds, batch_points, audit):
    """Read disjoint COPC nodes once, including ancestor-level points."""
    if not 1 <= batch_points <= 10_000_000:
        raise ValueError("COPC batch_points must be between 1 and 10000000")
    limit = min(batch_points, max(1, (256 * 1024**2) // reader.header.point_format.size))
    audit.update(hierarchy_points=0, node_count=0, decoded_points=0, max_batch_points=0)
    # Compact byte ranges detect duplicate/overlapping nodes without retaining
    # millions of OctreeNode objects and their child lists. Sixteen bytes/node.
    ranges = array("Q")
    batch, count, compressed_bytes = [], 0, 0

    def decode(nodes, expected):
        points = reader._fetch_and_decompress_points_of_nodes(nodes)
        if len(points) != expected:
            raise ValueError("COPC node decompression returned an incomplete batch")
        audit["decoded_points"] += expected
        audit["max_batch_points"] = max(audit["max_batch_points"], expected)
        return points

    for node in hierarchy_nodes(reader, bounds):
        if node.byte_size <= 0 or node.offset < 0:
            raise ValueError("Invalid COPC data range")
        if batch and (count + node.point_count > limit or len(batch) >= 4096 or
                      compressed_bytes + node.byte_size > 256 * 1024**2):
            yield decode(batch, count)
            batch, count, compressed_bytes = [], 0, 0
        ranges.extend((node.offset, node.offset + node.byte_size))
        audit["hierarchy_points"] += node.point_count
        audit["node_count"] += 1
        batch.append(node)
        count += node.point_count
        compressed_bytes += node.byte_size
    if batch:
        yield decode(batch, count)
    values = np.frombuffer(ranges, dtype=np.uint64).reshape(-1, 2)
    if len(values):
        ordered = values[np.argsort(values[:, 0])]
        if np.any(ordered[1:, 0] < ordered[:-1, 1]):
            raise ValueError("Overlapping or duplicate COPC nodes")


def selected_reader(path, required):
    with laspy.open(path) as reader:
        info = next((vlr for vlr in reader.header.vlrs
                     if vlr.user_id == "copc" and vlr.record_id == 1), None)
        if info is None or info.hierarchy_root_size <= 0:
            raise ValueError("Missing COPC info")
        selection, groups = decompression_groups(reader.header, required)
    return CopcReader.open(path, decompression_selection=selection), groups


def summarize_batch(states, points, instance, non_tree):
    ids = np.asarray(points[instance])
    valid = np.isfinite(ids) & (ids >= 0) & ~np.isin(ids, non_tree)
    selected = np.flatnonzero(valid)
    if not len(selected):
        return
    values = ids[selected]
    if np.any(values > 2**53) or np.any(values != np.floor(values)):
        raise ValueError("Instance IDs must be integers representable exactly as float64")
    order = np.argsort(values, kind="stable")
    selected = selected[order]
    values = values[order]
    starts = np.r_[0, np.flatnonzero(values[1:] != values[:-1]) + 1]
    counts = np.diff(np.r_[starts, len(values)])
    xyz = np.column_stack((points.X[selected], points.Y[selected], points.Z[selected]))
    minima = np.minimum.reduceat(xyz, starts)
    maxima = np.maximum.reduceat(xyz, starts)
    sums = np.add.reduceat(xyz[:, :2].astype(np.int64), starts)
    for i, start in enumerate(starts):
        key = str(int(values[start]))
        if key not in states:
            states[key] = {"id": int(values[start]), "count": 0,
                           "mins_raw": minima[i].tolist(), "maxs_raw": maxima[i].tolist(),
                           "sum_x_raw": 0, "sum_y_raw": 0}
        state = states[key]
        state["count"] += int(counts[i])
        state["sum_x_raw"] += int(sums[i, 0])
        state["sum_y_raw"] += int(sums[i, 1])
        state["mins_raw"] = [min(a, int(b)) for a, b in zip(state["mins_raw"], minima[i])]
        state["maxs_raw"] = [max(a, int(b)) for a, b in zip(state["maxs_raw"], maxima[i])]


def split_members(members, max_points):
    """Whole trees in groups of at most max_points (one oversized tree stays alone),
    split at the median point along the longer axis of the tree centroids."""
    if len(members) == 1 or sum(p["count"] for p in members) <= max_points:
        return [members]
    xs = [p["centroid"][0] for p in members]
    ys = [p["centroid"][1] for p in members]
    axis = 0 if max(xs) - min(xs) >= max(ys) - min(ys) else 1
    ordered = sorted(members, key=lambda p: (p["centroid"][axis], p["id"]))
    half, running, cut = sum(p["count"] for p in ordered) / 2, 0, 1
    for i, p in enumerate(ordered[:-1], start=1):
        running += p["count"]
        cut = i
        if running >= half:
            break
    return split_members(ordered[:cut], max_points) + split_members(ordered[cut:], max_points)


def plan_windows(instances, origin, size, buffer, max_points=None):
    """One window per core: its halo grown to the full extent of every owned tree.
    A core holding more than max_points is split into windows of whole trees, so
    the memory a window needs stays within the job's grant."""
    if not math.isfinite(size) or size <= 0 or not math.isfinite(buffer) or buffer < 0:
        raise ValueError("Invalid owned core size or buffer")
    cores = {}
    for instance in instances:
        core = tuple(math.floor((v - o) / size) for v, o in zip(instance["centroid"], origin))
        instance["owner_core"] = list(core)
        cores.setdefault(core, []).append(instance)
    windows = []
    for core, members in sorted(cores.items()):
        core_min = [o + c * size for o, c in zip(origin, core)]
        halo_min = [v - buffer for v in core_min]
        halo_max = [v + size + buffer for v in core_min]
        groups = split_members(members, max_points) if max_points else [members]
        for group in groups:
            group = sorted(group, key=lambda item: item["id"])
            if len(groups) == 1:
                minimum = [min(halo_min[axis], *(p["mins"][axis] for p in group)) for axis in (0, 1)]
                maximum = [max(halo_max[axis], *(p["maxs"][axis] for p in group)) for axis in (0, 1)]
            else:
                minimum = [min(p["mins"][axis] for p in group) for axis in (0, 1)]
                maximum = [max(p["maxs"][axis] for p in group) for axis in (0, 1)]
            windows.append({"owner_core": list(core), "ids": [p["id"] for p in group],
                            "counts": [p["count"] for p in group],
                            "mins": minimum, "maxs": maximum,
                            "expanded": minimum != halo_min or maximum != halo_max})
    return windows


def build_index(config):
    started = time.monotonic()
    path = config["input"]
    source_identity = identity(path)
    layouts = config["layouts"]
    required = {"X", "Y", "Z"} | {layout["instance"] for layout in layouts}
    reader, groups = selected_reader(path, required)
    states = {layout["slug"]: {} for layout in layouts}
    hull = np.empty((0, 2))
    audit = {}
    with reader:
        header = reader.header
        dimension_names = requested_dimensions(header, required | {name for layout in layouts for name in layout["columns"]})
        if config.get("index_workers", 0):
            from copc_parallel_index import index_summaries
            states, hull, audit = index_summaries(reader, config, dimension_names)
        else:
            for points in iter_batches(reader, None, config["batch_points"], audit):
                hull = update_hull(hull, np.asarray(points.x), np.asarray(points.y))
                for layout in layouts:
                    summarize_batch(states[layout["slug"]], points, dimension_names[layout["instance"]],
                                    np.atleast_1d(layout.get("non_tree_ids", [0])))
                print(f"copc_index={audit['decoded_points']}/{header.point_count}", file=sys.stderr, flush=True)
        if audit["decoded_points"] != header.point_count:
            raise ValueError("COPC hierarchy count disagrees with the LAS header")
        import hashlib
        audit["raw_summary_sha256"] = hashlib.sha256(json.dumps(states, sort_keys=True).encode()).hexdigest()
        for layout in layouts:
            instances = []
            for state in states[layout["slug"]].values():
                state["mins"] = (np.asarray(state.pop("mins_raw")) * header.scales + header.offsets).tolist()
                state["maxs"] = (np.asarray(state.pop("maxs_raw")) * header.scales + header.offsets).tolist()
                state["centroid"] = [state.pop(key) / state["count"] * header.scales[axis] + header.offsets[axis]
                                     for axis, key in enumerate(("sum_x_raw", "sum_y_raw"))]
                instances.append(state)
            layout["instances"] = sorted(instances, key=lambda item: item["id"])
            layout["windows"] = plan_windows(instances, header.mins[:2], config["core_size"],
                                                config["buffer"], config.get("max_window_points") or None)
        available = list(header.point_format.dimension_names)
        total = int(header.point_count)
    if identity(path) != source_identity:
        raise ValueError("Input COPC changed during indexing")
    return {"schema_version": 1, "source_identity": source_identity,
            "input_points": total, "available_dimensions": available,
            "decompressed_groups": groups, "layouts": layouts,
            "stand_hull": hull.tolist(), "audit": audit,
            "elapsed_seconds": time.monotonic() - started,
            "process_peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


_DECODER = None


def initialize_decoder(header, vlr, selection, columns, dimension_names, instance, ids):
    global _DECODER
    _DECODER = header, vlr, selection.to_lazrs(), columns, dimension_names, instance, ids


def decode_owned(payload, table, expected):
    """Decompress one batch of COPC nodes and keep the rows of owned IDs."""
    import lazrs
    header, vlr, selection, columns, dimension_names, instance, ids = _DECODER
    if len(payload) != sum(size for _, size in table) or expected != sum(n for n, _ in table):
        raise ValueError("Incomplete compressed COPC batch")
    decoded = np.zeros(expected * header.point_format.size, dtype=np.uint8)
    lazrs.decompress_points_with_chunk_table(payload, vlr, decoded, table, selection)
    record = laspy.PackedPointRecord.from_buffer(decoded, header.point_format)
    points = laspy.ScaleAwarePointRecord(record.array, record.point_format, header.scales, header.offsets)
    if len(points) != expected:
        raise ValueError("COPC node decompression returned an incomplete batch")
    return owned_rows(points, columns, dimension_names, instance, ids) + (expected,)


def owned_rows(points, columns, dimension_names, instance, ids):
    all_ids = np.asarray(points[dimension_names[instance]])
    selected = np.flatnonzero(np.isin(all_ids, ids))
    if not len(selected):
        return b"", {}
    kept_ids, counts = np.unique(all_ids[selected], return_counts=True)
    fields = []
    for column in columns:
        if column in ("X", "Y", "Z"):
            fields.append(np.asarray(getattr(points, column.lower()))[selected])
        else:
            fields.append(np.asarray(points[dimension_names[column]])[selected])
    payload = np.column_stack(fields).astype("<f8", copy=False).tobytes()
    return payload, {int(key): int(count) for key, count in zip(kept_ids, counts)}


def parallel_owned_batches(reader, bounds, config, columns, dimension_names, instance, ids, groups, audit):
    """Owned rows of a window, decoded by a process pool, delivered in node order.

    The main process reads compressed node bytes and hands them to decoders,
    which return only the owned rows; at most `workers` batches are in flight.
    """
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing
    from collections import deque
    from copc_parallel_index import batches, read_batch
    workers = int(config["decode_workers"])
    for name in ("RAYON_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[name] = "1"
    data = array("Q")
    for node in hierarchy_nodes(reader, bounds):
        if node.byte_size <= 0 or node.offset < 0:
            raise ValueError("Invalid COPC data range")
        data.extend((node.offset, node.byte_size, node.point_count))
    nodes = np.frombuffer(data, dtype=np.uint64).reshape(-1, 3)
    order = np.argsort(nodes[:, 0], kind="stable")
    nodes = nodes[order]
    if len(nodes) and np.any(nodes[1:, 0] < nodes[:-1, 0] + nodes[:-1, 1]):
        raise ValueError("Overlapping or duplicate COPC nodes")
    audit.update(hierarchy_points=int(nodes[:, 2].sum()) if len(nodes) else 0, node_count=len(nodes),
                 decoded_points=0, max_batch_points=0, decode_workers=workers)
    limit = min(config["batch_points"], max(1, (256 * 1024**2) // reader.header.point_format.size))
    initialization = (reader.header, reader.laszip_vlr.record_data,
                      decompression_groups(reader.header, set(columns))[0],
                      columns, dimension_names, instance, ids)
    pending = deque()
    with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("spawn"),
                             initializer=initialize_decoder, initargs=initialization) as pool:
        for batch, count, size in batches(nodes, limit):
            while len(pending) >= workers:
                yield pending.popleft().result()
            payload, table = read_batch(reader.source, batch, size)
            pending.append(pool.submit(decode_owned, payload, table, count))
            del payload
        while pending:
            yield pending.popleft().result()


def stream_window(config, index, slug, window_index, output):
    started = time.monotonic()
    path = config["input"]
    if identity(path) != index["source_identity"]:
        raise ValueError("Input COPC changed since indexing")
    layout = next(item for item in index["layouts"] if item["slug"] == slug)
    window = layout["windows"][window_index]
    ids = np.asarray(window["ids"])
    expected = dict(zip(map(int, window["ids"]), window["counts"]))
    observed = {key: 0 for key in expected}
    bounds = Bounds(np.array(window["mins"]), np.array(window["maxs"]))
    reader, groups = selected_reader(path, set(layout["columns"]))
    audit = {}
    with reader:
        if reader.header.point_count != index["input_points"]:
            raise ValueError("Input header changed since indexing")
        dimension_names = requested_dimensions(reader.header, layout["columns"])
        bounds = bounds.ensure_3d(reader.header.mins, reader.header.maxs)
        if np.all(bounds.mins <= reader.header.mins) and np.all(bounds.maxs >= reader.header.maxs):
            bounds = None
        if int(config.get("decode_workers", 1)) > 1:
            def owned():
                for payload, counts, decoded in parallel_owned_batches(
                        reader, bounds, config, layout["columns"], dimension_names, layout["instance"], ids, groups, audit):
                    audit["decoded_points"] += decoded
                    audit["max_batch_points"] = max(audit["max_batch_points"], decoded)
                    yield payload, counts
        else:
            def owned():
                for points in iter_batches(reader, bounds, config["batch_points"], audit):
                    yield owned_rows(points, layout["columns"], dimension_names, layout["instance"], ids)
        # The expanded window contains every owned ID's entire extent. Owned
        # points are kept even in nodes overlapping the window edge; nodes are
        # disjoint, so this cannot create duplicate points.
        for payload, counts in owned():
            if not counts:
                continue
            for key, count in counts.items():
                observed[key] += count
                if observed[key] > expected[key]:
                    raise ValueError(f"Duplicate/excess points for instance {key}")
            output.write(payload)
            print(f"copc_window={window_index} decoded={audit['decoded_points']} owned={sum(observed.values())}/{sum(expected.values())}", file=sys.stderr, flush=True)
        output.flush()
    if observed != expected:
        raise ValueError(f"Incomplete instances: expected {expected}, observed {observed}")
    if identity(path) != index["source_identity"]:
        raise ValueError("Input COPC changed during window read")
    return {"complete": True, "slug": slug, "window": window_index,
            "instance_ids": window["ids"], "expected_counts": window["counts"],
            "observed_counts": [observed[int(key)] for key in window["ids"]],
            "expanded": window["expanded"], "decompressed_groups": groups,
            **audit, "elapsed_seconds": time.monotonic() - started,
            "process_peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("index", "read", "dispatch"))
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--slug")
    parser.add_argument("--window", type=int)
    parser.add_argument("--audit", type=Path)
    parser.add_argument("--groups", type=Path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    if args.action == "index":
        result = build_index(config)
        args.index.write_text(json.dumps(result, allow_nan=False) + "\n")
    else:
        if args.slug is None or args.window is None or args.audit is None:
            parser.error("read requires --slug, --window and --audit")
        if args.action == "dispatch":
            if args.groups is None:
                parser.error("dispatch requires --groups")
            from copc_fanout import dispatch
            dispatch(config, json.loads(args.index.read_text()), args.slug, args.window,
                     json.loads(args.groups.read_text()), args.audit)
            return
        result = stream_window(config, json.loads(args.index.read_text()), args.slug,
                               args.window, sys.stdout.buffer)
        # The success record is written only after EOF and all count checks.
        args.audit.write_text(json.dumps(result, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
