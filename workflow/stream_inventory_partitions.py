#!/usr/bin/env python3
"""Read one merged COPC once and write bounded instance partition files.

The binary layout matches workflow/tool.R: little-endian float64 rows with
the configured columns in order. An XY spool retains all points for an exact
streaming convex hull in R, including points excluded from tree inventories.
"""

from __future__ import annotations

import argparse
from collections import OrderedDict
import json
from pathlib import Path
import resource
import shutil
import time

import laspy
import numpy as np


class PartitionFiles:
    def __init__(self, directory: Path, partitions: int, max_open: int = 256):
        self.directory = directory
        self.partitions = partitions
        self.max_open = max_open
        self.files: OrderedDict[int, object] = OrderedDict()

    def write(self, key: int, values: np.ndarray) -> None:
        if key in self.files:
            file = self.files[key]
            self.files.move_to_end(key)
        else:
            if len(self.files) >= self.max_open:
                self.files.popitem(last=False)[1].close()
            chunk_index, partition_index = divmod(key, self.partitions)
            path = self.directory / f"{chunk_index + 1:05d}_{partition_index + 1:05d}.bin"
            file = path.open("ab")
            self.files[key] = file
        file.write(values.astype("<f8", copy=False).tobytes(order="C"))

    def close(self) -> None:
        for file in self.files.values():
            file.close()
        self.files.clear()


def process(config: dict) -> dict:
    input_path = Path(config["input"])
    output_dir = Path(config["output_dir"])
    if not input_path.is_file() or output_dir.exists():
        raise ValueError("Input must exist and partition output directory must be new")
    partitions = int(config["partitions"])
    chunk_points = int(config.get("chunk_points", 1_000_000))
    grid = config["spatial_grid"]
    grid_x = float(grid["xmin"])
    grid_y = float(grid["ymin"])
    grid_size = float(grid["size"])
    grid_nx = int(grid["nx"])
    grid_ny = int(grid["ny"])
    if not 1 <= partitions <= 4096 or not 1 <= chunk_points <= 10_000_000:
        raise ValueError("Invalid partition count or point chunk size")
    if grid_size <= 0 or grid_nx < 1 or grid_ny < 1:
        raise ValueError("Invalid spatial grid")
    layouts = config["layouts"]
    if not layouts:
        raise ValueError("At least one segmentation layout is required")
    non_tree = np.asarray(config["non_tree_ids"], dtype=np.float64)
    output_dir.mkdir(parents=True)
    bucket_root = output_dir / ".spatial_buckets"
    bucket_root.mkdir()
    soft_limit = resource.getrlimit(resource.RLIMIT_NOFILE)[0]
    if soft_limit == resource.RLIM_INFINITY:
        soft_limit = 65536
    reserved_handles = min(64, soft_limit // 4)
    handles_per_layout = max(1, min(1024, (soft_limit - reserved_handles) // len(layouts)))
    handles = {}
    for layout in layouts:
        (output_dir / layout["slug"]).mkdir()
        bucket_dir = bucket_root / layout["slug"]
        bucket_dir.mkdir()
        handles[layout["slug"]] = PartitionFiles(
            bucket_dir, partitions, max_open=handles_per_layout
        )

    total_points = 0
    total_assigned = {}
    read_seconds = 0.0
    write_seconds = 0.0
    started = time.monotonic()
    xy_path = Path(config["xy_path"])
    if xy_path.exists():
        raise ValueError("XY spool already exists")
    try:
        with laspy.open(input_path, laz_backend=laspy.LazBackend.Lazrs) as reader:
            header = reader.header
            if not any(v.user_id == "copc" and v.record_id == 1 for v in header.vlrs):
                raise ValueError("Sequential inventory requires a COPC input")
            available = set(header.point_format.dimension_names)
            for layout in layouts:
                if layout["instance"] not in available:
                    raise ValueError(f"Missing instance dimension: {layout['instance']}")
            with xy_path.open("wb") as xy_file:
                while total_points < header.point_count:
                    read_started = time.monotonic()
                    points = reader.read_points(min(chunk_points, header.point_count - total_points))
                    read_seconds += time.monotonic() - read_started
                    if not len(points):
                        break
                    write_started = time.monotonic()
                    x = np.asarray(points.x, dtype=np.float64)
                    y = np.asarray(points.y, dtype=np.float64)
                    np.column_stack((x, y)).astype("<f8", copy=False).tofile(xy_file)
                    x_index = np.floor((x - grid_x) / grid_size).astype(np.int32)
                    y_index = np.floor((y - grid_y) / grid_size).astype(np.int32)
                    if np.any((x_index < 0) | (x_index >= grid_nx) |
                              (y_index < 0) | (y_index >= grid_ny)):
                        raise ValueError("COPC point lies outside the configured spatial grid")
                    spatial_index = y_index * grid_nx + x_index
                    for layout in layouts:
                        instance = layout["instance"]
                        ids = np.asarray(points[instance], dtype=np.float64)
                        valid = np.isfinite(ids) & (ids >= 0)
                        if non_tree.size:
                            valid &= ~np.isin(ids, non_tree)
                        if not np.any(valid):
                            continue
                        selected_ids = ids[valid]
                        partition_ids = (np.mod(np.abs(selected_ids), partitions).astype(np.int32) + 1)
                        keys = spatial_index[valid].astype(np.int64) * partitions + partition_ids - 1
                        order = np.argsort(keys, kind="stable")
                        ordered_keys = keys[order]
                        columns = []
                        for name in layout["columns"]:
                            if name == "X":
                                field = x
                            elif name == "Y":
                                field = y
                            elif name == "Z":
                                field = np.asarray(points.z, dtype=np.float64)
                            elif name in available:
                                field = np.asarray(points[name], dtype=np.float64)
                            else:
                                field = np.full(len(points), np.nan)
                            columns.append(field[valid][order])
                        values = np.column_stack(columns)
                        boundaries = np.r_[0, np.flatnonzero(np.diff(ordered_keys)) + 1, len(ordered_keys)]
                        writer = handles[layout["slug"]]
                        for start, end in zip(boundaries[:-1], boundaries[1:]):
                            writer.write(int(ordered_keys[start]), values[start:end])
                        slug = layout["slug"]
                        total_assigned[slug] = total_assigned.get(slug, 0) + len(selected_ids)
                    total_points += len(points)
                    write_seconds += time.monotonic() - write_started
                    print(f"inventory_scan={total_points}/{header.point_count}", flush=True)
            if total_points != header.point_count:
                raise RuntimeError(f"Read {total_points} of {header.point_count} COPC points")
            keep = {"X", "Y", "Z", "Classification"}
            for layout in layouts:
                keep.update(layout["columns"])
            loaded = [name for name in header.point_format.dimension_names if name in keep]
            dropped = [name for name in header.point_format.dimension_names if name not in keep]
    finally:
        for writer in handles.values():
            writer.close()

    bucket_bytes = 0
    for layout in layouts:
        slug = layout["slug"]
        bucket_dir = bucket_root / slug
        buckets: dict[int, list[tuple[int, Path]]] = {}
        for path in bucket_dir.glob("*.bin"):
            chunk_text, partition_text = path.stem.split("_")
            buckets.setdefault(int(partition_text), []).append((int(chunk_text), path))
            bucket_bytes += path.stat().st_size
        for partition, ordered_paths in buckets.items():
            destination = output_dir / slug / f"{partition:05d}.bin"
            with destination.open("wb") as output:
                for _, path in sorted(ordered_paths):
                    with path.open("rb") as source:
                        shutil.copyfileobj(source, output, length=8 * 1024 * 1024)
                    path.unlink()
    shutil.rmtree(bucket_root)

    return {
        "input": str(input_path), "input_points": total_points,
        "assigned_points": total_assigned, "xy_spool_bytes": xy_path.stat().st_size,
        "spatial_bucket_bytes": bucket_bytes,
        "loaded_dimensions": loaded, "dropped_dimensions": dropped,
        "input_read_seconds": read_seconds,
        "partition_write_seconds": write_seconds,
        "elapsed_seconds": time.monotonic() - started,
        "process_peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    args = parser.parse_args()
    result = process(json.loads(args.config.read_text()))
    args.metadata.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
