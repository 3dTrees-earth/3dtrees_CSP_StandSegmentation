#!/usr/bin/env python3
"""Read one merged COPC once and write bounded instance partition files.

The binary layout matches workflow/tool.R: little-endian float64 rows with
the configured columns in order. Only convex-hull vertices are retained for
the stand area, including points excluded from tree inventories.
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
from scipy.spatial import ConvexHull, QhullError

MAX_TRACKED_INSTANCE_IDS = 1_000_000
from inventory_policy import INSTANCE_POINT_LIMIT as STANDARD_INSTANCE_POINT_LIMIT


def update_hull(previous: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """The hull of a union equals the hull of its component hulls."""
    coordinates = np.column_stack((x, y))
    if len(previous):
        coordinates = np.concatenate((previous, coordinates))
    if len(coordinates) <= 2:
        return np.unique(coordinates, axis=0)
    try:
        return coordinates[ConvexHull(coordinates).vertices]
    except QhullError:
        # A collinear cloud has zero area; retain its extreme endpoints.
        ordered = np.lexsort((coordinates[:, 1], coordinates[:, 0]))
        return coordinates[ordered[[0, -1]]]


STANDARD_ATTRIBUTES = json.loads(Path(__file__).with_name("inventory_standard_attributes.json").read_text())


def requested_dimensions(header, required):
    """Resolve explicit R/LAS names and reject missing/vector attributes before reading."""
    available = set(header.point_format.dimension_names)
    result = {}
    for name in required:
        native = name if name in available else STANDARD_ATTRIBUTES.get(name, {}).get("python", name)
        if native not in available:
            raise ValueError(f"Missing requested LAS dimension: {name}")
        if header.point_format.dimension_by_name(native).num_elements != 1:
            raise ValueError(f"Requested LAS dimension must be scalar: {name}")
        result[name] = native
    return result


def decompression_groups(header, required: set[str]) -> tuple[object, list[str]]:
    """Select all compressed groups that contain a requested LAS dimension."""
    names = {name.lower() for name in requested_dimensions(header, required).values()}
    groups = {
        "XY_RETURNS_CHANNEL": {"x", "y", "return_number", "number_of_returns",
                               "scanner_channel", "scan_direction_flag", "edge_of_flight_line"},
        "Z": {"z"},
        "CLASSIFICATION": {"classification"},
        "FLAGS": {"synthetic", "key_point", "withheld", "overlap"},
        "INTENSITY": {"intensity"},
        "SCAN_ANGLE": {"scan_angle"},
        "USER_DATA": {"user_data"},
        "POINT_SOURCE_ID": {"point_source_id"},
        "GPS_TIME": {"gps_time"},
        "RGB": {"red", "green", "blue"},
        "NIR": {"nir"},
        "WAVEPACKET": {"wavepacket_index", "wavepacket_offset", "wavepacket_size",
                       "return_point_wave_location", "x_t", "y_t", "z_t"},
        "ALL_EXTRA_BYTES": {name.lower() for name in header.point_format.extra_dimension_names},
    }
    selected_names = [name for name, dimensions in groups.items() if names & dimensions]
    selection = laspy.DecompressionSelection(0)
    for name in selected_names:
        selection |= getattr(laspy.DecompressionSelection, name)
    return selection, selected_names


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
    instance_counts: dict[str, dict[float, int]] = {
        layout["slug"]: {} for layout in layouts
    }
    read_seconds = 0.0
    write_seconds = 0.0
    hull_seconds = 0.0
    started = time.monotonic()
    hull_path = Path(config["hull_path"])
    if hull_path.exists():
        raise ValueError("Hull output already exists")
    hull = np.empty((0, 2), dtype=np.float64)
    required = {"X", "Y", "Z"}
    for layout in layouts:
        required.update(layout["columns"])
    with laspy.open(input_path, laz_backend=laspy.LazBackend.Lazrs) as header_reader:
        header = header_reader.header
        point_format_id = header.point_format.id
        selection, selected_groups = decompression_groups(header, required)
    selective = point_format_id >= 6 and config.get("selective_decompression", True)
    reader_options = {"decompression_selection": selection} if selective else {}
    try:
        with laspy.open(input_path, laz_backend=laspy.LazBackend.Lazrs,
                        **reader_options) as reader:
            header = reader.header
            if not any(v.user_id == "copc" and v.record_id == 1 for v in header.vlrs):
                raise ValueError("Sequential inventory requires a COPC input")
            dimension_names = requested_dimensions(header, required)
            while total_points < header.point_count:
                read_started = time.monotonic()
                points = reader.read_points(min(chunk_points, header.point_count - total_points))
                read_seconds += time.monotonic() - read_started
                if not len(points):
                    break
                write_started = time.monotonic()
                x = np.asarray(points.x, dtype=np.float64)
                y = np.asarray(points.y, dtype=np.float64)
                hull_started = time.monotonic()
                hull = update_hull(hull, x, y)
                hull_seconds += time.monotonic() - hull_started
                fields = {"X": x, "Y": y, "Z": np.asarray(points.z, dtype=np.float64)}
                for name in required - fields.keys():
                    fields[name] = np.asarray(points[dimension_names[name]], dtype=np.float64)
                x_index = np.floor((x - grid_x) / grid_size).astype(np.int32)
                y_index = np.floor((y - grid_y) / grid_size).astype(np.int32)
                if np.any((x_index < 0) | (x_index >= grid_nx) |
                          (y_index < 0) | (y_index >= grid_ny)):
                    raise ValueError("COPC point lies outside the configured spatial grid")
                spatial_index = y_index * grid_nx + x_index
                for layout in layouts:
                    instance = layout["instance"]
                    ids = fields[instance]
                    valid = np.isfinite(ids) & (ids >= 0)
                    if non_tree.size:
                        valid &= ~np.isin(ids, non_tree)
                    if not np.any(valid):
                        continue
                    selected_ids = ids[valid]
                    unique_ids, counts = np.unique(selected_ids, return_counts=True)
                    totals = instance_counts[layout["slug"]]
                    for instance_id, count in zip(unique_ids, counts):
                        key = float(instance_id)
                        totals[key] = totals.get(key, 0) + int(count)
                    if sum(map(len, instance_counts.values())) > MAX_TRACKED_INSTANCE_IDS:
                        raise ValueError(
                            "More than 1,000,000 distinct instance IDs; "
                            "the bounded scan cannot retain exact preflight counts"
                        )
                    partition_ids = (np.mod(np.abs(selected_ids), partitions).astype(np.int32) + 1)
                    keys = spatial_index[valid].astype(np.int64) * partitions + partition_ids - 1
                    order = np.argsort(keys, kind="stable")
                    ordered_keys = keys[order]
                    selected_rows = np.flatnonzero(valid)[order]
                    values = np.column_stack(
                        tuple(fields[name][selected_rows] for name in layout["columns"])
                    )
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
            loaded = sorted(required)
            dropped = [name for name in header.point_format.dimension_names if name not in dimension_names.values()]
    finally:
        for writer in handles.values():
            writer.close()

    hull.astype("<f8", copy=False).tofile(hull_path)
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

    largest_instances = {}
    oversized_instances = {}
    for slug, totals in instance_counts.items():
        if totals:
            instance_id, count = max(totals.items(), key=lambda item: item[1])
            largest_instances[slug] = {"instance_id": instance_id, "points": count}
        oversized_instances[slug] = [
            {"instance_id": instance_id, "points": count}
            for instance_id, count in sorted(totals.items())
            if count > STANDARD_INSTANCE_POINT_LIMIT
        ]
    (output_dir / "partition_summary.json").write_text(json.dumps({
        "largest_instances": largest_instances,
        "oversized_instances": oversized_instances,
        "instance_counts": {slug: len(totals) for slug, totals in instance_counts.items()},
    }, indent=2) + "\n")

    return {
        "input": str(input_path), "input_points": total_points,
        "source_point_reads": 1,
        "selector": "python_sequential_selective" if selective else "python_sequential_all",
        "assigned_points": total_assigned, "xy_spool_bytes": 0,
        "largest_instances": largest_instances,
        "oversized_instances": oversized_instances,
        "hull_vertex_count": len(hull),
        "spatial_bucket_bytes": bucket_bytes,
        "loaded_dimensions": loaded, "dropped_dimensions": dropped,
        "decompressed_groups": selected_groups if selective else ["ALL"],
        "input_read_seconds": read_seconds,
        "hull_seconds": hull_seconds,
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
