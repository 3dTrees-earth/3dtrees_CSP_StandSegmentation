#!/usr/bin/env python3
"""Create a bounded, parallel low-surface candidate cloud for CSF/TIN."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
import json
import math
from pathlib import Path
import resource
import shutil
import time

import laspy
import numpy as np


def projection_records(header: laspy.LasHeader) -> list:
    """Return CRS VLRs without requiring the optional pyproj dependency."""
    records = list(header.vlrs)
    if header.evlrs is not None:
        records.extend(header.evlrs)
    return [record for record in records if record.user_id == "LASF_Projection"]


def projection_signature(header: laspy.LasHeader) -> tuple:
    return tuple(
        (record.user_id, record.record_id, record.record_data_bytes())
        for record in projection_records(header)
    )


def worker(
    input_path: str,
    parts_dir: str,
    worker_id: int,
    start: int,
    count: int,
    chunk_points: int,
    grid_min_x: float,
    grid_min_y: float,
    rows: int,
    cols: int,
    resolution: float,
) -> tuple[int, str, int, float]:
    started = time.monotonic()
    minimum_z = np.full((rows, cols), np.inf, dtype=np.float32)
    processed = 0

    with laspy.open(input_path, laz_backend=laspy.LazBackend.Lazrs) as reader:
        reader.seek(start)
        remaining = count
        while remaining:
            points = reader.read_points(min(chunk_points, remaining))
            if not len(points):
                break
            x = np.asarray(points.x)
            y = np.asarray(points.y)
            z = np.asarray(points.z, dtype=np.float32)
            col = np.floor((x - grid_min_x) / resolution).astype(np.int64)
            row = np.floor((y - grid_min_y) / resolution).astype(np.int64)
            np.clip(col, 0, cols - 1, out=col)
            np.clip(row, 0, rows - 1, out=row)
            cells = row * cols + col
            flat_z = minimum_z.ravel()
            np.minimum.at(flat_z, cells, z)
            processed += len(points)
            remaining -= len(points)

    if processed != count:
        raise RuntimeError(f"DTM worker {worker_id} read {processed} of {count} points")
    prefix = Path(parts_dir) / f"worker_{worker_id:03d}"
    z_path = f"{prefix}_z.npy"
    np.save(z_path, minimum_z)
    return (
        worker_id,
        z_path,
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        time.monotonic() - started,
    )


def write_candidates(
    output: Path,
    source_header: laspy.LasHeader,
    minimum_z: np.ndarray,
    grid_min_x: float,
    grid_min_y: float,
    resolution: float,
    row_block: int = 256,
) -> int:
    header = laspy.LasHeader(point_format=0, version="1.4")
    header.scales = np.array((0.001, 0.001, 0.001))
    header.offsets = np.array((grid_min_x, grid_min_y, 0.0))
    for record in projection_records(source_header):
        header.vlrs.append(deepcopy(record))
    written = 0
    with laspy.open(
        output,
        mode="w",
        header=header,
        do_compress=True,
        laz_backend=laspy.LazBackend.Lazrs,
    ) as writer:
        for start in range(0, minimum_z.shape[0], row_block):
            stop = min(start + row_block, minimum_z.shape[0])
            valid_rows, valid_cols = np.nonzero(np.isfinite(minimum_z[start:stop]))
            if not len(valid_rows):
                continue
            rows = valid_rows + start
            points = laspy.ScaleAwarePointRecord.zeros(len(rows), header=header)
            points.x = grid_min_x + (valid_cols.astype(np.float64) + 0.5) * resolution
            points.y = grid_min_y + (rows.astype(np.float64) + 0.5) * resolution
            points.z = minimum_z[rows, valid_cols]
            writer.write_points(points)
            written += len(points)
    return written


def process_file(
    input_path: Path,
    parts_dir: Path,
    workers: int,
    chunk_points: int,
    grid_min_x: float,
    grid_min_y: float,
    rows: int,
    cols: int,
    resolution: float,
) -> tuple[np.ndarray, int, int]:
    with laspy.open(input_path) as reader:
        point_count = int(reader.header.point_count)
    worker_count = min(workers, point_count)
    boundaries = np.linspace(0, point_count, worker_count + 1, dtype=np.int64)
    tasks = [
        (
            str(input_path),
            str(parts_dir),
            index,
            int(boundaries[index]),
            int(boundaries[index + 1] - boundaries[index]),
            chunk_points,
            grid_min_x,
            grid_min_y,
            rows,
            cols,
            resolution,
        )
        for index in range(worker_count)
    ]
    results = []
    with ProcessPoolExecutor(max_workers=worker_count) as executor:
        futures = [executor.submit(worker, *task) for task in tasks]
        completed = 0
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            worker_id = result[0]
            completed += tasks[worker_id][4]
            print(
                f"dtm_worker={worker_id} file={input_path.name} "
                f"completed={completed}/{point_count} elapsed_seconds={result[-1]:.3f}",
                flush=True,
            )

    minimum_z = np.full((rows, cols), np.inf, dtype=np.float32)
    peak_sum_kib = 0
    for _, z_path, peak_kib, _ in sorted(results):
        candidate_z = np.load(z_path, mmap_mode="r")
        np.minimum(minimum_z, candidate_z, out=minimum_z)
        peak_sum_kib += peak_kib
    return minimum_z, peak_sum_kib, point_count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", action="append", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--candidate-resolution", type=float, default=0.1)
    parser.add_argument("--output-resolution", type=float, default=0.2)
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--chunk-points", type=int, default=1_000_000)
    args = parser.parse_args()
    started = time.monotonic()

    if args.workers < 1:
        parser.error("--workers must be at least 1")
    if args.chunk_points < 1:
        parser.error("--chunk-points must be at least 1")
    if args.candidate_resolution <= 0 or args.output_resolution <= 0:
        parser.error("DTM resolutions must be greater than zero")
    if args.candidate_resolution > args.output_resolution:
        parser.error("--candidate-resolution may not exceed --output-resolution")
    for path in args.input:
        if not path.is_file():
            parser.error(f"input does not exist or is not a file: {path}")

    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise SystemExit(f"Output directory is not empty: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    parts_dir = args.output_dir / "parts"
    parts_dir.mkdir()

    headers = []
    for path in args.input:
        with laspy.open(path) as reader:
            if reader.header.point_count < 1:
                parser.error(f"input contains no points: {path}")
            headers.append(reader.header)
    source_crs = projection_signature(headers[0])
    for path, header in zip(args.input[1:], headers[1:]):
        if projection_signature(header) != source_crs:
            parser.error(f"input CRS differs from the first input: {path}")
    input_min_x = min(float(header.mins[0]) for header in headers)
    input_min_y = min(float(header.mins[1]) for header in headers)
    input_max_x = max(float(header.maxs[0]) for header in headers)
    input_max_y = max(float(header.maxs[1]) for header in headers)
    grid_min_x = math.floor(input_min_x / args.output_resolution) * args.output_resolution
    grid_min_y = math.floor(input_min_y / args.output_resolution) * args.output_resolution
    grid_max_x = math.ceil(input_max_x / args.output_resolution) * args.output_resolution
    grid_max_y = math.ceil(input_max_y / args.output_resolution) * args.output_resolution
    cols = max(1, round((grid_max_x - grid_min_x) / args.candidate_resolution))
    rows = max(1, round((grid_max_y - grid_min_y) / args.candidate_resolution))

    global_z = np.full((rows, cols), np.inf, dtype=np.float32)
    total_points = 0
    worker_peak_sum_kib = 0
    for file_index, input_path in enumerate(args.input):
        file_parts = parts_dir / f"file_{file_index:03d}"
        file_parts.mkdir()
        file_z, peak_kib, point_count = process_file(
            input_path,
            file_parts,
            args.workers,
            args.chunk_points,
            grid_min_x,
            grid_min_y,
            rows,
            cols,
            args.candidate_resolution,
        )
        np.minimum(global_z, file_z, out=global_z)
        total_points += point_count
        worker_peak_sum_kib = max(worker_peak_sum_kib, peak_kib)
        shutil.rmtree(file_parts)

    candidate_path = args.output_dir / "minimum_candidates.laz"
    candidate_points = write_candidates(
        candidate_path,
        headers[0],
        global_z,
        grid_min_x,
        grid_min_y,
        args.candidate_resolution,
    )
    elapsed_seconds = time.monotonic() - started
    self_usage = resource.getrusage(resource.RUSAGE_SELF)
    child_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu_user_seconds = self_usage.ru_utime + child_usage.ru_utime
    cpu_system_seconds = self_usage.ru_stime + child_usage.ru_stime
    metadata = {
        "input_files": [str(path) for path in args.input],
        "input_points": total_points,
        "candidate_points": candidate_points,
        "candidate_resolution_m": args.candidate_resolution,
        "candidate_xy_mode": "cell_center",
        "output_resolution_m": args.output_resolution,
        "min_x": grid_min_x,
        "min_y": grid_min_y,
        "max_x": grid_max_x,
        "max_y": grid_max_y,
        "rows": rows,
        "cols": cols,
        "workers": min(args.workers, total_points),
        "chunk_points": args.chunk_points,
        "worker_peak_rss_sum_upper_bound_kb": worker_peak_sum_kib,
        "parent_peak_rss_kb": self_usage.ru_maxrss,
        "cpu_user_seconds": cpu_user_seconds,
        "cpu_system_seconds": cpu_system_seconds,
        "average_cpu_cores": (cpu_user_seconds + cpu_system_seconds) / elapsed_seconds,
        "elapsed_seconds": elapsed_seconds,
    }
    (args.output_dir / "candidate_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(metadata, sort_keys=True), flush=True)
    shutil.rmtree(parts_dir)


if __name__ == "__main__":
    main()
