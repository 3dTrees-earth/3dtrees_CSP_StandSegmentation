#!/usr/bin/env python3
"""Build an exact large-instance audit for already completed inventory partitions.

Use this only when resuming partitions made before the scanner recorded all
oversized IDs. Layout syntax is SLUG:COLUMNS:ID_COLUMN (zero-based index).
The audit reads IDs in bounded batches and never loads a partition into RAM.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np

from inventory_policy import INSTANCE_POINT_LIMIT


def audit(partition_dir: Path, layouts: list[tuple[str, int, int]],
          limit: int = INSTANCE_POINT_LIMIT, batch_rows: int = 1_000_000) -> dict:
    if limit < 1 or batch_rows < 1:
        raise ValueError("Limit and batch size must be positive")
    largest = {}
    oversized = {}
    instance_counts = {}
    partition_counts = {}
    for slug, column_count, id_column in layouts:
        if column_count < 4 or not 0 <= id_column < column_count:
            raise ValueError(f"Invalid layout for {slug}")
        files = sorted((partition_dir / slug).glob("*.bin"))
        if not files:
            raise ValueError(f"No partition files for {slug}")
        counts: dict[float, int] = defaultdict(int)
        owner: dict[float, Path] = {}
        total_rows = 0
        for path in files:
            size = path.stat().st_size
            if size % (8 * column_count):
                raise ValueError(f"Invalid partition byte length: {path}")
            rows = size // (8 * column_count)
            total_rows += rows
            if rows == 0:
                continue
            records = np.memmap(path, dtype="<f8", mode="r",
                                shape=(rows, column_count))
            for start in range(0, rows, batch_rows):
                ids = np.asarray(records[start:start + batch_rows, id_column])
                if not np.all(np.isfinite(ids) & (ids >= 0)):
                    raise ValueError(f"Invalid instance ID in {path}")
                unique, frequencies = np.unique(ids, return_counts=True)
                for instance_id, frequency in zip(unique, frequencies):
                    key = float(instance_id)
                    previous = owner.setdefault(key, path)
                    if previous != path:
                        raise ValueError(
                            f"Instance {key} occurs in {previous} and {path}; "
                            "hash partitions must keep an instance in one file"
                        )
                    counts[key] += int(frequency)
            del records
        if not counts:
            raise ValueError(f"No instances for {slug}")
        max_id, max_points = max(counts.items(), key=lambda item: item[1])
        largest[slug] = {"instance_id": max_id, "points": max_points}
        oversized[slug] = [
            {"instance_id": instance_id, "points": points}
            for instance_id, points in sorted(counts.items()) if points > limit
        ]
        instance_counts[slug] = len(counts)
        partition_counts[slug] = {"files": len(files), "rows": total_rows,
                                  "bytes": sum(path.stat().st_size for path in files)}
    return {"largest_instances": largest, "oversized_instances": oversized,
            "instance_counts": instance_counts,
            "partition_counts": partition_counts,
            "instance_point_limit": limit}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partition-dir", required=True, type=Path)
    parser.add_argument("--layout", action="append", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--limit", type=int, default=INSTANCE_POINT_LIMIT)
    args = parser.parse_args()
    layouts = []
    for value in args.layout:
        try:
            slug, columns, id_column = value.split(":")
            layouts.append((slug, int(columns), int(id_column)))
        except ValueError as error:
            parser.error(f"Invalid --layout {value!r}: {error}")
    result = audit(args.partition_dir, layouts, args.limit)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"output": str(args.output),
                      "oversized_instance_counts": {
                          slug: len(items) for slug, items in
                          result["oversized_instances"].items()},
                      "partition_counts": result["partition_counts"]}))


if __name__ == "__main__":
    main()
