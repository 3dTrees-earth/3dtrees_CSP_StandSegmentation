#!/usr/bin/env python3
"""Estimate a bounded CSP inventory job request from COPC and partition metadata.

The estimator does not read point payloads. Its recommendations are conservative job
requests, not measured peak memory or a guarantee against instance skew.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import laspy


def _round_up(value: float, multiple: int) -> int:
    return int(math.ceil(value / multiple) * multiple)


def _power_of_two_at_least(value: int) -> int:
    return 1 << (value - 1).bit_length()


def request_resources(
    input_path: Path,
    *,
    dtm_supplied: bool,
    enable_csp: bool = False,
    profile: str = "stable",
    max_cpus: int = 20,
    max_memory_gb: int = 100,
    inventory_columns: int = 8,
    partition_dir: Path | None = None,
) -> dict:
    """Return Docker/Galaxy resources and CLI tuning for one merged COPC file.

    Spatial chunk sizing uses global XY density; dense local clusters or one
    unusually large instance can exceed the predicted working set. Callers
    must enforce the returned hard caps and provide the requested scratch.
    """
    if not 1 <= max_cpus <= 20 or not 1 <= max_memory_gb <= 100:
        raise ValueError("Approved maximum is 20 CPUs and 100 GB RAM")
    if profile not in {"stable", "throughput"}:
        raise ValueError("profile must be stable or throughput")
    if inventory_columns < 4:
        raise ValueError("inventory_columns must be at least four")
    input_path = Path(input_path)
    with laspy.open(input_path) as reader:
        header = reader.header
        if not any(vlr.user_id == "copc" and vlr.record_id == 1 for vlr in header.vlrs):
            raise ValueError("Input must be a COPC file with a COPC info VLR")
        points = int(header.point_count)
        xmin, ymin = map(float, header.mins[:2])
        xmax, ymax = map(float, header.maxs[:2])
        point_bytes = int(header.point_format.size)
    if points < 1 or xmax <= xmin or ymax <= ymin:
        raise ValueError("COPC header has invalid point count or XY bounds")

    input_bytes = input_path.stat().st_size
    input_gb = input_bytes / 1e9
    area_m2 = (xmax - xmin) * (ymax - ymin)
    density = points / area_m2

    if enable_csp:
        # CSP currently reads and preserves the complete cloud in memory.
        estimated_working_gb = 4 + points * max(180, point_bytes * 3) / 1e9
        memory_gb = max(8, _round_up(estimated_working_gb * 1.35, 4))
        cpu = 1
        read_mode = "full_cloud"
        partitions = None
        chunk_size_m = None
        scratch_gb = max(8, math.ceil(input_gb * 3))
    else:
        # The measured 100 MB run used about 1 GB anonymous memory. Larger
        # inventories need headroom for the largest instance partition and
        # for a spatial read; the point-count term is intentionally cautious.
        memory_gb = max(4, _round_up(6 + points * 0.012 / 1e6, 4))
        if points <= 10_000_000 and dtm_supplied:
            # A complete 7.8-million-point run passed a hard 4 GB cap.
            memory_gb = 4
        if points >= 120_000_000:
            memory_gb = max(memory_gb, 16)
        if points >= 500_000_000:
            memory_gb = max(memory_gb, 24)
        if points > 2_000_000_000:
            memory_gb = max(memory_gb, 100)
        if not dtm_supplied:
            memory_gb = max(memory_gb, 8 if points <= 10_000_000 else 16)
        dtm_cpu = 1 if dtm_supplied else (4 if points < 50_000_000 else
                                         8 if points < 500_000_000 else 10)
        cpu = dtm_cpu
        inventory_workers = 1
        if profile == "throughput" and points > 10_000_000:
            # One CPU wins for the 7.8M-point run. More than one inventory
            # worker uses partition-local random seeds, so callers select
            # this profile explicitly when turnaround matters more than
            # reproducing the single-worker measurement stream.
            if points <= 120_000_000:
                inventory_workers, memory_gb = 4, max(memory_gb, 16)
            elif points <= 600_000_000:
                inventory_workers, memory_gb = 4, max(memory_gb, 40)
            elif points <= 2_000_000_000:
                # The 925.7M-point complete file finished with this request.
                inventory_workers, memory_gb = 20, max(memory_gb, 100)
            else:
                # The larger tier has scan/partition evidence only; keep the
                # request capped and label it as unvalidated for inventory.
                inventory_workers, memory_gb = 20, max(memory_gb, 100)
            cpu = max(cpu, inventory_workers)
        read_mode = "sequential_copc_and_instance_partitions"
        # The 925.7M-point sample showed severe superlinear work at roughly
        # 7M points per partition. Keep each partition near 1.5–3M points.
        partitions = min(4096, _power_of_two_at_least(max(64, math.ceil(points / 3_000_000))))
        # The COPC scan uses fixed point-count batches; this width determines
        # only the spatial ordering of disk buckets. Keep dense local tiles
        # bounded despite a misleadingly sparse whole-cloud extent.
        target_chunk_points = (2_000_000 if points <= 10_000_000 else
                               3_000_000 if points <= 120_000_000 else 12_000_000)
        chunk_size_m = max(20, min(80, round(math.sqrt(target_chunk_points / density) / 10) * 10))
        partition_bytes = points * inventory_columns * 8
        # The ordered one-pass writer temporarily holds spatial buckets and
        # final partitions together, plus the all-point XY hull spool.
        scratch_gb = max(4, math.ceil(
            (2 * partition_bytes + 16 * points + input_bytes * 0.5) * 1.2 / 1e9
        ))

    largest_partition_bytes = None
    partition_file_count = None
    estimated_peak_anon_gb = None
    if partition_dir is not None:
        if enable_csp:
            raise ValueError("Partition preflight applies only to inventory-only runs")
        partition_dir = Path(partition_dir)
        if (partition_dir / ".spatial_buckets").exists():
            raise ValueError("Partition preflight requires completed partitions")
        sizes = [path.stat().st_size for path in partition_dir.glob("*/*.bin")]
        if not sizes:
            raise ValueError("Partition preflight found no finished partition files")
        largest_partition_bytes = max(sizes)
        partition_file_count = len(sizes)
        # Measured 5.7 GB and 11.6 GB runs peaked at roughly 17 times the
        # largest partition per active worker, plus a small parent allowance.
        # This is a risk screen, not a worst-case memory proof.
        worker_gb = 17 * largest_partition_bytes / 2**30
        permitted_workers = int((max_memory_gb - 8) // worker_gb) if worker_gb else max_cpus
        if permitted_workers < 1:
            raise ValueError(
                f"Largest partition is {largest_partition_bytes / 2**30:.2f} GiB; "
                f"one inventory worker is estimated above the {max_memory_gb} GB RAM cap. "
                "Reduce the largest instance or use a bounded inventory algorithm; "
                "more hash partitions cannot split one instance."
            )
        inventory_workers = min(inventory_workers, permitted_workers)
        cpu = max(dtm_cpu, inventory_workers)
        estimated_peak_anon_gb = round(8 + worker_gb * inventory_workers, 1)
        memory_gb = max(memory_gb, _round_up(estimated_peak_anon_gb, 4))

    if memory_gb > max_memory_gb:
        raise ValueError(
            f"Estimated {memory_gb} GB exceeds the {max_memory_gb} GB RAM cap; "
            "use a smaller spatial input or change the tool's memory architecture"
        )
    if cpu > max_cpus:
        cpu = max_cpus
    if not enable_csp:
        inventory_workers = min(inventory_workers, cpu)

    return {
        "input": str(input_path), "input_bytes": input_bytes,
        "input_points": points, "xy_density_points_per_m2": round(density, 2),
        "dtm_supplied": dtm_supplied, "csp_enabled": enable_csp,
        "profile": profile,
        "cpus": cpu, "memory_gb": memory_gb, "scratch_gb": scratch_gb,
        "read_mode": read_mode, "read_chunk_size_m": chunk_size_m,
        "inventory_read_mode": "sequential" if not enable_csp else None,
        "inventory_workers": inventory_workers if not enable_csp else None,
        "inventory_partitions": partitions,
        "preflight_partition_file_count": partition_file_count,
        "preflight_largest_partition_bytes": largest_partition_bytes,
        "preflight_estimated_peak_anon_gb": estimated_peak_anon_gb,
        "dtm_workers": 0 if dtm_supplied or enable_csp else cpu,
        "confidence": (
            "partition_preflight_risk_screen" if partition_dir is not None
            else "calibrated_size_range" if points <= 925_666_260 and not enable_csp
            else "above_end_to_end_test_range" if points > 925_666_260 and not enable_csp
            else "unvalidated_full_cloud_csp"
        ),
        "limitations": [
            "Global density cannot detect the densest spatial tile or largest instance.",
            "Recalculate for each new COPC; enforce Docker/Galaxy CPU, RAM and scratch limits.",
            "Multiple inventory workers use per-partition random seeds; derived DBH values may differ.",
            "Partition preflight uses an empirical memory multiplier; it cannot bound a giant instance.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Local merged COPC file")
    parser.add_argument("--dtm-supplied", action="store_true")
    parser.add_argument("--enable-csp", action="store_true")
    parser.add_argument("--profile", choices=("stable", "throughput"), default="stable")
    parser.add_argument("--max-cpus", type=int, default=20)
    parser.add_argument("--max-memory-gb", type=int, default=100)
    parser.add_argument("--inventory-columns", type=int, default=8,
                        help="Total binary fields across requested inventories (FM uses eight)")
    parser.add_argument("--partition-dir", type=Path,
                        help="Finished instance partition root for a memory-skew preflight")
    args = parser.parse_args()
    try:
        result = request_resources(
            args.input, dtm_supplied=args.dtm_supplied,
            enable_csp=args.enable_csp, profile=args.profile,
            max_cpus=args.max_cpus,
            max_memory_gb=args.max_memory_gb,
            inventory_columns=args.inventory_columns,
            partition_dir=args.partition_dir,
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
