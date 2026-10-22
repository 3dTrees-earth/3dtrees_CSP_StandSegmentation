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


from inventory_policy import INSTANCE_POINT_LIMIT, SAMPLE_LIMIT


def _instance_label(value: object) -> str:
    number = float(value)
    return str(int(number)) if number.is_integer() else str(number)


def _largest_instance(summary_path: Path) -> tuple[str, int] | None:
    """Read either the scanner summary or an older partition audit."""
    summary = json.loads(summary_path.read_text())
    entries = summary.get("largest_instances", {})
    if isinstance(entries, dict):
        candidates = [(f"{slug}:{_instance_label(item['instance_id'])}", int(item["points"]))
                      for slug, item in entries.items()]
    else:
        candidates = [(str(item[0]), int(item[1])) for item in entries]
    if not candidates and "largest_partition_instances" in summary:
        candidates = [(str(item[0]), int(item[1])) for item in
                      summary["largest_partition_instances"].get("largest_instances", [])]
    return max(candidates, key=lambda item: item[1]) if candidates else None


def _oversized_instance_ids(summary_path: Path) -> set[str] | None:
    """Return audited large IDs; None means the audit predates streaming fallback."""
    summary = json.loads(summary_path.read_text())
    entries = summary.get("oversized_instances")
    if not isinstance(entries, dict):
        return None
    return {f"{slug}:{_instance_label(item['instance_id'])}" for slug, items in entries.items()
            for item in items}


def _round_up(value: float, multiple: int) -> int:
    return int(math.ceil(value / multiple) * multiple)


def _power_of_two_at_least(value: int) -> int:
    return 1 << (value - 1).bit_length()


def resource_cli_args(request):
    """Executable resource settings; reporting-only zero counts stay out of argv.

    Input/output paths, segmentation specs and DTM path remain caller inputs.
    Never derive CLI flags by blindly translating every reporting field.
    """
    args = []
    def add(flag, key):
        args.extend((flag, str(request[key])))
    mode = request["inventory_read_mode"]
    if request["csp_enabled"]:
        args.extend(("--enable-csp", "true"))
    else:
        add("--inventory-read-mode", "inventory_read_mode")
        add("--inventory-workers", "inventory_workers")
        add("--read-chunk-size", "read_chunk_size_m")
        if mode == "copc-owned":
            for flag, key in (
                ("--inventory-core-size", "inventory_core_size_m"),
                ("--inventory-buffer", "inventory_buffer_m"),
                ("--inventory-batch-points", "inventory_batch_points"),
                ("--inventory-retained-points", "inventory_retained_points"),
                ("--inventory-decode-threads", "inventory_decode_threads"),
                ("--inventory-prefetch-mb", "inventory_prefetch_mb"),
            ):
                add(flag, key)
        else:
            add("--inventory-partitions", "inventory_partitions")
    if request["dtm_workers"] > 0:
        add("--dtm-workers", "dtm_workers")
    return args


def owned_resources(input_path, points, point_bytes, density, dtm_supplied,
                    profile, max_cpus, max_memory_gb, inventory_columns, index_path=None):
    """Request resources from bounded buffers, optionally refined by owned windows.

    Four copies of retained numeric columns allow for R tables, LAS and geometry
    temporaries. This is a conservative model, not measured peak memory. Exact
    convex-hull vertex sets can still grow with input geometry.
    """
    index = None
    if index_path is not None:
        index = json.loads(Path(index_path).read_text())
        stat = input_path.stat()
        if (index.get("schema_version") != 1 or index.get("input_points") != points or
                index.get("source_identity") != {"bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}):
            raise ValueError("Owned instance index does not match the input COPC")
    # One decoder shared by all instance workers, even inside a single core.
    decode_threads = min(4, max(1, max_cpus - 1))
    workers = 1 if profile == "stable" else max(1, max_cpus - decode_threads)
    windows = None
    if index is not None:
        windows = [window for layout in index["layouts"] for window in layout["windows"]]
        if not windows:
            raise ValueError("Owned index contains no inventory windows")
        workers = min(workers, max(len(window["ids"]) for window in windows))
    base_gb = 2 if dtm_supplied else 4
    batch_points = 2_000_000
    # Native decode + float64 routing buffers are shared; each R worker needs
    # runtime/geometry headroom. Retention is a total window budget, not copied
    # for every worker: an instance belongs to exactly one consumer.
    batch_gb = min(batch_points * point_bytes, 256 * 1024**2) * 6 / 1e9
    worker_gb = 0.5
    bytes_per_retained_point = 8 * (inventory_columns + 1) * 4
    minimum_retention_gb = INSTANCE_POINT_LIMIT * bytes_per_retained_point / 1e9
    while workers > 1 and base_gb + batch_gb + workers * worker_gb + minimum_retention_gb > max_memory_gb:
        workers -= 1
    # A bounded queue lets the decoder advance beyond a slow spatial cluster.
    # With an index, target one largest instance (capped at 8 GiB and 1/8 RAM).
    prefetch_mb = 256
    if windows is not None and workers > 1:
        largest_bytes = max(max(w["counts"]) for w in windows) * inventory_columns * 8
        while prefetch_mb * 1024**2 < largest_bytes and prefetch_mb < 8192:
            prefetch_mb *= 2
    while prefetch_mb > 1 and (prefetch_mb * 1024**2 / 1e9 > max_memory_gb / 8 or
            base_gb + batch_gb + workers * worker_gb + minimum_retention_gb +
            prefetch_mb * 1024**2 / 1e9 > max_memory_gb):
        prefetch_mb //= 2
    queue_gb = prefetch_mb * 1024**2 / 1e9 if workers > 1 else 0
    budget = max_memory_gb - base_gb - batch_gb - workers * worker_gb - queue_gb
    retained_points = min(40_000_000, int(budget * 1e9 / bytes_per_retained_point))
    if retained_points < INSTANCE_POINT_LIMIT:
        raise ValueError("RAM cap cannot fit one standard owned instance plus the COPC batch")
    held_points = min(points, retained_points)
    if windows is not None:
        held_points = min(held_points, max(sum(n if n <= INSTANCE_POINT_LIMIT else SAMPLE_LIMIT
                                             for n in window["counts"]) for window in windows))
    estimate = base_gb + batch_gb + workers * worker_gb + queue_gb + held_points * bytes_per_retained_point / 1e9
    memory = min(max_memory_gb, max(4, _round_up(estimate, 4)))
    # A 1-CPU cgroup time-shares the producer and one consumer.
    cpu = min(max_cpus, workers + decode_threads)
    result = {
        "input": str(input_path), "input_bytes": input_path.stat().st_size,
        "input_points": points, "xy_density_points_per_m2": round(density, 2),
        "dtm_supplied": dtm_supplied, "csp_enabled": False, "profile": profile,
        "cpus": cpu, "memory_gb": memory, "scratch_gb": 4,
        "read_mode": "copc_owned_windows", "inventory_read_mode": "copc-owned",
        "inventory_core_size_m": 50, "inventory_buffer_m": 20,
        "inventory_batch_points": batch_points, "inventory_retained_points": retained_points,
        "inventory_decode_threads": decode_threads, "inventory_workers": workers,
        "inventory_prefetch_mb": prefetch_mb,
        "inventory_partitions": 0, "read_chunk_size_m": 50,
        "dtm_workers": 0 if dtm_supplied else min(4, cpu),
        "confidence": "owned_window_memory_model" if index is None else "owned_instance_index_memory_model",
        "limitations": [
            "Provisional buffer-based estimate; enforce the returned CPU/RAM caps and measure cgroup peak.",
            "One complete index pass plus spatial window rereads; each window has one decoder shared by disjoint ID workers.",
            "Worker and decoder counts share the CPU budget; a one-CPU cap time-shares both processes.",
            "Scratch excludes the existing input. DTM scratch for very large spatial extents needs separate sizing.",
            "Instances above two million points retain complete counts/extrema/hulls but sample medians and omit DBH.",
            "Exact hull vertex counts can grow with geometry; the hard RAM limit remains required.",
            "The 20 m buffer expands to cover owned instance extents; point counts must match the global index.",
        ],
    }
    result["cli_args"] = resource_cli_args(result)
    return result


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
    instance_preflight: Path | None = None,
    calibration: Path | None = None,
    inventory_read_mode: str = "sequential",
    owned_index: Path | None = None,
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
    if inventory_read_mode not in {"sequential", "copc-owned"}:
        raise ValueError("Resource estimates support sequential or copc-owned reads")
    if owned_index is not None and inventory_read_mode != "copc-owned":
        raise ValueError("--owned-index requires --inventory-read-mode copc-owned")
    if inventory_read_mode == "copc-owned" and (enable_csp or partition_dir is not None or calibration is not None):
        raise ValueError("Owned COPC inventory cannot use CSP or partition calibration")
    if instance_preflight is not None and partition_dir is None:
        raise ValueError("--instance-preflight requires --partition-dir")
    if calibration is not None and partition_dir is None:
        raise ValueError("--calibration requires --partition-dir")
    if instance_preflight is not None and not Path(instance_preflight).is_file():
        raise ValueError(f"Instance preflight file does not exist: {instance_preflight}")
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
    calibration_data = None
    if calibration is not None:
        calibration_data = json.loads(Path(calibration).read_text())
        try:
            matches = (int(calibration_data.get("input_bytes", -1)) == input_bytes and
                       int(calibration_data.get("input_points", -1)) == points and
                       int(calibration_data.get("cgroup_peak_bytes", 0)) > 0 and
                       1 <= int(calibration_data.get("inventory_workers", 0)) <= 20)
        except (TypeError, ValueError, AttributeError):
            matches = False
        if not matches:
            raise ValueError("Calibration must record the same COPC and a valid full inventory run")
    input_gb = input_bytes / 1e9
    area_m2 = (xmax - xmin) * (ymax - ymin)
    density = points / area_m2

    if inventory_read_mode == "copc-owned":
        return owned_resources(input_path, points, point_bytes, density, dtm_supplied,
                               profile, max_cpus, max_memory_gb, inventory_columns, owned_index)

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
            memory_gb = max(memory_gb, 60)
        if not dtm_supplied:
            memory_gb = max(memory_gb, 8 if points <= 10_000_000 else 16)
        # COPC DTM candidate reads are sequential: generic LAZ random seeks
        # can yield zero-filled records for a valid COPC hierarchy.
        dtm_cpu = 1 if dtm_supplied or points >= 50_000_000 else 4
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
        # final partitions together. Inventory later holds final partitions
        # plus one normalized spool with Zref and source-order fields.
        normalized_spool_bytes = points * (inventory_columns + 2) * 8
        scratch_peak_bytes = max(2 * partition_bytes,
                                 partition_bytes + normalized_spool_bytes)
        scratch_gb = max(4, math.ceil(
            (scratch_peak_bytes + input_bytes * 0.5) * 1.2 / 1e9
        ))

    largest_partition_bytes = None
    partition_file_count = None
    estimated_peak_anon_gb = None
    largest_instance_id = None
    largest_instance_points = None
    oversized_instance_count = None
    calibrated_peak_gb = None
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
        summary_path = instance_preflight or partition_dir / "partition_summary.json"
        largest = _largest_instance(summary_path) if summary_path.is_file() else None
        oversized_ids = (_oversized_instance_ids(summary_path)
                         if summary_path.is_file() else None)
        if calibration_data is not None and not oversized_ids:
            raise ValueError("Calibration requires an oversized-instance audit")
        oversized_instance_count = (len(oversized_ids)
                                    if oversized_ids is not None else None)
        if largest is not None:
            largest_instance_id, largest_instance_points = largest
            if (largest_instance_points > INSTANCE_POINT_LIMIT and
                    (oversized_ids is None or largest_instance_id not in oversized_ids)):
                raise ValueError(
                    f"Instance {largest_instance_id} has {largest_instance_points:,} points, "
                    f"above the standard inventory limit of {INSTANCE_POINT_LIMIT:,} "
                    "without a complete oversized-instance audit; rescan or supply that audit"
                )
            if oversized_ids:
                # The 45.3 GB full run peaked at 13.7 GB with 20 workers.
                # Use a 2 GB/worker audit request until an exact same-input
                # cgroup measurement is supplied for a tighter request.
                if calibration_data is not None:
                    calibrated_peak_gb = calibration_data["cgroup_peak_bytes"] / 1e9
                    worker_gb = max(0.8, calibrated_peak_gb /
                                    calibration_data["inventory_workers"])
                    parent_gb = 4
                else:
                    worker_gb = 2
                    parent_gb = 8
            else:
                # Standard instances still materialize a capped LAS.
                worker_gb = max(3, 2 + 8 * largest_instance_points / INSTANCE_POINT_LIMIT)
                parent_gb = 1 if points <= 10_000_000 else 3 if points <= 120_000_000 else 8
        else:
            # Without an instance audit, retain the previous empirical risk
            # screen; partition bytes alone cannot prove a giant tree is safe.
            worker_gb = 17 * largest_partition_bytes / 2**30
            parent_gb = 8
        permitted_workers = int((max_memory_gb - parent_gb) // worker_gb) if worker_gb else max_cpus
        if permitted_workers < 1:
            raise ValueError(
                f"Largest partition is {largest_partition_bytes / 2**30:.2f} GiB; "
                f"one inventory worker is estimated above the {max_memory_gb} GB RAM cap. "
                "Reduce the largest instance or use a bounded inventory algorithm; "
                "more hash partitions cannot split one instance."
            )
        inventory_workers = min(inventory_workers, permitted_workers)
        cpu = max(dtm_cpu, inventory_workers)
        estimated_peak_anon_gb = round(parent_gb + worker_gb * inventory_workers, 1)
        if oversized_ids:
            memory_gb = max(16, _round_up(estimated_peak_anon_gb, 4))
            if calibrated_peak_gb is not None:
                memory_gb = max(memory_gb, _round_up(calibrated_peak_gb * 1.6, 4))
        else:
            memory_gb = max(memory_gb, _round_up(estimated_peak_anon_gb, 4))
        # Completed partitions replace the header's worst-case point-count
        # scratch estimate. During scanning, buckets and final files coexist;
        # during inventory, only active workers have normalized spools.
        partition_bytes_actual = sum(sizes)
        active_spool_bytes = 1.5 * sum(sorted(sizes, reverse=True)[:inventory_workers])
        scratch_peak_bytes = max(2 * partition_bytes_actual,
                                 partition_bytes_actual + active_spool_bytes)
        scratch_gb = max(4, math.ceil(
            (scratch_peak_bytes + input_bytes * 0.5) * 1.2 / 1e9
        ))

    if memory_gb > max_memory_gb:
        raise ValueError(
            f"Estimated {memory_gb} GB exceeds the {max_memory_gb} GB RAM cap; "
            "use a smaller spatial input or change the tool's memory architecture"
        )
    if cpu > max_cpus:
        cpu = max_cpus
    if not enable_csp:
        inventory_workers = min(inventory_workers, cpu)

    result = {
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
        "preflight_largest_instance_id": largest_instance_id,
        "preflight_largest_instance_points": largest_instance_points,
        "preflight_oversized_instance_count": oversized_instance_count,
        "preflight_estimated_peak_anon_gb": estimated_peak_anon_gb,
        "calibrated_cgroup_peak_gb": (round(calibrated_peak_gb, 2)
                                       if calibrated_peak_gb is not None else None),
        "dtm_workers": 0 if dtm_supplied or enable_csp else min(dtm_cpu, cpu),
        "confidence": (
            "same_input_full_run_calibration" if calibrated_peak_gb is not None
            else "bounded_instance_preflight" if largest_instance_points is not None
            else "partition_preflight_risk_screen" if partition_dir is not None
            else "calibrated_size_range" if points <= 925_666_260 and not enable_csp
            else "above_end_to_end_test_range" if points > 925_666_260 and not enable_csp
            else "unvalidated_full_cloud_csp"
        ),
        "limitations": [
            "Global density cannot detect the densest spatial tile or largest instance.",
            "Recalculate for each new COPC; enforce Docker/Galaxy CPU, RAM and scratch limits.",
            "Multiple inventory workers use per-partition random seeds; derived DBH values may differ.",
            "Generated DTM candidates from COPC use one verified sequential point read.",
            "Without an instance-count audit, partition preflight cannot detect a giant instance.",
            "Instances above two million points use bounded summaries when the scan supplies an oversized-instance audit; their DBH is unavailable and selected medians are sampled.",
            "DBH geometry is also capped at 500,000 normalized points per group during inventory.",
        ],
    }
    result["cli_args"] = resource_cli_args(result)
    return result


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
    parser.add_argument("--inventory-read-mode", choices=("sequential", "copc-owned"), default="sequential")
    parser.add_argument("--owned-index", type=Path, help="Optional complete owned-COPC index for resource sizing")
    parser.add_argument("--partition-dir", type=Path,
                        help="Finished instance partition root for a memory-skew preflight")
    parser.add_argument("--instance-preflight", type=Path,
                        help="Optional scanner summary or existing partition audit JSON")
    parser.add_argument("--calibration", type=Path,
                        help="Same-input full-run cgroup peak JSON for a tighter RAM request")
    args = parser.parse_args()
    try:
        result = request_resources(
            args.input, dtm_supplied=args.dtm_supplied,
            enable_csp=args.enable_csp, profile=args.profile,
            max_cpus=args.max_cpus,
            max_memory_gb=args.max_memory_gb,
            inventory_columns=args.inventory_columns,
            partition_dir=args.partition_dir,
            instance_preflight=args.instance_preflight,
            calibration=args.calibration,
            inventory_read_mode=args.inventory_read_mode, owned_index=args.owned_index,
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
