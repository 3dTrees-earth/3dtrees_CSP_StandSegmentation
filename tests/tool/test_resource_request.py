"""Header-only resource request checks using the committed COPC fixture."""

from pathlib import Path
import json
import sys
import tempfile
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow"))
from resource_request import request_resources, owned_resources  # noqa: E402


FIXTURE = ROOT / "tests" / "tool" / "fixtures" / "multi-instance.copc.laz"


class ResourceRequestTests(unittest.TestCase):
    def test_recommendations_round_trip_through_real_cli(self):
        requests = []
        for mode in ("sequential", "copc-owned"):
            for supplied in (True, False):
                for cpus in (1, 5, 20):
                    requests.append(request_resources(FIXTURE, dtm_supplied=supplied,
                                    inventory_read_mode=mode, max_cpus=cpus, profile="throughput"))
        requests.append(request_resources(FIXTURE, dtm_supplied=True, enable_csp=True))
        script = r"""
args <- commandArgs(trailingOnly = TRUE)
source(file.path(args[[1L]], "workflow", "cli.R"))
requests <- jsonlite::fromJSON(args[[2L]], simplifyVector = FALSE)
for (request in requests) {
  base <- c("--input", request$input, "--output-dir", "unused-output",
            "--segmentation-spec", "PredInstance_FM")
  if (request$dtm_supplied) base <- c(base, "--dtm", "supplied.tif")
  config <- parse_cli_args(c(base, unlist(request$cli_args)))
  if (!request$csp_enabled) stopifnot(config$inventory_workers == request$inventory_workers)
  if (identical(request$inventory_read_mode, "copc-owned")) {
    stopifnot(config$inventory_core_size == request$inventory_core_size_m,
              config$inventory_buffer == request$inventory_buffer_m,
              config$inventory_batch_points == request$inventory_batch_points,
              config$inventory_retained_points == request$inventory_retained_points,
              config$inventory_decode_threads == request$inventory_decode_threads,
              config$inventory_prefetch_mb == request$inventory_prefetch_mb,
              !"--inventory-partitions" %in% unlist(request$cli_args))
  }
  if (request$dtm_workers > 0) stopifnot(config$dtm_workers == request$dtm_workers)
}
cat("all resource requests accepted by CLI\n")
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"requests.json"
            path.write_text(json.dumps(requests))
            result = subprocess.run(["Rscript", "-", str(ROOT), str(path)], input=script,
                                    text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_owned_requests_bound_buffers_and_never_request_point_partitions(self):
        for cpus, ram in ((1, 4), (5, 16), (20, 100)):
            result = owned_resources(FIXTURE, 4_257_196_560, 67, 1_700_000,
                                     True, "throughput", cpus, ram, 8)
            self.assertLessEqual(result["cpus"], cpus)
            self.assertLessEqual(result["memory_gb"], ram)
            self.assertEqual(result["inventory_read_mode"], "copc-owned")
            self.assertEqual(result["inventory_partitions"], 0)
            self.assertEqual(result["scratch_gb"], 4)
            self.assertGreaterEqual(result["inventory_retained_points"], 2_000_000)

    def test_owned_single_core_can_request_multiple_workers(self):
        with tempfile.TemporaryDirectory() as temporary:
            index = Path(temporary) / "index.json"
            stat = FIXTURE.stat()
            index.write_text(json.dumps({
                "schema_version": 1, "input_points": 4000000,
                "source_identity": {"bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns},
                "layouts": [{"windows": [{"ids": list(range(1, 17)), "counts": [250000] * 16}]}],
            }))
            result = owned_resources(FIXTURE, 4000000, 67, 1600, True,
                                     "throughput", 20, 100, 8, index)
            self.assertEqual(result["inventory_workers"], 16)
            self.assertEqual(result["inventory_decode_threads"], 4)
            self.assertEqual(result["cpus"], 20)
            self.assertLess(result["memory_gb"], 32)

    def test_owned_index_must_match_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            index = Path(temporary) / "index.json"
            index.write_text(json.dumps({"schema_version": 1, "input_points": 1,
                                         "source_identity": {"bytes": 1, "mtime_ns": 1}}))
            with self.assertRaisesRegex(ValueError, "does not match"):
                request_resources(FIXTURE, dtm_supplied=True,
                                  inventory_read_mode="copc-owned", owned_index=index)

    def test_supplied_dtm_uses_one_cpu_and_stays_within_caps(self):
        result = request_resources(FIXTURE, dtm_supplied=True)
        self.assertEqual(result["cpus"], 1)
        self.assertGreater(result["input_points"], 0)
        self.assertGreaterEqual(result["memory_gb"], 4)
        self.assertLessEqual(result["memory_gb"], 100)
        self.assertGreaterEqual(result["inventory_partitions"], 64)
        self.assertGreater(result["scratch_gb"], 0)

    def test_generated_dtm_adds_workers_and_memory_headroom(self):
        supplied = request_resources(FIXTURE, dtm_supplied=True)
        generated = request_resources(FIXTURE, dtm_supplied=False)
        self.assertEqual(generated["cpus"], 4)
        self.assertEqual(generated["dtm_workers"], 4)
        self.assertGreaterEqual(generated["memory_gb"], supplied["memory_gb"])

    def test_respects_smaller_cpu_cap(self):
        result = request_resources(FIXTURE, dtm_supplied=False, max_cpus=2)
        self.assertEqual(result["cpus"], 2)
        self.assertEqual(result["dtm_workers"], 2)

    def test_rejects_unsupported_cap(self):
        with self.assertRaisesRegex(ValueError, "Approved maximum"):
            request_resources(FIXTURE, dtm_supplied=True, max_memory_gb=101)

    def test_partition_preflight_rejects_one_worker_above_cap(self):
        with tempfile.TemporaryDirectory() as temporary:
            partition_dir = Path(temporary) / "partitions"
            fm_dir = partition_dir / "FM"
            fm_dir.mkdir(parents=True)
            with (fm_dir / "00001.bin").open("wb") as output:
                output.truncate(10 * 2**30)
            with self.assertRaisesRegex(ValueError, "bounded inventory algorithm"):
                request_resources(FIXTURE, dtm_supplied=True,
                                  partition_dir=partition_dir, max_memory_gb=100)

    def test_instance_audit_rejects_giant_tree_before_worker_request(self):
        with tempfile.TemporaryDirectory() as temporary:
            partition_dir = Path(temporary) / "partitions"
            fm_dir = partition_dir / "FM"
            fm_dir.mkdir(parents=True)
            (fm_dir / "00001.bin").write_bytes(bytes(64))
            summary = partition_dir / "partition_summary.json"
            summary.write_text(json.dumps({"largest_instances": {
                "FM": {"instance_id": 116, "points": 109_750_559}
            }}))
            with self.assertRaisesRegex(ValueError, "Instance FM:116.*without a complete oversized-instance audit"):
                request_resources(FIXTURE, dtm_supplied=True, partition_dir=partition_dir)
            summary.write_text(json.dumps({
                "largest_instances": {"FM": {"instance_id": 116, "points": 109_750_559}},
                "oversized_instances": {"FM": [
                    {"instance_id": 116, "points": 109_750_559}]},
            }))
            large_request = request_resources(FIXTURE, dtm_supplied=True,
                                              partition_dir=partition_dir)
            self.assertEqual(large_request["preflight_oversized_instance_count"], 1)
            self.assertLessEqual(large_request["memory_gb"], 100)
            calibration = Path(temporary) / "calibration.json"
            calibration.write_text(json.dumps({
                "input_bytes": FIXTURE.stat().st_size,
                "input_points": large_request["input_points"],
                "cgroup_peak_bytes": 1_000_000_000,
                "inventory_workers": 1,
            }))
            calibrated = request_resources(FIXTURE, dtm_supplied=True,
                                           partition_dir=partition_dir,
                                           calibration=calibration)
            self.assertEqual(calibrated["confidence"], "same_input_full_run_calibration")
            self.assertEqual(calibrated["calibrated_cgroup_peak_gb"], 1.0)
            calibration.write_text(json.dumps({
                "input_bytes": 1, "input_points": 1,
                "cgroup_peak_bytes": 1_000_000_000, "inventory_workers": 1,
            }))
            with self.assertRaisesRegex(ValueError, "same COPC"):
                request_resources(FIXTURE, dtm_supplied=True,
                                  partition_dir=partition_dir,
                                  calibration=calibration)
            summary.write_text(json.dumps({"largest_instances": {
                "FM": {"instance_id": 2, "points": 400_000}
            }}))
            request = request_resources(FIXTURE, dtm_supplied=True,
                                        partition_dir=partition_dir)
            self.assertEqual(request["preflight_largest_instance_points"], 400_000)
            self.assertEqual(request["confidence"], "bounded_instance_preflight")


if __name__ == "__main__":
    unittest.main()
