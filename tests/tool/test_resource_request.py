"""Header-only resource request checks using the committed COPC fixture."""

from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow"))
from resource_request import request_resources  # noqa: E402


FIXTURE = ROOT / "tests" / "tool" / "fixtures" / "multi-instance.copc.laz"


class ResourceRequestTests(unittest.TestCase):
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
                                  partition_dir=partition_dir)


if __name__ == "__main__":
    unittest.main()
