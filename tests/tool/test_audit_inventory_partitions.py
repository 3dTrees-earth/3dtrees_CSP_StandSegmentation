"""Checks for resuming already completed binary inventory partitions."""

from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow"))
from audit_inventory_partitions import audit  # noqa: E402


class AuditTests(unittest.TestCase):
    def test_counts_large_ids_across_batches(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            directory = root / "FM"
            directory.mkdir()
            rows = np.array([[1, 2, 3, 116], [1, 2, 3, 7],
                             [1, 2, 3, 116], [1, 2, 3, 116]], dtype="<f8")
            rows.tofile(directory / "00001.bin")
            result = audit(root, [("FM", 4, 3)], limit=2, batch_rows=2)
            self.assertEqual(result["largest_instances"]["FM"]["points"], 3)
            self.assertEqual(result["oversized_instances"]["FM"],
                             [{"instance_id": 116.0, "points": 3}])
            self.assertEqual(result["partition_counts"]["FM"]["rows"], 4)

    def test_rejects_instance_split_between_partitions(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            directory = root / "FM"
            directory.mkdir()
            row = np.array([[1, 2, 3, 116]], dtype="<f8")
            row.tofile(directory / "00001.bin")
            row.tofile(directory / "00002.bin")
            with self.assertRaisesRegex(ValueError, "hash partitions"):
                audit(root, [("FM", 4, 3)], limit=2)


if __name__ == "__main__":
    unittest.main()
