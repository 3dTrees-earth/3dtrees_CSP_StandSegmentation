"""64-bit COPC metadata used before any R LAS header parsing."""

from pathlib import Path
import sys
import unittest

import laspy

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow"))
from copc_header import read_header  # noqa: E402

FIXTURE = ROOT / "tests" / "tool" / "fixtures" / "multi-instance.copc.laz"


class CopcHeaderTests(unittest.TestCase):
    def test_reads_copc_count_crs_and_extra_dimensions(self):
        result = read_header(FIXTURE)
        with laspy.open(FIXTURE) as reader:
            self.assertEqual(result["point_count"], reader.header.point_count)
        self.assertIn("PredInstance_FM", result["available_dimensions"])
        self.assertIn("PROJ", result["crs_wkt"].upper())
        self.assertLess(result["xmin"], result["xmax"])


if __name__ == "__main__":
    unittest.main()
