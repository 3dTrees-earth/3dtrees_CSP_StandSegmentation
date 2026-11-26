"""64-bit COPC metadata used before any R LAS header parsing."""

from pathlib import Path
import struct
import sys
from tempfile import TemporaryDirectory
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
        self.assertLess(result["zmin"], result["zmax"])

    @staticmethod
    def write_unlabelled_las(path):
        cloud = laspy.create(point_format=6, file_version="1.4")
        cloud.x, cloud.y, cloud.z = [0, 2], [0, 3], [4, 5]
        cloud.vlrs.append(laspy.VLR("LASF_Projection", 2112, record_data=b""))
        cloud.write(path)

    def test_empty_wkt_is_unlabelled_and_z_bounds_are_reported(self):
        with TemporaryDirectory() as work:
            path = Path(work) / "empty-wkt.las"
            self.write_unlabelled_las(path)
            header = read_header(path)
        self.assertIsNone(header["crs_wkt"])
        self.assertEqual((header["zmin"], header["zmax"]), (4.0, 5.0))

    def test_invalid_bounds_are_refused(self):
        # The LAS public header stores max X at 179 and max Z at 211.
        with TemporaryDirectory() as work:
            path = Path(work) / "invalid-bounds.las"
            for offset, value in [(179, float("inf")), (211, float("nan")), (211, 3.0)]:
                with self.subTest(offset=offset, value=value):
                    self.write_unlabelled_las(path)
                    with path.open("r+b") as stream:
                        stream.seek(offset)
                        stream.write(struct.pack("<d", value))
                    with self.assertRaisesRegex(ValueError, "invalid point count or bounds"):
                        read_header(path)


if __name__ == "__main__":
    unittest.main()
