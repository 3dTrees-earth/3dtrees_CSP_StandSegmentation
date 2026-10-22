"""DTM candidate scanning must preserve real COPC XYZ values."""

from pathlib import Path
import math
import sys
import unittest

import laspy
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow"))
from stream_dtm_candidates import process_copc_sequential  # noqa: E402

FIXTURE = ROOT / "tests" / "tool" / "fixtures" / "multi-instance.copc.laz"


class SequentialCopcDtmTests(unittest.TestCase):
    def test_cell_minima_match_all_fixture_points(self):
        with laspy.open(FIXTURE) as reader:
            header = reader.header
            points = reader.read()
        resolution = 1.0
        xmin = math.floor(header.mins[0])
        ymin = math.floor(header.mins[1])
        cols = math.ceil(header.maxs[0]) - xmin
        rows = math.ceil(header.maxs[1]) - ymin
        actual, _, count = process_copc_sequential(
            FIXTURE, 5000, xmin, ymin, rows, cols, resolution)
        expected = np.full((rows, cols), np.inf, dtype=np.float32)
        x = np.asarray(points.x)
        y = np.asarray(points.y)
        z = np.asarray(points.z, dtype=np.float32)
        col = np.clip(np.floor((x - xmin) / resolution).astype(int), 0, cols - 1)
        row = np.clip(np.floor((y - ymin) / resolution).astype(int), 0, rows - 1)
        np.minimum.at(expected.ravel(), row * cols + col, z)
        self.assertEqual(count, header.point_count)
        np.testing.assert_array_equal(actual, expected)


if __name__ == "__main__":
    unittest.main()
