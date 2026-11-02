"""Local COPC scanner regression tests; no production data or services."""

from pathlib import Path
import sys
import tempfile
import unittest

import laspy
import numpy as np
from scipy.spatial import ConvexHull
import json
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow"))
from stream_inventory_partitions import process  # noqa: E402


FIXTURE = ROOT / "tests" / "tool" / "fixtures" / "multi-instance.copc.laz"


class ScannerTests(unittest.TestCase):
    def test_many_distinct_ids_fail_before_count_map_grows(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = {
                "input": str(FIXTURE), "output_dir": str(root / "partitions"),
                "hull_path": str(root / "hull.bin"), "partitions": 8,
                "chunk_points": 5000,
                "spatial_grid": {"xmin": -50, "ymin": -70, "size": 20,
                                 "nx": 2, "ny": 2},
                "layouts": [{"slug": "SAT", "instance": "PredInstance_SAT",
                             "columns": ["X", "Y", "Z", "PredInstance_SAT"]}],
                "non_tree_ids": [0],
            }
            with patch("stream_inventory_partitions.MAX_TRACKED_INSTANCE_IDS", 1):
                with self.assertRaisesRegex(ValueError, "distinct instance IDs"):
                    process(config)

    def test_selected_groups_preserve_partitions_and_exact_hull(self):
        with tempfile.TemporaryDirectory() as temporary:
            outputs = []
            for selective in (False, True):
                root = Path(temporary) / str(selective)
                root.mkdir()
                config = {
                    "input": str(FIXTURE),
                    "output_dir": str(root / "partitions"),
                    "hull_path": str(root / "hull.bin"),
                    "partitions": 8,
                    "chunk_points": 5000,
                    "spatial_grid": {"xmin": -50, "ymin": -70, "size": 20,
                                     "nx": 2, "ny": 2},
                    "layouts": [
                        {"slug": "FM", "instance": "PredInstance_FM",
                         "columns": ["X", "Y", "Z", "PredInstance_FM",
                                     "species_id_FM", "species_prob_FM",
                                     "PredScore_FM", "PredSemantic_FM"]},
                        {"slug": "SAT", "instance": "PredInstance_SAT",
                         "columns": ["X", "Y", "Z", "PredInstance_SAT"]},
                    ],
                    "non_tree_ids": [0],
                    "selective_decompression": selective,
                }
                outputs.append((root, process(config)))

            baseline, selected = outputs
            self.assertEqual(selected[1]["decompressed_groups"],
                             ["XY_RETURNS_CHANNEL", "Z", "ALL_EXTRA_BYTES"])
            self.assertEqual(selected[1]["source_point_reads"], 1)
            self.assertEqual(selected[1]["selector"], "python_sequential_selective")
            self.assertNotIn("classification", selected[1]["loaded_dimensions"])
            self.assertEqual(selected[1]["xy_spool_bytes"], 0)
            self.assertEqual(selected[1]["largest_instances"],
                             baseline[1]["largest_instances"])
            summary = json.loads((selected[0] / "partitions" /
                                  "partition_summary.json").read_text())
            self.assertEqual(summary["largest_instances"],
                             selected[1]["largest_instances"])
            self.assertEqual(summary["oversized_instances"],
                             selected[1]["oversized_instances"])
            for source in (baseline[0] / "partitions").glob("*/*.bin"):
                target = selected[0] / "partitions" / source.relative_to(
                    baseline[0] / "partitions")
                self.assertEqual(source.read_bytes(), target.read_bytes())

            vertices = np.fromfile(selected[0] / "hull.bin", dtype="<f8").reshape(-1, 2)
            with laspy.open(FIXTURE) as reader:
                points = reader.read()
            full_xy = np.column_stack((points.x, points.y))
            full_hull = ConvexHull(full_xy)
            self.assertAlmostEqual(ConvexHull(vertices).volume, full_hull.volume, places=9)

    def test_scanner_reports_all_large_ids(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = {
                "input": str(FIXTURE), "output_dir": str(root / "partitions"),
                "hull_path": str(root / "hull.bin"), "partitions": 8,
                "chunk_points": 5000,
                "spatial_grid": {"xmin": -50, "ymin": -70, "size": 20,
                                 "nx": 2, "ny": 2},
                "layouts": [{"slug": "SAT", "instance": "PredInstance_SAT",
                             "columns": ["X", "Y", "Z", "PredInstance_SAT"]}],
                "non_tree_ids": [0],
            }
            with patch("stream_inventory_partitions.STANDARD_INSTANCE_POINT_LIMIT", 1):
                result = process(config)
            self.assertEqual(len(result["oversized_instances"]["SAT"]), 2)
            self.assertGreater(result["oversized_instances"]["SAT"][0]["points"], 1)


if __name__ == "__main__":
    unittest.main()
