"""Complete-instance COPC reads, including IDs crossing cores and halos."""
from io import BytesIO
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import laspy
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow"))
import copc_owned_inventory as owned

FIXTURE = ROOT / "tests/tool/fixtures/multi-instance.copc.laz"


def configuration():
    return dict(input=str(FIXTURE), batch_points=25000, core_size=1, buffer=0,
                retained_points=40000000, non_tree_ids=[0], layouts=[
                    dict(slug="FM", instance="PredInstance_FM", columns=[
                        "X", "Y", "Z", "PredInstance_FM", "species_id_FM",
                        "species_prob_FM", "PredScore_FM", "PredSemantic_FM"]),
                    dict(slug="SAT", instance="PredInstance_SAT", columns=[
                        "X", "Y", "Z", "PredInstance_SAT", "species_id_FM", "species_prob_FM"]),
                    dict(slug="ULS", instance="PredInstance_ULS", columns=[
                        "X", "Y", "Z", "PredInstance_ULS"]),
                ])


class OwnedTests(unittest.TestCase):
    def test_crossing_instances_return_every_point_once_with_exact_attributes(self):
        config = configuration()
        index = owned.build_index(config)
        source = laspy.read(FIXTURE)
        self.assertEqual(index["input_points"], len(source))
        self.assertLessEqual(index["audit"]["max_batch_points"], config["batch_points"])
        for layout in index["layouts"]:
            seen = []
            actual = []
            for number, window in enumerate(layout["windows"]):
                seen.extend(window["ids"])
                self.assertTrue(window["expanded"])
                output = BytesIO()
                audit = owned.stream_window(config, index, layout["slug"], number, output)
                self.assertTrue(audit["complete"])
                self.assertEqual(audit["expected_counts"], audit["observed_counts"])
                actual.append(np.frombuffer(output.getvalue(), dtype="<f8").reshape(-1, len(layout["columns"])))
            self.assertEqual(len(seen), len(set(seen)))
            valid = np.asarray(source[layout["instance"]]) > 0
            expected = np.column_stack([
                np.asarray(getattr(source, name.lower()) if name in ("X", "Y", "Z") else source[name])[valid]
                for name in layout["columns"]])
            actual = np.concatenate(actual)
            def ordered(rows):
                return rows[np.lexsort(rows.T)]
            np.testing.assert_array_equal(ordered(actual), ordered(expected))
            self.assertEqual(sum(item["count"] for item in layout["instances"]), int(valid.sum()))

    def test_missing_nodes_fail_completeness(self):
        config = configuration()
        config["core_size"] = 50
        config["buffer"] = 20
        index = owned.build_index(config)
        load = owned.hierarchy_nodes
        def omit(*args, **kwargs):
            nodes = list(load(*args, **kwargs))
            return iter(nodes[:-1])
        with patch.object(owned, "hierarchy_nodes", omit):
            with self.assertRaisesRegex(ValueError, "Incomplete instances"):
                owned.stream_window(config, index, "SAT", 0, BytesIO())

    def test_changed_input_is_rejected_before_read(self):
        config = configuration()
        index = owned.build_index(config)
        index["source_identity"]["mtime_ns"] += 1
        with self.assertRaisesRegex(ValueError, "changed since indexing"):
            owned.stream_window(config, index, "FM", 0, BytesIO())

    def test_single_large_node_fails_before_decompression(self):
        config = configuration()
        config["batch_points"] = 1
        with self.assertRaisesRegex(ValueError, "above batch limit"):
            owned.build_index(config)

    def test_core_boundary_ownership_and_retention_budget(self):
        instances = [dict(id=i, count=1500000, centroid=[x, 5],
                          mins=[x - 50, 0, 0], maxs=[x + 50, 10, 20])
                     for i, x in enumerate((0, 49, 50, 100), start=1)]
        windows = owned.plan_windows(instances, [0, 0], 50, 20, 2100000)
        self.assertEqual([item["owner_core"] for item in instances],
                         [[0, 0], [0, 0], [1, 0], [2, 0]])
        self.assertEqual(len(windows), 4)
        self.assertEqual(sorted(i for window in windows for i in window["ids"]), [1, 2, 3, 4])
        self.assertTrue(all(window["expanded"] for window in windows))



class ExplicitAttributeTests(unittest.TestCase):
    def test_minimal_and_explicit_columns_match_original_values(self):
        source = laspy.read(FIXTURE)
        for optional in ([], ['PredScore_FM', 'unused_1', 'Intensity']):
            config = configuration()
            config['layouts'] = config['layouts'][:1]
            config['core_size'] = 100
            config['buffer'] = 20
            columns = ['X', 'Y', 'Z', 'PredInstance_FM'] + optional
            config['layouts'][0]['columns'] = columns
            index = owned.build_index(config)
            result = BytesIO()
            audit = owned.stream_window(config, index, 'FM', 0, result)
            actual = np.frombuffer(result.getvalue(), dtype='<f8').reshape(-1, len(columns))
            requested = owned.requested_dimensions(source.header, columns)
            keep = np.asarray(source['PredInstance_FM']) > 0
            expected = np.column_stack([np.asarray(getattr(source, name.lower()) if name in 'XYZ'
                                                  else source[requested[name]])[keep] for name in columns])
            order = lambda rows: rows[np.lexsort((rows[:, 2], rows[:, 1], rows[:, 0]))]
            np.testing.assert_allclose(order(actual), order(expected), rtol=0, atol=1e-9)
            self.assertNotIn('CLASSIFICATION', audit['decompressed_groups'])
            self.assertEqual('INTENSITY' in audit['decompressed_groups'], bool(optional))

    def test_explicit_missing_attribute_is_rejected_before_scan(self):
        config = configuration()
        config['layouts'][0]['columns'].append('does_not_exist')
        with patch.object(owned, 'iter_batches', side_effect=AssertionError('scan started')):
            with self.assertRaisesRegex(ValueError, 'Missing requested.*does_not_exist'):
                owned.build_index(config)


if __name__ == "__main__":
    unittest.main()
