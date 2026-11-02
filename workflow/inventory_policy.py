"""Shared, versioned measurement limits used by R consumers and Python planners."""
import json
from pathlib import Path

POLICY = json.loads(Path(__file__).with_suffix(".json").read_text())
if POLICY["schema_version"] != 2 or any(type(v) is not int or v <= 0 for v in POLICY.values()):
    raise ValueError("Invalid inventory policy")
INSTANCE_POINT_LIMIT = POLICY["instance_point_limit"]
DBH_BAND_POINT_BUDGET = POLICY["dbh_band_point_budget"]
SAMPLE_LIMIT = POLICY["large_instance_sample_limit"]
