#!/usr/bin/env python3
"""Regenerate the small COPC integration fixture from inst/extdata/beech.las.

Requires laspy, numpy, pyproj, and PDAL with writers.copc. The generated COPC
is committed so the container test does not need a fixture-generation toolchain.
"""

import argparse
from pathlib import Path
import subprocess
import tempfile

import laspy
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--plain-wkt", type=Path, help="write plain LAS/LAZ with this WKT instead of COPC")
    args = parser.parse_args()
    source = laspy.read(args.source)
    sampled = source[::8]
    point_count = len(sampled.points)
    header = laspy.LasHeader(point_format=6, version="1.4")
    header.scales = source.header.scales.copy()
    header.offsets = source.header.offsets.copy()
    if args.plain_wkt:
        header.vlrs.append(laspy.vlrs.known.WktCoordinateSystemVlr(args.plain_wkt.read_text().strip()))
    else:
        from pyproj import CRS

        header.add_crs(CRS.from_epsg(32632))
    cloud = laspy.LasData(header)
    cloud.x = sampled.x
    cloud.y = sampled.y
    cloud.z = sampled.z
    cloud.intensity = sampled.intensity
    cloud.classification = sampled.classification
    cloud.return_number = np.ones(point_count, dtype=np.uint8)
    cloud.number_of_returns = np.ones(point_count, dtype=np.uint8)

    x = np.asarray(cloud.x)
    y = np.asarray(cloud.y)
    sat = np.where(x < np.median(x), 1, 2).astype(np.int32)
    fm = np.digitize(y, np.quantile(y, [1 / 3, 2 / 3])).astype(np.int32) + 11
    uls = np.where(x + y < np.median(x + y), 21, 22).astype(np.int32)

    dimensions = [
        ("PredInstance_SAT", "int32", sat),
        ("PredInstance_FM", "int32", fm),
        ("species_id_FM", "int32", (fm - 10).astype(np.int32)),
        ("species_prob_FM", "float32", np.full(point_count, 0.8, dtype=np.float32)),
        ("PredScore_FM", "float32", np.full(point_count, 0.9, dtype=np.float32)),
        ("PredSemantic_FM", "int8", np.where(np.arange(point_count) % 3 == 0, 1, 2).astype(np.int8)),
    ]
    dimensions.extend((f"unused_{index}", "int8", np.zeros(point_count, dtype=np.int8)) for index in range(1, 5))
    dimensions.append(("PredInstance_ULS", "int32", uls))
    for name, type_name, values in dimensions:
        cloud.add_extra_dim(laspy.ExtraBytesParams(name=name, type=type_name))
        cloud[name] = values

    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.plain_wkt:
        cloud.write(args.output)
        return
    with tempfile.TemporaryDirectory() as directory:
        intermediate = Path(directory) / "multi-instance.las"
        cloud.write(intermediate)
        subprocess.run(
            ["pdal", "translate", str(intermediate), str(args.output),
             "-w", "writers.copc", "--writers.copc.extra_dims=all"],
            check=True,
        )


if __name__ == "__main__":
    main()
