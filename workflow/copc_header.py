#!/usr/bin/env python3
"""Read a merged LAS/LAZ's 64-bit count, dimensions, bounds and WKT CRS (if any)."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import laspy


def read_header(path: Path) -> dict:
    with laspy.open(path, laz_backend=laspy.LazBackend.Lazrs) as reader:
        header = reader.header
        records = list(header.vlrs) + list(header.evlrs or [])
        wkts = [vlr.record_data_bytes().rstrip(b"\x00").decode("utf-8")
                for vlr in records if vlr.user_id == "LASF_Projection"
                and vlr.record_id == 2112]
        wkts = [wkt for wkt in wkts if wkt.strip()]
        if len(set(wkts)) > 1:
            raise ValueError("Point cloud has conflicting WKT CRS records")
        count = int(header.point_count)
        bounds = [float(header.mins[0]), float(header.maxs[0]),
                  float(header.mins[1]), float(header.maxs[1])]
        zmin, zmax = float(header.mins[2]), float(header.maxs[2])
        if (count < 1 or not all(map(math.isfinite, [*bounds, zmin, zmax]))
                or not bounds[0] < bounds[1] or not bounds[2] < bounds[3] or zmin > zmax):
            raise ValueError("Point cloud has invalid point count or bounds")
        return {
            "input": str(path), "point_count": count,
            "xmin": bounds[0], "xmax": bounds[1],
            "ymin": bounds[2], "ymax": bounds[3],
            "zmin": zmin, "zmax": zmax,
            "crs_wkt": wkts[0] if wkts else None,
            "available_dimensions": list(header.point_format.dimension_names),
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(read_header(args.input), indent=2) + "\n")


if __name__ == "__main__":
    main()
