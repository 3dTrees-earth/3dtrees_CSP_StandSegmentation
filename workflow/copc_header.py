#!/usr/bin/env python3
"""Read a merged COPC's 64-bit count, dimensions, bounds and WKT header."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import laspy


def read_header(path: Path) -> dict:
    with laspy.open(path, laz_backend=laspy.LazBackend.Lazrs) as reader:
        header = reader.header
        if not any(vlr.user_id == "copc" and vlr.record_id == 1
                   for vlr in header.vlrs):
            raise ValueError("Sequential inventory requires a COPC input")
        records = list(header.vlrs) + list(header.evlrs or [])
        wkts = [vlr.record_data_bytes().rstrip(b"\x00").decode("utf-8")
                for vlr in records if vlr.user_id == "LASF_Projection"
                and vlr.record_id == 2112]
        if not wkts or len(set(wkts)) != 1:
            raise ValueError("COPC has no unambiguous WKT CRS record")
        count = int(header.point_count)
        bounds = [float(header.mins[0]), float(header.maxs[0]),
                  float(header.mins[1]), float(header.maxs[1])]
        if count < 1 or not bounds[0] < bounds[1] or not bounds[2] < bounds[3]:
            raise ValueError("COPC has invalid point count or XY bounds")
        return {
            "input": str(path), "point_count": count,
            "xmin": bounds[0], "xmax": bounds[1],
            "ymin": bounds[2], "ymax": bounds[3],
            "crs_wkt": wkts[0],
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
