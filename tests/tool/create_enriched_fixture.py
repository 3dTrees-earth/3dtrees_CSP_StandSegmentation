#!/usr/bin/env python3
"""Create a deterministic multi-attribute LAZ fixture from mikro_segmented.laz."""

from pathlib import Path
import sys

import laspy
import numpy as np


def add_dimension(las: laspy.LasData, name: str, type_: str, description: str) -> None:
    las.add_extra_dim(laspy.ExtraBytesParams(name=name, type=type_, description=description))


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: create_enriched_fixture.py INPUT.laz OUTPUT.laz")
    input_path, output_path = map(Path, sys.argv[1:])
    las = laspy.read(input_path)
    source_ids = np.asarray(las.PredInstance, dtype=np.int32)
    valid = source_ids > 0

    add_dimension(las, "PredInstance_FM", "int32", "Synthetic FM instance")
    add_dimension(las, "species_id_FM", "int32", "Synthetic species ID")
    add_dimension(las, "species_prob_FM", "float32", "Synthetic species probability")
    add_dimension(las, "PredScore_FM", "float32", "Synthetic FM score")
    add_dimension(las, "PredSemantic_FM", "int8", "Synthetic FM semantic class")

    las.PredInstance_FM = source_ids
    species = np.where(valid, source_ids % 3 + 1, -1).astype(np.int32)
    probability = np.where(valid, 0.65 + (source_ids % 4) * 0.05, -1).astype(np.float32)
    score = np.where(valid, 0.70 + (source_ids % 5) * 0.04, -1).astype(np.float32)
    semantic = np.where(valid, np.where(np.arange(len(las.points)) % 3 == 0, 1, 2), 0).astype(np.int8)

    first_id = next((value for value in np.unique(source_ids) if value > 0), None)
    if first_id is not None:
        indices = np.flatnonzero(source_ids == first_id)
        if len(indices) >= 3:
            species[indices[-1]] = species[indices[0]] + 10
            probability[indices[-1]] = 0.95
            score[indices[-1]] = 0.15

    las.species_id_FM = species
    las.species_prob_FM = probability
    las.PredScore_FM = score
    las.PredSemantic_FM = semantic
    las.write(output_path)


if __name__ == "__main__":
    main()
