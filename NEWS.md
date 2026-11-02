# Unreleased (3Dtrees fork)

## Breaking changes

-   Inventory: `pred_score_fm` and `pred_score_mixed` are replaced by `pred_score_fm_values`, every distinct `PredScore_FM` of the tree with its point count (JSON), so an instance merged from several FM instances shows all their scores. Large streamed trees are counted over all points, not a sample.
-   CLI: `--segmentation-spec INSTANCE[,SPECIES,SPECIES_PROB]` is replaced by three options: `--instance-attribute` (repeat per segmentation), and `--species-attribute` and `--species-probability-attribute`, which belong to the `--instance-attribute` before them.
-   CLI (`exec/run.R`): one inventory table in two files, `inventory.csv` and `inventory.parquet` (same rows and columns, sorted by `(instance_dimension, instance_id)`, Parquet in 2,048-row groups, snappy), plus `stand_summary.csv` with one row per segmentation. Each tree row carries its DBH evidence: `dbh_fits` and `dbh_slices` as JSON arrays and `crown_hull` as WKT, after a `diagnostics` column. The per-segmentation files, the `species`, `validation` and `issues` tables, the `dbh/` folder, `inventory.gpkg` and the `--output-formats` option (TSV/JSON/GeoPackage copies) are gone; species composition and validation counts follow from the inventory and the stand summary.
-   CLI: `--inventory-read-mode` defaults to `auto` (it was `spatial`). One COPC file is read with `copc-owned` and 8 index workers, one plain LAS/LAZ file with the `sequential` scan; several files and `--enable-csp true` keep `spatial`. The spatial reader ran out of memory on large clouds (40 GB on a 457 M-point COPC crop, 60 GB on a 668 M-point LAZ). Pass `--inventory-read-mode spatial` for the previous reader.
-   CLI: `--inventory-partitions` defaults to about 3 M points per partition (a power of two from 64 to 4096, from the header point count); it was a fixed 64, which put 10 M points into each partition of a 668 M-point cloud.
-   CLI: the `sequential` reader accepts plain (non-COPC) LAS/LAZ, and reads the CRS with lidR when the file has no WKT record. Its `read_mode` in `run_metadata.json` is now `sequential_scan_and_instance_partitions`.

## Bug fixes

-   CLI: a generated DTM no longer fails with `[crop] extents do not overlap` when a tile's only points lie in its buffer (seen on a 668 M-point LAZ); such a tile is skipped like an empty one.
-   Container: the pinned lidR no longer fails with "Internal error in spatial index: point out of range" when a tree's extreme coordinate lies on a 0.1 m multiple that floating point rounds up (4021833.4 m: `floor(x / 0.1) * 0.1` is 4021833.4000000004). The grid bounds now always contain the points (`patches/lidR-grid-bounds.patch`). It stopped the inventory of a 668 M-point cloud after 70 min.
-   Container: the pinned rlas prints its read progress bar only to a terminal. A filtered (tiled) read reported NaN progress, printed as `-2147483648%`, and the bar grew job logs by megabytes (`patches/rlas-quiet-progress.patch`).

-   `forest_inventory()`: a rejected stem fit (spline position more than 5 m off) returned four values instead of five, so its quality flag (3) was reported as a DBH of 3 m. It now returns DBH = NA with quality flag 3.
-   `forest_inventory()`: a slice whose planarity or verticality contains NaN (duplicate points in a neighbourhood) no longer aborts the whole tree; NaN values are ignored in the quantile filter.

# CspStandSegmentation 0.2.1

## New Features and Enhancements

### Performance Improvements
-   Replaced R distance calculation function with a C++ implementation for faster computation in fds

# CspStandSegmentation 0.2.0

## New Features and Enhancements

### Stem Point Filtering for DBH Estimation

-   Added semantic segmentation support to the forest inventory function, allowing users to filter points by classification (e.g., stem-only points) for more accurate DBH estimation.

### Custom TreeID Column Support

-   **Flexible tree identification**: The `forest_inventory()` function now supports custom TreeID column names through the `tree_id_col` parameter
-   Added `non_tree_id` parameter to specify which values represent non-tree elements (can be a vector)

### Performance Improvements

-   Significantly optimized circle fitting algorithm with vectorized distance calculations
-   Implemented weighted point sampling based on local point density for better convergence
-   Reduced iterations to 500 and added early stopping criteria (30 angular segments)
-   Completely rewrote core `forest_inventory()` logic to use `data.table` instead of loops and foreach parallelization
-   Improved readability with helper functions: `.fit_circle()`, `.fit_circles()`, `.spline_predict()`
-   Eliminated parallelization overhead for single-core efficiency

### Stability and Quality Enhancements

-   **Improved height normalization**: Better handling of Z-coordinates when reference height (`Zref`) is available
-   **Enhanced quality control**: Added validation checks for unrealistic circle diameters and center positions
-   **Robust data types**: Fixed data.frame structure issues by converting merged results properly
-   **Better error handling**: Added validation for stem segmentation parameters and semantic column existence

### Bug Fixes

-   Fixed spline prediction height (now uses 1.3 + 0.5\*increment for better accuracy)
-   Corrected data.frame structure in inventory output
-   Fixed cosmetic issues and improved code formatting
-   Resolved issues with NA value handling in RANSAC circle fitting
