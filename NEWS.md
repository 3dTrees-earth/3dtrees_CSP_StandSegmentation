# Unreleased (3Dtrees fork)

## Breaking changes

-   CLI (`exec/run.R`): every table is now written as Parquet and the DBH crown hulls as GeoParquet 1.1; `--output-formats` defaults to `parquet` (it was `csv,tsv,json,gpkg`). Parquet is always written; `csv`, `tsv`, `json` and `gpkg` add copies of the same tables. Pass `--output-formats parquet,csv,tsv,json,gpkg` for the previous file set. Tables are sorted by `(instance_dimension, instance_id)` and Parquet row groups hold at most 16,384 rows.

## Bug fixes

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
