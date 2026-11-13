# CspStandSegmentation 1.0.0 (3Dtrees fork)

## Breaking changes

-   Segmentations: `--score-attribute` and `--non-tree-ids` belong to the `--instance-attribute` before them, like the species options; non-tree IDs are no longer one list for all segmentations (default 0 each). Any segmentation with a score attribute gets `pred_score_values` (every distinct score with its point count); this replaces the FM-only `pred_score_fm_values`. The `PredSemantic_FM` wood and leaf counts (`wood_point_count`, `leaf_point_count`, `wood_share`, `leaf_share`) are removed.
-   Inventory: every tree is collected the same way and measured once complete. DBH, height, crown and position use one point per tree and `--metrics-grid` cell (default 0.01 m); `point_count`, species, probability and FM fields count every input point; new column `metrics_point_count`. The separate summary for trees above two million points, `inventory_method`, `diagnostics`, `--inventory-retained-points` and the COPC window retention budget are gone (one window per core). `workflow/resource_request.py` and the partition audit helper are removed.
-   Outputs: `inventory.csv` and `inventory.parquet` (same table) and `stand_summary.csv`, plus the DTM. The crown hull is a GeoJSON Polygon geometry in the `crown_hull` column, next to `dbh_fits` and `dbh_slices`; `stand_summary.csv` names the CRS in a `crs` column. `manifest.json`, `run_metadata.json`, `resource_summary.json` and `instance_index.json` are written only with `--debug-outputs true`.
-   CSP: each CSP tree row carries the seed it grew from in `seed_x`, `seed_y` and `seed_height_m` (height above the terrain), right after `z`; other segmentations have no seed columns. `effective_seeds.tsv` is written only with `--debug-outputs true`.
-   Inventory: no point limits or memory guards. Every stem-band point is measured; `--dbh-band-points`, `dbh_band_source_points` and `dbh_band_voxel_m` are gone. Streamed large trees use exact medians instead of a 100,000-point sample. Size checks on instance ids, COPC nodes and hierarchy, window retention and RAM no longer abort a run.
-   CLI: `--segmentation-spec INSTANCE[,SPECIES,SPECIES_PROB]` is replaced by three options: `--instance-attribute` (repeat per segmentation), and `--species-attribute` and `--species-probability-attribute`, which belong to the `--instance-attribute` before them.
-   CLI (`exec/run.R`): one inventory table in two files, `inventory.csv` and `inventory.parquet` (same rows and columns, sorted by `(instance_dimension, instance_id)`, Parquet in 2,048-row groups, snappy), plus `stand_summary.csv` with one row per segmentation. Each tree row carries its DBH evidence: `dbh_fits` and `dbh_slices` as JSON arrays and `crown_hull` as WKT, after a `diagnostics` column. The per-segmentation files, the `species`, `validation` and `issues` tables, the `dbh/` folder, `inventory.gpkg` and the `--output-formats` option (TSV/JSON/GeoPackage copies) are gone; species composition and validation counts follow from the inventory and the stand summary.
-   CLI: `--inventory-read-mode` defaults to `auto` (it was `spatial`). One COPC file is read with `copc-owned` and 8 index workers, one plain LAS/LAZ file with the `sequential` scan; several files and `--enable-csp true` keep `spatial`. The spatial reader ran out of memory on large clouds (40 GB on a 457 M-point COPC crop, 60 GB on a 668 M-point LAZ). Pass `--inventory-read-mode spatial` for the previous reader.
-   CLI: `--inventory-partitions` defaults to about 3 M points per partition (a power of two from 64 to 4096, from the header point count); it was a fixed 64, which put 10 M points into each partition of a 668 M-point cloud.
-   CLI: the `sequential` reader accepts plain (non-COPC) LAS/LAZ, and reads the CRS with lidR when the file has no WKT record. Its `read_mode` in `run_metadata.json` is now `sequential_scan_and_instance_partitions`.

## Performance

-   Inventory memory: the 1 cm grid points are kept as the LAS integers (12 instead of 32 bytes per point), batches are collected in a list instead of copying a tree's grid each time, and heights above the terrain are looked up by DTM cell index instead of `terra::extract`, which needed about 1 GB per million points. Same results byte for byte; the 457 M-point crop needs 17.1 instead of 23.5 GiB and 13.3 instead of 15.3 min, the 2 B-point dataset 1887 26 instead of 30 GiB and 30 instead of 35 min at 10 cores.

-   Inventory speed: a pool of COPC decoders reads each window (reading kept ~3.6 of 10 cores busy before, ~6 now); a tree with at least 2 M grid points computes its geometry and seeded fits on cores handed back by workers that have finished; each worker measures its largest trees first; windows of all segmentations run side by side while their memory, estimated from the index (32 bytes per owned point), fits `GALAXY_MEMORY_MB`, and a window needing more than half of it runs alone with all workers. A window that would need more than 80 % of `GALAXY_MEMORY_MB` is split into windows of whole trees, so peak memory follows the grant instead of the largest 50 m cell. The seeded DBH fits get only the points upstream measures (the stem band, the height and x/y extremes and the crown hull above 0.5 m, on the coordinates lidR quantizes) instead of a copy of the whole tree: the largest 1924 tree (22 M grid points) needs 5.5 instead of 13.5 GiB, or 30 GiB with five parallel fits, and is measured in 43 instead of 77 s. With a fixed 32 GB grant at 10 cores, 1924 (4.3 B points) peaks at 25 GiB, 1887 at 18 GiB and 2623 at 8 GiB. Same results byte for byte; with the DTM supplied at 10 cores 1887 takes 17 instead of 21 min, 1924 61 instead of 76 min, 2518 52 instead of 74 min.
-   Container: release tags also publish a Singularity image for Galaxy, built from the pushed Docker digest: `oras://ghcr.io/3dtrees-earth/3dtrees_csp_standsegmentation:<version>-sif`. Docker stays for local use.
-   Container: rlas is patched to report LAS headers of more than 2^31 points (`patches/rlas-64bit-point-count.patch`, as in the Raster tool), so lidR reads the header of large merged clouds.

## Bug fixes

-   DTM: a point cloud with a 3D projected CRS (for example "Promoted to 3D from EPSG:32616") got its DTM CRS written to a `.aux.xml` side file, and the run failed its output check. DTMs are now written with the 2D CRS, and GDAL side files are switched off.

-   CSP: `segmented_csp.laz` keeps every dimension of the input; it held only the dimensions the run read (instance, species and requested attributes).
-   CSP: points are normalized like in the inventory readers: the same DTM cell lookup and the written Float32 DTM. lidR assigned points lying exactly on a DTM cell edge to the neighbouring cell, and the in-memory DTM kept double precision, so with CSP on, the same FM or SAT tree could get a slightly different stem band and position.
-   DTM: `dtm_full.tif` and `dtm_aoi.tif` are written with DEFLATE instead of LZW, so readers without LZW support (such as Galaxy's image assertions) can open them.
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
