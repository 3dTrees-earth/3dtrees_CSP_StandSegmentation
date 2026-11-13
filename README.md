## CspStandSegmentation is an R-package for the segmentation of single trees from forest point clouds scanned with terrestrial, mobile or unmanned LiDAR systems <img src="https://github.com/JulFrey/CspStandSegmentation/blob/main/inst/figures/csp_logo.png" align="right" width = 300/>

Authors: Julian Frey and Zoe Schindler, University of Freiburg, Chair of Forest Growth and Dendroecology


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17294732.svg)](https://doi.org/10.5281/zenodo.17294732)  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## 3Dtrees command-line workflow

This fork retains the upstream segmentation implementation and adds a headless
CLI and container for a future 3Dtrees Galaxy tool. It accepts a supplied
GeoTIFF DTM or generates one and can
inventory any number of existing instance dimensions. CSP segmentation is
optional; its output preserves the original points and explicitly selected
attributes and adds `PredInstance_CSP`.

The default inventory read is **XYZ + the instance dimension**. Optional
species/probability fields are read only when named with `--species-attribute` and `--species-probability-attribute`.
Each segmentation can also name its prediction score attribute with
`--score-attribute` (for example `PredScore_FM`) and its own non-tree IDs with
`--non-tree-ids` (default 0); like the species options, both belong to the
`--instance-attribute` before them. `--inventory-attributes` accepts a comma-separated list of numeric,
scalar LAS attributes and can be repeated; requested missing attributes fail
validation. Arbitrary selected attributes are available to the inventory path
and preserved in CSP point-cloud output; selection alone does not define new
per-tree aggregations or output columns. Standard attributes use R/LAS names,
for example `Intensity`, `Classification`, `ReturnNumber`, or `gpstime`.

Classification is read automatically only by a DTM-generation stage that needs
it, not by the inventory reader. Unrequested extra fields are removed before
points reach the R inventory calculations or point transport.

COPC/LAZ decompression operates on compressed groups: extra attributes share
one `ALL_EXTRA_BYTES` group, so selecting an instance dimension may still
require decoding unrelated extra bytes inside the bounded native decoder
buffer. Only selected columns are copied into the working point tables. This
is not a claim of per-attribute physical decompression. LASlib's spatial reader
can address the first nine extra-byte records individually; later fields require
reading the extra-byte block and immediately projecting it. The container's
pinned `rlas` patch ensures those later fields are available. Metadata reports
retained dimensions separately from decompressed groups.

For example, the minimal existing-instance inventory is:

```bash
Rscript exec/run.R --input merged.copc.laz --output-dir results \
  --instance-attribute PredInstance_FM --dtm terrain.tif
```

The default `--inventory-read-mode auto` keeps memory bounded whatever the
cloud size: one COPC file is read with `copc-owned` (8 index workers unless
`--inventory-index-workers` says otherwise), and one plain LAS/LAZ file with
the `sequential` scan, so it needs no COPC conversion. Several input files and
`--enable-csp true` use the `spatial` reader, which loads whole tiles and needs
memory in proportion to the points per tile. `run_metadata.json` records the
reader that was used as `read_mode`.

Include the selected column count (XYZ + instance + optional
fields) in resource estimation through `--inventory-columns` when it exceeds
the estimator's conservative eight-column default.

```bash
Rscript exec/run.R \
  --input merged.copc.laz \
  --output-dir results \
  --instance-attribute PredInstance_SAT \
    --species-attribute species_id_SAT --species-probability-attribute species_prob_SAT \
  --instance-attribute PredInstance_FM \
    --species-attribute species_id_FM --species-probability-attribute species_prob_FM \
  --enable-csp false \
  --dtm existing_dtm.tif \
  --dtm-resolution 0.2 \
  --random-seed 42
```

`--instance-attribute` is repeatable, one per segmentation.
`--species-attribute` and `--species-probability-attribute` belong to the
`--instance-attribute` before them and are optional, but come as a pair. A
segmentation can borrow another model's species, for example RayCloudTools
trees with FM species:
`--instance-attribute PredInstance_RCT --species-attribute species_id_FM
--species-probability-attribute species_prob_FM`. Without species attributes,
the inventory omits species columns. Run `Rscript exec/run.R --help` for all controls.

Common controls are the merged COPC or LAS/LAZ input, repeatable segmentation specs with their species, score and non-tree IDs,
optional CSP and seed source, optional native-CRS AOI GeoJSON, an optional DTM,
DTM resolution (default 0.2 m when generated), and random seed.
The source point cloud must declare a projected metre CRS. The CLI can recover
the CRS from a valid extended WKT record when a COPC file has a truncated
ordinary WKT record, and rejects conflicting valid records. A supplied DTM must
use the same CRS or an equivalent horizontal CRS and cover the full cloud
extent. When the GeoTIFF omits a vertical datum, its heights must already use
the point cloud's vertical datum. The supplied DTM is copied to
`dtm_full.tif` in the output. With no `--dtm`, the tool generates that file
(DEFLATE-compressed, no predictor). Both inventory paths look up the terrain
under each point the same way, in the cell that contains it.
Missing DTM cells beneath inventoried points stop the run instead of silently
dismissing those points.
Fine-tuning controls retain upstream defaults,
including a 0.3 m CSP voxel and one routing worker. CSP geometry features are
computed only when a non-zero geometry weight requires them.

### Outputs

**Results:** `inventory.csv` and `inventory.parquet` (the same table, one row
per tree, with the DBH runs, slices and crown hull embedded as JSON),
`stand_summary.csv` and the DTM. Run records are written only with
`--debug-outputs true`.

| File | One row per | Written |
| --- | --- | --- |
| `inventory.csv`, `inventory.parquet` | tree and segmentation | always |
| `stand_summary.csv` | segmentation (instance dimension) | always |
| `dtm_full.tif`, `dtm_aoi.tif` | DTM (supplied or generated) | always; AOI raster only with `--aoi-json` |
| `segmented_csp.laz` | point (every input dimension plus `PredInstance_CSP`) | only with `--enable-csp true` |
| `effective_seeds.tsv` | CSP seed | only with `--enable-csp true --debug-outputs true`; the seeds are also in the inventory |
| `manifest.json` | artifact keys and formats, fields and units, coordinate reference, segmentations, DTM | only with `--debug-outputs true` |
| `run_metadata.json`, `resource_summary.json` | run settings, timings and memory | only with `--debug-outputs true` |
| `instance_index.json` | instance index | only with `--debug-outputs true` and the `copc-owned` reader |

**Row order.** Rows are sorted by `(instance_dimension, instance_id)`; strings
sort in C (byte) order, so `PredInstance_FM` comes before `PredInstance_SAT`.
`inventory.parquet` is written in row groups of at most 2,048 rows with min/max
statistics, snappy-compressed (no zstd), so a reader fetches one tree by reading
only the row group that covers it, and only the columns it needs.

**Coordinates.** `x`, `y`, `z`, the slice circles and the crown hulls are in
the source CRS, named in the `crs` column of `stand_summary.csv` (for example
`EPSG:2056`); `longitude`/`latitude` are WGS84.

**Types.** In `inventory.csv` (`data.table::fwrite`, RFC 4180 quoting), `NA` is
an empty field, booleans are `TRUE`/`FALSE` and doubles are printed to 15
significant digits. The JSON evidence cells contain commas and quotes, so read
the CSV with a reader that follows RFC 4180 quoting (`utils::read.csv`, pandas,
Excel); `data.table::fread` keeps the doubled quotes.

#### Inventory columns (`inventory.csv`, `inventory.parquet`)

Units are also listed per field in `manifest.json`.

| Column | Unit | Meaning |
| --- | --- | --- |
| `segmentation_source` | | Short segmentation name (`SAT`, `FM`, `CSP`, ...): the instance dimension without `PredInstance_` |
| `instance_dimension` | | LAS dimension holding the tree IDs; the Potree attribute to filter by |
| `instance_id` | | Tree ID in that dimension |
| `x`, `y`, `z` | m (source CRS) | Tree position from the DBH measurement; median x/y and lowest z of the tree's grid points when no DBH was measured |
| `seed_x`, `seed_y`, `seed_height_m` | m (source CRS; height above terrain) | CSP rows only: the seed the tree grew from; the same X, Y and Z a `--seed-file` takes |
| `longitude`, `latitude` | degree (WGS84) | `x`, `y` transformed to EPSG:4326 |
| `dbh_m` | m | Diameter at breast height (median of the DBH runs); null when unavailable |
| `height_m` | m | Tree height above the DTM |
| `convex_hull_area_m2` | m² | Crown convex-hull area |
| `point_count` | | Input points in the tree (every point, full resolution) |
| `metrics_point_count` | | Grid points the geometry was measured on (`--metrics-grid`) |
| `species_id`, `species_prob` | | Modal species ID and median probability of the points with that ID (only with species dimensions) |
| `species_conflict_fraction` | | Share of the tree's points with another species ID |
| `scientific_name`, `common_name`, `taxonomy_version`, `model_version` | | From `--species-lookup`, when supplied |
| `pred_score_values` | | Every distinct value of the segmentation's `--score-attribute` in the tree with its point count as JSON (`[{"score":0.83,"point_count":12000}, ...]`, most points first; counted over all points). An instance merged from several predicted instances shows each one's score. Only for segmentations with a score attribute |
| `dbh_quality_code`, `dbh_quality_label` | | 1 `too_few_stem_slice_points`, 2 `fallback_circle_estimate`, 3 `rejected_or_unstable_fit`, 4 `spline_fit`, 9 no measurement; the label is `dbh_unavailable` whenever `dbh_m` is null |
| `dbh_runs`, `dbh_runs_fitted`, `dbh_spread_m` | –, –, m | DBH runs made, runs that produced a DBH, and max − min of their DBH values |
| `dbh_band_points` | – | Grid points between `--slice-min` and `--slice-max`, all used for DBH |
| `dbh_fits` | JSON | Every DBH run: array of `{run, seed, dbh_m, quality_code, reported}`; `reported` is true for the run whose DBH the row reports, `dbh_m` null when the run fitted none |
| `dbh_slices` | JSON | The stem slices of the reported run: array of `{slice_id, z_lo_m, z_hi_m, ground_z, n_points, n_fit_points, circle_x, circle_y, radius_m, status, spline_x, spline_y, spline_r}` |
| `crown_hull` | GeoJSON | The crown as a Polygon geometry (`{"type":"Polygon","coordinates":[[[x,y],...]]}`, source CRS, mm): convex hull of the grid points above 0.5 m, or of all points when fewer than three are; empty for trees with fewer than three points |

The table has the union of the per-segmentation columns; a column that a
segmentation does not produce is empty for its rows. The last three columns are
the evidence the measurements were computed from; `--dbh-evidence false`
leaves them empty and changes no value.

In `dbh_slices`, `z_lo_m`/`z_hi_m` are above the DTM, `ground_z` and
`circle_x`/`circle_y` are absolute source-CRS coordinates (m), `n_fit_points`
are the points handed to RANSAC, `radius_m` is the raw fit, and
`spline_x`/`spline_y`/`spline_r` are the smoothing splines at the slice's lower
edge, the abscissa the 1.3 m prediction uses (null when no spline was fitted).
`status` names the upstream rule that set a slice aside: `fitted`, `outlier_2sd`
(radius more than two standard deviations from the median), `growth_1_3x`
(radius above 1.3 times the smallest radius below it), `radius_below_0_04m`,
`too_few_points` (fewer than three points) or `no_fit` (RANSAC found no circle).
Upstream applies the growth and minimum-radius rules from the second slice on, so the lowest slice is never set aside by them.

#### Stand summary (`stand_summary.csv`)

One row per segmentation: `segmentation_source`, `instance_dimension`,
`area_source` (`aoi` or `point_cloud_convex_hull`), `area_m2` (m²),
`tree_count`, `dbh_measured_tree_count`, `dbh_coverage_fraction`,
`trees_per_ha` (trees/ha), `basal_area_m2_per_ha` (m²/ha; null unless every
tree has a DBH), `mean_dbh_m`, `median_dbh_m`, `dbh_p05_m`, `dbh_p25_m`,
`dbh_p75_m`, `dbh_p95_m`, `mean_height_m`, `median_height_m`, `height_p05_m`,
`height_p25_m`, `height_p75_m`, `height_p95_m` (m), and `crs`, the source CRS of
all coordinates (`EPSG:<code>`, or WKT when the CRS has no EPSG code).

If species predictions are present, `--species-lookup lookup.tsv` can add names
and taxonomy/model versions. The lookup needs `segmentation_source`,
`species_id`, `scientific_name`, `taxonomy_version`, and `model_version` columns;
`common_name` is optional. Without a lookup the inventory has species IDs
only; no names are invented.

Processing is fail-atomic: the final output directory is published only after
all requested products succeed.

Existing-instance inventory keeps upstream's 500 RANSAC iterations but uses a
quiet, vectorized implementation of the same Pratt circle equations. It also
projects only the point attributes required by each requested segmentation and
skips the full preservation copy when CSP output is disabled. On the local
5.9-million-point GFZ benchmark these changes reduced tool time from 59.4 to
41.4 seconds and process peak RSS from 4.85 to 3.19 GiB.

With selective reading and the standardization `rlas` patch, the same GFZ file
successfully loaded 14 extra-byte attributes and inventoried SAT and FM together
in 53.7 seconds at 2.14 GiB process peak RSS. Sampled CPU averaged 100.1%,
confirming that the one-thread default consumes approximately one core. The run
produced 64 tree rows plus DTM, stand, and species products without creating a
point-cloud output.

Inventory-only execution uses two bounded passes. DTM generation defaults to
300 m tiles with a 5 m buffer and up to 10 workers. Below 50 million points the
tiles read the source directly. Larger inputs are scanned once in parallel by
point range, retaining the minimum Z at each deterministic 0.1 m cell centre;
CSF and TIN rasterization then run on that reduced surface in parallel tiles.
Each spatial worker uses one lidR thread, avoiding nested oversubscription. The
result raster retains the aligned input extent, with unsupported edge areas as
NoData.

The second pass collects every tree the same way, whatever its size, and
measures it once it is complete:

- **Points.** `point_count`, the species vote (`species_id`, the median
  `species_prob` of its points, `species_conflict_fraction`) and the FM fields
  count every input point.
- **Geometry.** DBH, height, crown area and position use one point per tree and
  `--metrics-grid` cell (default 0.01 m, world-aligned, absolute Z). The point
  kept is the smallest by X, Y, Z, so the result does not depend on file,
  window or worker order. `metrics_point_count` is the number of grid points;
  `--metrics-grid 0` measures every point.
- **DBH.** Upstream `forest_inventory()` runs `--dbh-runs` times (default 5)
  with seeds derived from the base seed, segmentation and tree ID; the run
  with the median DBH is reported.

The readers only decide how a tree's points arrive:

- `copc-owned` (the `auto` choice for one COPC file) indexes the COPC once
  (counts, centroids, extents), then reads one window per 50 m core
  (`--inventory-core-size`): the core plus a 20 m halo (`--inventory-buffer`),
  grown to the full extent of every tree whose centroid is in the core. One
  Python reader decodes the window's COPC nodes and pipes each owned tree's
  points to one of `--inventory-workers` R workers; counts are checked against
  the index before anything is published. No point files are written.
- `sequential` (one plain LAS/LAZ) and `spatial` (several files, and with CSP)
  write the selected columns once to instance partitions (ID modulo the
  partition count), then collect the trees of each partition in parallel.

There are no point limits or memory guards: memory follows the grid clouds
of the trees in flight, so enforce a container RAM cap.

CSP segmentation keeps the full-cloud path, because upstream routing needs the
whole cloud, and then measures every segmentation with the same tree
collection.

Build and run the container based on Rocker geospatial 4.6.1 with pinned
current `rlas` and `lidR` source revisions:

```bash
docker build -t 3dtrees-csp .
docker run --rm -v "$PWD:/work" -w /work 3dtrees-csp \
  Rscript /opt/CspStandSegmentation/exec/run.R --help
```

The checked-in `tests/tool/fixtures/multi-instance.copc.laz` contains 29,011
points and three synthetic instance dimensions. `PredInstance_ULS` is deliberately
the eleventh extra-byte dimension. Regenerate it in an environment with laspy,
NumPy, pyproj, and PDAL `writers.copc` support:

```bash
python tests/tool/generate_multi_instance_copc.py \
  inst/extdata/beech.las tests/tool/fixtures/multi-instance.copc.laz
```

The container test exercises those three instance dimensions, generated and
supplied DTMs, and both inventory-only and CSP-enabled execution:

```bash
docker run --rm -v "$PWD:/work" -w /work 3dtrees-csp \
  Rscript tests/tool/test_end_to_end.R
```

`tests/tool/test_table_outputs.R` checks the result files on 35,000
synthetic trees:

```bash
docker run --rm -v "$PWD:/work" -w /work 3dtrees-csp \
  Rscript tests/tool/test_table_outputs.R
```

### Inventory implementation boundaries

- `workflow/instance_inventory.R` owns normalization, the per-tree collection
  (`new_tree`, `add_tree_points`, `finish_tree`), the metrics grid, the DBH
  runs and the instance seeds.
- `workflow/partition_inventory.R` and `workflow/owned_inventory.R` deliver
  points to the trees and check completeness.
- `workflow/inventory_readers.R` prepares the selected reader and exposes input
  statistics, the stand hull and a measurement operation.
- `workflow/inventory_run.R` handles DTM preparation, exports and run metadata.
  `workflow/outputs.R` builds the tables and manifest; `workflow/table_formats.R`
  writes CSV, Parquet and GeoJSON.

## Upstream installation history
The upstream CRAN package was archived in June 2026 with `lidR`. The container
above installs pinned source revisions of `rlas` and `lidR`.

The older upstream instructions below are retained for reference and may no
longer work with current CRAN:
```R
install.packages("CspStandSegmentation")
```

### Latest version from GitHub:
If you are working on Windows operating systems, you will need to install Rtools prior to installation: https://cran.r-project.org/bin/windows/Rtools/>. On Mac, Xcode is required. 

```R
install.packages(c('devtools', 'Rcpp', 'lidR', 'dbscan', 'igraph', 'foreach', 'doParallel','magrittr', 'data.table'))

devtools::install_github('https://github.com/JulFrey/CspStandSegmentation')

# Check if it is working
library(CspStandSegmentation)
example("csp_cost_segmentation", run.dontrun = TRUE)

```

## Usage
The package is firmly based on the `lidR` package and uses the las file structure. Smaller point clouds can be directly segmented using the ```csp_cost_segmentation``` function. This requires a set of tree positions (map) as starting points, which can be derived using the ```find_base_coordinates_raster``` function, which might require parameter optimization. Theoretically, tree positions might also come from field measurements or manual assignments.:

```R
# read example data
file = system.file("extdata", "beech.las", package="CspStandSegmentation")
tls = lidR::readTLSLAS(file)

# find tree positions as starting points for segmentation
map <- CspStandSegmentation::find_base_coordinates_raster(tls)

# segment trees
segmented <- tls |>
  CspStandSegmentation::add_geometry(n_cores = parallel::detectCores()/2) |>
  CspStandSegmentation::csp_cost_segmentation(map, 1, N_cores = parallel::detectCores()/2)

# show results
lidR::plot(segmented, color = "TreeID")

# create inventory
inventory <- CspStandSegmentation::forest_inventory(segmented)
head(inventory)
lidR::plot(segmented, color = "TreeID") |> CspStandSegmentation::plot_inventory(inventory)
```

For large areas, the package can be used within the lidR LAScatalogue engine to cope with memory limitations. The following example shows how this can be done. The single tiles of segmented trees are saved in a folder in this example and merged afterwards. 

```R
# packages
library(lidR)
library(CspStandSegmentation)

# parameters
las_file <- "your_file.laz"
base_dir <- "~/your_project_folder/" # with trailing /
cores <- parallel::detectCores()/2 # number od cpu cores 
res <- 0.3 # voxel resolution for segmentation
chunk_size <- 50 # size of one tile in m excl. buffer
chunk_buffer <- 10 # buffer around tile in m

# main

# create dir for segmentation tiles
if(!dir.exists(paste0(base_dir,"segmentation_tiles/"))) {
  dir.create(paste0(base_dir,"segmentation_tiles/"))
}

uls = lidR::readTLSLAScatalog(paste0(base_dir,las_file), select = "XYZ0", chunk_size = chunk_size, chunk_buffer = chunk_buffer)
plot(uls, chunk_pattern = TRUE)
# plot(dtm,add = TRUE)
# sf::as_Spatial(sf::st_as_sf(map, coords = c("X", "Y"))) |> plot(add = TRUE)

opt_output_files(uls) <- paste0(base_dir,"segmentation_tiles/{ID}")
segmented <- catalog_apply(uls, function(cluster) {
  
  las <- suppressWarnings(readLAS(cluster)) # read files
  if (is.empty(las) ) return(NULL) # stop if empty
  message(str(cluster))
  # find tree positions as starting points for segmentation
  map <- CspStandSegmentation::find_base_coordinates_raster(las)
  
  # add the tile ID*100,000 to the TreeID to ensure unique IDs across all tiles
  map$TreeID <- map$TreeID + as.numeric(basename(cluster@save)) * 100000
  
  # only use seed positions within the tile+buffer and save the tile bbox to only return tree pos within the tile (excl. buffer)
  inv <- map
  invb <- map
  # the bbox includes the buffer, the bbbox excludes the buffer 
  bbox <- cluster@bbox
  bbbox <- cluster@bbbox
  inv <- inv[inv$X < bbox[1,2] & inv$X > bbox[1,1] & inv$Y < bbox[2,2] & inv$Y > bbox[2,1],]
  invb <- invb[invb$X < bbbox[1,2] & invb$X > bbbox[1,1] & invb$Y < bbbox[2,2] & invb$Y > bbbox[2,1],]
  if (nrow(inv) == 0) return(NULL) # stop if no tree pos in tile found
  if (is.empty(las) ) return(NULL) # stop if empty
  
  # Assign all points to trees
  las <- las |> add_geometry(n_cores = cores) |> csp_cost_segmentation(invb,res, N_cores = cores, V_w = 0.5)
  # las <- las |> csp_cost_segmentation(map,res, N_cores = cores, V_w = 0.5) # this is a faster version which does not make use of the geometric feature weights
  if (is.empty(las)) return(NULL)
  
  las <- las |>  filter_poi(TreeID %in% c(0,inv$TreeID)) # only return trees within the tile
  if (is.empty(las)) return(NULL) # stop if empty
  
  # remove unneccesary attributes for further processing 
  las <- las |> remove_lasattribute('x_vox') |> 
    remove_lasattribute('y_vox') |> 
    remove_lasattribute('z_vox') |> 
    remove_lasattribute('buffer') |>
    remove_lasattribute('Linearity') |>
    remove_lasattribute('Sphericity') |>
    remove_lasattribute('Verticality')
  
  # validate las
  las <- las  |>  las_quantize()  |> las_update()
  if (is.empty(las)) return(NULL)
  return(las)
}, .options = list(automerge = TRUE))

# merge segmented trees
segmented <- readTLSLAScatalog(paste0(base_dir,"segmentation_tiles/"), select = "xyz0", chunk_buffer = 0)
opt_merge(segmented) <- TRUE
opt_output_files(segmented) <- paste0("")
segmented <- catalog_apply(segmented, function(cluster) {
  las <- suppressWarnings(readLAS(cluster)) # read files
  if (is.empty(las) ) return(NULL) # stop if empty
  return(las)
}, .options = list(automerge = TRUE))

# write results in a single file
writeLAS(segmented, paste0(base_dir,"segmented.las"))
```

## Citation
If you publish work related to CspStandSegmentation please cite the following article:

Larysch E, Frey J, Schindler Z, Sprengel L, Hillenmeyer K, Kohnle U, Seifert T, Spiecker H (2025) Quantifying and mapping the ready-to-use veneer volume of European beech trees based on terrestrial laser scanning data. European Journal of Forest Research. https://doi.org/10.1007/s10342-025-01796-z

BibTex:
```
@article{larysch_2025,
	title = {Quantifying and mapping the ready-to-use veneer volume of European beech trees based on terrestrial laser scanning data},
	issn = {1612-4669},
	url = {https://link.springer.com/epdf/10.1007/s10342-025-01796-z},
	doi = {10.1007/s10342-025-01796-z},
	abstract = {Using 3D point clouds obtained with terrestrial laser scanning ({TLS}), we automatically and non-destructively quantified and mapped the estimated veneer wood volume of standing trees in differently structured beech stands. To mitigate climate change, we need to utilise wood for long-term carbon storage in products like construction wood and for substituting building materials based on fossil fuels. As the supply of wood from Norway spruce decreases, alternative species like beech must be considered for construction purposes. We present an approach to quantify and map the volume available for veneer production in beech forests. Our method is based on point clouds derived from {TLS}. We studied three forest plots, each with two different treatments (moderate vs. heavy thinning), resulting in varying stand basal areas ranging from 25 m2 to 36 m2 per hectare. We fitted different configurations of veneer rolls into point clouds of tree stems, choosing the configuration that yielded the highest volume of veneer wood. Our automatic optimisation algorithm ensured no misplaced veneer rolls. At the tree level, veneer wood volume was higher in intensely thinned stands. At the stand level, overall veneer volume was higher in moderately thinned stands, whereas the overall veneer share was higher in the heavily thinned stands. The veneer volume of a tree depended on diameter at breast height, crown base height, taper and curvature depth. Our approach detects all trees in a forest potentially ready for veneer production and shows the direct volumetric outcome under bark. This enables the planning of tree selection for harvest based on adaptable requirements for the veneer production.},
	journaltitle = {European Journal of Forest Research},
	author = {Larysch, Elena and Frey, Julian and Schindler, Zoe and Sprengel, Lars and Hillenmeyer, Katharina and Kohnle, Ulrich and Seifert, Thomas and Spiecker, Heinrich},
	urldate = {2025-06-25},
	date = {2025},
	langid = {english},
	file = {Full Text PDF:O\:\\Research\\Projects\\Confobi_IWW\\Literatur\\lit_database\\storage\\R7Q8BFU5\\Larysch et al. - 2025 - Quantifying and mapping the ready-to-use veneer volume of European beech trees based on terrestrial.pdf:application/pdf},
}
```
