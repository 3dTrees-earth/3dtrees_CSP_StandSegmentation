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
species/probability fields are read only when listed in `--segmentation-spec`.
Use `--inventory-attributes PredScore_FM,PredSemantic_FM` to opt into FM score
and semantic summaries. This flag accepts a comma-separated list of numeric,
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
  --segmentation-spec PredInstance_FM --inventory-read-mode copc-owned \
  --dtm terrain.tif
```

Add `--inventory-attributes PredScore_FM,PredSemantic_FM` for those optional
measurements. Include the selected column count (XYZ + instance + optional
fields) in resource estimation through `--inventory-columns` when it exceeds
the estimator's conservative eight-column default.

```bash
Rscript exec/run.R \
  --input merged.copc.laz \
  --output-dir results \
  --segmentation-spec PredInstance_SAT,species_id_SAT,species_prob_SAT \
  --segmentation-spec PredInstance_FM,species_id_FM,species_prob_FM \
  --enable-csp false \
  --dtm existing_dtm.tif \
  --output-formats csv,tsv,json,gpkg \
  --dtm-resolution 0.2 \
  --random-seed 42
```

`--segmentation-spec INSTANCE[,SPECIES,SPECIES_PROB]` is repeatable. Species
dimensions are optional per segmentation but must be supplied as a pair. When
none are supplied, the combined inventory omits species columns and the species
table is empty. Run `Rscript exec/run.R --help` for all controls.

Common controls are the merged COPC input, repeatable segmentation specs, non-tree IDs,
optional CSP and seed source, optional native-CRS AOI GeoJSON, an optional DTM,
selected output formats, DTM resolution (default 0.2 m when generated), and random seed.
The source point cloud must declare a projected metre CRS. The CLI can recover
the CRS from a valid extended WKT record when a COPC file has a truncated
ordinary WKT record, and rejects conflicting valid records. A supplied DTM must
use the same CRS or an equivalent horizontal CRS and cover the full cloud
extent. When the GeoTIFF omits a vertical datum, its heights must already use
the point cloud's vertical datum. The supplied DTM is copied to
`dtm_full.tif` in the output. With no `--dtm`, the tool generates that file.
Missing DTM cells beneath inventoried points stop the run instead of silently
dismissing those points.
Fine-tuning controls retain upstream defaults,
including a 0.3 m CSP voxel and one routing worker. CSP geometry features are
computed only when a non-zero geometry weight requires them.

All four table formats are written by default; `--output-formats` can select any
combination. Outputs include:

- `dtm_full.tif` and optional `dtm_aoi.tif`;
- one CSV, TSV, and/or JSON per instance dimension plus combined Inventory,
  Species, Stand Summary, Validation, and Issues tables in each selected format;
- `inventory.gpkg` when selected, with one 3D point layer per segmentation and
  those five attribute tables;
- `manifest.json` with version, artifact keys, fields and units, source/WGS84
  coordinate metadata, segmentation dimensions, and DTM metadata;
- `effective_seeds.tsv` and `segmented_csp.laz` only when CSP is enabled;
- `run_metadata.json` and `resource_summary.json`.

JSON table exports retain double-precision coordinates so the selected formats
describe the same tree positions.

Every run also writes, as CSV and Parquet, the evidence each DBH and crown area
was computed from; `--dbh-evidence false` switches it off. The fit and every
other output stay byte-identical either way; nothing here feeds back into a
measurement.

| File | One row per | Columns |
| --- | --- | --- |
| `dbh/dbh_slices.*` | tree × stem slice of the reported (median) run | `instance_dimension`, `instance_id`, `run`, `slice_id`, `z_lo_m`, `z_hi_m` (above the DTM), `ground_z` (m, source CRS), `n_points`, `n_fit_points` (points handed to RANSAC), `circle_x`, `circle_y`, `radius_m` (raw fit), `status`, `spline_x`, `spline_y`, `spline_r` (the smoothing splines at the slice's lower edge, the abscissa the 1.3 m prediction uses; empty when no spline was fitted) |
| `dbh/dbh_runs.*` | tree × DBH run | `run`, `seed`, `dbh_m`, `quality_code`, `reported` |
| `dbh/crown_hulls.*` | tree | `vertex_count` and the convex hull polygon (`wkt` in the CSV, GeoParquet 1.1 WKB geometry in the source CRS) of the points above 0.5 m, or of all points when fewer than three are |

`status` names the upstream rule that set a slice aside: `fitted`, `outlier_2sd`
(radius more than two standard deviations from the median), `growth_1_3x`
(radius above 1.3 times the smallest radius below it), `radius_below_0_04m`,
`too_few_points` (fewer than three points) or `no_fit` (RANSAC found no circle).
Upstream applies the growth and minimum-radius rules from the second slice on, so the lowest slice is never set aside by them.

Inventory rows include source `x`, `y`, `z` and WGS84 `longitude`, `latitude`.
DBH measurement status uses `dbh_quality_code` and `dbh_quality_label`.
If species predictions are present, `--species-lookup lookup.tsv` can add names
and taxonomy/model versions. The lookup needs `segmentation_source`,
`species_id`, `scientific_name`, `taxonomy_version`, and `model_version` columns;
`common_name` is optional. Missing lookup metadata is reported in the Issues
table without inventing species names.

Inventory rows include point count, position, height, DBH, crown convex-hull
area, and an explicit measurement-quality field. ForestMamba inventories also
include median `PredScore_FM`, a mixed-score flag, and wood/leaf point counts and
shares when those source dimensions are explicitly selected. Processing is fail-atomic: the final
output directory is published only after all requested products succeed.

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

The second pass reads only selected inventory fields into disk-backed hash
partitions by instance ID. Each tree remains complete even when its points span
spatial chunks. Inventory samples the DTM in 100,000-point batches, then spools
normalized points by instance and runs the existing measurements in groups of
at most two million points. The scanner records the largest instance and the
tool rejects any instance above two million points before inventory workers
start; the spatial fallback enforces the same limit while spooling. Group
geometry is also capped at 500,000 DBH-height points. Group boundaries can
change the neighbor set used for DBH geometry, so DBH values
may differ from a whole-partition run when a partition spans multiple groups.
The scan rejects more than one million distinct instance IDs across layouts
instead of letting its exact preflight-count map grow without a limit.
CSP continues to use the full-cloud path to preserve upstream global voxel
routing and the optional point-cloud output.

For a single merged COPC input, `--inventory-read-mode sequential` uses a
bounded Python scan of the compressed point stream. It writes spatial buckets
and joins them in the same tile order as the R reader, preserving the serial
inventory output while avoiding repeated COPC spatial reads. The default
`spatial` mode remains available for LAS/LAZ directories and comparison runs.
The sequential scan needs scratch for both bucket and final partition files.
It maintains exact convex-hull vertices without an XY spool. LAS 1.4 point
formats 6–10 use selective LAZ decompression for requested standard groups;
extra-byte dimensions share one compressed group and cannot be selected
individually. Older formats use the full decompression path.

### Complete-instance COPC windows without point partitions

Use `--inventory-read-mode copc-owned` to read COPC nodes in bounded batches
and send selected points directly to R through a pipe. The first pass records
all instance counts, global XY centroids and extents. Each ID belongs to one
non-overlapping core (default 50 m, origin at the input minimum). Its initial
20 m read buffer expands to include every owned ID's full extent. Both Python
and R check observed counts against the global index before finalizing a
window. All COPC hierarchy levels are included, without LOD reduction. Input
changes, overlapping nodes, oversized decode nodes and incomplete reads fail
without publishing the staged inventory.

```bash
Rscript exec/run.R --input merged.copc.laz --output-dir inventory \
  --segmentation-spec PredInstance_FM,species_id_FM,species_prob_FM \
  --dtm terrain.tif --inventory-read-mode copc-owned \
  --inventory-core-size 50 --inventory-buffer 20 \
  --inventory-batch-points 2000000 --inventory-retained-points 40000000 \
  --inventory-workers 16 --inventory-decode-threads 4
```

Omit `--dtm` to use the existing automatic terrain generation. DTM generation
precedes the index. The initial 20 m buffer is not a maximum tree radius.
Owned instances are selected using global centroids, never centroids of clipped
tiles. Neighbouring points are discarded as each COPC batch is decoded; an
entire buffered cloud is never materialized. A node batch has at most the
configured point count and 256 MiB of native point records. A single node above
that budget is rejected before decompression.

A core may be split into ID batches using the retention budget: small instances
cost their full count; instances above two million points cost up to 100,000
sample points. Complete small instances use the existing forest inventory.
Large instances retain exact all-point counts, extrema, hull and class totals;
medians use the existing bounded stride sample, and DBH remains unavailable.
The 500,000-point DBH geometry guard also remains active. Every ID is measured
once. Within a window, complete IDs are assigned to workers by balancing their
point counts. One Python COPC reader routes each bounded batch to those workers
through backpressured kernel pipes. Increasing worker count does not add COPC
reads. Windows run in sequence, with the configured worker budget available to
each; even a single dense core can use several R workers. Named FIFO entries are
removed after each window and never store point payload on disk.

Each instance keeps its source point order and receives an `instance_identity_v1`
seed derived from the base seed, segmentation dimension and exact instance ID.
Worker count, window order, core size and retention budget do not change that
seed. This replaces the earlier window-dependent seed policy; stochastic DBH
values can change once when upgrading. The partition reader keeps its existing
seed policy, so its stochastic DBH/positions need not be byte identical. Sampled medians can differ from the partition reader because COPC
traversal order differs. This mode applies to existing segmentation inventory only.

This path writes no all-point buckets, partitions or per-instance point files.
It exports `instance_index.json` and per-window completeness audits in
`run_metadata.json`, alongside the existing table formats and DTM. The index
can refine resource estimates. Exact hulls and index entries consume additional
memory, so enforce a container RAM cap even with bounded decode/retention
budgets. `inventory-workers + inventory-decode-threads` may not exceed 20.
For example, 16 R workers share four decoder threads. The point retention limit
applies to the entire window, divided among workers, with additional per-worker
runtime/geometry overhead. Each worker checks its complete counts; a failed
consumer or reader aborts the window and prevents publication. Per-window
metadata records actual worker count, disjoint ID assignments and point counts.

For `copc-owned` parallel inventory, `--inventory-prefetch-mb` bounds the shared
in-memory delivery queue (default 256 MiB). The reader can decode subsequent
spatial batches while R workers consume earlier batches. The resource estimator
includes this queue and can increase it using the complete instance index, up to
8 GiB and one eighth of the requested RAM cap. It never writes point payloads to
disk. This overlaps work but cannot guarantee full CPU utilization for spatially
clustered or uneven instances.


The cost is one full index pass plus spatial window rereads. Overlapping windows
can decode points repeatedly; only owned points are delivered to inventory.
Metadata records decoded points, including repeats. Benchmark this tradeoff
against the partition reader for the input's density and ID distribution.

```bash
python3 workflow/resource_request.py merged.copc.laz --dtm-supplied \
  --inventory-read-mode copc-owned
python3 workflow/resource_request.py merged.copc.laz --dtm-supplied \
  --inventory-read-mode copc-owned --owned-index inventory/instance_index.json
```

The new resource estimates use buffer sizes and optional instance counts;
they are provisional until calibrated against full-run cgroup measurements.

A local shared-reader benchmark on a 45.295 GB COPC (4.257 billion points,
154 complete FM instances) used 16 R workers, four decoder threads and an
8 GiB queue. Measurement took 1,845.8 s versus 4,367.7 s with one worker (2.37x
faster); cgroup peak was 21.48 GiB under a 100 GiB cap, average CPU use 4.04
cores, and no point payload files were written. Adding the separately measured
1,590.9 s index to the resumed measurement/export process gives 57.33 minutes,
versus 101.42 minutes for the earlier serial CLI. This is a sum of measured
phases with a reused verified index, not a fresh end-to-end or cold-cache run.
The indexed resource estimate requests 32 GB RAM with headroom. Twenty CPUs
is the tested allowance, not a proven optimum for this large file.

All 154 IDs, 2,546,177,540 assigned points, 23 artifacts, and DTM checks passed.
Non-DBH inventory metrics matched the earlier serial result. DBH changed with
the per-instance seed policy: four finite fits versus five, one newly unavailable
fit, and up to 0.217887 m difference among mutually finite values. This is not a
claim of DBH equivalence. New serial/parallel fixture and 100 MB runs are byte
identical across worker counts using the new seed policy.

On the same 100 MB COPC with a verified reused index, single warm-cache runs
at 2/4/8/12/20 CPU caps took 74.30/47.42/34.58/26.59/24.09 seconds, respectively.
The corresponding R-worker/decoder pairs were 1/1, 2/2, 4/4, 8/4 and 16/4;
peak RAM was 1.23/1.92/2.64/4.30/7.11 GiB. All outputs validated and inventories
were byte-identical. Twelve CPUs was a useful compromise for this small input:
40% fewer requested CPUs than 20, for 10.4% more elapsed time. This does not
establish the large-input CPU optimum.

Request resources from the COPC header before submission:

```bash
python3 workflow/resource_request.py merged.copc.laz --dtm-supplied
python3 workflow/resource_request.py merged.copc.laz --dtm-supplied --profile throughput
python3 workflow/resource_request.py merged.copc.laz --dtm-supplied --profile throughput \
  --partition-dir /path/to/completed/partitions
```

The JSON returns CPUs, RAM, scratch, spatial chunk width, partition count,
inventory read mode, and worker count within the requested 20-CPU/100-GB
limits. The returned `cli_args` array is the executable resource configuration:
append it to `Rscript exec/run.R` together with your input/output paths,
segmentation specs, and supplied DTM path (if applicable). Pass it as an argument
array, not as a shell string. Enforce `cpus` and `memory_gb` in the container or
scheduler. Other fields are reporting/estimation data: for example, zero
partitions in COPC mode is a statistic, not a valid `--inventory-partitions`
argument. The resource-to-CLI integration test parses recommendations for both
readers, supplied/generated DTM, and multiple CPU caps using the real R parser.
The default `stable` profile requests one inventory worker. For the partition
reader this preserves the serial random stream; `copc-owned` uses the same
per-instance seed policy at every worker count. The opt-in `throughput` profile requests more workers and
RAM when the COPC point count warrants it. The estimate reads only the header;
it does not run the inventory or guarantee that a highly skewed tree fits.
Header-wide density cannot reveal unusually dense local tiles or one giant
instance; apply the returned hard memory and scratch limits and inspect the
run metadata. For inputs above the completed end-to-end size, use a bounded
scan/partition preflight and pass its completed partition directory back with
`--partition-dir`. New scans save `partition_summary.json` in that directory;
the estimator uses its largest-instance count to enforce the two-million-point
limit before recommending workers. An older external audit can be supplied with
`--instance-preflight AUDIT.json`. Without an instance audit, the estimator
retains the older partition-size risk screen. Completed partitions also let it
replace the header-only scratch estimate with actual partition bytes and the
largest spools that can coexist across active workers. The default single
inventory worker retains the upstream random stream. `--inventory-workers`
can parallelize independent partitions, but
values above one use per-partition seeds; derived DBH measurements can differ
from a single-worker run even when tree IDs and point counts agree.

On a segmented dataset 2924 COPC, the 1.015 GB inventory completed in
327.7 seconds at 4 CPUs/16 GiB, versus 209.6 seconds at 8 CPUs/32 GiB and
185.6 seconds at 20 CPUs/100 GiB. The complete 11.593 GB / 925.7-million-point
COPC finished in 7,990.6 seconds with 20 CPUs/100 GiB and 512 partitions;
its cgroup peak was 88.95 GiB with no OOM. These are observed requests, not
minimum safe resources for other clouds. Larger inputs require an exact
oversized-instance audit before the resource function can request the bounded
streaming path.
An existing 45.3 GB COPC completed the scan/partition stage in 4,279 seconds
at 1 CPU/16 GiB. Its largest 6.54-GiB partition contains one 109.8-million-
point instance. The scanner now records every ID above two million points;
those IDs use bounded streaming summaries rather than a whole-instance LAS.
Point counts, extrema, species frequencies, FM semantic counts, and convex
hulls remain exact. X/Y, species probability, and FM score medians use a
deterministic sample of at most 100,000 points. DBH is marked unavailable for
these IDs, and the exported issues table identifies each fallback. Standard
instances whose DBH slice exceeds the geometry cap use the same fallback.
The stand summary reports DBH coverage and leaves total basal area empty when
any tree lacks DBH.

The complete 45.295 GB / 4.257-billion-point inventory was then run from its
audited, completed partitions with a generated 1 m DTM from the same COPC.
All 154 FM IDs and 2,546,177,540 assigned points were exported as CSV, TSV,
JSON, and GeoPackage with a validated manifest and DTM. Under a hard 20 CPU /
24 GB cap, largest-first dynamic partition scheduling took 372.4 seconds for
inventory and peaked at 13.44 GiB cgroup memory. The same inputs took 503.8
seconds with fixed scheduling and 438.4 seconds with filename-order dynamic
scheduling. The inventory CSV was byte-identical across all three runs.
Twenty CPUs and 24 GB is the measured same-input request; before a full-run
calibration, the resource estimator requests 20 CPUs and 48 GB from the
oversized-instance audit. A fresh ordered scan needs about 419 GB scratch by
the estimator; the completed-partition recovery avoided that rescan. This run
did not build an image or exercise the future Galaxy wrapper.

COPC DTM candidate generation now reads sequentially because random LAS seeks
on this valid 4.257-billion-point COPC returned zero-filled records. The 45 GB
run reused a previously generated, full-extent DTM after checking its CRS,
coverage, finite cells, and source provenance. Automatic DTM generation in a
new image still needs its own full CLI run.

A 94.7 MB format-7 COPC scan on 25 September 2026 took 5.44 s
with selective decompression versus 6.99 s with full decompression; all 60
instance partitions were byte-identical. On those same partitions and DTM, the
bounded inventory CSV matched the previous whole-partition inventory byte for
byte for 172 trees. Its inventory phase took 85.9 s versus 71.0 s, while cgroup
peak memory fell from 0.89 to 0.72 GiB at one CPU. A 1.015 GB selective scan
finished in 55.9 s at 0.44 GiB process peak RSS; the earlier scan took 70.6 s.
A 5.713 GB selective scan finished in 326.1 s at 0.45 GiB process peak RSS;
the earlier scan took 403.5 s. The scanner emitted no all-point XY spool in
these runs. The complete 11.593 GB COPC scan finished in 665.9 s at 0.50 GiB
process peak RSS versus 817.4 s previously. Its largest 2.24-million-point
partition produced a byte-identical 72-tree inventory with the bounded path;
cgroup peak fell from 2.84 to 1.44 GiB at one CPU, and wall time rose from
37.3 to 40.6 s. The new bounded code has not repeated the complete 11.593 GB
inventory/export phase in a newly built image.

An offline 15,915,252-point merged COPC run exercised two requested instance
dimensions under a 10-CPU/50-GB Docker limit. With eight 20 m spatial chunks,
the supplied-DTM run completed in 283.4 seconds and the generated-DTM run in
461.6 seconds. Both emitted all four formats and exactly the same 96 FM and 98
FM2 instance IDs and point counts: 11,049,398 and 11,048,497 assigned points,
respectively. The generated 0.2 m DTM had 49,289 finite cells. Process peak RSS
was about 1.45 GiB in both runs; the generated-DTM run recorded a conservative
9.16 GiB sum of worker peak RSS values. The per-instance inventory pass remains
serial and dominated the measured runtime.

On dataset 2056 (1,141,911,324 points, 3.7 GB compact LAZ), the memory-focused
streaming candidate stage completed in 56.8 seconds with 10 workers, averaged
9.51 CPU cores, and retained 1,267,661 candidates. Parent peak RSS was 0.78 GiB
and the conservative sum of all worker peaks was 6.72 GiB under a 50 GB Docker
limit. The 300 m + 5 m CSF/TIN stage completed in 76.4 seconds; its conservative
aggregate worker peak was 2.98 GiB. An earlier exact-coordinate prototype
scaled from 543.4 seconds with one worker to 306.6 seconds with two, 206.6
seconds with four, and 127.9 seconds with ten. Replacing per-worker XY grids
with deterministic cell centres produced the final 56.8-second result and cut
the conservative 10-worker peak sum from 16.6 GiB to 6.72 GiB.

On the 5.9-million-point GFZ reference, reducing to deterministic 0.1 m cell
centres before CSF/TIN changed the DTM relative to direct full-cloud CSF/TIN by
about 9.1 cm RMSE (4.0 cm median absolute difference). Automatic mode therefore
keeps the direct spatial method for clouds below the 50-million-point threshold.

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

The container test exercises those three instance dimensions, every output
format alone and together, generated and
supplied DTMs, and both inventory-only and CSP-enabled execution:

```bash
docker run --rm -v "$PWD:/work" -w /work 3dtrees-csp \
  Rscript tests/tool/test_end_to_end.R
```




### Inventory implementation boundaries

- `workflow/inventory_policy.json` is the canonical source of measurement and
  retention limits for the R consumer, Python planners and resource estimator.
- `workflow/instance_inventory.R` owns normalization, ordinary measurements,
  geometry guards, large-instance summaries, and stable instance seeds.
- `workflow/partition_inventory.R` and `workflow/owned_inventory.R` own their
  respective point transports and completeness checks.
- `workflow/inventory_readers.R` prepares the selected reader and exposes input
  statistics, the stand hull and a measurement operation. It owns mode-specific
  metadata and artifacts; the common runner never examines partition paths.
- `workflow/inventory_run.R` handles shared DTM preparation, exports and run
  metadata. Legacy spatial grids are allocated only by readers or DTM methods
  that actually need them. A supplied-DTM COPC run does not allocate that grid.

The September 2026 large-file benchmarks above describe the earlier image and
seed policy. They are retained as historical measurements, not as a rerun of
this structural refactor.

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
