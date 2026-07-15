## CspStandSegmentation is an R-package for the segmentation of single trees from forest point clouds scanned with terrestrial, mobile or unmanned LiDAR systems <img src="https://github.com/JulFrey/CspStandSegmentation/blob/main/inst/figures/csp_logo.png" align="right" width = 300/>

Authors: Julian Frey and Zoe Schindler, University of Freiburg, Chair of Forest Growth and Dendroecology


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17294732.svg)](https://doi.org/10.5281/zenodo.17294732)  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## 3Dtrees command-line workflow

This fork retains the upstream segmentation implementation and adds a headless
CLI and container for the 3Dtrees Galaxy tool. It always creates a DTM and can
inventory any number of existing instance dimensions. CSP segmentation is
optional; when enabled, the original points and dimensions are preserved and
the output gains exactly one `PredInstance_CSP` dimension.

Inventory-only runs avoid loading dimensions that cannot affect the DTM or the
requested inventories. They retain XYZ, Classification, requested
instance/species fields, and optional ForestMamba score/semantic fields. CSP
runs continue to load every dimension because the emitted point cloud must
preserve them. LASlib can select only the first nine extra-byte records by
position; when a requested field occurs later, the reader safely falls back to
all extra bytes, then immediately projects the in-memory cloud back to the
required fields. The container uses the same pinned `rlas` patch as 3Dtrees
standardization, so that fallback correctly loads every declared extra byte
rather than stopping after nine.

```bash
Rscript exec/run.R \
  --input input.laz \
  --output-dir results \
  --segmentation-spec PredInstance_SAT,species_id_SAT,species_prob_SAT \
  --segmentation-spec PredInstance_FM,species_id_FM,species_prob_FM \
  --enable-csp false \
  --dtm-resolution 0.2 \
  --random-seed 42
```

`--segmentation-spec INSTANCE[,SPECIES,SPECIES_PROB]` is repeatable. Species
dimensions are optional per segmentation but must be supplied as a pair. When
none are supplied, the combined inventory omits species columns and no species
composition file is created. Run `Rscript exec/run.R --help` for all controls.

Common controls are the input, repeatable segmentation specs, non-tree IDs,
optional CSP and seed source, optional native-CRS AOI GeoJSON, DTM resolution
(default 0.2 m), and random seed. Fine-tuning controls retain upstream defaults,
including a 0.3 m CSP voxel and one routing worker. CSP geometry features are
computed only when a non-zero geometry weight requires them.

Outputs include:

- `dtm_full.tif` and optional `dtm_aoi.tif`;
- one TSV per instance dimension plus `inventory_combined.tsv`;
- `stand_summary.tsv` and optional `species_composition.tsv`;
- `effective_seeds.tsv` and `segmented_csp.laz` only when CSP is enabled;
- `run_metadata.json` and `resource_summary.json`.

Inventory rows include point count, position, height, DBH, crown convex-hull
area, and an explicit measurement-quality field. ForestMamba inventories also
include median `PredScore_FM`, a mixed-score flag, and wood/leaf point counts and
shares when those source dimensions exist. Processing is fail-atomic: the final
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
spatial chunks. CSP continues to use the full-cloud path to preserve upstream
global voxel routing and the optional point-cloud output.

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

Build and run the pinned container with:

```bash
docker build -t 3dtrees-csp .
docker run --rm -v "$PWD:/work" -w /work 3dtrees-csp \
  Rscript /opt/CspStandSegmentation/exec/run.R --help
```




## Installation
We are on CRAN 🎉🥳 you can install the packe in R the simple way:
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
