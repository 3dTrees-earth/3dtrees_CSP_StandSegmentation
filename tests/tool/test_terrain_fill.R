#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = FALSE)
test_file <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[[1L]])
repository_root <- normalizePath(file.path(dirname(test_file), "..", ".."))
source(file.path(repository_root, "workflow", "cli.R"))
source(file.path(repository_root, "workflow", "tool.R"))
source(file.path(repository_root, "workflow", "outputs.R"))

expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); "" }, error = conditionMessage)
  stopifnot(grepl(pattern, error, fixed = TRUE))
}

dtm <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 10, crs = "")
grid <- rep(NA_real_, 25)
grid[c(8L, 12L, 14L)] <- c(100, 2.123456789, 5)
terra::values(dtm) <- grid
fill <- new_terrain_fill(dtm, list())
stopifnot(
  identical(terrain_under_cells(dtm, grid, c(12, 13, 13, 14), fill),
            c(grid[[12L]], rep(mean(grid[c(12L, 14L)]), 2L), grid[[14L]])),
  identical(terrain_fill_distances(fill), c(`13` = 1)),
  identical(merge_terrain_fills(list(a = c(`13` = 1), b = c(`13` = 1, `20` = 2))), c(`13` = 1, `20` = 2))
)
# Neither a newly estimated cell nor a finite centre just beyond the metre cap
# can support a second estimate. This also exercises the non-square grid.
expect_error(terrain_under_cells(dtm, grid, 18, new_terrain_fill(dtm, list(dtm_fill_max_distance = 2))),
             "DTM NoData hole at cell 18")
stopifnot(terrain_under_cells(dtm, grid, 18, fill) == mean(grid[c(12L, 14L)]),
          terrain_fill_distances(fill)[["18"]] == sqrt(5))

work <- tempfile("inventory-fill-")
dir.create(work)
input <- file.path(repository_root, "tests/tool/fixtures/multi-instance.copc.laz")
cloud <- suppressWarnings(lidR::readLAS(input, select = "*"))
dtm <- terra::rast(xmin = floor(min(cloud$X)) - 7, xmax = ceiling(max(cloud$X)) + 7,
                   ymin = floor(min(cloud$Y)) - 7, ymax = ceiling(max(cloud$Y)) + 7,
                   resolution = 1, crs = sf::st_crs(cloud)$wkt)
cells <- terra::cellFromXY(dtm, cbind(cloud$X, cloud$Y))
# Choose an occupied cell with at least one tree elsewhere. Assign its centre
# a distinct elevation so a failed lookup or unwanted change is observable.
valid_ids <- cloud$PredInstance_FM > 0
cell <- cells[which(valid_ids)[[1L]]]
touched <- unique(cloud$PredInstance_FM[cells == cell & valid_ids])
stopifnot(!1L %in% cells, length(setdiff(unique(cloud$PredInstance_FM[valid_ids]), touched)) > 0L)
grid <- rep(2, terra::ncell(dtm))
grid[[1L]] <- NA_real_
grid[[cell]] <- 2.3
write_fixture <- function(label, elevations) {
  path <- file.path(work, paste0(label, ".tif"))
  terra::values(dtm) <- elevations
  terra::writeRaster(dtm, path, datatype = "FLT8S")
  path
}
original <- write_fixture("original", grid)
grid[[cell]] <- NA_real_
holed <- write_fixture("single-hole", grid)
near <- colSums((t(terra::xyFromCell(dtm, seq_along(grid))) - terra::xyFromCell(dtm, cell)[1L, ])^2)
grid[near <= 25] <- NA_real_
large_hole <- write_fixture("large-hole", grid)

# Crop an AOI containing the point cloud and the hole. Verify the same AOI
# raster that would be produced directly from the unmodified holed input.
xmin <- min(cloud$X) - 1
xmax <- max(cloud$X) + 1
ymin <- min(cloud$Y) - 1
ymax <- max(cloud$Y) + 1
polygon <- rbind(c(xmin, ymin), c(xmax, ymin), c(xmax, ymax), c(xmin, ymax), c(xmin, ymin))
aoi <- file.path(work, "aoi.geojson")
crs <- sf::st_crs(cloud)
sf::st_write(sf::st_sfc(sf::st_polygon(list(polygon)), crs = crs), aoi, quiet = TRUE)
expected_aoi <- file.path(work, "expected-aoi.tif")
write_aoi_dtm(terra::rast(holed), read_aoi(aoi, crs), expected_aoi, local = FALSE)
hash <- function(path) unname(tools::md5sum(path))

run_case <- function(label, mode, path, flags = character(), failure = NULL) {
  output <- file.path(work, label)
  config <- parse_cli_args(c(
    "--input", input, "--output-dir", output, "--instance-attribute", "PredInstance_FM",
    "--dtm", path, "--inventory-workers", "2", "--aoi-json", aoi,
    if (mode == "csp") c("--enable-csp", "true") else c("--inventory-read-mode", mode),
    flags
  ))
  if (!is.null(failure)) {
    expect_error(suppressWarnings(run_tool(config)), failure)
    stopifnot(!dir.exists(output))
    return(invisible(NULL))
  }
  log <- capture.output(run_tool(config), type = "message")
  filled <- grep("^WARNING: DTM height lookup", log, value = TRUE)
  if (identical(path, holed)) {
    metadata <- jsonlite::fromJSON(file.path(output, "run_metadata.json"))
    stopifnot(length(filled) == 1L, grepl("1 filled NoData cells", filled, fixed = TRUE),
              grepl("largest fill distance: 1 m", filled, fixed = TRUE),
              metadata$dtm_filled_cells == 1L, metadata$dtm_max_fill_distance_m == 1,
              metadata$dtm_fill_cap_m == 5, hash(file.path(output, "dtm_aoi.tif")) == hash(expected_aoi))
  } else {
    stopifnot(!length(filled), !file.exists(file.path(output, "run_metadata.json")))
  }
  stopifnot(hash(file.path(output, "dtm_full.tif")) == hash(path))
  output
}

read_rows <- function(path) utils::read.csv(file.path(path, "inventory.csv"))
for (reader in c("copc-owned", "sequential", "spatial", "csp")) {
  baseline <- run_case(paste0(reader, "-original"), reader, original)
  strict <- run_case(paste0(reader, "-strict-original"), reader, original, c("--dtm-fill-radius", "0"))
  tables <- c("inventory.csv", "inventory.parquet", "stand_summary.csv")
  stopifnot(identical(hash(file.path(baseline, tables)), hash(file.path(strict, tables))))
  patched <- run_case(paste0(reader, "-filled"), reader, holed)
  before <- read_rows(baseline)
  after <- read_rows(patched)
  keep <- before$instance_dimension == "PredInstance_FM" & !before$instance_id %in% touched
  unchanged <- after$instance_dimension == "PredInstance_FM" & !after$instance_id %in% touched
  stopifnot(identical(before[keep, ], after[unchanged, ]))
  cat(sprintf("%s: no-hole parity, one used fill, untouched tree rows and published DTMs checked\n", reader))
}
run_case("strict-holed", "copc-owned", holed, c("--dtm-fill-radius", "0"),
         "DTM has missing or non-finite cells under inventoried points")
run_case("too-wide", "sequential", large_hole, failure = "extends more than 5 m")
unlink(work, recursive = TRUE)
cat("Bounded terrain fill, strict refusal and cap checks passed\n")
