#!/usr/bin/env Rscript

# DBH evidence records what the measurement already computed: the inventory
# row is identical with and without it, the reported run's slices match the
# reported DBH, and every run and the crown hull are written once per tree.

arguments <- commandArgs(trailingOnly = FALSE)
test_file <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[[1L]])
root <- normalizePath(file.path(dirname(test_file), "..", ".."), mustWork = TRUE)
source(file.path(root, "workflow", "cli.R"))
source(file.path(root, "workflow", "tool.R"))

set.seed(1)
config <- list(slice_min = 0.3, slice_max = 4, slice_increment = 0.2,
               slice_width = 0.1, max_dbh = 1, geometry_threads = 1L,
               non_tree_ids = 0, random_seed = 42L,
               dbh_band_points = inventory_dbh_band_point_budget, dbh_runs = 5L, geometry_k = 10L)
spec <- list(source = "SAT", instance = "PredInstance_SAT")
terrain <- 2
radius <- 0.2
grid <- expand.grid(angle = seq(0, 2 * pi, length.out = 121)[-121], z = seq(0, 6, by = 0.02))
crown <- 4000L
points <- data.table::data.table(
  X = c(10 + (radius + stats::rnorm(nrow(grid), sd = 0.002)) * cos(grid$angle), 10 + stats::runif(crown, -2, 2)),
  Y = c(10 + (radius + stats::rnorm(nrow(grid), sd = 0.002)) * sin(grid$angle), 10 + stats::runif(crown, -2, 2)),
  Z = c(terrain + grid$z, terrain + stats::runif(crown, 6, 12)),
  PredInstance_SAT = 1
)
dtm <- terra::rast(xmin = 0, xmax = 20, ymin = 0, ymax = 20, resolution = 1, crs = "EPSG:32632")
terra::values(dtm) <- terrain
normalized <- normalize_inventory_points(data.table::copy(points), dtm)
crs <- sf::st_crs(32632)

set.seed(7)
plain <- measure_inventory_points(data.table::copy(normalized), spec, config, NULL, crs)

stage <- tempfile("dbh-evidence-")
dir.create(stage)
with_evidence <- modifyList(config, list(dbh_evidence_dir = file.path(stage, ".dbh-evidence")))
dir.create(with_evidence$dbh_evidence_dir)
set.seed(7)
recorded <- measure_inventory_points(data.table::copy(normalized), spec, with_evidence, NULL, crs)
stopifnot(identical(plain, recorded))

artifacts <- write_dbh_evidence(stage, list(SAT = recorded), with_evidence, crs)
stopifnot(length(artifacts) == 6L, !dir.exists(with_evidence$dbh_evidence_dir))
slices <- data.table::fread(file.path(stage, "dbh", "dbh_slices.csv"))
runs <- data.table::fread(file.path(stage, "dbh", "dbh_runs.csv"))
hulls <- data.table::fread(file.path(stage, "dbh", "crown_hulls.csv"))

stopifnot(nrow(runs) == recorded$dbh_runs, sum(runs$reported) == 1L,
          isTRUE(all.equal(runs[reported == TRUE, dbh_m], recorded$dbh_m)))
stopifnot(nrow(slices) == 19L, all(slices$run == runs[reported == TRUE, run]),
          all(slices$status %in% dbh_slice_statuses),
          all(abs(slices[status == "fitted", radius_m] - radius) < 0.02),
          all(abs(slices$ground_z - terrain) < 1e-9))
# The reported DBH is the radius spline at 1.3 m (+ half an increment), near the fitted slices.
stopifnot(abs(recorded$dbh_m / 2 - stats::median(slices$spline_r, na.rm = TRUE)) < 0.01)
stopifnot(nrow(hulls) == 1L, grepl("^POLYGON \\(\\(", hulls$wkt))
hull <- sf::st_as_sfc(hulls$wkt, crs = crs)
stopifnot(abs(as.numeric(sf::st_area(hull)) - recorded$convex_hull_area_m2) < 1e-3)
geo <- arrow::read_parquet(file.path(stage, "dbh", "crown_hulls.parquet"), as_data_frame = FALSE)
stopifnot(!is.null(geo$metadata$geo), "geometry" %in% names(geo))
cat("dbh evidence: ok\n")
