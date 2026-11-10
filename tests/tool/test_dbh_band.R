#!/usr/bin/env Rscript

# Every tree gets a DBH from all its stem-band points: small trees as before
# and streamed large trees from the band retained during streaming.

arguments <- commandArgs(trailingOnly = FALSE)
test_file <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[[1L]])
root <- normalizePath(file.path(dirname(test_file), "..", ".."), mustWork = TRUE)
source(file.path(root, "workflow", "cli.R"))
source(file.path(root, "workflow", "tool.R"))

set.seed(1)
config <- list(slice_min = 0.3, slice_max = 4, slice_increment = 0.2,
               slice_width = 0.1, max_dbh = 1, geometry_threads = 1L,
               non_tree_ids = 0, random_seed = 42L,
               dbh_runs = 5L, geometry_k = 10L, metrics_grid = 0)
spec <- list(source = "SAT", instance = "PredInstance_SAT")
terrain <- 2
radius <- 0.2

# A vertical stem of 40 cm diameter from 0 to 6 m, 1.5 degree / 1 cm
# sampling with 2 mm noise, plus a crown of loose points above 6 m.
stem <- function(id, cx, cy) {
  grid <- expand.grid(angle = seq(0, 2 * pi, length.out = 241)[-241], z = seq(0, 6, by = 0.01))
  noise <- stats::rnorm(nrow(grid), sd = 0.002)
  crown <- 4000L
  data.table::data.table(
    X = c(cx + (radius + noise) * cos(grid$angle), cx + stats::runif(crown, -2, 2)),
    Y = c(cy + (radius + noise) * sin(grid$angle), cy + stats::runif(crown, -2, 2)),
    Z = c(terrain + grid$z, terrain + stats::runif(crown, 6, 12)),
    PredInstance_SAT = id
  )
}
points <- stem(1, 10, 10)
dtm <- terra::rast(xmin = 0, xmax = 20, ymin = 0, ymax = 20, resolution = 1, crs = "EPSG:32632")
terra::values(dtm) <- terrain
normalized <- normalize_inventory_points(data.table::copy(points), dtm)
band_total <- sum(in_dbh_band(normalized, config))

dbh_close <- function(value) is.finite(value) && abs(value - 2 * radius) < 0.02

crs <- sf::st_crs(32632)
measure <- function(points, config) inventory_from_points(points, spec, config, NULL, crs)

# 1. One batch.
full <- measure(data.table::copy(normalized), config)
stopifnot(nrow(full) == 1L, dbh_close(full$dbh_m),
          full$dbh_band_points == band_total, full$point_count == nrow(points),
          full$dbh_runs == 5L, full$dbh_runs_fitted >= 3L, is.finite(full$dbh_spread_m))

# Same tree in another point order (another file writer, worker or window)
# gives the identical measurement: points are sorted, seeds follow the tree.
set.seed(123)
shuffled <- measure(data.table::copy(normalized)[sample(.N)], config)
stopifnot(identical(shuffled$dbh_m, full$dbh_m), identical(shuffled$x, full$x),
          identical(shuffled$y, full$y), identical(shuffled$dbh_quality_code, full$dbh_quality_code))

# 2. Many batches give exactly the same row as one.
stream <- function(points, config) {
  tree <- new_tree()
  for (start in seq.int(1L, nrow(points), by = 50000L)) {
    tree <- add_tree_points(tree, points[start:min(start + 49999L, nrow(points))], spec, config)
  }
  finish_tree(tree, 1, spec, config, NULL, crs)
}
same_tree <- function(a, b) {
  columns <- c("x", "y", "z", "dbh_m", "height_m", "convex_hull_area_m2", "point_count",
               "metrics_point_count", "dbh_quality_code", "dbh_runs_fitted", "dbh_band_points")
  isTRUE(all.equal(as.list(a[, ..columns]), as.list(b[, ..columns]), check.attributes = FALSE))
}
stopifnot(same_tree(stream(normalized, config), full))

# 3. --metrics-grid 0.01: geometry on one point per 1 cm cell, the count on
# every point; the same result in any point order and in batches.
grid_config <- modifyList(config, list(metrics_grid = 0.01))
grid_full <- measure(data.table::copy(normalized), grid_config)
cells <- nrow(unique(grid_keys(normalized, 0.01)))
stopifnot(dbh_close(grid_full$dbh_m), grid_full$point_count == nrow(points),
          grid_full$metrics_point_count == cells, cells < nrow(points))
set.seed(7)
stopifnot(same_tree(measure(data.table::copy(normalized)[sample(.N)], grid_config), grid_full),
          same_tree(stream(normalized, grid_config), grid_full))

# Exact medians from value counts: even and odd counts, across batches, and
# species probabilities per species with missing species ignored.
counts <- add_value_counts(add_value_counts(NULL, c(1, 2, 2)), c(NaN, 4))
stopifnot(median_of_counts(counts) == 2, median_of_counts(add_value_counts(counts, 5)) == 2,
          median_of_counts(add_value_counts(counts, 7)) == stats::median(c(1, 2, 2, 4, 7)),
          is.na(median_of_counts(NULL)))
species_spec <- list(source = "SAT", instance = "PredInstance_SAT", species = "sp", species_prob = "prob")
labelled <- data.table::copy(normalized)[1:1000]
labelled[, `:=`(sp = rep(c(14, NA, 3, 14), 250L), prob = rep(c(0.2, 0.9, 0.5, 0.4), 250L))]
voted <- inventory_from_points(labelled, species_spec, config, NULL, crs)
stopifnot(voted$species_id == 14, abs(voted$species_prob - 0.3) < 1e-12,
          abs(voted$species_conflict_fraction - 1 / 3) < 1e-12)

# 4. A tree with no stem points in the band keeps an explicit unavailable DBH.
empty <- measure(normalized[Z > 6], config)
stopifnot(is.na(empty$dbh_m), empty$dbh_quality_label == "dbh_unavailable", empty$dbh_band_points == 0)

cat("stem-band DBH tests passed\n")
