#!/usr/bin/env Rscript

# Every tree gets a DBH from its stem band: small trees as before, trees
# whose band exceeds the budget from a voxel-thinned band, and streamed
# large trees from the band retained during streaming.

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

# The thinning keeps original points, starts at 1 mm and doubles until it fits.
unthinned <- dbh_band_keep(normalized[in_dbh_band(normalized, config)], band_total)
stopifnot(all(unthinned$keep), unthinned$voxel == 0)
thinned <- dbh_band_keep(normalized[in_dbh_band(normalized, config)], 20000)
stopifnot(sum(thinned$keep) <= 20000, thinned$voxel >= inventory_dbh_band_initial_voxel_m,
          isTRUE(all.equal(log2(thinned$voxel / inventory_dbh_band_initial_voxel_m) %% 1, 0)))
again <- dbh_band_keep(normalized[in_dbh_band(normalized, config)], 20000)
stopifnot(identical(thinned, again))

dbh_close <- function(value) is.finite(value) && abs(value - 2 * radius) < 0.02

# 1. Standard path: band within budget, unchanged measurement.
set.seed(instance_inventory_seed(42L, spec$instance, "1"))
full <- measure_inventory_points(data.table::copy(normalized), spec, config, NULL, sf::st_crs(32632))
stopifnot(nrow(full) == 1L, dbh_close(full$dbh_m), full$dbh_band_voxel_m == 0,
          full$dbh_band_points == band_total, full$point_count == nrow(points))

stopifnot(full$dbh_runs == 5L, full$dbh_runs_fitted >= 3L, is.finite(full$dbh_spread_m))

# Same tree in another point order (another file writer, worker or window)
# gives the identical measurement: points are sorted, seeds follow the tree.
set.seed(123)
shuffled <- measure_inventory_points(data.table::copy(normalized)[sample(.N)], spec, config, NULL, sf::st_crs(32632))
stopifnot(identical(shuffled$dbh_m, full$dbh_m), identical(shuffled$x, full$x),
          identical(shuffled$y, full$y), identical(shuffled$dbh_quality_code, full$dbh_quality_code))

# 2. Band above the budget: thinned band, exact point count, same DBH.
small_budget <- modifyList(config, list(dbh_band_points = 20000))
set.seed(instance_inventory_seed(42L, spec$instance, "1"))
reduced <- measure_inventory_points(data.table::copy(normalized), spec, small_budget, NULL, sf::st_crs(32632))
stopifnot(nrow(reduced) == 1L, dbh_close(reduced$dbh_m), reduced$dbh_band_voxel_m > 0,
          reduced$dbh_band_points <= 20000, reduced$dbh_band_source_points == band_total,
          reduced$point_count == nrow(points),
          abs(reduced$dbh_m - full$dbh_m) < 0.01,
          abs(reduced$height_m - full$height_m) < 0.002)

# 3. Streamed large tree: band retained across batches, DBH measured.
state <- new_large_instance_summary(nrow(normalized))
for (start in seq.int(1L, nrow(normalized), by = 50000L)) {
  end <- min(start + 49999L, nrow(normalized))
  state <- update_large_instance_summary(state, normalized[start:end], spec, small_budget)
}
stopifnot(nrow(state$band) <= 20000, state$band_source_points == band_total)
# Incremental thinning across batches equals thinning the whole band at once.
whole <- normalized[in_dbh_band(normalized, config), .(X, Y, Z, Zref)]
once <- dbh_band_keep(whole, 20000)
stopifnot(state$band_voxel == once$voxel, isTRUE(all.equal(state$band, whole[once$keep], check.attributes = FALSE)))
set.seed(instance_inventory_seed(42L, spec$instance, "1"))
large <- finish_large_instance_summary(state, "1", spec, NULL, sf::st_crs(32632), small_budget)
stopifnot(dbh_close(large$dbh_m), large$dbh_quality_code %in% c(2L, 4L),
          large$point_count == nrow(points), large$inventory_method == "streaming_large_instance",
          abs(large$x - 10) < 0.05, abs(large$y - 10) < 0.05,
          abs(large$height_m - full$height_m) < 0.002,  # LAS objects quantize to 1 mm
          large$dbh_band_source_points == band_total, large$dbh_band_voxel_m > 0,
          large$dbh_runs == 5L, large$dbh_runs_fitted >= 3L)

# 4. A tree with no stem points in the band keeps an explicit unavailable DBH.
crown_only <- normalized[Z > 6]
state <- update_large_instance_summary(new_large_instance_summary(nrow(crown_only)), crown_only, spec, config)
empty <- finish_large_instance_summary(state, "1", spec, NULL, sf::st_crs(32632), config)
stopifnot(is.na(empty$dbh_m), empty$dbh_quality_label == "dbh_unavailable", empty$dbh_band_points == 0)

cat("stem-band DBH tests passed\n")
