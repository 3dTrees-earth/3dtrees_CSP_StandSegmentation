#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = FALSE)
test_file <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[[1L]])
root <- normalizePath(file.path(dirname(test_file), "..", ".."), mustWork = TRUE)
source(file.path(root, "workflow", "cli.R"))
source(file.path(root, "workflow", "tool.R"))

work <- tempfile("csp-bounded-")
dir.create(work)
on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)
dtm <- terra::rast(xmin = 0, xmax = 20, ymin = 0, ymax = 20,
                   resolution = 1, crs = "EPSG:32632")
terra::values(dtm) <- 2
columns <- c("X", "Y", "Z", "PredInstance_SAT")
points <- data.table::data.table(
  X = c(1, 2, 3, 4), Y = c(1, 2, 3, 4),
  Z = c(3, 5, 4, 6), PredInstance_SAT = c(1, 2, 1, 2)
)
test_config <- list(slice_min = 0.3, slice_max = 4, slice_increment = 0.2,
                    slice_width = 0.1, max_dbh = 1, geometry_threads = 1L,
                    non_tree_ids = 0, random_seed = 42L, dbh_runs = 5L, geometry_k = 10L, metrics_grid = 0)
partition <- file.path(work, "partition.bin")
append_partition(partition, points)
sat_spec <- list(source = "SAT", instance = "PredInstance_SAT")
trees <- collect_partition_trees(partition, columns, dtm, sat_spec, test_config, batch_points = 1L)
grid_1 <- tree_grid(trees[["1"]], 0)
stopifnot(identical(sort(names(trees)), c("1", "2")),
          trees[["1"]]$count == 2, identical(grid_1$Z, c(1, 2)),
          identical(grid_1$Zref, c(3, 4)), identical(tree_grid(trees[["2"]], 0)$X, c(2, 4)))
# With the input's scale and offset the grid is stored as LAS integers and
# decodes to the same coordinates and heights.
framed <- collect_partition_trees(partition, columns, dtm, sat_spec,
                                  c(test_config, list(coordinate_frame = list(
                                    scale = c(X = 0.001, Y = 0.001, Z = 0.001), offset = c(X = 0, Y = 0, Z = 0)))),
                                  batch_points = 1L)
stopifnot(is.integer(framed[["1"]]$chunks[[1L]]$X),
          identical(tree_grid(framed[["1"]], 0)[, .(X, Y, Z, Zref)], grid_1[, .(X, Y, Z, Zref)]))

fm_spec <- list(source = "FM", instance = "PredInstance_FM", species = "species_id_FM",
                species_prob = "species_prob_FM", score = "PredScore_FM", non_tree_ids = 0)
fm_columns <- c("X", "Y", "Z", "PredInstance_FM", "species_id_FM",
                "species_prob_FM", "PredScore_FM")
fm_points <- data.table::data.table(
  X = c(1, 2, 3, 4), Y = c(1, 2, 3, 4), Z = c(3, 5, 4, 6),
  PredInstance_FM = c(1, 2, 1, 2), species_id_FM = c(7, 8, 7, 8),
  species_prob_FM = c(0.8, 0.7, 0.9, 0.6),
  PredScore_FM = c(0.4, 0.5, 0.6, 0.7)
)
fm_partition <- file.path(work, "fm.bin")
append_partition(fm_partition, fm_points)
fm_trees <- collect_partition_trees(fm_partition, fm_columns, dtm, fm_spec, test_config, batch_points = 2L)
row <- finish_tree(fm_trees[["1"]], 1, fm_spec, test_config, NULL, NULL)
stopifnot(row$point_count == 2, row$x == 2, row$y == 2,
          row$z == 3.5, row$height_m == 1,  # mean stem-band Zref, extent
          row$species_id == 7,
          isTRUE(all.equal(row$species_prob, 0.85)),
          !any(c("wood_point_count", "leaf_point_count") %in% names(row)),
          # Every score of a merged instance, counted over all batches.
          identical(row$pred_score_values,
                    '[{"score":0.4,"point_count":1},{"score":0.6,"point_count":1}]'),
          is.na(row$dbh_m), row$dbh_quality_code == 1L)  # too few stem points

# Exercise the full inventory handoff on the committed COPC fixture. The
# per-instance spool must still produce the expected tree IDs and counts.
input <- file.path(root, "tests", "tool", "fixtures", "multi-instance.copc.laz")
header <- lidR::readLASheader(input)
bounds <- header@PHB
provided <- file.path(work, "provided.tif")
fixture_dtm <- terra::rast(
  xmin = floor(bounds[["Min X"]]) - 2, xmax = ceiling(bounds[["Max X"]]) + 2,
  ymin = floor(bounds[["Min Y"]]) - 2, ymax = ceiling(bounds[["Max Y"]]) + 2,
  resolution = 1, crs = sf::st_crs(32632)$wkt
)
terra::values(fixture_dtm) <- 2
terra::writeRaster(fixture_dtm, provided)
output <- file.path(work, "result")
log <- file.path(work, "run.log")
status <- suppressWarnings(system2("Rscript", c(
  file.path(root, "exec", "run.R"), "--input", input,
  "--output-dir", output, "--dtm", provided,
  "--instance-attribute", "PredInstance_SAT",
  "--read-chunk-size", "10",
  "--inventory-partitions", "8", "--inventory-workers", "1",
  "--inventory-read-mode", "spatial"
), stdout = log, stderr = log))
if (!identical(status, 0L)) stop(paste(readLines(log, warn = FALSE), collapse = "\n"))
inventory <- data.table::as.data.table(utils::read.csv(file.path(output, "inventory.csv")))
stopifnot(identical(inventory$instance_id, c(1L, 2L)),
          sum(inventory$point_count) > 0,
          all(is.finite(inventory$x)), all(is.finite(inventory$y)),
          # Without --debug-outputs only the results are written.
          setequal(list.files(output), c("dtm_full.tif", "inventory.csv", "inventory.parquet", "stand_summary.csv")))

generated_output <- file.path(work, "generated")
generated_log <- file.path(work, "generated.log")
generated_status <- suppressWarnings(system2("Rscript", c(
  file.path(root, "exec", "run.R"), "--input", input,
  "--output-dir", generated_output,
  "--instance-attribute", "PredInstance_SAT",
  "--read-chunk-size", "10",
  "--inventory-partitions", "8", "--inventory-workers", "1",
  "--inventory-read-mode", "spatial", "--dtm-workers", "1", "--debug-outputs", "true"
), stdout = generated_log, stderr = generated_log))
if (!identical(generated_status, 0L)) {
  stop(paste(readLines(generated_log, warn = FALSE), collapse = "\n"))
}
generated_inventory <- data.table::as.data.table(utils::read.csv(file.path(generated_output, "inventory.csv")))
generated_manifest <- jsonlite::fromJSON(file.path(generated_output, "manifest.json"))
stopifnot(file.exists(file.path(generated_output, "dtm_full.tif")),
          !identical(generated_manifest$dtm$source, "supplied"),
          identical(generated_inventory$instance_id, inventory$instance_id),
          identical(generated_inventory$point_count, inventory$point_count))
cat("partition inventory test passed\n")

# Without a score attribute there is no score column; non-tree IDs are per segmentation.
unscored <- normalize_inventory_points(data.table::copy(fm_points), dtm)
rows <- inventory_from_points(unscored, modifyList(fm_spec, list(score = NULL, non_tree_ids = 2)), test_config, NULL, NULL)
stopifnot(!"pred_score_values" %in% names(rows), identical(rows$instance_id, 1))
