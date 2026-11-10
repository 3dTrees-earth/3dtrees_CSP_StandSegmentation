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
stopifnot(identical(sort(names(trees)), c("1", "2")),
          trees[["1"]]$count == 2, identical(trees[["1"]]$grid$Z, c(1, 2)),
          identical(trees[["1"]]$grid$Zref, c(3, 4)), identical(trees[["2"]]$grid$X, c(2, 4)))

fm_spec <- list(source = "FM", instance = "PredInstance_FM",
                species = "species_id_FM", species_prob = "species_prob_FM")
fm_columns <- c("X", "Y", "Z", "PredInstance_FM", "species_id_FM",
                "species_prob_FM", "PredScore_FM", "PredSemantic_FM")
fm_points <- data.table::data.table(
  X = c(1, 2, 3, 4), Y = c(1, 2, 3, 4), Z = c(3, 5, 4, 6),
  PredInstance_FM = c(1, 2, 1, 2), species_id_FM = c(7, 8, 7, 8),
  species_prob_FM = c(0.8, 0.7, 0.9, 0.6),
  PredScore_FM = c(0.4, 0.5, 0.6, 0.7), PredSemantic_FM = c(1, 2, 2, 1)
)
fm_partition <- file.path(work, "fm.bin")
append_partition(fm_partition, fm_points)
fm_trees <- collect_partition_trees(fm_partition, fm_columns, dtm, fm_spec, test_config, batch_points = 2L)
row <- finish_tree(fm_trees[["1"]], 1, fm_spec, test_config, NULL, NULL)
stopifnot(row$point_count == 2, row$x == 2, row$y == 2,
          row$z == 3.5, row$height_m == 1,  # mean stem-band Zref, extent
          row$species_id == 7,
          isTRUE(all.equal(row$species_prob, 0.85)),
          row$wood_point_count == 1, row$leaf_point_count == 1,
          # Every FM score of a merged instance, counted over all batches.
          identical(row$pred_score_fm_values,
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

without_semantic <- data.table::copy(fm_points)
without_semantic[, PredSemantic_FM := NULL]
without_semantic <- normalize_inventory_points(without_semantic, dtm)
row <- inventory_from_points(without_semantic, fm_spec, test_config, NULL, NULL)[instance_id == 1]
stopifnot(is.na(row$wood_point_count), is.na(row$leaf_point_count))
