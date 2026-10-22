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
partition <- file.path(work, "partition.bin")
append_partition(partition, points)
spool <- spool_partition_instances(partition, columns, "PredInstance_SAT", dtm,
                                   file.path(work, "spool"), batch_points = 2L)
stopifnot(length(spool$paths) == 2L, identical(unlist(spool$counts), c(`1` = 2L, `2` = 2L)))
first <- read_binary_partition(spool$paths[["1"]], spool$columns)
second <- read_binary_partition(spool$paths[["2"]], spool$columns)
stopifnot(identical(first$X, c(1, 3)), identical(second$X, c(2, 4)),
          identical(first$Z, c(1, 2)), identical(second$Z, c(3, 4)),
          identical(first$Zref, c(3, 4)),
          identical(first$source_order__, c(1, 3)),
          identical(second$source_order__, c(2, 4)))

failure <- tryCatch({
  spool_partition_instances(partition, columns, "PredInstance_SAT", dtm,
                            file.path(work, "limited"), point_limit = 1L,
                            batch_points = 2L)
  NULL
}, error = conditionMessage)
if (!is.null(failure)) message(failure)
stopifnot(!is.null(failure), grepl("bounded inventory limit", failure))
giant <- tryCatch({
  inventory_from_partitions(
    list(largest_instances = list(FM = list(instance_id = 116,
                                             points = 109750559))),
    dtm, list(), NULL, NULL
  )
  NULL
}, error = conditionMessage)
stopifnot(!is.null(giant), grepl("no inventory workers started", giant))

test_config <- list(slice_min = 0.3, slice_max = 4, slice_increment = 0.2,
                    slice_width = 0.1, max_dbh = 1, geometry_threads = 1L,
                    non_tree_ids = 0, random_seed = 42L, dbh_band_points = inventory_dbh_band_point_budget, dbh_runs = 5L, geometry_k = 10L)
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
streamed <- spool_partition_instances(
  fm_partition, fm_columns, "PredInstance_FM", dtm,
  file.path(work, "streaming"), batch_points = 2L,
  large_counts = c(`1` = 2), spec = fm_spec, config = test_config
)
stopifnot(identical(names(streamed$paths), "2"),
          streamed$large_states[["1"]]$count == 2)
large_row <- finish_large_instance_summary(streamed$large_states[["1"]],
                                           "1", fm_spec, NULL, NULL, test_config)
stopifnot(large_row$point_count == 2, large_row$x == 2, large_row$y == 2,
          large_row$z == 3, large_row$height_m == 1,
          large_row$species_id == 7,
          isTRUE(all.equal(large_row$species_prob, 0.85)),
          large_row$wood_point_count == 1, large_row$leaf_point_count == 1,
          is.na(large_row$dbh_m), large_row$dbh_quality_code == 9L,
          large_row$inventory_method == "streaming_large_instance")

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
  "--segmentation-spec", "PredInstance_SAT",
  "--output-formats", "csv,tsv,json,gpkg", "--read-chunk-size", "10",
  "--inventory-partitions", "8", "--inventory-workers", "1",
  "--inventory-read-mode", "spatial"
), stdout = log, stderr = log))
if (!identical(status, 0L)) stop(paste(readLines(log, warn = FALSE), collapse = "\n"))
inventory <- data.table::fread(file.path(output, "inventories", "SAT.csv"))
stopifnot(identical(inventory$instance_id, c(1L, 2L)),
          sum(inventory$point_count) > 0,
          all(is.finite(inventory$x)), all(is.finite(inventory$y)),
          all(file.exists(file.path(output, "inventories",
                                    paste0("SAT.", c("csv", "tsv", "json"))))),
          file.exists(file.path(output, "inventory.gpkg")),
          file.exists(file.path(output, "manifest.json")))

generated_output <- file.path(work, "generated")
generated_log <- file.path(work, "generated.log")
generated_status <- suppressWarnings(system2("Rscript", c(
  file.path(root, "exec", "run.R"), "--input", input,
  "--output-dir", generated_output,
  "--segmentation-spec", "PredInstance_SAT",
  "--output-formats", "csv", "--read-chunk-size", "10",
  "--inventory-partitions", "8", "--inventory-workers", "1",
  "--inventory-read-mode", "spatial", "--dtm-workers", "1"
), stdout = generated_log, stderr = generated_log))
if (!identical(generated_status, 0L)) {
  stop(paste(readLines(generated_log, warn = FALSE), collapse = "\n"))
}
generated_inventory <- data.table::fread(file.path(generated_output, "inventories", "SAT.csv"))
generated_manifest <- jsonlite::fromJSON(file.path(generated_output, "manifest.json"))
stopifnot(file.exists(file.path(generated_output, "dtm_full.tif")),
          !identical(generated_manifest$dtm$source, "supplied"),
          identical(generated_inventory$instance_id, inventory$instance_id),
          identical(generated_inventory$point_count, inventory$point_count))
cat("bounded inventory spool test passed\n")

without_semantic <- data.table::copy(fm_points)
without_semantic[, PredSemantic_FM := NULL]
without_semantic <- normalize_inventory_points(without_semantic, dtm)
state <- update_large_instance_summary(new_large_instance_summary(nrow(without_semantic)), without_semantic, fm_spec, test_config)
row <- finish_large_instance_summary(state, '1', fm_spec, NULL, NULL, test_config)
stopifnot(is.na(row$wood_point_count), is.na(row$leaf_point_count))
