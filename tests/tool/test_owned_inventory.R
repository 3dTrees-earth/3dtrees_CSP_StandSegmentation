#!/usr/bin/env Rscript
# End-to-end ownership, complete-instance counts, exports and streaming fallback.
args <- commandArgs(trailingOnly = FALSE)
test_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[[1L]])
repository_root <- normalizePath(file.path(dirname(test_file), "..", ".."))
source(file.path(repository_root, "workflow", "cli.R"))
source(file.path(repository_root, "workflow", "tool.R"))
source(file.path(repository_root, "workflow", "owned_inventory.R"))
source(file.path(repository_root, "workflow", "outputs.R"))
work <- tempfile("csp-owned-test-")
dir.create(work)
input <- file.path(repository_root, "tests/tool/fixtures/multi-instance.copc.laz")
cloud <- suppressWarnings(lidR::readLAS(input, select = "*"))
dtm <- terra::rast(xmin = floor(min(cloud$X)) - 1, xmax = ceiling(max(cloud$X)) + 1,
                   ymin = floor(min(cloud$Y)) - 1, ymax = ceiling(max(cloud$Y)) + 1,
                   resolution = 1, crs = sf::st_crs(cloud)$wkt)
terra::values(dtm) <- 2
dtm_path <- file.path(work, "dtm.tif")
terra::writeRaster(dtm, dtm_path)
config <- parse_cli_args(c("--input", input, "--output-dir", file.path(work, "owned"),
                          "--segmentation-spec", "PredInstance_FM,species_id_FM,species_prob_FM",
                          "--segmentation-spec", "PredInstance_SAT,species_id_FM,species_prob_FM",
                          "--segmentation-spec", "PredInstance_ULS", "--dtm", dtm_path,
                          "--inventory-attributes", "PredScore_FM,PredSemantic_FM",
                          "--inventory-read-mode", "copc-owned", "--inventory-core-size", "1",
                          "--inventory-buffer", "0", "--inventory-batch-points", "25000",
                          "--inventory-workers", "2", "--inventory-decode-threads", "1"))
# A supplied-DTM owned read must never allocate the legacy bounding-box grid.
legacy_spatial_chunks <- spatial_chunks
spatial_chunks <- function(...) stop("Owned reader requested a legacy spatial grid")
run_tool(config)
spatial_chunks <- legacy_spatial_chunks
metadata <- jsonlite::fromJSON(file.path(config$output_dir, "run_metadata.json"), simplifyVector = FALSE)
stopifnot(metadata$read_mode == "copc_owned_windows", metadata$partition_bytes == 0,
          metadata$xy_spool_bytes == 0, metadata$spatial_bucket_bytes == 0,
          metadata$owned_copc$index_audit$decoded_points == nrow(cloud@data))
for (spec in config$segmentation_specs) {
  slug <- spec$source
  inventory <- data.table::fread(file.path(config$output_dir, "inventories", paste0(slug, ".csv")))
  expected <- cloud@data[get(spec$instance) > 0, .(count = .N), by = c(spec$instance)]
  data.table::setorderv(expected, spec$instance)
  stopifnot(identical(as.numeric(inventory$point_count), as.numeric(expected$count)),
            identical(as.numeric(inventory$instance_id), as.numeric(expected[[spec$instance]])),
            !anyDuplicated(inventory$instance_id))
  for (audit in metadata$owned_copc$windows[[slug]]) {
    stopifnot(audit$complete, audit$expanded,
              identical(audit$observed_counts, audit$expected_counts),
              identical(audit$r_observed_counts, audit$expected_counts))
  }
}
manifest <- jsonlite::fromJSON(file.path(config$output_dir, "manifest.json"), simplifyVector = FALSE)
stopifnot(setequal(vapply(manifest$artifacts, `[[`, character(1), "key"),
                   setdiff(list.files(config$output_dir, recursive = TRUE), "manifest.json")))
stopifnot(all(c("csv", "tsv", "json", "gpkg") %in% vapply(manifest$artifacts, `[[`, character(1), "format")))

# Compare deterministic metrics against the existing partition implementation.
baseline_config <- config
baseline_config$output_dir <- file.path(work, "partition-baseline")
baseline_config$inventory_read_mode <- "sequential"
baseline_config$inventory_workers <- 1L
run_tool(baseline_config)
for (spec in config$segmentation_specs) {
  columns <- c("instance_id", "point_count", "height_m", "convex_hull_area_m2")
  if (spec$source == "FM") columns <- c(columns, "species_id", "species_prob", "wood_point_count", "leaf_point_count")
  read <- function(path) data.table::fread(file.path(path, "inventories", paste0(spec$source, ".csv")))[, ..columns]
  stopifnot(isTRUE(all.equal(read(config$output_dir), read(baseline_config$output_dir),
                            tolerance = 1e-7, check.attributes = FALSE)))
}

# A single spatial core must still use multiple ID workers. Per-instance
# order and RNG are invariant to the number/scheduling of workers.
parallel_config <- config
parallel_config$inventory_core_size <- 100
parallel_config$inventory_buffer <- 20
parallel_config$inventory_workers <- 3L
parallel_config$inventory_index_workers <- 8L
parallel_config$output_dir <- file.path(work, "one-core-parallel")
run_tool(parallel_config)
serial_config <- parallel_config
serial_config$inventory_workers <- 1L
serial_config$inventory_index_workers <- 0L
serial_config$output_dir <- file.path(work, "one-core-serial")
run_tool(serial_config)
parallel_metadata <- jsonlite::fromJSON(file.path(parallel_config$output_dir, "run_metadata.json"), simplifyVector = FALSE)
serial_metadata <- jsonlite::fromJSON(file.path(serial_config$output_dir, "run_metadata.json"), simplifyVector = FALSE)
stopifnot(parallel_metadata$owned_copc$index_audit$index_workers == 8L,
          parallel_metadata$owned_copc$index_audit$peak_inflight_batches <= 8L)
stopifnot(parallel_metadata$owned_copc$decoded_points_total == serial_metadata$owned_copc$decoded_points_total)
for (spec in config$segmentation_specs) {
  audits <- parallel_metadata$owned_copc$windows[[spec$source]]
  stopifnot(length(audits) == 1L, audits[[1L]]$worker_count > 1L,
            identical(audits[[1L]]$expected_counts, audits[[1L]]$r_observed_counts))
  read <- function(path) readLines(file.path(path, "inventories", paste0(spec$source, ".csv")))
  stopifnot(identical(read(parallel_config$output_dir), read(serial_config$output_dir)))
}

# Repartition the identical points into different windows. The measurement
# seed belongs to an instance, not to a resource-dependent window number.
regroup_dir <- file.path(work, "regrouped")
dir.create(regroup_dir)
regroup_stream <- index_owned_inventory(input, config$segmentation_specs, parallel_config, regroup_dir)
# Preserve nanosecond source identities as Python integers; an R JSON rewrite
# would round these 64-bit values and correctly fail the source-change guard.
regroup_script <- file.path(regroup_dir, "regroup.py")
writeLines(c(
  "import json, sys",
  "from pathlib import Path",
  "path = Path(sys.argv[1])",
  "index = json.loads(path.read_text())",
  "for layout in index['layouts']:",
  "    original = layout['windows'][0]",
  "    layout['windows'] = [dict(original, ids=[key], counts=[count]) for key, count in zip(original['ids'], original['counts'])]",
  "path.write_text(json.dumps(index))"
), regroup_script)
stopifnot(system2("python3", c(shQuote(regroup_script), shQuote(regroup_stream$owned_index_path))) == 0L)
regroup_stream$owned_index <- jsonlite::fromJSON(regroup_stream$owned_index_path, simplifyVector = FALSE)
regrouped <- inventory_from_owned_windows(regroup_stream, terra::rast(dtm_path),
                                           parallel_config, NULL, sf::st_crs(cloud))
regroup_output <- file.path(work, "regroup-output")
dir.create(regroup_output)
invisible(write_outputs(regrouped, list(area_m2 = 1, area_source = "fixture"), regroup_output,
              parallel_config, sf::st_crs(cloud)))
for (spec in config$segmentation_specs) {
  read <- function(path) readLines(file.path(path, "inventories", paste0(spec$source, ".csv")))
  stopifnot(identical(read(serial_config$output_dir), read(regroup_output)))
}
stopifnot(instance_inventory_seed(42, "PredInstance_FM", "2") == 2091750084L,
          instance_inventory_seed(42, "PredInstance_FM", "2") != instance_inventory_seed(42, "PredInstance_SAT", "2"),
          instance_inventory_seed(42, "PredInstance_FM", "2") != instance_inventory_seed(43, "PredInstance_FM", "2"))

# Same full-point count contract for the bounded large-instance path. Lower
# the threshold for this fixture; no synthetic giant cloud or disk spool.
config$inventory_workers <- 3L
config$inventory_core_size <- 100
config$inventory_buffer <- 20
stream <- index_owned_inventory(input, config$segmentation_specs, config, work)
inventory_instance_point_limit <- 100L
result <- inventory_from_owned_windows(stream, terra::rast(dtm_path), config, NULL, sf::st_crs(cloud))
for (slug in names(result)) {
  baseline <- data.table::fread(file.path(config$output_dir, "inventories", paste0(slug, ".csv")))
  rows <- result[[slug]]
  stopifnot(identical(as.numeric(rows$point_count), as.numeric(baseline$point_count)),
            all(rows$inventory_method == "streaming_large_instance"), all(is.na(rows$dbh_m)))
  spec <- stream$layouts[[slug]]$spec
  for (id in rows$instance_id) {
    points <- cloud@data[get(spec$instance) == id]
    row <- rows[instance_id == id]
    stopifnot(row$x == median(points$X), row$y == median(points$Y),
              abs(row$height_m - diff(range(points$Z))) < 1e-8)
  }
}

# The consumer independently rejects disagreement with the global count.
bad_plan <- stream$owned_index$layouts[[1L]]$windows[[1L]]
bad_plan$counts[[1L]] <- bad_plan$counts[[1L]] + 1
error <- tryCatch({
  read_owned_window(stream, stream$layouts[[1L]], bad_plan, 1L,
                    terra::rast(dtm_path), config, NULL, sf::st_crs(cloud))
  NULL
}, error = conditionMessage)
stopifnot(!is.null(error), grepl("Incomplete owned instances", error))
# Consumer failure after a valid read must still abort publication and clean
# up worker processes/FIFOs. A complete Python audit alone is not sufficient.
failure_stream <- stream
failure_stream$work_dir <- file.path(work, "failed-consumer")
dir.create(failure_stream$work_dir)
bad_dtm <- terra::rast(dtm_path)
terra::values(bad_dtm) <- NA_real_
bad_dtm_path <- file.path(work, "bad-dtm.tif")
terra::writeRaster(bad_dtm, bad_dtm_path)
error <- tryCatch({
  suppressWarnings(read_parallel_owned_window(failure_stream, stream$layouts[[1L]],
    stream$owned_index$layouts[[1L]]$windows[[1L]], 1L, bad_dtm_path,
    config, NULL, sf::st_crs(cloud)))
  NULL
}, error = conditionMessage)
stopifnot(!is.null(error), grepl("worker failed|dispatcher failed", error),
          !length(list.files(failure_stream$work_dir, pattern = "owned-pipes", all.files = TRUE)))
stopifnot(!length(list.files(work, pattern = "\\.bin$", recursive = TRUE)))
unlink(work, recursive = TRUE)
cat("complete-instance owned COPC integration tests passed\n")
