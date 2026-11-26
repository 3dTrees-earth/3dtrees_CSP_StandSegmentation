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
# One segmentation's rows of inventory.csv; without the evidence columns for
# outputs written without DBH evidence.
read_inventory <- function(path, slug, evidence = TRUE) {
  rows <- data.table::as.data.table(utils::read.csv(file.path(path, "inventory.csv")))
  rows <- rows[segmentation_source == slug]
  if (!evidence) rows[, c("dbh_fits", "dbh_slices", "crown_hull") := NULL]
  rows
}
config <- parse_cli_args(c("--input", input, "--output-dir", file.path(work, "owned"), "--debug-outputs", "true",
                          "--instance-attribute", "PredInstance_FM", "--species-attribute", "species_id_FM", "--species-probability-attribute", "species_prob_FM",
                          "--score-attribute", "PredScore_FM",
                          "--instance-attribute", "PredInstance_SAT", "--species-attribute", "species_id_FM", "--species-probability-attribute", "species_prob_FM",
                          "--instance-attribute", "PredInstance_ULS", "--dtm", dtm_path,
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
  inventory <- read_inventory(config$output_dir, slug)
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
stopifnot(all(c("parquet", "csv") %in% vapply(manifest$artifacts, `[[`, character(1), "format")))

# Compare deterministic metrics against the existing partition implementation.
baseline_config <- config
baseline_config$output_dir <- file.path(work, "partition-baseline")
baseline_config$inventory_read_mode <- "sequential"
baseline_config$inventory_workers <- 1L
run_tool(baseline_config)
for (spec in config$segmentation_specs) {
  columns <- c("instance_id", "point_count", "height_m", "convex_hull_area_m2")
  if (spec$source == "FM") columns <- c(columns, "species_id", "species_prob", "pred_score_values")
  read <- function(path) read_inventory(path, spec$source)[, ..columns]
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
  read <- function(path) read_inventory(path, spec$source, evidence = FALSE)
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
  read <- function(path) read_inventory(path, spec$source, evidence = FALSE)
  stopifnot(identical(read(serial_config$output_dir), read(regroup_output)))
}
stopifnot(instance_inventory_seed(42, "PredInstance_FM", "2") == 2091750084L,
          instance_inventory_seed(42, "PredInstance_FM", "2") != instance_inventory_seed(42, "PredInstance_SAT", "2"),
          instance_inventory_seed(42, "PredInstance_FM", "2") != instance_inventory_seed(43, "PredInstance_FM", "2"))

# Other cores and worker counts give the same rows: every tree is collected
# the same way and seeded from its identity.
config$inventory_workers <- 3L
config$inventory_core_size <- 100
config$inventory_buffer <- 20
# Called without run_tool, so resolve the auto defaults (index workers) here.
config <- resolve_inventory_read_mode(config, input)
stream <- index_owned_inventory(input, config$segmentation_specs, config, work)
result <- inventory_from_owned_windows(stream, terra::rast(dtm_path), config, NULL, sf::st_crs(cloud))
for (slug in names(result)) {
  baseline <- read_inventory(config$output_dir, slug)
  rows <- result[[slug]]
  stopifnot(identical(as.numeric(rows$point_count), as.numeric(baseline$point_count)))
  for (column in c("x", "y", "z", "dbh_m", "height_m", "convex_hull_area_m2", "metrics_point_count",
                   "dbh_quality_code", "dbh_band_points")) {
    stopifnot(isTRUE(all.equal(as.numeric(rows[[column]]), as.numeric(baseline[[column]]), tolerance = 1e-9)))
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
tree_point <- which(cloud$PredInstance_FM > 0)[[1L]]
bad_cell <- terra::cellFromXY(bad_dtm, cbind(cloud$X[tree_point], cloud$Y[tree_point]))
bad_values <- terra::values(bad_dtm, mat = FALSE)
bad_values[bad_cell] <- NA_real_
terra::values(bad_dtm) <- bad_values
bad_dtm_path <- file.path(work, "bad-dtm.tif")
terra::writeRaster(bad_dtm, bad_dtm_path)
error <- tryCatch({
  suppressWarnings(read_parallel_owned_window(failure_stream, stream$layouts[[1L]],
    stream$owned_index$layouts[[1L]]$windows[[1L]], 1L, bad_dtm_path,
    modifyList(config, list(dtm_fill_radius = 0L)), NULL, sf::st_crs(cloud)))
  NULL
}, error = conditionMessage)
nodata_message <- "DTM has missing or non-finite cells under inventoried points"
stopifnot(!is.null(error), grepl("instance worker reported an R error", error, fixed = TRUE),
          grepl(nodata_message, error, fixed = TRUE),
          !length(list.files(failure_stream$work_dir, pattern = "owned-pipes", all.files = TRUE)))

failed_config <- config
failed_config$dtm <- bad_dtm_path
failed_config$dtm_fill_radius <- 0L
failed_config$inventory_workers <- 4L  # fork concurrent windows as well as instance workers
failed_config$output_dir <- file.path(work, "failed-nodata")
prior_memory <- Sys.getenv("GALAXY_MEMORY_MB", unset = NA_character_)
Sys.setenv(GALAXY_MEMORY_MB = "8192")
error <- tryCatch({ run_tool(failed_config); "" }, error = conditionMessage, finally = {
  if (is.na(prior_memory)) Sys.unsetenv("GALAXY_MEMORY_MB")
  else Sys.setenv(GALAXY_MEMORY_MB = prior_memory)
})
stopifnot(grepl("window worker reported an R error", error, fixed = TRUE),
          grepl(nodata_message, error, fixed = TRUE), !dir.exists(failed_config$output_dir),
          length(regmatches(error, gregexpr("; no inventory published", error, fixed = TRUE))[[1L]]) == 1L)

# Exercise the real parent and FIFO workers with a dispatcher that writes a
# partial record before failing. Worker EOF errors must preserve the cause.
dispatcher_probe <- file.path(work, "fail-dispatch.py")
writeLines(c(
  "import json, os, signal, struct, sys, time",
  "from contextlib import ExitStack",
  "from pathlib import Path",
  "mode, group_path = sys.argv[1:]",
  "with ExitStack() as stack:",
  "    groups = json.loads(Path(group_path).read_text())",
  "    for group in groups:",
  "        pipe = stack.enter_context(open(group['pipe'], 'wb', buffering=0))",
  "        pipe.write(struct.pack('<d', 0))",
  "    if mode == 'kill':",
  "        print('dispatcher probe before signal', file=sys.stderr, flush=True)",
  "        os.kill(os.getpid(), signal.SIGKILL)",
  "time.sleep(0.5)",
  "raise RuntimeError('injected COPC decode error')"
), dispatcher_probe)
original_reader_command <- owned_reader_command
for (mode in c("raise", "kill")) {
  owned_reader_command <- function(config, arguments) {
    if (arguments[[1L]] != "dispatch") return(original_reader_command(config, arguments))
    groups <- arguments[[match("--groups", arguments) + 1L]]
    paste("python3", shQuote(dispatcher_probe), shQuote(mode), shQuote(groups))
  }
  error <- tryCatch({
    read_parallel_owned_window(failure_stream, stream$layouts[[1L]],
      stream$owned_index$layouts[[1L]]$windows[[1L]], 1L, dtm_path,
      config, NULL, sf::st_crs(cloud))
    ""
  }, error = conditionMessage)
  expected <- if (mode == "raise") {
    c("exit status 1", "Traceback", "RuntimeError: injected COPC decode error")
  } else {
    c("exit status 137", "signal 9", "out-of-memory", "dispatcher probe before signal")
  }
  stopifnot(grepl("Owned COPC dispatcher failed", error, fixed = TRUE),
            all(vapply(expected, grepl, logical(1), x = error, fixed = TRUE)))
}
owned_reader_command <- original_reader_command

killed <- parallel::mcparallel(tools::pskill(Sys.getpid(), signal = 9L))
killed_result <- suppressWarnings(parallel::mccollect(killed))
error <- tryCatch({ check_owned_worker(killed_result[[as.character(killed$pid)]], "window"); "" },
                  error = conditionMessage)
stopifnot(grepl("exited without a result", error, fixed = TRUE))

sloped <- terra::rast(dtm_path)
slope <- 2 + 0.02 * terra::colFromCell(sloped, seq_len(terra::ncell(sloped)))
slope[[bad_cell]] <- NA_real_
terra::values(sloped) <- slope
sloped_path <- file.path(work, "sloped.tif")
terra::writeRaster(sloped, sloped_path)
bilinear <- parallel_config
bilinear$segmentation_specs <- bilinear$segmentation_specs[1L]
bilinear$dtm <- sloped_path
bilinear$dtm_lookup <- "bilinear"
bilinear$dtm_fill_radius <- 1L
bilinear$output_dir <- file.path(work, "bilinear-owned")
logs <- capture.output(run_tool(bilinear), type = "message")
stopifnot("CSP workflow: bilinear DTM lookup" %in% logs)
manifest <- jsonlite::fromJSON(file.path(bilinear$output_dir, "manifest.json"))
metadata <- jsonlite::fromJSON(file.path(bilinear$output_dir, "run_metadata.json"))
stopifnot(identical(manifest$dtm$terrain_lookup, "bilinear"),
          identical(manifest$dtm$fill_radius_cells, 1L), manifest$dtm$filled_cells == 1L,
          manifest$dtm$fill_cap_m == 1, manifest$dtm$max_fill_distance_m == 1,
          identical(metadata$terrain_lookup, "bilinear"),
          identical(metadata$dtm_fill_radius_cells, 1L), metadata$dtm_filled_cells == 1L,
          metadata$dtm_fill_cap_m == 1, metadata$dtm_max_fill_distance_m == 1)
owned_bilinear <- read_inventory(bilinear$output_dir, "FM")
bilinear$output_dir <- file.path(work, "bilinear-sequential")
bilinear$inventory_read_mode <- "sequential"
bilinear$inventory_index_workers <- 0L
run_tool(bilinear)
columns <- c("instance_id", "point_count", "x", "y", "z", "height_m", "convex_hull_area_m2")
stopifnot(isTRUE(all.equal(owned_bilinear[, ..columns],
                          read_inventory(bilinear$output_dir, "FM")[, ..columns],
                          tolerance = 1e-7)))
stopifnot(!length(list.files(work, pattern = "\\.bin$", recursive = TRUE)))
unlink(work, recursive = TRUE)
cat("complete-instance owned COPC integration tests passed\n")
