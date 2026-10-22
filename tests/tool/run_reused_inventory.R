#!/usr/bin/env Rscript

# Local recovery benchmark for a completed partition scan and generated DTM.
# Inputs: COPC PARTITIONS AUDIT HULL DTM OUTPUT_DIR WORKERS
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop("Usage: run_reused_inventory.R COPC PARTITIONS AUDIT HULL DTM OUTPUT_DIR WORKERS")
}
arguments <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[[1L]])
repository_root <- normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
source(file.path(repository_root, "workflow", "cli.R"))
source(file.path(repository_root, "workflow", "tool.R"))
source(file.path(repository_root, "workflow", "outputs.R"))

run <- function() {
  started <- Sys.time()
  input <- normalizePath(args[[1L]], mustWork = TRUE)
  partition_root <- normalizePath(args[[2L]], mustWork = TRUE)
  audit <- jsonlite::fromJSON(args[[3L]], simplifyVector = TRUE)
  hull_file <- normalizePath(args[[4L]], mustWork = TRUE)
  supplied_dtm <- normalizePath(args[[5L]], mustWork = TRUE)
  output_dir <- args[[6L]]
  workers <- as.integer(args[[7L]])
  if (is.na(workers) || workers < 1L || workers > 20L) abort("Workers must be 1–20")
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    abort("Output directory is not empty")
  }
  header <- jsonlite::fromJSON(file.path(dirname(supplied_dtm), "copc_header.json"),
                               simplifyVector = TRUE)
  bounds <- list(xmin = header$xmin, xmax = header$xmax,
                 ymin = header$ymin, ymax = header$ymax,
                 crs = sf::st_crs(header$crs_wkt))
  validate_inventory_crs(bounds$crs)
  input_points <- as.numeric(header$point_count)
  original_scan <- jsonlite::fromJSON(file.path(dirname(partition_root), "metadata.json"),
                                      simplifyVector = TRUE)
  if (input_points != original_scan$input_points ||
      audit$partition_counts$FM$rows != original_scan$assigned_points$FM) {
    abort("COPC header, partition audit and original scan point counts disagree")
  }
  coords <- readBin(hull_file, numeric(), n = file.info(hull_file)$size / 8,
                    size = 8L, endian = "little")
  if (length(coords) < 6L || length(coords) %% 2L) abort("Invalid hull artifact")
  xy <- matrix(coords, ncol = 2L, byrow = TRUE)
  area <- list(area_m2 = polygon_area(xy[, 1L], xy[, 2L]),
               area_source = "recovered_all_point_xy_hull")
  spec <- list(source = "FM", instance = "PredInstance_FM",
               species = "species_id_FM", species_prob = "species_prob_FM")
  layout <- list(spec = spec,
                 columns = c("X", "Y", "Z", "PredInstance_FM", "species_id_FM",
                             "species_prob_FM", "PredScore_FM", "PredSemantic_FM"),
                 directory = file.path(partition_root, "FM"))
  stream <- list(layouts = list(FM = layout), point_count = input_points,
                 largest_instances = audit$largest_instances,
                 oversized_instances = audit$oversized_instances)
  spool_parent <- tempfile("csp-reused-spool-", tmpdir = dirname(output_dir))
  dir.create(spool_parent)
  on.exit(unlink(spool_parent, recursive = TRUE, force = TRUE), add = TRUE)
  config <- list(
    segmentation_specs = list(spec), output_formats = c("csv", "tsv", "json", "gpkg"),
    species_lookup = NULL, inventory_workers = workers,
    inventory_spool_dir = spool_parent,
    non_tree_ids = 0, random_seed = 42L, geometry_threads = 1L,
    slice_min = 0.3, slice_max = 4, slice_increment = 0.2,
    slice_width = 0.1, max_dbh = 1,
    dbh_band_points = inventory_dbh_band_point_budget, dbh_runs = 5L, geometry_k = 10L
  )
  stage_dir <- tempfile("csp-reused-stage-", tmpdir = dirname(output_dir))
  dir.create(stage_dir)
  on.exit(unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)
  dtm_result <- load_supplied_dtm(supplied_dtm, bounds$crs, bounds, stage_dir)
  dtm_generation <- jsonlite::fromJSON(
    file.path(dirname(supplied_dtm), "dtm_generation.json"),
    simplifyVector = TRUE)
  dtm_result$method <- dtm_generation$method
  dtm_result$strategy <- dtm_generation$strategy
  dtm_result$candidate_metadata <- if (!is.null(dtm_generation$candidate)) {
    jsonlite::fromJSON(
      file.path(dirname(dtm_generation$candidate), "candidate_metadata.json"),
      simplifyVector = TRUE)
  } else NULL
  inventory_started <- Sys.time()
  inventories <- inventory_from_partitions(stream, dtm_result$dtm, config,
                                           NULL, bounds$crs)
  inventory_seconds <- as.numeric(difftime(Sys.time(), inventory_started, units = "secs"))
  if (nrow(inventories$FM) != audit$instance_counts$FM ||
      sum(inventories$FM$point_count) != audit$partition_counts$FM$rows ||
      anyDuplicated(inventories$FM$instance_id)) {
    abort("Inventory output IDs or point totals do not match the partition audit")
  }
  exports <- write_outputs(inventories, area, stage_dir, config, bounds$crs)
  jsonlite::write_json(list(
    tool = "3Dtrees: CSP StandSegmentation", input_files = input,
    input_bytes = file.info(input)$size, input_points = input_points,
    read_mode = "reused_completed_instance_partitions",
    partition_root = partition_root,
    partition_bytes = audit$partition_counts$FM$bytes,
    largest_instances = audit$largest_instances,
    oversized_instances = audit$oversized_instances,
    inventory_instance_point_limit = inventory_instance_point_limit,
    inventory_dbh_band_points = config$dbh_band_points,
    inventory_large_instance_sample_limit = inventory_large_instance_sample_limit,
    inventory_workers = workers, inventory_elapsed_seconds = inventory_seconds,
    csp_enabled = FALSE, output_formats = config$output_formats,
    area_source = area$area_source, area_m2 = area$area_m2,
    ground_method = dtm_result$method,
    dtm_candidate_metadata = dtm_result$candidate_metadata
  ), file.path(stage_dir, "run_metadata.json"), auto_unbox = TRUE,
  pretty = TRUE, na = "null")
  jsonlite::write_json(list(
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    process_peak_rss_kb = process_peak_rss_kb(),
    inventory_workers = workers,
    note = "Process-local VmHWM; cgroup peak printed by this recovery benchmark."
  ), file.path(stage_dir, "resource_summary.json"), auto_unbox = TRUE,
  pretty = TRUE, na = "null")
  write_inventory_manifest(exports, config, bounds$crs, dtm_result, stage_dir)
  publish_stage(stage_dir, output_dir)
  cat(sprintf("inventory_elapsed_seconds=%.3f trees=%d points=%.0f\n",
              inventory_seconds, nrow(inventories$FM),
              sum(inventories$FM$point_count)))
  if (file.exists("/sys/fs/cgroup/memory.peak")) {
    cat(sprintf("cgroup_memory_peak_bytes=%s\n",
                readLines("/sys/fs/cgroup/memory.peak", n = 1L)))
  }
}

run()
