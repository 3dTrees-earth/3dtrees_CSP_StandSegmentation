run_chunked_inventory <- function(config, files, stage_dir, output_parent, started) {
  work_dir <- tempfile("csp-chunks-", tmpdir = scratch_parent(config, output_parent))
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
  bounds <- if (config$inventory_read_mode %in% c("sequential", "copc-owned") &&
                length(files) == 1L) {
    point_cloud_bounds_python(files[[1L]], work_dir)
  } else point_cloud_bounds(files)
  validate_inventory_crs(bounds$crs)
  available <- unique(c(
    "X", "Y", "Z",
    if (is.null(bounds$available_dimensions)) names(inventory_standard_attributes),
    if (!is.null(bounds$available_dimensions)) inventory_canonical_names(bounds$available_dimensions) else
      unlist(lapply(bounds$headers, extra_byte_names), use.names = FALSE)
  ))
  for (spec in config$segmentation_specs) {
    required <- inventory_spec_dimensions(spec)
    missing <- setdiff(required[!is.na(required)], available)
    if (length(missing)) abort(sprintf("Requested dimensions are missing for %s: %s", spec$instance, paste(missing, collapse = ", ")))
  }
  las_crs <- bounds$crs
  aoi <- read_aoi(config$aoi_json, las_crs)
  message("CSP workflow: preparing DTM")
  dtm_started <- Sys.time()
  dtm_result <- make_inventory_dtm(files, bounds, config, work_dir, stage_dir)
  if (!is.null(aoi)) {
    aoi_dtm <- terra::mask(terra::crop(dtm_result$dtm, terra::vect(aoi)), terra::vect(aoi))
    terra::writeRaster(aoi_dtm, file.path(stage_dir, "dtm_aoi.tif"), overwrite = TRUE)
  }
  dtm_elapsed_seconds <- as.numeric(difftime(Sys.time(), dtm_started, units = "secs"))
  message(sprintf("CSP workflow: DTM ready in %.1f s; preparing inventory reader", dtm_elapsed_seconds))
  prepared <- prepare_inventory_reader(files, bounds, config, work_dir, stage_dir)
  area <- if (!is.null(aoi)) {
    list(area_m2 = as.numeric(sf::st_area(aoi)), area_source = "aoi")
  } else {
    list(area_m2 = polygon_area(prepared$hull$x, prepared$hull$y), area_source = "point_cloud_convex_hull")
  }
  inventory_started <- Sys.time()
  measured <- prepared$measure(dtm_result$dtm, aoi)
  inventories <- measured$inventories
  inventory_elapsed_seconds <- as.numeric(difftime(Sys.time(), inventory_started, units = "secs"))
  message(sprintf("CSP workflow: inventory ready in %.1f s; exporting", inventory_elapsed_seconds))
  export_started <- Sys.time()
  exports <- write_outputs(inventories, area, stage_dir, config, las_crs)
  export_elapsed_seconds <- as.numeric(difftime(Sys.time(), export_started, units = "secs"))

  metadata <- c(list(
    tool = "3Dtrees: CSP StandSegmentation",
    package_version = as.character(utils::packageVersion("CspStandSegmentation")),
    input_files = unname(files), input_bytes = unname(sum(file.info(files)$size)),
    segmentation_dimensions = vapply(config$segmentation_specs, `[[`, character(1), "instance"),
    species_dimensions_supplied = any(vapply(config$segmentation_specs, function(spec) !is.null(spec$species), logical(1))),
    csp_enabled = FALSE, csp_seed_count = 0L, csp_unassigned_point_count = 0L,
    dtm_elapsed_seconds = dtm_elapsed_seconds,
    inventory_instance_point_limit = inventory_instance_point_limit,
    inventory_dbh_band_points = config$dbh_band_points,
    inventory_large_instance_sample_limit = inventory_large_instance_sample_limit,
    inventory_elapsed_seconds = inventory_elapsed_seconds,
    export_elapsed_seconds = export_elapsed_seconds,
    dtm_workers = dtm_result$workers,
    dtm_valid_chunk_count = dtm_result$valid_chunk_count,
    dtm_strategy = dtm_result$strategy,
    dtm_candidate_resolution_m = if (is.null(dtm_result$candidate_metadata)) NULL else dtm_result$candidate_metadata$candidate_resolution_m,
    dtm_candidate_xy_mode = if (is.null(dtm_result$candidate_metadata)) NULL else dtm_result$candidate_metadata$candidate_xy_mode,
    dtm_candidate_points = if (is.null(dtm_result$candidate_metadata)) NULL else dtm_result$candidate_metadata$candidate_points,
    dtm_candidate_elapsed_seconds = if (is.null(dtm_result$candidate_metadata)) NULL else dtm_result$candidate_metadata$elapsed_seconds,
    dtm_candidate_average_cpu_cores = if (is.null(dtm_result$candidate_metadata)) NULL else dtm_result$candidate_metadata$average_cpu_cores,
    dtm_candidate_worker_peak_rss_sum_upper_bound_kb = if (is.null(dtm_result$candidate_metadata)) NULL else dtm_result$candidate_metadata$worker_peak_rss_sum_upper_bound_kb,
    dtm_candidate_parent_peak_rss_kb = if (is.null(dtm_result$candidate_metadata)) NULL else dtm_result$candidate_metadata$parent_peak_rss_kb,
    dtm_worker_peak_rss_sum_kb = dtm_result$worker_peak_rss_sum_kb,
    dtm_worker_peak_rss_max_kb = dtm_result$worker_peak_rss_max_kb,
    inventory_workers = config$inventory_workers,
    output_formats = config$output_formats,
    ground_method = dtm_result$method, dtm_resolution_m = as.numeric(terra::res(dtm_result$dtm)),
    area_source = area$area_source, area_m2 = area$area_m2,
    random_seed = config$random_seed, routing_workers = config$routing_workers,
    geometry_threads = config$geometry_threads, voxel_resolution_m = config$voxel_resolution
  ), prepared$metadata, measured$metadata)
  jsonlite::write_json(metadata, file.path(stage_dir, "run_metadata.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
  completed <- Sys.time()
  resources <- list(
    elapsed_seconds = as.numeric(difftime(completed, started, units = "secs")),
    process_peak_rss_kb = process_peak_rss_kb(), input_bytes = unname(sum(file.info(files)$size)),
    input_points = prepared$point_count, routing_workers = config$routing_workers,
    inventory_workers = config$inventory_workers,
    geometry_threads = config$geometry_threads,
    note = "Process-local VmHWM only; container/cgroup telemetry is recorded by benchmark runs."
  )
  jsonlite::write_json(resources, file.path(stage_dir, "resource_summary.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
  write_inventory_manifest(exports, config, las_crs, dtm_result, stage_dir)
  resources
}
