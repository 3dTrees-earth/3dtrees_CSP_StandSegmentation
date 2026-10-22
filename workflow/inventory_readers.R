# Reader-specific preparation stays behind this boundary. The runner consumes
# a hull, input statistics and one measurement operation, never partition paths.
prepare_inventory_reader <- function(files, bounds, config, work_dir, stage_dir) {
  prepare <- switch(config$inventory_read_mode,
                    "copc-owned" = prepare_owned_reader,
                    "sequential" = prepare_partition_reader,
                    "spatial" = prepare_partition_reader,
                    abort("Unsupported inventory reader"))
  prepare(files, bounds, config, work_dir, stage_dir)
}

reader_input_metadata <- function(stream) {
  list(input_points = stream$point_count, input_dimensions = stream$loaded_dimensions,
       input_read_selectors = unname(stream$selectors),
       input_source_point_reads = stream$source_point_reads,
       input_dimensions_dropped = stream$dropped_dimensions,
       input_decompressed_groups = stream$decompressed_groups,
       largest_instances = stream$largest_instances, oversized_instances = stream$oversized_instances)
}

prepare_owned_reader <- function(files, bounds, config, work_dir, stage_dir) {
  stream <- index_owned_inventory(files, config$segmentation_specs, config, work_dir)
  if (!file.copy(stream$owned_index_path, file.path(stage_dir, "instance_index.json"))) {
    abort("Could not preserve the complete-instance index")
  }
  metadata <- c(reader_input_metadata(stream), list(
    read_mode = "copc_owned_windows", partition_elapsed_seconds = 0,
    partition_bytes = 0, largest_partition_bytes = 0, inventory_partitions = 0,
    partition_write_seconds = 0, xy_spool_bytes = 0, spatial_bucket_bytes = 0,
    inventory_memory_strategy = "owned_copc_pipe_complete_instances",
    inventory_seed_policy = "instance_identity_v1"))
  measure <- function(dtm, aoi) {
    inventories <- inventory_from_owned_windows(stream, dtm, config, aoi, bounds$crs)
    audits <- attr(inventories, "copc_audits")
    owned <- list(
      core_size_m = config$inventory_core_size, initial_buffer_m = config$inventory_buffer,
      batch_point_limit = config$inventory_batch_points, retained_point_budget = config$inventory_retained_points,
      shared_decode_threads = config$inventory_decode_threads,
      index_elapsed_seconds = stream$owned_index$elapsed_seconds,
      index_audit = stream$owned_index$audit, windows = audits,
      instance_index = lapply(stream$owned_index$layouts, function(x) list(source = x$slug, instances = x$instances)),
      decoded_points_total = stream$owned_index$audit$decoded_points +
        sum(vapply(unlist(audits, recursive = FALSE), `[[`, numeric(1), "decoded_points")))
    list(inventories = inventories, metadata = list(owned_copc = owned))
  }
  list(point_count = stream$point_count, hull = stream$hull, metadata = metadata, measure = measure)
}

prepare_partition_reader <- function(files, bounds, config, work_dir, stage_dir) {
  chunks <- spatial_chunks(bounds, config$read_chunk_size)
  sequential <- identical(config$inventory_read_mode, "sequential")
  started <- Sys.time()
  stream <- if (sequential) {
    stream_inventory_partitions_python(files, chunks, config$segmentation_specs, config, work_dir)
  } else {
    stream_inventory_partitions(files, chunks, config$segmentation_specs, config, work_dir, bounds$crs)
  }
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  paths <- unlist(lapply(stream$layouts, function(layout) {
    list.files(layout$directory, pattern = "\\.bin$", full.names = TRUE)
  }), use.names = FALSE)
  sizes <- file.info(paths)$size
  metadata <- c(reader_input_metadata(stream), list(
    read_mode = if (sequential) "sequential_copc_and_instance_partitions" else "spatial_chunks_and_instance_partitions",
    read_chunk_size_m = config$read_chunk_size, chunk_buffer_m = config$chunk_buffer,
    spatial_chunk_count = if (sequential) 0L else length(chunks),
    spatial_bucket_count = if (sequential) length(chunks) else 0L,
    partition_elapsed_seconds = elapsed, spatial_read_seconds = stream$spatial_read_seconds,
    partition_write_seconds = stream$partition_write_seconds, hull_elapsed_seconds = stream$hull_elapsed_seconds,
    xy_spool_bytes = stream$xy_spool_bytes, spatial_bucket_bytes = stream$spatial_bucket_bytes,
    sequential_scan_elapsed_seconds = stream$sequential_elapsed_seconds,
    sequential_scan_peak_rss_kb = stream$sequential_peak_rss_kb,
    partition_bytes = sum(sizes), largest_partition_bytes = if (length(sizes)) max(sizes) else 0,
    inventory_partitions = config$inventory_partitions,
    inventory_memory_strategy = "batched_dtm_ordered_instances_and_streaming_large_summaries",
    inventory_seed_policy = if (config$inventory_workers > 1L) "base_seed_plus_partition_index" else "sequential_stream"))
  measure <- function(dtm, aoi) {
    list(inventories = inventory_from_partitions(stream, dtm, config, aoi, bounds$crs), metadata = list())
  }
  list(point_count = stream$point_count, hull = stream$hull, metadata = metadata, measure = measure)
}
