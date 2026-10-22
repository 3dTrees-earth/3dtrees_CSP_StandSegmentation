# Disk-backed inventory reader; measurements live in instance_inventory.R.

append_partition <- function(path, values) {
  connection <- file(path, open = "ab")
  on.exit(close(connection))
  writeBin(as.double(t(as.matrix(values))), connection, size = 8L)
}

partition_layouts <- function(specs, work_dir, create_directories = TRUE) {
  layouts <- inventory_layouts(specs)
  for (slug in names(layouts)) {
    directory <- file.path(work_dir, "partitions", slug)
    if (create_directories) dir.create(directory, recursive = TRUE)
    layouts[[slug]]$directory <- directory
  }
  layouts
}

stream_inventory_partitions <- function(files, chunks, specs, config, work_dir, crs) {
  selectors <- vapply(files, read_selector, character(1), specs = specs)
  layouts <- partition_layouts(specs, work_dir)

  point_count <- 0
  hull <- NULL
  loaded_dimensions <- character()
  dropped_dimensions <- character()
  spatial_read_seconds <- 0
  partition_write_seconds <- 0
  for (chunk_index in seq_along(chunks)) {
    core <- chunks[[chunk_index]]
    read_started <- Sys.time()
    las <- read_spatial_chunk(files, selectors, core, specs, crs)
    spatial_read_seconds <- spatial_read_seconds + as.numeric(difftime(Sys.time(), read_started, units = "secs"))
    if (is.null(las)) next
    points <- las@data[
      X >= core[["xmin"]] & X < core[["xmax"]] &
        Y >= core[["ymin"]] & Y < core[["ymax"]]
    ]
    if (!nrow(points)) next
    point_count <- point_count + nrow(points)
    hull <- update_streaming_hull(hull, points$X, points$Y)
    loaded_dimensions <- union(loaded_dimensions, names(points))
    dropped_dimensions <- union(dropped_dimensions, attr(las, "dropped_dimensions"))

    write_started <- Sys.time()
    for (layout in layouts) {
      spec <- layout$spec
      instance <- spec$instance
      layout_columns <- layout$columns
      selected <- points[
        !is.na(get(instance)) & is.finite(get(instance)) &
          get(instance) >= 0 & !(get(instance) %in% config$non_tree_ids),
        ..layout_columns
      ]
      if (!nrow(selected)) next
      data.table::set(
        selected,
        j = "partition__",
        value = as.integer(abs(selected[[instance]]) %% config$inventory_partitions) + 1L
      )
      for (partition in unique(selected$partition__)) {
        path <- file.path(layout$directory, sprintf("%05d.bin", partition))
        partition_columns <- setdiff(names(selected), "partition__")
        append_partition(path, selected[partition__ == partition, ..partition_columns])
      }
    }
    partition_write_seconds <- partition_write_seconds + as.numeric(difftime(Sys.time(), write_started, units = "secs"))
    message(sprintf("CSP workflow: partitioned chunk %d/%d (%d points)",
                    chunk_index, length(chunks), nrow(points)))
    rm(las, points)
    gc(FALSE)
  }
  list(
    layouts = layouts,
    point_count = point_count,
    hull = hull,
    loaded_dimensions = loaded_dimensions,
    dropped_dimensions = dropped_dimensions,
    selectors = selectors,
    spatial_read_seconds = spatial_read_seconds,
    partition_write_seconds = partition_write_seconds
  )
}

stream_inventory_partitions_python <- function(files, chunks, specs, config, work_dir) {
  if (length(files) != 1L) abort("Sequential inventory requires exactly one merged COPC file")
  layouts <- partition_layouts(specs, work_dir, create_directories = FALSE)
  script <- file.path(repository_root, "workflow", "stream_inventory_partitions.py")
  hull_path <- file.path(work_dir, "stand_hull_xy.bin")
  config_path <- file.path(work_dir, "sequential_config.json")
  metadata_path <- file.path(work_dir, "sequential_metadata.json")
  log_path <- file.path(work_dir, "sequential.log")
  helper_config <- list(
    input = files[[1L]], output_dir = file.path(work_dir, "partitions"),
    hull_path = hull_path, partitions = config$inventory_partitions,
    chunk_points = 1000000L, non_tree_ids = as.numeric(config$non_tree_ids),
    spatial_grid = list(
      xmin = min(vapply(chunks, `[[`, numeric(1), "xmin")),
      ymin = min(vapply(chunks, `[[`, numeric(1), "ymin")),
      size = config$read_chunk_size,
      nx = length(unique(vapply(chunks, `[[`, numeric(1), "xmin"))),
      ny = length(unique(vapply(chunks, `[[`, numeric(1), "ymin")))
    ),
    layouts = lapply(names(layouts), function(slug) list(
      slug = slug, instance = layouts[[slug]]$spec$instance,
      columns = layouts[[slug]]$columns
    ))
  )
  jsonlite::write_json(helper_config, config_path, auto_unbox = TRUE, pretty = TRUE, na = "null")
  status <- system2("python3", c(script, "--config", config_path, "--metadata", metadata_path),
                    stdout = log_path, stderr = log_path)
  if (!identical(status, 0L) || !file.exists(metadata_path)) {
    detail <- if (file.exists(log_path)) paste(tail(readLines(log_path, warn = FALSE), 20L), collapse = "\n") else ""
    abort(sprintf("Sequential COPC inventory scan failed%s",
                  if (nzchar(detail)) paste0(":\n", detail) else ""))
  }
  metadata <- jsonlite::fromJSON(metadata_path, simplifyVector = TRUE)
  hull_started <- Sys.time()
  values <- readBin(hull_path, numeric(), n = file.info(hull_path)$size / 8,
                    size = 8L, endian = "little")
  if (length(values) %% 2L != 0L || !length(values)) abort("Invalid sequential hull vertices")
  xy <- matrix(values, ncol = 2L, byrow = TRUE)
  hull <- update_streaming_hull(NULL, xy[, 1L], xy[, 2L])
  hull_elapsed_seconds <- as.numeric(difftime(Sys.time(), hull_started, units = "secs"))
  unlink(hull_path)
  list(
    layouts = layouts, point_count = metadata$input_points, hull = hull,
    loaded_dimensions = metadata$loaded_dimensions,
    dropped_dimensions = metadata$dropped_dimensions,
    decompressed_groups = metadata$decompressed_groups,
    largest_instances = metadata$largest_instances,
    oversized_instances = metadata$oversized_instances,
    selectors = stats::setNames(metadata$selector, basename(files)),
    source_point_reads = metadata$source_point_reads,
    spatial_read_seconds = metadata$input_read_seconds,
    partition_write_seconds = metadata$partition_write_seconds,
    hull_elapsed_seconds = hull_elapsed_seconds + metadata$hull_seconds,
    xy_spool_bytes = metadata$xy_spool_bytes,
    spatial_bucket_bytes = metadata$spatial_bucket_bytes,
    sequential_elapsed_seconds = metadata$elapsed_seconds,
    sequential_peak_rss_kb = metadata$process_peak_rss_kb
  )
}

read_binary_partition <- function(path, columns) {
  count <- file.info(path)$size / 8
  values <- readBin(path, numeric(), n = count, size = 8L)
  if (length(values) %% length(columns) != 0L) abort(sprintf("Invalid inventory partition: %s", path))
  data.table::as.data.table(matrix(values, ncol = length(columns), byrow = TRUE, dimnames = list(NULL, columns)))
}

large_instance_counts <- function(stream, slug) {
  entries <- stream$oversized_instances[[slug]]
  if (is.null(entries) || !length(entries)) return(numeric())
  if (is.data.frame(entries)) {
    ids <- entries$instance_id
    counts <- entries$points
  } else {
    ids <- vapply(entries, `[[`, numeric(1), "instance_id")
    counts <- vapply(entries, `[[`, numeric(1), "points")
  }
  result <- stats::setNames(as.numeric(counts), as.character(ids))
  if (anyDuplicated(names(result)) || any(!is.finite(result)) ||
      any(result <= inventory_instance_point_limit)) {
    abort(sprintf("Invalid oversized instance audit for %s", slug))
  }
  result
}

spool_partition_instances <- function(path, columns, instance, dtm, directory,
                                      point_limit = inventory_instance_point_limit,
                                      batch_points = inventory_consumer_batch_points, large_counts = numeric(),
                                      spec = NULL) {
  row_bytes <- 8L * length(columns)
  bytes <- file.info(path)$size
  if (!is.finite(bytes) || bytes %% row_bytes != 0L) {
    abort(sprintf("Invalid inventory partition: %s", path))
  }
  dir.create(directory, recursive = TRUE)
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  instance_paths <- list()
  instance_counts <- list()
  large_states <- lapply(large_counts, new_large_instance_summary)
  normalized_columns <- c(columns, "Zref", "source_order__")
  rows_seen <- 0
  repeat {
    values <- readBin(connection, numeric(), n = as.integer(batch_points * length(columns)),
                      size = 8L, endian = "little")
    if (!length(values)) break
    if (length(values) %% length(columns) != 0L) abort(sprintf("Truncated inventory partition: %s", path))
    points <- data.table::as.data.table(matrix(
      values, ncol = length(columns), byrow = TRUE,
      dimnames = list(NULL, columns)
    ))
    data.table::set(points, j = "source_order__",
                    value = rows_seen + seq_len(nrow(points)))
    rows_seen <- rows_seen + nrow(points)
    points <- normalize_inventory_points(points, dtm)
    points <- points[
      is.finite(Z) & is.finite(get(instance)) & get(instance) >= 0
    ]
    if (!nrow(points)) next
    for (id in unique(points[[instance]])) {
      key <- as.character(id)
      selected <- points[get(instance) == id, ..normalized_columns]
      if (key %in% names(large_states)) {
        large_states[[key]] <- update_large_instance_summary(
          large_states[[key]], selected, spec)
        next
      }
      count <- (if (is.null(instance_counts[[key]])) 0L else instance_counts[[key]]) + nrow(selected)
      if (count > point_limit) {
        abort(sprintf(
          "Instance %s in %s has at least %d points, above the bounded inventory limit of %d; no partial output was written",
          key, basename(path), count, point_limit
        ))
      }
      if (is.null(instance_paths[[key]])) {
        instance_paths[[key]] <- file.path(directory, sprintf("%08d.bin", length(instance_paths) + 1L))
      }
      append_partition(instance_paths[[key]], selected)
      instance_counts[[key]] <- count
    }
  }
  list(paths = instance_paths, columns = normalized_columns,
       counts = instance_counts, large_states = large_states)
}

inventory_from_partitions <- function(stream, dtm, config, aoi, crs) {
  if (!is.null(stream$largest_instances)) {
    for (slug in names(stream$largest_instances)) {
      largest <- stream$largest_instances[[slug]]
      if (largest$points > inventory_instance_point_limit &&
          !as.character(largest$instance_id) %in%
          names(large_instance_counts(stream, slug))) {
        abort(sprintf(
          "Instance %s (%s) has %s points, above the bounded inventory limit of %d and lacks an oversized-instance audit; no inventory workers started",
          largest$instance_id, slug, format(largest$points, scientific = FALSE),
          inventory_instance_point_limit
        ))
      }
    }
  }
  inventories <- list()
  for (slug in names(stream$layouts)) {
    layout <- stream$layouts[[slug]]
    large_counts <- large_instance_counts(stream, slug)
    paths <- sort(list.files(layout$directory, pattern = "\\.bin$", full.names = TRUE))
    # Start the largest instances first so they do not become the final
    # single-worker tail. Retain stable seed indices from filename order.
    path_order <- order(file.info(paths)$size, decreasing = TRUE)
    seed_indices <- seq_along(paths)[path_order]
    paths <- paths[path_order]
    workers <- min(config$inventory_workers, length(paths))
    use_parallel <- workers > 1L
    if (use_parallel) {
      dtm_source <- terra::sources(dtm)
      if (length(dtm_source) != 1L || !file.exists(dtm_source)) {
        abort("Parallel inventory requires a file-backed DTM")
      }
    }
    process_partition <- function(index) {
      path <- paths[[index]]
      local_dtm <- dtm
      if (use_parallel) {
        # Parallel scheduling has no shared random stream. Keep its per-partition
        # stream repeatable; the one-worker default retains upstream seed order.
        set.seed(as.integer((as.double(config$random_seed) + seed_indices[[index]]) %%
                            .Machine$integer.max))
        # Each worker opens its own raster handle after the fork.
        local_dtm <- terra::rast(dtm_source)
        lidR::set_lidr_threads(1L)
      }
      spool_parent <- if (is.null(config$inventory_spool_dir)) dirname(path) else
        config$inventory_spool_dir
      dir.create(spool_parent, recursive = TRUE, showWarnings = FALSE)
      spool_dir <- tempfile("csp-instances-", tmpdir = spool_parent)
      on.exit(unlink(spool_dir, recursive = TRUE, force = TRUE), add = TRUE)
      spool <- spool_partition_instances(path, layout$columns, layout$spec$instance,
                                         local_dtm, spool_dir,
                                         large_counts = large_counts, spec = layout$spec)
      if (!length(spool$paths) && !length(spool$large_states)) {
        return(list(empty = TRUE, inventory = NULL))
      }
      groups <- list()
      current <- character()
      current_points <- 0L
      for (id in names(spool$paths)) {
        count <- spool$counts[[id]]
        if (length(current) && current_points + count > inventory_instance_point_limit) {
          groups[[length(groups) + 1L]] <- current
          current <- character()
          current_points <- 0L
        }
        current <- c(current, id)
        current_points <- current_points + count
      }
      if (length(current)) groups[[length(groups) + 1L]] <- current
      results <- lapply(groups, function(ids) {
        points <- data.table::rbindlist(lapply(spool$paths[ids], read_binary_partition,
                                               columns = spool$columns), use.names = TRUE)
        data.table::setorder(points, source_order__)
        points[, source_order__ := NULL]
        result <- measure_inventory_points(points, layout$spec, config, aoi, crs)
        rm(points)
        gc(FALSE)
        result
      })
      results <- c(results, lapply(names(spool$large_states), function(id) {
        state <- spool$large_states[[id]]
        if (!state$count) return(NULL)
        finish_large_instance_summary(state, id, layout$spec, aoi, crs)
      }))
      result <- data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
      gc(FALSE)
      list(empty = FALSE, inventory = result)
    }
    partials <- if (use_parallel) {
      parallel::mclapply(seq_along(paths), process_partition, mc.cores = workers,
                         mc.preschedule = FALSE, mc.set.seed = FALSE)
    } else {
      lapply(seq_along(paths), process_partition)
    }
    failed <- which(vapply(partials, function(result) is.null(result) || inherits(result, "try-error"), logical(1)))
    if (length(failed)) {
      failed_result <- partials[[failed[[1L]]]]
      detail <- if (inherits(failed_result, "try-error")) {
        paste0(": ", as.character(failed_result))
      } else {
        ": worker returned no result"
      }
      abort(paste0(sprintf("Inventory worker failed for %s partition %s", layout$spec$instance,
                           basename(paths[[failed[[1L]]]])), detail))
    }
    partials <- lapply(Filter(function(result) !result$empty, partials), `[[`, "inventory")
    if (!length(partials)) abort(sprintf("No valid inventory partitions for %s", layout$spec$instance))
    inventory <- data.table::rbindlist(partials, use.names = TRUE, fill = TRUE)
    data.table::setorder(inventory, instance_id)
    attr(inventory, "segmentation_source") <- layout$spec$source
    attr(inventory, "instance_dimension") <- layout$spec$instance
    inventories[[slug]] <- inventory
  }
  inventories
}
