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
          get(instance) >= 0 & !(get(instance) %in% spec$non_tree_ids),
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
  if (length(files) != 1L) abort("Sequential inventory requires exactly one merged LAS/LAZ file")
  layouts <- partition_layouts(specs, work_dir, create_directories = FALSE)
  script <- file.path(repository_root, "workflow", "stream_inventory_partitions.py")
  hull_path <- file.path(work_dir, "stand_hull_xy.bin")
  config_path <- file.path(work_dir, "sequential_config.json")
  metadata_path <- file.path(work_dir, "sequential_metadata.json")
  log_path <- file.path(work_dir, "sequential.log")
  helper_config <- list(
    input = files[[1L]], output_dir = file.path(work_dir, "partitions"),
    hull_path = hull_path, partitions = config$inventory_partitions,
    chunk_points = 1000000L,
    spatial_grid = list(
      xmin = min(vapply(chunks, `[[`, numeric(1), "xmin")),
      ymin = min(vapply(chunks, `[[`, numeric(1), "ymin")),
      size = config$read_chunk_size,
      nx = length(unique(vapply(chunks, `[[`, numeric(1), "xmin"))),
      ny = length(unique(vapply(chunks, `[[`, numeric(1), "ymin")))
    ),
    layouts = lapply(names(layouts), function(slug) list(
      slug = slug, instance = layouts[[slug]]$spec$instance,
      non_tree_ids = as.list(as.numeric(layouts[[slug]]$spec$non_tree_ids)),
      columns = layouts[[slug]]$columns
    ))
  )
  jsonlite::write_json(helper_config, config_path, auto_unbox = TRUE, pretty = TRUE, na = "null")
  status <- system2("python3", c(script, "--config", config_path, "--metadata", metadata_path),
                    stdout = log_path, stderr = log_path)
  if (!identical(status, 0L) || !file.exists(metadata_path)) {
    detail <- if (file.exists(log_path)) paste(tail(readLines(log_path, warn = FALSE), 20L), collapse = "\n") else ""
    abort(sprintf("Sequential inventory scan failed%s",
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

# Reads one instance partition in batches and collects every tree in it.
# A tree always lands in exactly one partition (ID modulo the partition count).
collect_partition_trees <- function(path, columns, dtm, spec, config,
                                    batch_points = inventory_batch_points) {
  instance <- spec$instance
  if (file.info(path)$size %% (8L * length(columns)) != 0L) abort(sprintf("Invalid inventory partition: %s", path))
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  trees <- list()
  repeat {
    values <- readBin(connection, numeric(), n = as.integer(batch_points * length(columns)),
                      size = 8L, endian = "little")
    if (!length(values)) break
    if (length(values) %% length(columns) != 0L) abort(sprintf("Truncated inventory partition: %s", path))
    points <- data.table::as.data.table(matrix(values, ncol = length(columns), byrow = TRUE,
                                               dimnames = list(NULL, columns)))
    points <- normalize_inventory_points(points, dtm)
    points <- points[is.finite(Z) & is.finite(get(instance)) & get(instance) >= 0 &
                       !(get(instance) %in% spec$non_tree_ids)]
    for (id in unique(points[[instance]])) {
      key <- sprintf("%.0f", id)
      tree <- if (is.null(trees[[key]])) new_tree(tree_frame(config, dtm)) else trees[[key]]
      trees[[key]] <- add_tree_points(tree, points[get(instance) == id], spec, config)
    }
  }
  trees
}

inventory_from_partitions <- function(stream, dtm, config, aoi, crs) {
  inventories <- list()
  for (slug in names(stream$layouts)) {
    layout <- stream$layouts[[slug]]
    # Largest partitions first so they do not become the single-worker tail.
    paths <- list.files(layout$directory, pattern = "\\.bin$", full.names = TRUE)
    paths <- paths[order(file.info(paths)$size, decreasing = TRUE)]
    workers <- min(config$inventory_workers, length(paths))
    dtm_source <- terra::sources(dtm)
    if (workers > 1L && (length(dtm_source) != 1L || !file.exists(dtm_source))) {
      abort("Parallel inventory requires a file-backed DTM")
    }
    # Every measurement seeds itself from the tree's identity, so the worker
    # count and partition order do not change any value.
    process_partition <- function(path) {
      local_dtm <- if (workers > 1L) terra::rast(dtm_source) else dtm
      if (workers > 1L) lidR::set_lidr_threads(1L)
      trees <- collect_partition_trees(path, layout$columns, local_dtm, layout$spec, config)
      rows <- lapply(names(trees), function(key) {
        finish_tree(trees[[key]], as.numeric(key), layout$spec, config, aoi, crs)
      })
      data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
    }
    partials <- if (workers > 1L) {
      parallel::mclapply(paths, process_partition, mc.cores = workers,
                         mc.preschedule = FALSE, mc.set.seed = FALSE)
    } else {
      lapply(paths, process_partition)
    }
    failed <- which(vapply(partials, function(result) is.null(result) || inherits(result, "try-error"), logical(1)))
    if (length(failed)) {
      detail <- if (inherits(partials[[failed[[1L]]]], "try-error")) as.character(partials[[failed[[1L]]]]) else
        "worker returned no result"
      abort(sprintf("Inventory worker failed for %s partition %s: %s", layout$spec$instance,
                    basename(paths[[failed[[1L]]]]), detail))
    }
    inventory <- inventory_table(partials, layout$spec)
    if (!nrow(inventory) && is.null(aoi)) abort(sprintf("No valid inventory partitions for %s", layout$spec$instance))
    inventories[[slug]] <- inventory
  }
  inventories
}
