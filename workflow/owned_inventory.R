# COPC window adapter. Python owns index/query completeness; R owns DTM
# normalization and measurements. Point rows travel through a bounded pipe.

owned_reader_command <- function(config, arguments) {
  paste(c("env", sprintf("RAYON_NUM_THREADS=%d", config$inventory_decode_threads),
          "OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1", "python3",
          shQuote(file.path(repository_root, "workflow", "copc_owned_inventory.py")),
          vapply(arguments, shQuote, character(1))), collapse = " ")
}

index_owned_inventory <- function(files, specs, config, work_dir) {
  if (length(files) != 1L) abort("Owned COPC inventory requires exactly one merged COPC file")
  layouts <- inventory_layouts(specs)
  config_path <- file.path(work_dir, "owned_config.json")
  index_path <- file.path(work_dir, "owned_index.json")
  log_path <- file.path(work_dir, "owned_index.log")
  helper <- list(input = files[[1L]], batch_points = config$inventory_batch_points,
                 index_workers = config$inventory_index_workers,
                 max_window_points = window_point_cap(),
                 decode_workers = config$inventory_decode_threads,
                 index_queue_mb = config$inventory_index_queue_mb,
                 core_size = config$inventory_core_size, buffer = config$inventory_buffer,
                 prefetch_mb = config$inventory_prefetch_mb,
                 layouts = lapply(names(layouts), function(slug) list(
                   slug = slug, instance = layouts[[slug]]$spec$instance,
                   non_tree_ids = as.list(as.numeric(layouts[[slug]]$spec$non_tree_ids)),
                   columns = as.list(layouts[[slug]]$columns))))
  jsonlite::write_json(helper, config_path, auto_unbox = TRUE, digits = NA)
  command <- owned_reader_command(config, c("index", "--config", config_path, "--index", index_path))
  status <- system(paste(command, ">", shQuote(log_path), "2>&1"))
  if (status != 0L || !file.exists(index_path)) {
    abort(paste("COPC instance index failed:", paste(tail(readLines(log_path, warn = FALSE), 20L), collapse = "\n")))
  }
  index <- jsonlite::fromJSON(index_path, simplifyVector = FALSE)
  xy <- do.call(rbind, index$stand_hull)
  xy <- matrix(as.numeric(xy), ncol = 2L)
  hull <- update_streaming_hull(NULL, xy[, 1L], xy[, 2L])
  for (layout in index$layouts) {
    if (!length(layout$instances)) abort(sprintf("No valid instances for %s", layout$slug))
  }
  loaded <- unique(unlist(lapply(layouts, `[[`, "columns")))
  list(layouts = layouts, point_count = index$input_points, hull = hull,
       loaded_dimensions = loaded, dropped_dimensions = setdiff(inventory_canonical_names(unlist(index$available_dimensions)), loaded),
       decompressed_groups = unlist(index$decompressed_groups),
       selectors = "python_copc_owned", source_point_reads = NULL,
       owned_index = index, owned_config = config_path, owned_index_path = index_path,
       work_dir = work_dir)
}

read_owned_window <- function(stream, layout, plan, window_index, dtm, config, aoi, crs,
                              connection = NULL, audit_path = NULL) {
  spec <- layout$spec
  columns <- layout$columns
  ids <- sprintf("%.0f", as.numeric(unlist(plan$ids)))
  expected <- stats::setNames(as.numeric(unlist(plan$counts)), ids)
  observed <- stats::setNames(numeric(length(ids)), ids)
  fill <- new_terrain_fill(dtm, config)
  frame <- tree_frame(config, dtm, fill)
  trees <- stats::setNames(replicate(length(ids), new_tree(frame), simplify = FALSE), ids)
  if (is.null(audit_path)) audit_path <- file.path(stream$work_dir, sprintf("owned_%s_%06d.json", spec$source, window_index))
  log_path <- paste0(audit_path, ".log")
  command <- owned_reader_command(config, c(
    "read", "--config", stream$owned_config, "--index", stream$owned_index_path,
    "--slug", spec$source, "--window", as.character(window_index - 1L), "--audit", audit_path))
  if (is.null(connection)) connection <- pipe(paste(command, "2>", shQuote(log_path)), open = "rb")
  open <- TRUE
  on.exit(if (open) close(connection), add = TRUE)
  repeat {
    values <- readBin(connection, numeric(), n = inventory_batch_points * length(columns), size = 8L, endian = "little")
    if (!length(values)) break
    if (length(values) %% length(columns)) abort("Truncated owned COPC point stream")
    points <- data.table::as.data.table(matrix(values, ncol = length(columns), byrow = TRUE,
                                               dimnames = list(NULL, columns)))
    points <- normalize_inventory_points(points, dtm, config$dtm_lookup, fill)
    for (id in unique(points[[spec$instance]])) {
      key <- sprintf("%.0f", id)
      if (!key %in% ids) abort("COPC reader emitted an instance owned by another core")
      keep <- points[[spec$instance]] == id
      selected <- points[keep]
      observed[[key]] <- observed[[key]] + nrow(selected)
      if (observed[[key]] > expected[[key]]) abort(sprintf("Excess points for owned instance %s", key))
      trees[[key]] <- add_tree_points(trees[[key]], selected, spec, config)
    }
  }
  status <- close(connection)
  open <- FALSE
  if ((!is.null(status) && status != 0L) || !file.exists(audit_path)) {
    detail <- paste(tail(readLines(log_path, warn = FALSE), 20L), collapse = "\n")
    abort(sprintf("Owned COPC window failed: %s", detail))
  }
  audit <- jsonlite::fromJSON(audit_path, simplifyVector = FALSE)
  positions <- match(as.numeric(ids), as.numeric(unlist(audit$instance_ids)))
  if (!isTRUE(audit$complete) || anyNA(positions) || any(observed != expected) ||
      !identical(as.numeric(unlist(audit$observed_counts))[positions], unname(observed))) {
    abort("Incomplete owned instances; no inventory output published")
  }
  # Only after both reader and consumer agree may measurements be finalized;
  # largest first, so the longest trees start while other workers still run.
  results <- lapply(ids[order(-expected[ids])], function(key) {
    tree <- trees[[key]]
    trees[[key]] <<- NULL
    finish_tree(tree, as.numeric(key), spec, config, aoi, crs)
  })
  audit$r_observed_counts <- as.list(unname(observed))
  audit$r_peak_rss_kb <- process_peak_rss_kb()
  list(inventory = data.table::rbindlist(results, use.names = TRUE, fill = TRUE), audit = audit,
       instance_ids = ids, observed_counts = unname(observed), terrain_fills = terrain_fill_distances(fill))
}

# Largest instances first: balance point work, while never splitting one ID.
owned_worker_groups <- function(plan, workers) {
  ids <- as.numeric(unlist(plan$ids))
  counts <- as.numeric(unlist(plan$counts))
  workers <- min(workers, length(ids))
  groups <- replicate(workers, integer(), simplify = FALSE)
  loads <- numeric(workers)
  for (i in order(-counts, ids)) {
    worker <- which.min(loads)
    groups[[worker]] <- c(groups[[worker]], i)
    loads[[worker]] <- loads[[worker]] + counts[[i]]
  }
  lapply(groups, function(indices) {
    indices <- indices[order(ids[indices])]
    result <- plan
    result$ids <- as.list(ids[indices])
    result$counts <- as.list(counts[indices])
    result
  })
}

abort_owned_worker <- function(worker, detail = NULL, dispatcher = NULL) {
  failure <- if (is.null(detail)) {
    "exited without a result (it may have been killed, for example by an out-of-memory limit)"
  } else paste("reported an R error:", detail)
  diagnostic <- sprintf("Owned %s worker %s", worker, failure)
  suffix <- "; no inventory published"
  if (!grepl(suffix, diagnostic, fixed = TRUE)) diagnostic <- paste0(diagnostic, suffix)
  if (!is.null(dispatcher)) diagnostic <- paste(diagnostic, dispatcher, sep = "\n")
  abort(diagnostic)
}

check_owned_worker <- function(result, worker) {
  if (is.null(result)) abort_owned_worker(worker)
  if (inherits(result, "try-error")) {
    error <- attr(result, "condition")
    abort_owned_worker(worker, if (is.null(error)) as.character(result) else conditionMessage(error))
  }
}

read_parallel_owned_window <- function(stream, layout, plan, window_index, dtm_path, config, aoi, crs) {
  groups <- owned_worker_groups(plan, config$inventory_workers)
  if (length(groups) == 1L) {
    dtm <- terra::rast(dtm_path, guessCRS = FALSE)
    result <- read_owned_window(stream, layout, plan, window_index, dtm, config, aoi, crs)
    result$audit$worker_count <- 1L
    return(result)
  }
  pipe_dir <- tempfile("owned-pipes-", tmpdir = fifo_parent())
  if (!dir.create(pipe_dir)) abort("Cannot create owned worker transport directory")
  on.exit(unlink(pipe_dir, recursive = TRUE), add = TRUE)
  config$core_tokens <- file.path(pipe_dir, "cores")
  dir.create(config$core_tokens)
  paths <- file.path(pipe_dir, sprintf("worker-%02d.fifo", seq_along(groups)))
  error_paths <- file.path(pipe_dir, sprintf("worker-%02d.error", seq_along(groups)))
  if (system2("mkfifo", c("-m", "600", vapply(paths, shQuote, character(1)))) != 0L) {
    abort("Cannot create owned worker memory pipes")
  }
  audit_path <- file.path(stream$work_dir, sprintf("owned_%s_%06d.json", layout$spec$source, window_index))
  log_path <- paste0(audit_path, ".log")
  group_path <- file.path(pipe_dir, "groups.json")
  jsonlite::write_json(lapply(seq_along(groups), function(i) list(
    ids = groups[[i]]$ids, pipe = paths[[i]])), group_path, auto_unbox = TRUE, digits = NA)
  jobs <- list()
  collected <- FALSE
  on.exit({
    if (!collected && length(jobs)) {
      for (job in jobs) tools::pskill(job$pid, signal = 9L)
      parallel::mccollect(jobs)
    }
  }, add = TRUE)
  for (i in seq_along(groups)) {
    jobs[[i]] <- parallel::mcparallel(withCallingHandlers({
      lidR::set_lidr_threads(1L)
      data.table::setDTthreads(1L)
      # fifo() exposes short OS reads (often 8 KiB), which both fragments
      # records and makes terrain/summary updates run on tiny chunks. A raw
      # buffered file connection uses fread to fill readBin's requested batch.
      connection <- file(paths[[i]], open = "rb", raw = TRUE)
      result <- read_owned_window(stream, layout, groups[[i]], window_index,
                                  terra::rast(dtm_path, guessCRS = FALSE), config, aoi, crs,
                                  connection = connection, audit_path = audit_path)
      hand_back_core(config$core_tokens)
      result
    }, error = function(error) {
      # Save the R error before unwinding closes the FIFO and the dispatcher
      # can see a broken pipe. The worker still returns its normal try-error.
      try(writeLines(conditionMessage(error), error_paths[[i]]), silent = TRUE)
    }), mc.set.seed = FALSE)
  }
  command <- owned_reader_command(config, c(
    "dispatch", "--config", stream$owned_config, "--index", stream$owned_index_path,
    "--slug", layout$spec$source, "--window", as.character(window_index - 1L),
    "--groups", group_path, "--audit", audit_path))
  status <- system(paste(command, ">", shQuote(log_path), "2>&1"))
  if (status != 0L) {
    log_tail <- tail(readLines(log_path, warn = FALSE), 20L)
    reason <- sprintf("exit status %d", status)
    # The shell reports a signalled dispatcher as 128 plus its signal number.
    if (status == 137L) {
      reason <- paste0(reason, ", killed by signal 9 (possible out-of-memory kill)")
    } else if (status > 128L && status <= 192L) {
      reason <- sprintf("%s, terminated by signal %d", reason, status - 128L)
    }
    dispatcher <- paste(c(sprintf("Owned COPC dispatcher failed (%s):", reason), log_tail),
                        collapse = "\n")
    # A broken pipe points to a stopped consumer. Other dispatcher failures
    # cause worker EOF errors, which must not replace the dispatcher cause.
    if (any(grepl("^BrokenPipeError:", log_tail))) {
      diagnostics <- error_paths[file.exists(error_paths)]
      detail <- if (length(diagnostics)) {
        paste(readLines(diagnostics[[1L]], warn = FALSE), collapse = "\n")
      }
      abort_owned_worker("instance", detail, dispatcher)
    }
    abort(dispatcher)
  }
  results <- parallel::mccollect(jobs)
  collected <- TRUE
  for (result in results) check_owned_worker(result, "instance")
  if (length(results) != length(groups)) abort_owned_worker("instance")
  ids <- unlist(lapply(results, `[[`, "instance_ids"), use.names = FALSE)
  counts <- unlist(lapply(results, `[[`, "observed_counts"), use.names = FALSE)
  expected_ids <- sprintf("%.0f", as.numeric(unlist(plan$ids)))
  expected <- as.numeric(unlist(plan$counts))
  if (anyDuplicated(ids) || !setequal(ids, expected_ids) ||
      any(counts[match(expected_ids, ids)] != expected)) abort("Incomplete or duplicate parallel instance ownership")
  audit <- jsonlite::fromJSON(audit_path, simplifyVector = FALSE)
  audit$r_observed_counts <- as.list(counts[match(expected_ids, ids)])
  audit$workers <- lapply(results, function(result) list(
    instance_ids = as.list(as.numeric(result$instance_ids)),
    observed_counts = as.list(result$observed_counts),
    r_peak_rss_kb = result$audit$r_peak_rss_kb))
  list(inventory = data.table::rbindlist(lapply(results, `[[`, "inventory"), use.names = TRUE, fill = TRUE),
       audit = audit, terrain_fills = merge_terrain_fills(lapply(results, `[[`, "terrain_fills")))
}

# Bytes of memory per owned point held at once (grid trees, batches and
# measurement, including parallel fits of large trees): 1924, 1927 and 1887
# peaked at 14-28 GiB per billion points.
owned_bytes_per_point <- 32
owned_window_overhead_bytes <- 1024^3

# Memory the job may use: Galaxy's allocation, else 80 % of what is available.
memory_budget_bytes <- function() {
  galaxy <- suppressWarnings(as.numeric(Sys.getenv("GALAXY_MEMORY_MB", "")))
  if (is.finite(galaxy) && galaxy > 0) return(galaxy * 1024^2)
  info <- tryCatch(readLines("/proc/meminfo"), error = function(error) character())
  available <- suppressWarnings(as.numeric(sub("^MemAvailable:\\s+([0-9]+) kB$", "\\1",
                                                grep("^MemAvailable:", info, value = TRUE))))
  if (length(available) == 1L && is.finite(available)) 0.8 * available * 1024 else Inf
}

# Owned points one window may hold so that it fits 80 % of the memory budget
# on its own; larger windows are split into whole-tree windows by the index.
window_point_cap <- function() {
  budget <- memory_budget_bytes()
  if (!is.finite(budget)) return(0)
  max(1e6, floor((0.8 * budget - owned_window_overhead_bytes) / owned_bytes_per_point))
}

# Every window of every segmentation, largest first. A window that needs more
# than half of the memory runs alone with all workers; smaller ones run side
# by side with two workers each while their estimated memory fits.
owned_window_tasks <- function(stream, config) {
  tasks <- list()
  for (indexed in stream$owned_index$layouts) {
    for (i in seq_along(indexed$windows)) {
      points <- sum(as.numeric(unlist(indexed$windows[[i]]$counts)))
      tasks[[length(tasks) + 1L]] <- list(slug = indexed$slug, index = i, total = length(indexed$windows),
                                         window = indexed$windows[[i]],
                                         bytes = points * owned_bytes_per_point + owned_window_overhead_bytes)
    }
  }
  tasks[order(-vapply(tasks, `[[`, numeric(1), "bytes"))]
}

run_owned_window <- function(task, stream, dtm_path, config, aoi, crs, workers) {
  lidR::set_lidr_threads(1L)
  data.table::setDTthreads(1L)
  config$inventory_workers <- workers
  config$inventory_decode_threads <- min(config$inventory_decode_threads, max(1L, workers))
  message(sprintf("CSP workflow: owned %s window %d/%d (%d complete IDs, %d workers)",
                  task$slug, task$index, task$total, length(task$window$ids), workers))
  result <- read_parallel_owned_window(stream, stream$layouts[[task$slug]], task$window, task$index,
                                       dtm_path, config, aoi, crs)
  gc(FALSE)
  result
}

inventory_from_owned_windows <- function(stream, dtm, config, aoi, crs) {
  dtm_path <- terra::sources(dtm)
  if (length(dtm_path) != 1L || !file.exists(dtm_path)) abort("Owned inventory requires a file-backed DTM")
  tasks <- owned_window_tasks(stream, config)
  budget <- memory_budget_bytes()
  slots <- config$inventory_workers
  small_workers <- min(2L, slots)
  results <- vector("list", length(tasks))
  running <- list()
  # Workers still using the shared scratch directory must finish before the
  # caller removes it on error.
  on.exit(if (length(running)) parallel::mccollect(lapply(running, `[[`, "job")), add = TRUE)
  next_task <- 1L
  collect <- function(wait) {
    done <- parallel::mccollect(lapply(running, `[[`, "job"), wait = wait, timeout = if (wait) 0 else 5)
    for (pid in names(done)) {
      entry <- Filter(function(item) identical(as.character(item$job$pid), pid), running)[[1L]]
      running <<- Filter(function(item) !identical(as.character(item$job$pid), pid), running)
      check_owned_worker(done[[pid]], "window")
      results[[entry$task]] <<- done[[pid]]
    }
  }
  while (next_task <= length(tasks) || length(running)) {
    used_bytes <- sum(vapply(running, `[[`, numeric(1), "bytes"))
    used_slots <- sum(vapply(running, `[[`, numeric(1), "workers"))
    if (next_task <= length(tasks)) {
      task <- tasks[[next_task]]
      alone <- task$bytes > budget / 2 || slots <= small_workers
      workers <- if (alone) slots else small_workers
      fits <- if (alone) !length(running) else used_slots + workers <= slots && used_bytes + task$bytes <= budget
      if (fits) {
        if (alone) {
          results[[next_task]] <- run_owned_window(task, stream, dtm_path, config, aoi, crs, workers)
        } else {
          job <- parallel::mcparallel(run_owned_window(task, stream, dtm_path, config, aoi, crs, workers),
                                      mc.set.seed = FALSE)
          running[[length(running) + 1L]] <- list(job = job, task = next_task, bytes = task$bytes, workers = workers)
        }
        next_task <- next_task + 1L
        next
      }
    }
    collect(FALSE)
  }
  inventories <- audits <- list()
  for (indexed in stream$owned_index$layouts) {
    slug <- indexed$slug
    layout <- stream$layouts[[slug]]
    mine <- results[vapply(tasks, function(task) identical(task$slug, slug), logical(1))]
    inventory <- data.table::rbindlist(lapply(mine, `[[`, "inventory"), use.names = TRUE, fill = TRUE)
    if (anyDuplicated(inventory$instance_id)) abort("An owned instance was measured more than once")
    if (is.null(aoi) && nrow(inventory) != length(indexed$instances)) abort("Missing complete instance outputs")
    data.table::setorder(inventory, instance_id)
    attr(inventory, "segmentation_source") <- layout$spec$source
    attr(inventory, "instance_dimension") <- layout$spec$instance
    inventories[[slug]] <- inventory
    audits[[slug]] <- lapply(mine, `[[`, "audit")
  }
  attr(inventories, "terrain_fills") <- merge_terrain_fills(lapply(results, `[[`, "terrain_fills"))
  attr(inventories, "copc_audits") <- audits
  inventories
}
