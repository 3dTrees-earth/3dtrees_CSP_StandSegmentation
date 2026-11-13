# Canonical instance measurements. Readers supply normalized rows in source order.

# Points read per batch from a COPC pipe or partition file.
inventory_batch_points <- 1000000L

# DTM cell values, read once per file-backed DTM. terra::extract needs about
# 1 GB per million points; the cell index gives the same values for a fraction.
terrain_cache <- new.env()
terrain_values <- function(dtm) {
  key <- paste(terra::sources(dtm), collapse = "|")
  if (!nzchar(key)) return(terra::values(dtm, mat = FALSE))
  if (!identical(terrain_cache$key, key)) {
    terrain_cache$values <- terra::values(dtm, mat = FALSE)
    terrain_cache$key <- key
  }
  terrain_cache$values
}

normalize_inventory_points <- function(points, dtm) {
  terrain <- terrain_values(dtm)[terra::cellFromXY(dtm, cbind(points$X, points$Y))]
  if (any(!is.finite(terrain))) abort("DTM has missing or non-finite cells under inventoried points")
  data.table::set(points, j = "Zref", value = points$Z)
  data.table::set(points, j = "Z", value = points$Z - terrain)
  points
}

# Version 1: UTF-8 dimension + exact decimal ID, polynomial hash modulo a
# 31-bit prime. Double arithmetic stays exact (the multiplier is only 131).
# Identity, not traversal/window placement or worker assignment, owns the RNG.
instance_inventory_seed <- function(base_seed, dimension, id) {
  key <- paste0(dimension, ":", sprintf("%.0f", as.numeric(id)))
  seed <- as.double(base_seed)
  for (byte in as.integer(charToRaw(enc2utf8(key)))) seed <- (seed * 131 + byte) %% 2147483647
  as.integer(seed)
}

# The pipe stores float64 values; restore LAS standard-field storage types
# before lidR validation. Extra attributes remain numeric measurement inputs.
inventory_las <- function(points) {
  data <- as.data.frame(points)
  for (name in intersect(names(data), names(inventory_standard_attributes))) {
    storage <- inventory_standard_attributes[[name]]$r_storage
    data[[name]] <- switch(storage, integer = as.integer(data[[name]]),
                           logical = as.logical(data[[name]]), double = as.double(data[[name]]))
  }
  lidR::LAS(data)
}

in_dbh_band <- function(points, config) {
  is.finite(points$Z) & points$Z > config$slice_min & points$Z < config$slice_max
}

# Grid cell of every point on one world-aligned grid (absolute Z, so the cell
# does not depend on the terrain model or on tile and window borders).
grid_keys <- function(points, size) {
  z <- if ("Zref" %in% names(points)) points$Zref else points$Z
  data.table::data.table(kx = floor(points$X / size), ky = floor(points$Y / size), kz = floor(z / size))
}

# One point of a tree per grid cell: the smallest by X, Y, Z, so the choice
# does not depend on file, window or worker order. Coordinates are kept.
grid_points <- function(points, size) {
  if (!size) return(points)
  z <- if ("Zref" %in% names(points)) points$Zref else points$Z
  ord <- order(points$X, points$Y, z)
  points[ord][!duplicated(grid_keys(points, size)[ord])]
}

quality_label <- function(code, dbh) {
  labels <- c(
    `1` = "too_few_stem_slice_points",
    `2` = "fallback_circle_estimate",
    `3` = "rejected_or_unstable_fit",
    `4` = "spline_fit"
  )
  result <- unname(labels[as.character(code)])
  result[is.na(result)] <- "measurement_unavailable"
  result[is.na(dbh)] <- "dbh_unavailable"
  result
}

# Every distinct FM score with its point count, most points first, as a JSON
# array: an instance merged from several FM instances keeps each one's score.
score_values <- function(counts) {
  if (!length(counts)) return(NA_character_)
  counts <- counts[order(-counts, as.numeric(names(counts)))]
  sprintf("[%s]", paste(sprintf('{"score":%s,"point_count":%.0f}', names(counts), counts), collapse = ","))
}

filter_inventory_to_aoi <- function(inventory, aoi, crs) {
  if (is.null(aoi) || !nrow(inventory)) return(inventory)
  points <- sf::st_as_sf(inventory, coords = c("x", "y"), crs = crs, remove = FALSE)
  inventory[lengths(sf::st_intersects(points, aoi)) > 0L, ]
}

# One upstream forest_inventory() call; returns one row per tree or NULL.
# Upstream simplifies a one-tree result to an 8x1 frame keyed by row names.
# With DBH evidence on, the slice circles ride along as attribute "dbh_slices".
upstream_inventory <- function(las, instance, config) {
  measured <- tryCatch(
    CspStandSegmentation::forest_inventory(
      las,
      slice_min = config$slice_min,
      slice_max = config$slice_max,
      increment = config$slice_increment,
      width = config$slice_width,
      max_dbh = config$max_dbh,
      n_cores = config$geometry_threads,
      tree_id_col = instance,
      non_tree_id = if (is.null(config$non_tree_ids)) 0 else config$non_tree_ids,
      return_evidence = evidence_enabled(config)
    ),
    error = function(error) {
      warning(sprintf("Inventory measurements failed for %s: %s", instance, conditionMessage(error)))
      NULL
    }
  )
  if (is.null(measured) || !nrow(measured)) return(NULL)
  slices <- attr(measured, "dbh_slices")
  if (!instance %in% names(measured) && ncol(measured) == 1L) {
    values <- stats::setNames(as.numeric(measured[[1L]]), rownames(measured))
    measured <- data.table::as.data.table(as.list(values))
  }
  if (!instance %in% names(measured)) return(NULL)
  measured <- data.table::as.data.table(measured)
  data.table::setnames(measured, instance, "instance_id")
  if (!is.null(slices)) data.table::setattr(measured, "dbh_slices", slices)
  measured
}

# Upstream fits stem circles with RANSAC on a small planar/vertical subset of
# each slice, so one run depends on the random draw and on point order. The
# upstream algorithm is kept unchanged; reproducibility and stability come
# from around it: points are sorted, stem-band geometry is computed once,
# the fit runs config$dbh_runs times with seeds derived from the tree's
# identity, and the run with the median DBH is reported (its position and
# quality code too). Height and crown area do not depend on the draw.
# Trees this large (grid points) take half a minute and more per fit on one
# core. Their geometry and seeded fits also use cores that finished workers
# have handed back (config$core_tokens), so a few huge trees do not leave the
# other cores idle at the end of a run. No extra process starts while all
# workers are busy, so memory stays as without it. Each fit keeps its own seed.
large_tree_points <- 2e6

# Spare cores are files in a shared directory; renaming one claims it.
claim_cores <- function(directory, wanted) {
  if (is.null(directory) || wanted < 1L) return(character())
  claimed <- character()
  for (free in list.files(directory, pattern = "^free-", full.names = TRUE)) {
    taken <- sub("/free-", sprintf("/taken-%d-", Sys.getpid()), free, fixed = TRUE)
    if (file.rename(free, taken)) claimed <- c(claimed, taken)
    if (length(claimed) >= wanted) break
  }
  claimed
}

release_cores <- function(claimed) {
  for (taken in claimed) file.rename(taken, sub("/taken-[0-9]+-", "/free-", taken))
}

# A finished worker hands its core to workers still measuring large trees.
hand_back_core <- function(directory) {
  if (!is.null(directory)) file.create(file.path(directory, sprintf("free-%d", Sys.getpid())))
}

# Upstream fits the DBH to the stem band only; of the other points it uses the
# height range, the x/y extent and the convex hull above 0.5 m. Those points,
# in their original order, give the same measurements as the whole tree, and
# every seeded fit then holds a fraction of the tree instead of a copy of it.
# Upstream sees the coordinates lidR::LAS() quantizes (0.001, offsets from the
# raw minima), so the selection uses the same values and keeps the raw minima.
upstream_rows <- function(points, config) {
  quantize <- utils::getFromNamespace("fast_quantization", "lidR")
  quantized <- function(values, offset) {
    values <- values + 0
    quantize(values, 0.001, offset)
    values
  }
  x <- quantized(points$X, floor(min(points$X)))
  y <- quantized(points$Y, floor(min(points$Y)))
  z <- quantized(points$Z, 0)
  keep <- (z > config$slice_min & z < config$slice_max) |
    (points$Z > config$slice_min & points$Z < config$slice_max)
  keep[c(which.min(points$X), which.min(points$Y), which.min(points$Z),
         which.min(x), which.max(x), which.min(y), which.max(y), which.min(z), which.max(z))] <- TRUE
  crown <- which(z > 0.5)
  if (length(crown)) keep[crown[grDevices::chull(x[crown], y[crown])]] <- TRUE
  which(keep)
}

measure_tree <- function(points, instance, id, config) {
  data.table::setorderv(points, c("X", "Y", "Z"))
  tree_points <- nrow(points)
  measured <- points[upstream_rows(points, config)]
  band <- which(measured$Z > config$slice_min & measured$Z < config$slice_max)
  runs <- if (length(band) >= 3L) config$dbh_runs else 1L
  claimed <- if (tree_points >= large_tree_points) claim_cores(config$core_tokens, runs - 1L) else character()
  on.exit(release_cores(claimed), add = TRUE)
  cores <- 1L + length(claimed)
  if (length(band) > config$geometry_k) {
    geometry <- CspStandSegmentation::add_geometry(
      suppressMessages(suppressWarnings(lidR::LAS(as.data.frame(measured[band, .(X, Y, Z)])))),
      k = config$geometry_k, n_cores = max(config$geometry_threads, cores))@data
    for (feature in c("Planarity", "Linearity", "Verticality")) {
      values <- rep(NA_real_, nrow(measured))
      values[band] <- geometry[[feature]]
      data.table::set(measured, j = feature, value = values)
    }
  }
  las <- suppressMessages(suppressWarnings(lidR::LAS(as.data.frame(measured))))
  seed <- instance_inventory_seed(config$random_seed, instance, id)
  fit <- function(run) {
    set.seed((seed + run - 1L) %% .Machine$integer.max)
    suppressMessages(upstream_inventory(las, instance, config))
  }
  fits <- if (cores > 1L) {
    parallel::mclapply(seq_len(runs), fit, mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE)
  } else lapply(seq_len(runs), fit)
  if (any(vapply(fits, inherits, logical(1), "try-error"))) abort(sprintf("DBH fit process failed for %s %s", instance, id))
  slices <- lapply(fits, function(fit) if (is.null(fit)) NULL else attr(fit, "dbh_slices"))
  fits <- data.table::rbindlist(fits, use.names = TRUE, fill = TRUE)
  if (!nrow(fits)) {
    record_tree_evidence(config, instance, id, seed, fits, slices, 0L, points)
    return(NULL)
  }
  fitted <- which(is.finite(fits$DBH))
  chosen <- if (length(fitted) * 2L > nrow(fits)) {
    fitted[order(fits$DBH[fitted])][(length(fitted) + 1L) %/% 2L]
  } else {
    # Most runs failed: report the most frequent failure code, no DBH.
    codes <- fits$quality_flag
    if (all(is.na(codes))) 1L else which(codes == as.numeric(names(which.max(table(codes)))))[[1L]]
  }
  result <- fits[chosen]
  if (!length(fitted) * 2L > nrow(fits)) result[, DBH := NA_real_]
  result[, `:=`(Height = fits$Height[[1L]], ConvexHullArea = fits$ConvexHullArea[[1L]],
                dbh_runs = nrow(fits), dbh_runs_fitted = length(fitted),
                dbh_spread_m = if (length(fitted)) diff(range(fits$DBH[fitted])) else NA_real_)]
  record_tree_evidence(config, instance, id, seed, fits, slices, chosen, points)
  result
}

inventory_columns <- c(
  "segmentation_source", "instance_dimension", "instance_id", "x", "y", "z",
  "dbh_m", "height_m", "convex_hull_area_m2", "point_count", "metrics_point_count",
  "species_id", "species_prob", "species_conflict_fraction",
  "pred_score_values", "dbh_quality_code", "dbh_quality_label",
  "dbh_runs", "dbh_runs_fitted", "dbh_spread_m", "dbh_band_points"
)

# One tree's row from its --metrics-grid cloud (the instance column set):
# upstream's measurements where it produced them, otherwise the cloud's median
# position, lowest point, extent and convex hull.
measure_tree_geometry <- function(points, spec, id, config) {
  z_original <- if ("Zref" %in% names(points)) points$Zref else points$Z
  fit <- measure_tree(data.table::copy(points), spec$instance, id, config)
  if (is.null(fit)) {
    fit <- data.table::data.table(X = NA_real_, Y = NA_real_, Z = NA_real_, DBH = NA_real_,
                                  quality_flag = NA_real_, Height = NA_real_, ConvexHullArea = NA_real_,
                                  dbh_runs = NA_integer_, dbh_runs_fitted = NA_integer_, dbh_spread_m = NA_real_)
  }
  code <- data.table::fcoalesce(as.integer(fit$quality_flag), 9L)
  data.table::data.table(
    segmentation_source = spec$source, instance_dimension = spec$instance, instance_id = as.numeric(id),
    x = data.table::fcoalesce(as.numeric(fit$X), stats::median(points$X)),
    y = data.table::fcoalesce(as.numeric(fit$Y), stats::median(points$Y)),
    z = data.table::fcoalesce(as.numeric(fit$Z), min(z_original)),
    dbh_m = as.numeric(fit$DBH),
    height_m = data.table::fcoalesce(as.numeric(fit$Height), max(points$Z) - min(points$Z)),
    convex_hull_area_m2 = data.table::fcoalesce(as.numeric(fit$ConvexHullArea), polygon_area(points$X, points$Y)),
    metrics_point_count = nrow(points),
    dbh_quality_code = code, dbh_quality_label = quality_label(code, fit$DBH),
    dbh_runs = as.integer(fit$dbh_runs), dbh_runs_fitted = as.integer(fit$dbh_runs_fitted),
    dbh_spread_m = as.numeric(fit$dbh_spread_m),
    dbh_band_points = as.numeric(sum(in_dbh_band(points, config)))
  )
}

add_value_counts <- function(counts, values) {
  batch <- data.table::data.table(value = values)[is.finite(value), .(n = as.numeric(.N)), by = value]
  data.table::rbindlist(list(counts, batch))[, .(n = sum(n)), by = value]
}

median_of_counts <- function(counts) {
  if (is.null(counts) || !nrow(counts)) return(NA_real_)
  counts <- counts[order(value)]
  cumulative <- cumsum(counts$n)
  total <- cumulative[[nrow(counts)]]
  lower <- counts$value[[which(cumulative >= (total + 1) %/% 2)[[1L]]]]
  upper <- counts$value[[which(cumulative >= total %/% 2 + 1)[[1L]]]]
  (lower + upper) / 2
}


# Every tree is collected the same way, in one batch or many: its
# --metrics-grid cloud for the geometry, and exact counts of every input point
# for point_count, species, probability and the FM fields.
# Grid points are kept as the LAS integers (4 bytes per axis instead of 8 for
# X, Y, Z and Zref) when the run has one coordinate frame; heights above the
# terrain are recomputed from the same DTM lookup when the tree is measured.
tree_frame <- function(config, dtm) {
  if (is.null(config$coordinate_frame) || is.null(dtm)) return(NULL)
  c(config$coordinate_frame, list(dtm = dtm))
}

encode_grid <- function(points, frame) {
  if (is.null(frame)) return(points)
  z <- if ("Zref" %in% names(points)) points$Zref else points$Z
  axes <- list(X = points$X, Y = points$Y, Z = z)
  encoded <- Map(function(value, scale, offset) as.integer(round((value - offset) / scale)),
                 axes, frame$scale, frame$offset)
  for (axis in names(axes)) {
    if (!identical(as.numeric(encoded[[axis]]) * frame$scale[[axis]] + frame$offset[[axis]], axes[[axis]])) {
      abort(sprintf("%s coordinates are not on the point cloud's scale and offset grid", axis))
    }
  }
  data.table::as.data.table(encoded)
}

decode_grid <- function(grid, frame) {
  for (axis in c("X", "Y", "Z")) {
    data.table::set(grid, j = axis, value = as.numeric(grid[[axis]]) * frame$scale[[axis]] + frame$offset[[axis]])
  }
  grid
}

new_tree <- function(frame = NULL) {
  list(count = 0, chunks = list(), rows = 0, compacted = 0, frame = frame, species = NULL,
       probability = NULL, scores = NULL)
}

# All grid points collected so far, reduced to one per cell; with a frame
# the heights above the terrain are recomputed for the measurement.
compact_grid <- function(tree, size) {
  grid <- data.table::rbindlist(tree$chunks)
  if (is.null(tree$frame)) return(grid_points(grid, size))
  encode_grid(grid_points(decode_grid(grid, tree$frame), size), tree$frame)
}

tree_grid <- function(tree, size) {
  grid <- compact_grid(tree, size)
  if (is.null(tree$frame)) grid else normalize_inventory_points(decode_grid(grid, tree$frame), tree$frame$dtm)
}

add_tree_points <- function(tree, points, spec, config) {
  tree$count <- tree$count + nrow(points)
  coordinates <- intersect(c("X", "Y", "Z", "Zref"), names(points))
  reduced <- grid_points(points[, ..coordinates], config$metrics_grid)
  tree$chunks[[length(tree$chunks) + 1L]] <- encode_grid(reduced, tree$frame)
  tree$rows <- tree$rows + nrow(reduced)
  # Re-reduce once the cloud has doubled: memory stays near the grid cloud
  # and the total work stays proportional to the points read.
  if (config$metrics_grid && tree$rows > 2 * max(tree$compacted, 1e6)) {
    tree$chunks <- list(compact_grid(tree, config$metrics_grid))
    tree$rows <- tree$compacted <- nrow(tree$chunks[[1L]])
  }
  if (!is.null(spec$species)) {
    species <- points[[spec$species]]
    valid <- is.finite(species) & species >= 0
    tree$species <- add_value_counts(tree$species, species[valid])
    if (!is.null(spec$species_prob)) {
      probability <- points[[spec$species_prob]]
      counted <- valid & is.finite(probability) & probability >= 0
      batch <- data.table::data.table(species = species[counted], value = probability[counted])
      tree$probability <- data.table::rbindlist(list(tree$probability, batch[, .(n = as.numeric(.N)), by = .(species, value)]))[
        , .(n = sum(n)), by = .(species, value)]
    }
  }
  if (!is.null(spec$score)) {
    scores <- points[[spec$score]]
    tree$scores <- add_value_counts(tree$scores, scores[is.finite(scores) & scores >= 0])
  }
  tree
}

# The most frequent species (the smallest id on ties), the median probability
# of its points, and the share of points voting for another species.
species_vote <- function(tree) {
  if (is.null(tree$species) || !nrow(tree$species)) return(list(NA_real_, NA_real_, NA_real_))
  chosen <- min(tree$species[n == max(n), value])
  votes <- tree$species[value == chosen, n]
  probability <- if (is.null(tree$probability)) NULL else tree$probability[species == chosen, .(value, n)]
  list(chosen, median_of_counts(probability), 1 - votes / sum(tree$species$n))
}

# Every distinct prediction score of the tree with its point count (JSON).
tree_scores <- function(tree) {
  if (is.null(tree$scores) || !nrow(tree$scores)) return(NA_character_)
  score_values(stats::setNames(tree$scores$n, as.character(tree$scores$value)))
}

finish_tree <- function(tree, id, spec, config, aoi, crs) {
  grid <- tree_grid(tree, config$metrics_grid)
  data.table::set(grid, j = spec$instance, value = as.numeric(id))
  config$non_tree_ids <- spec$non_tree_ids
  row <- measure_tree_geometry(grid, spec, id, config)
  row[, point_count := tree$count]
  if (!is.null(spec$species)) {
    row[, c("species_id", "species_prob", "species_conflict_fraction") := species_vote(tree)]
  }
  if (!is.null(spec$score)) row[, pred_score_values := tree_scores(tree)]
  data.table::setcolorder(row, intersect(inventory_columns, names(row)))
  filter_inventory_to_aoi(row, aoi, crs)
}

# Rows of one segmentation, sorted by instance and labelled with it.
inventory_table <- function(rows, spec) {
  inventory <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (nrow(inventory)) data.table::setorder(inventory, instance_id)
  attr(inventory, "segmentation_source") <- spec$source
  attr(inventory, "instance_dimension") <- spec$instance
  inventory
}

# Inventory of trees whose points are all in memory (the CSP path and tests).
inventory_from_points <- function(points, spec, config, aoi, crs, dtm = NULL) {
  instance <- spec$instance
  points <- points[is.finite(get(instance)) & get(instance) >= 0 & !(get(instance) %in% spec$non_tree_ids)]
  if (!nrow(points)) abort(sprintf("No valid instances remain for requested dimension %s", instance))
  rows <- lapply(split(points, by = instance), function(tree) {
    finish_tree(add_tree_points(new_tree(tree_frame(config, dtm)), tree, spec, config),
                tree[[instance]][[1L]], spec, config, aoi, crs)
  })
  inventory_table(rows, spec)
}

inventory_layouts <- function(specs) {
  layouts <- list()
  for (spec in specs) {
    slug <- gsub("[^A-Za-z0-9_]+", "_", spec$source)
    columns <- unique(c("X", "Y", "Z", inventory_spec_dimensions(spec)))
    layouts[[slug]] <- list(spec = spec, columns = columns[!is.na(columns) & nzchar(columns)])
  }
  layouts
}
