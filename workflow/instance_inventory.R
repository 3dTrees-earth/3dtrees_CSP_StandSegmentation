# Canonical instance measurements. Readers supply normalized rows in source order.

# Points read per batch from a COPC pipe or partition file.
inventory_batch_points <- 1000000L

normalize_inventory_points <- function(points, dtm) {
  terrain <- terra::extract(dtm, data.frame(X = points$X, Y = points$Y), ID = FALSE)[[1L]]
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
fm_score_values <- function(counts) {
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
      non_tree_id = config$non_tree_ids,
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
measure_tree <- function(points, instance, id, config) {
  data.table::setorderv(points, c("X", "Y", "Z"))
  band <- which(points$Z > config$slice_min & points$Z < config$slice_max)
  runs <- if (length(band) >= 3L) config$dbh_runs else 1L
  if (length(band) > config$geometry_k) {
    geometry <- CspStandSegmentation::add_geometry(
      suppressMessages(suppressWarnings(lidR::LAS(as.data.frame(points[band, .(X, Y, Z)])))),
      k = config$geometry_k, n_cores = config$geometry_threads)@data
    for (feature in c("Planarity", "Linearity", "Verticality")) {
      values <- rep(NA_real_, nrow(points))
      values[band] <- geometry[[feature]]
      data.table::set(points, j = feature, value = values)
    }
  }
  las <- suppressMessages(suppressWarnings(lidR::LAS(as.data.frame(points))))
  seed <- instance_inventory_seed(config$random_seed, instance, id)
  fits <- lapply(seq_len(runs), function(run) {
    set.seed((seed + run - 1L) %% .Machine$integer.max)
    suppressMessages(upstream_inventory(las, instance, config))
  })
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
  "pred_score_fm_values", "wood_point_count", "leaf_point_count",
  "wood_share", "leaf_share", "dbh_quality_code", "dbh_quality_label",
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

is_fm_spec <- function(spec) identical(spec$instance, "PredInstance_FM") || identical(spec$source, "FM")

# Every tree is collected the same way, in one batch or many: its
# --metrics-grid cloud for the geometry, and exact counts of every input point
# for point_count, species, probability and the FM fields.
new_tree <- function() {
  list(count = 0, grid = NULL, compacted = 0, species = NULL, probability = NULL,
       scores = NULL, semantic = FALSE, wood = 0, leaf = 0)
}

add_tree_points <- function(tree, points, spec, config) {
  tree$count <- tree$count + nrow(points)
  coordinates <- intersect(c("X", "Y", "Z", "Zref"), names(points))
  tree$grid <- data.table::rbindlist(list(tree$grid, grid_points(points[, ..coordinates], config$metrics_grid)))
  # Re-reduce once the cloud has doubled: memory stays near the grid cloud
  # and the total work stays proportional to the points read.
  if (config$metrics_grid && nrow(tree$grid) > 2 * max(tree$compacted, 1e6)) {
    tree$grid <- grid_points(tree$grid, config$metrics_grid)
    tree$compacted <- nrow(tree$grid)
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
  if (is_fm_spec(spec) && "PredScore_FM" %in% names(points)) {
    scores <- points$PredScore_FM
    tree$scores <- add_value_counts(tree$scores, scores[scores >= 0])
  }
  if (is_fm_spec(spec) && "PredSemantic_FM" %in% names(points)) {
    tree$semantic <- TRUE
    tree$wood <- tree$wood + sum(points$PredSemantic_FM == 1, na.rm = TRUE)
    tree$leaf <- tree$leaf + sum(points$PredSemantic_FM == 2, na.rm = TRUE)
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

fm_fields <- function(tree) {
  scores <- if (is.null(tree$scores) || !nrow(tree$scores)) numeric() else
    stats::setNames(tree$scores$n, as.character(tree$scores$value))
  classified <- tree$wood + tree$leaf
  list(fm_score_values(scores),
       if (tree$semantic) tree$wood else NA_real_, if (tree$semantic) tree$leaf else NA_real_,
       if (tree$semantic && classified > 0) tree$wood / classified else NA_real_,
       if (tree$semantic && classified > 0) tree$leaf / classified else NA_real_)
}

finish_tree <- function(tree, id, spec, config, aoi, crs) {
  grid <- grid_points(tree$grid, config$metrics_grid)
  data.table::set(grid, j = spec$instance, value = as.numeric(id))
  row <- measure_tree_geometry(grid, spec, id, config)
  row[, point_count := tree$count]
  if (!is.null(spec$species)) {
    row[, c("species_id", "species_prob", "species_conflict_fraction") := species_vote(tree)]
  }
  if (is_fm_spec(spec)) {
    row[, c("pred_score_fm_values", "wood_point_count", "leaf_point_count", "wood_share", "leaf_share") :=
          fm_fields(tree)]
  }
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
inventory_from_points <- function(points, spec, config, aoi, crs) {
  instance <- spec$instance
  points <- points[is.finite(get(instance)) & get(instance) >= 0 & !(get(instance) %in% config$non_tree_ids)]
  if (!nrow(points)) abort(sprintf("No valid instances remain for requested dimension %s", instance))
  rows <- lapply(split(points, by = instance), function(tree) {
    finish_tree(add_tree_points(new_tree(), tree, spec, config), tree[[instance]][[1L]], spec, config, aoi, crs)
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
