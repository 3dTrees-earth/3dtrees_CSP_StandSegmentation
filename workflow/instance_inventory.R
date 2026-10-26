# Canonical instance measurements. Readers supply normalized rows in source order.
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

# Keep every stem-band point while the band fits the budget. Above it, keep
# the first point (source order) per voxel, starting at the policy's initial
# size and doubling until the band fits; once thinned, later batches use the
# same grid. Original coordinates are kept, nothing is averaged.
voxel_keys <- function(points, size) {
  z <- if ("Zref" %in% names(points)) points$Zref else points$Z
  data.table::data.table(kx = floor(points$X / size), ky = floor(points$Y / size), kz = floor(z / size))
}

dbh_band_keep <- function(points, budget, voxel = 0) {
  if (voxel <= 0 && nrow(points) <= budget) return(list(keep = rep(TRUE, nrow(points)), voxel = 0))
  size <- if (voxel > 0) voxel else inventory_dbh_band_initial_voxel_m
  repeat {
    keep <- !duplicated(voxel_keys(points, size))
    if (sum(keep) <= budget) return(list(keep = keep, voxel = size))
    size <- size * 2
  }
}

with_dbh_band <- function(result, band_points, source_points, voxel) {
  result[, `:=`(dbh_band_points = as.numeric(band_points),
                dbh_band_source_points = as.numeric(source_points),
                dbh_band_voxel_m = voxel)]
  result
}

measure_inventory_points <- function(points, spec, config, aoi, crs) {
  # Preserve the partition reader's existing grouped calculation and RNG order.
  # Split only when the aggregate stem band exceeds the per-tree budget.
  budget <- config$dbh_band_points
  measure <- function(selected) {
    band <- in_dbh_band(selected, config)
    if (sum(band) <= budget) return(inventory_for_spec(inventory_las(selected), spec, config, aoi, crs))
    thinned <- dbh_band_keep(selected[band], budget)
    keep <- rep(TRUE, nrow(selected))
    keep[which(band)[!thinned$keep]] <- FALSE
    result <- inventory_for_spec(inventory_las(selected[keep]), spec, config, aoi, crs)
    result[, point_count := nrow(selected)]
    with_dbh_band(result, sum(thinned$keep), sum(band), thinned$voxel)
  }
  if (sum(in_dbh_band(points, config)) <= budget) return(measure(points))
  groups <- split(points, by = spec$instance, keep.by = TRUE)
  data.table::rbindlist(lapply(groups, measure), use.names = TRUE, fill = TRUE)
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

modal_value <- function(values) {
  values <- values[!is.na(values) & is.finite(values) & values >= 0]
  if (!length(values)) return(NA_real_)
  counts <- table(values)
  as.numeric(sort(names(counts)[counts == max(counts)])[1L])
}

species_by_instance <- function(points, spec) {
  id <- spec$instance
  species <- spec$species
  probability <- spec$species_prob
  points[, {
    chosen <- modal_value(get(species))
    valid_species <- get(species)[!is.na(get(species)) & is.finite(get(species)) & get(species) >= 0]
    matching_prob <- get(probability)[get(species) == chosen & is.finite(get(probability)) & get(probability) >= 0]
    list(
      species_id = chosen,
      species_prob = if (length(matching_prob)) stats::median(matching_prob) else NA_real_,
      species_conflict_fraction = if (length(valid_species)) 1 - sum(valid_species == chosen) / length(valid_species) else NA_real_
    )
  }, by = id]
}

fm_by_instance <- function(points, instance) {
  has_score <- "PredScore_FM" %in% names(points)
  has_semantic <- "PredSemantic_FM" %in% names(points)
  points[, {
    scores <- if (has_score) get("PredScore_FM") else numeric()
    scores <- scores[is.finite(scores) & scores >= 0]
    semantics <- if (has_semantic) get("PredSemantic_FM") else numeric()
    wood <- if (has_semantic) sum(semantics == 1, na.rm = TRUE) else NA_integer_
    leaf <- if (has_semantic) sum(semantics == 2, na.rm = TRUE) else NA_integer_
    classified <- if (has_semantic) wood + leaf else 0L
    list(
      pred_score_fm = if (length(scores)) stats::median(scores) else NA_real_,
      pred_score_mixed = if (length(scores)) data.table::uniqueN(scores) > 1L else NA,
      wood_point_count = wood,
      leaf_point_count = leaf,
      wood_share = if (classified > 0L) wood / classified else NA_real_,
      leaf_share = if (classified > 0L) leaf / classified else NA_real_
    )
  }, by = instance]
}

filter_inventory_to_aoi <- function(inventory, aoi, crs) {
  if (is.null(aoi) || !nrow(inventory)) return(inventory)
  points <- sf::st_as_sf(inventory, coords = c("x", "y"), crs = crs, remove = FALSE)
  inventory[lengths(sf::st_intersects(points, aoi)) > 0L, ]
}

# One upstream forest_inventory() call; returns one row per tree or NULL.
# Upstream simplifies a one-tree result to an 8x1 frame keyed by row names.
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
      non_tree_id = config$non_tree_ids
    ),
    error = function(error) {
      warning(sprintf("Inventory measurements failed for %s: %s", instance, conditionMessage(error)))
      NULL
    }
  )
  if (is.null(measured) || !nrow(measured)) return(NULL)
  if (!instance %in% names(measured) && ncol(measured) == 1L) {
    values <- stats::setNames(as.numeric(measured[[1L]]), rownames(measured))
    measured <- data.table::as.data.table(as.list(values))
  }
  if (!instance %in% names(measured)) return(NULL)
  measured <- data.table::as.data.table(measured)
  data.table::setnames(measured, instance, "instance_id")
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
  fits <- data.table::rbindlist(fits, use.names = TRUE, fill = TRUE)
  if (!nrow(fits)) return(NULL)
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
  result
}

measure_trees <- function(las, instance, config) {
  columns <- intersect(c("X", "Y", "Z", "Zref", instance), names(las@data))
  trees <- split(las@data[!is.na(get(instance)) & !(get(instance) %in% config$non_tree_ids), ..columns],
                 by = instance, keep.by = TRUE)
  measured <- data.table::rbindlist(lapply(trees, function(points) {
    measure_tree(data.table::copy(points), instance, points[[instance]][[1L]], config)
  }), use.names = TRUE, fill = TRUE)
  if (!nrow(measured)) {
    measured <- data.table::data.table(
      instance_id = numeric(), X = numeric(), Y = numeric(), Z = numeric(),
      DBH = numeric(), quality_flag = integer(), Height = numeric(), ConvexHullArea = numeric(),
      dbh_runs = integer(), dbh_runs_fitted = integer(), dbh_spread_m = numeric())
  }
  measured
}

inventory_for_spec <- function(las, spec, config, aoi, crs) {
  instance <- spec$instance
  available <- names(las@data)
  point_columns <- unique(c(
    "X", "Y", "Z", intersect("Zref", available),
    instance, spec$species, spec$species_prob,
    if (identical(spec$source, "FM")) intersect(c("PredScore_FM", "PredSemantic_FM"), available)
  ))
  points <- las@data[
    !is.na(get(instance)) & is.finite(get(instance)) &
      get(instance) >= 0 & !(get(instance) %in% config$non_tree_ids),
    ..point_columns
  ]
  if (!nrow(points)) abort(sprintf("No valid instances remain for requested dimension %s", instance))

  z_original <- if ("Zref" %in% names(points)) "Zref" else "Z"
  base <- points[, list(
    point_count = .N,
    fallback_x = stats::median(X),
    fallback_y = stats::median(Y),
    fallback_z = min(get(z_original), na.rm = TRUE),
    fallback_height = max(Z, na.rm = TRUE) - min(Z, na.rm = TRUE),
    fallback_hull = polygon_area(X, Y),
    dbh_band_points = as.numeric(sum(Z > config$slice_min & Z < config$slice_max, na.rm = TRUE))
  ), by = instance]
  data.table::setnames(base, instance, "instance_id")

  measured <- measure_trees(las, instance, config)

  result <- merge(base, measured, by = "instance_id", all.x = TRUE, sort = FALSE)
  result[, `:=`(
    segmentation_source = spec$source,
    instance_dimension = instance,
    x = data.table::fcoalesce(X, fallback_x),
    y = data.table::fcoalesce(Y, fallback_y),
    z = data.table::fcoalesce(Z, fallback_z),
    dbh_m = DBH,
    height_m = data.table::fcoalesce(Height, fallback_height),
    convex_hull_area_m2 = data.table::fcoalesce(ConvexHullArea, fallback_hull),
    dbh_quality_code = data.table::fcoalesce(as.integer(quality_flag), 9L),
    dbh_band_source_points = dbh_band_points,
    dbh_band_voxel_m = 0
  )]
  result[, dbh_quality_label := quality_label(dbh_quality_code, dbh_m)]
  result[, inventory_method := "standard"]

  if (!is.null(spec$species)) {
    species <- species_by_instance(points, spec)
    data.table::setnames(species, instance, "instance_id")
    result <- merge(result, species, by = "instance_id", all.x = TRUE, sort = FALSE)
  }

  if (identical(instance, "PredInstance_FM") || identical(spec$source, "FM")) {
    fm <- fm_by_instance(points, instance)
    data.table::setnames(fm, instance, "instance_id")
    result <- merge(result, fm, by = "instance_id", all.x = TRUE, sort = FALSE)
  }

  keep <- c(
    "segmentation_source", "instance_dimension", "instance_id", "x", "y", "z",
    "dbh_m", "height_m", "convex_hull_area_m2", "point_count",
    "species_id", "species_prob", "species_conflict_fraction",
    "pred_score_fm", "pred_score_mixed", "wood_point_count", "leaf_point_count",
    "wood_share", "leaf_share", "dbh_quality_code", "dbh_quality_label",
    "dbh_runs", "dbh_runs_fitted", "dbh_spread_m",
    "dbh_band_points", "dbh_band_source_points", "dbh_band_voxel_m",
    "inventory_method"
  )
  result <- result[, intersect(keep, names(result)), with = FALSE]
  result <- filter_inventory_to_aoi(result, aoi, crs)
  attr(result, "segmentation_source") <- spec$source
  attr(result, "instance_dimension") <- instance
  result
}

update_streaming_hull <- function(hull, x, y) {
  candidate <- data.table::data.table(x = x, y = y)
  if (nrow(candidate) > 3L) candidate <- candidate[chull(x, y)]
  if (!is.null(hull)) candidate <- data.table::rbindlist(list(hull, candidate))
  if (nrow(candidate) > 3L) candidate <- candidate[chull(x, y)]
  candidate
}

new_large_instance_summary <- function(expected_points) {
  list(
    expected_points = as.numeric(expected_points), count = 0,
    min_zref = Inf, min_z = Inf, max_z = -Inf, hull = NULL,
    sample = NULL, species_counts = numeric(), valid_species_count = 0,
    wood_count = 0, leaf_count = 0, has_semantic = FALSE, score_first = NA_real_,
    score_mixed = FALSE, crown_hull = NULL,
    band = NULL, band_keys = NULL, band_source_points = 0, band_voxel = 0
  )
}

update_large_instance_summary <- function(state, points, spec, config) {
  previous_count <- state$count
  state$count <- previous_count + nrow(points)
  state$min_zref <- min(state$min_zref, points$Zref)
  state$min_z <- min(state$min_z, points$Z)
  state$max_z <- max(state$max_z, points$Z)
  state$hull <- update_streaming_hull(state$hull, points$X, points$Y)
  # Same crown definition as the full measurement: points above 0.5 m.
  crown <- points$Z > 0.5
  if (any(crown)) state$crown_hull <- update_streaming_hull(state$crown_hull, points$X[crown], points$Y[crown])

  # Stem band for DBH, in source order, bounded by the per-tree budget.
  band <- in_dbh_band(points, config)
  if (any(band)) {
    state$band_source_points <- state$band_source_points + sum(band)
    incoming <- points[band, .(X, Y, Z, Zref)]
    if (state$band_voxel > 0) {
      # Already thinned: drop incoming points whose voxel is taken by an
      # earlier point. Same result as re-thinning the whole band, without
      # rescanning the retained points on every batch.
      keys <- voxel_keys(incoming, state$band_voxel)
      keep <- !duplicated(keys)
      keep[keys[state$band_keys, on = c("kx", "ky", "kz"), which = TRUE, nomatch = 0L]] <- FALSE
      incoming <- incoming[keep]
      state$band_keys <- data.table::rbindlist(list(state$band_keys, keys[keep]))
    }
    state$band <- data.table::rbindlist(list(state$band, incoming))
    if (nrow(state$band) > config$dbh_band_points) {
      thinned <- dbh_band_keep(state$band, config$dbh_band_points, state$band_voxel)
      state$band <- state$band[thinned$keep]
      state$band_voxel <- thinned$voxel
      state$band_keys <- voxel_keys(state$band, state$band_voxel)
    }
  }

  # A fixed stride gives a deterministic, bounded sample across the full
  # source order. Only medians use it; counts, extrema and hull remain exact.
  stride <- max(1, ceiling(state$expected_points / inventory_large_instance_sample_limit))
  positions <- which(((previous_count + seq_len(nrow(points)) - 1) %% stride) == 0)
  if (length(positions)) {
    sample_columns <- intersect(c("X", "Y", spec$species, spec$species_prob,
                                  "PredScore_FM"), names(points))
    selected <- points[positions, ..sample_columns]
    state$sample <- if (is.null(state$sample)) selected else
      data.table::rbindlist(list(state$sample, selected), use.names = TRUE)
  }

  if (!is.null(spec$species) && spec$species %in% names(points)) {
    values <- points[[spec$species]]
    values <- values[is.finite(values) & values >= 0]
    state$valid_species_count <- state$valid_species_count + length(values)
    if (length(values)) {
      counts <- table(values)
      for (name in names(counts)) {
        if (!name %in% names(state$species_counts)) state$species_counts[name] <- 0
        state$species_counts[name] <- state$species_counts[name] + as.numeric(counts[[name]])
      }
    }
  }
  if (identical(spec$source, "FM")) {
    if ("PredSemantic_FM" %in% names(points)) {
      state$has_semantic <- TRUE
      state$wood_count <- state$wood_count + sum(points$PredSemantic_FM == 1, na.rm = TRUE)
      state$leaf_count <- state$leaf_count + sum(points$PredSemantic_FM == 2, na.rm = TRUE)
    }
    if ("PredScore_FM" %in% names(points)) {
      scores <- points$PredScore_FM
      scores <- scores[is.finite(scores) & scores >= 0]
      if (length(scores)) {
        if (is.na(state$score_first)) state$score_first <- scores[[1L]]
        state$score_mixed <- state$score_mixed || any(scores != state$score_first)
      }
    }
  }
  state
}

# DBH, its quality and the DBH-height position come from the retained stem
# band; counts, height, crown and species stay exact streaming summaries.
large_instance_dbh <- function(state, id, spec, config, crs) {
  if (is.null(state$band) || nrow(state$band) < 3L) return(NULL)
  band <- data.table::copy(state$band)
  data.table::set(band, j = spec$instance, value = as.numeric(id))
  band_spec <- list(instance = spec$instance, source = spec$source)
  measured <- inventory_for_spec(inventory_las(band), band_spec, config, NULL, crs)
  if (!nrow(measured) || measured$dbh_quality_code[[1L]] == 9L) return(NULL)
  measured[1L]
}

finish_large_instance_summary <- function(state, id, spec, aoi, crs, config,
                                          method = "streaming_large_instance") {
  if (state$count != state$expected_points || is.null(state$sample) ||
      !nrow(state$sample)) {
    abort(sprintf("Large instance %s has incomplete partition data: %.0f of %.0f points",
                  id, state$count, state$expected_points))
  }
  sample <- state$sample
  species_id <- NA_real_
  species_prob <- NA_real_
  conflict <- NA_real_
  if (length(state$species_counts)) {
    winners <- names(state$species_counts)[
      state$species_counts == max(state$species_counts)]
    species_id <- min(as.numeric(winners))
    conflict <- 1 - state$species_counts[[as.character(species_id)]] /
      state$valid_species_count
    if (!is.null(spec$species_prob) && spec$species_prob %in% names(sample)) {
      probabilities <- sample[[spec$species_prob]][
        sample[[spec$species]] == species_id]
      probabilities <- probabilities[is.finite(probabilities) & probabilities >= 0]
      if (length(probabilities)) species_prob <- stats::median(probabilities)
    }
  }
  score <- NA_real_
  score_mixed <- NA
  wood_count <- NA_real_
  leaf_count <- NA_real_
  wood_share <- NA_real_
  leaf_share <- NA_real_
  if (identical(spec$source, "FM")) {
    scores <- sample$PredScore_FM
    scores <- scores[is.finite(scores) & scores >= 0]
    if (length(scores)) score <- stats::median(scores)
    score_mixed <- if (is.na(state$score_first)) NA else state$score_mixed
    if (state$has_semantic) {
      wood_count <- state$wood_count
      leaf_count <- state$leaf_count
      classified <- wood_count + leaf_count
      if (classified > 0) {
        wood_share <- wood_count / classified
        leaf_share <- leaf_count / classified
      }
    }
  }
  dbh <- large_instance_dbh(state, id, spec, config, crs)
  # The stem position replaces the sample median only when a DBH was fitted.
  located <- !is.null(dbh) && is.finite(dbh$dbh_m)
  crown <- if (is.null(state$crown_hull)) state$hull else state$crown_hull
  result <- data.table::data.table(
    segmentation_source = spec$source, instance_dimension = spec$instance,
    instance_id = as.numeric(id),
    x = if (located) dbh$x else stats::median(sample$X),
    y = if (located) dbh$y else stats::median(sample$Y),
    z = if (located) dbh$z else state$min_zref,
    dbh_m = if (located) dbh$dbh_m else NA_real_,
    height_m = state$max_z - state$min_z,
    convex_hull_area_m2 = polygon_area(crown$x, crown$y),
    point_count = state$count, species_id = species_id,
    species_prob = species_prob, species_conflict_fraction = conflict,
    pred_score_fm = score, pred_score_mixed = score_mixed,
    wood_point_count = wood_count, leaf_point_count = leaf_count,
    wood_share = wood_share, leaf_share = leaf_share,
    dbh_quality_code = if (is.null(dbh)) 9L else dbh$dbh_quality_code,
    dbh_quality_label = if (is.null(dbh)) "dbh_unavailable" else dbh$dbh_quality_label,
    dbh_runs = if (is.null(dbh)) NA_integer_ else dbh$dbh_runs,
    dbh_runs_fitted = if (is.null(dbh)) NA_integer_ else dbh$dbh_runs_fitted,
    dbh_spread_m = if (is.null(dbh)) NA_real_ else dbh$dbh_spread_m,
    dbh_band_points = if (is.null(state$band)) 0 else as.numeric(nrow(state$band)),
    dbh_band_source_points = state$band_source_points,
    dbh_band_voxel_m = state$band_voxel,
    inventory_method = method
  )
  filter_inventory_to_aoi(result, aoi, crs)
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
