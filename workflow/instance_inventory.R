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

measure_inventory_points <- function(points, spec, config, aoi, crs) {
  # Preserve the partition reader's existing grouped calculation and RNG order.
  # Split only when the aggregate geometry guard requires per-ID decisions.
  geometry_points <- function(x) sum(x$Z > config$slice_min & x$Z < config$slice_max, na.rm = TRUE)
  measure <- function(selected) {
    if (geometry_points(selected) <= inventory_dbh_geometry_limit) {
      return(inventory_for_spec(inventory_las(selected), spec, config, aoi, crs))
    }
    id <- as.character(selected[[spec$instance]][[1L]])
    state <- new_large_instance_summary(nrow(selected))
    for (start in seq.int(1L, nrow(selected), by = inventory_consumer_batch_points)) {
      end <- min(start + inventory_consumer_batch_points - 1L, nrow(selected))
      state <- update_large_instance_summary(state, selected[start:end], spec)
    }
    finish_large_instance_summary(state, id, spec, aoi, crs, method = "streaming_geometry_limit")
  }
  if (geometry_points(points) <= inventory_dbh_geometry_limit) return(measure(points))
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
    fallback_hull = polygon_area(X, Y)
  ), by = instance]
  data.table::setnames(base, instance, "instance_id")

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

  # Upstream simplifies a one-tree inventory to an 8x1 matrix/data frame whose
  # field names are row names. Restore the same one-row schema returned for
  # multi-tree calls before applying the normal merge path.
  if (!is.null(measured) && nrow(measured) && !instance %in% names(measured) && ncol(measured) == 1L) {
    fields <- rownames(measured)
    required_fields <- c(instance, "X", "Y", "Z", "DBH", "quality_flag", "Height", "ConvexHullArea")
    if (all(required_fields %in% fields)) {
      values <- as.numeric(measured[[1L]])
      names(values) <- fields
      measured <- data.table::as.data.table(as.list(values))
    }
  }

  if (is.null(measured) || !nrow(measured) || !instance %in% names(measured)) {
    if (!is.null(measured) && nrow(measured) && !instance %in% names(measured)) {
      warning(sprintf(
        "Inventory measurements for %s omitted the instance column; using fallback measurements for this partition",
        instance
      ))
    }
    measured <- data.table::data.table(
      instance_id = numeric(), X = numeric(), Y = numeric(), Z = numeric(),
      DBH = numeric(), quality_flag = integer(), Height = numeric(), ConvexHullArea = numeric()
    )
  } else {
    measured <- data.table::as.data.table(measured)
    data.table::setnames(measured, instance, "instance_id")
  }

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
    dbh_quality_code = data.table::fcoalesce(as.integer(quality_flag), 9L)
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
    score_mixed = FALSE
  )
}

update_large_instance_summary <- function(state, points, spec) {
  previous_count <- state$count
  state$count <- previous_count + nrow(points)
  state$min_zref <- min(state$min_zref, points$Zref)
  state$min_z <- min(state$min_z, points$Z)
  state$max_z <- max(state$max_z, points$Z)
  state$hull <- update_streaming_hull(state$hull, points$X, points$Y)

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

finish_large_instance_summary <- function(state, id, spec, aoi, crs,
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
  result <- data.table::data.table(
    segmentation_source = spec$source, instance_dimension = spec$instance,
    instance_id = as.numeric(id), x = stats::median(sample$X),
    y = stats::median(sample$Y), z = state$min_zref, dbh_m = NA_real_,
    height_m = state$max_z - state$min_z,
    convex_hull_area_m2 = polygon_area(state$hull$x, state$hull$y),
    point_count = state$count, species_id = species_id,
    species_prob = species_prob, species_conflict_fraction = conflict,
    pred_score_fm = score, pred_score_mixed = score_mixed,
    wood_point_count = wood_count, leaf_point_count = leaf_count,
    wood_share = wood_share, leaf_share = leaf_share,
    dbh_quality_code = 9L, dbh_quality_label = "dbh_unavailable",
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
