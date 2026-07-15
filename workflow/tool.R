resolve_input_files <- function(path) {
  if (file.exists(path) && !dir.exists(path)) return(normalizePath(path, mustWork = TRUE))
  if (!dir.exists(path)) abort(sprintf("Input does not exist: %s", path))
  files <- list.files(path, pattern = "\\.(las|laz)$", full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) abort(sprintf("Input directory contains no LAS/LAZ files: %s", path))
  sort(normalizePath(files, mustWork = TRUE))
}

inventory_read_dimensions <- function(specs) {
  dimensions <- unlist(lapply(specs, function(spec) {
    result <- c(spec$instance, spec$species, spec$species_prob)
    if (identical(spec$source, "FM")) {
      result <- c(result, "PredScore_FM", "PredSemantic_FM")
    }
    result
  }), use.names = FALSE)
  unique(dimensions[!is.na(dimensions) & nzchar(dimensions)])
}

extra_byte_names <- function(header) {
  extra_bytes <- header@VLR$Extra_Bytes[["Extra Bytes Description"]]
  if (is.null(extra_bytes)) character() else names(extra_bytes)
}

build_read_selector <- function(dimensions, extra_dimensions) {
  standard_codes <- c(
    gpstime = "t", ScanAngle = "a", Intensity = "i",
    NumberOfReturns = "n", ReturnNumber = "r", Classification = "c",
    Synthetic_flag = "s", Keypoint_flag = "k", Withheld_flag = "w",
    Overlap_flag = "o", UserData = "u", PointSourceID = "p",
    EdgeOfFlightline = "e", ScanDirectionFlag = "d", R = "R", G = "G",
    B = "B", NIR = "N", ScannerChannel = "C", Waveform = "W"
  )
  standard <- unname(standard_codes[intersect(dimensions, names(standard_codes))])
  positions <- sort(match(intersect(dimensions, extra_dimensions), extra_dimensions))

  # LASlib can address only the first nine extra-byte records individually.
  # If a requested dimension is later, load all extra bytes rather than risk
  # silently omitting a requested inventory field.
  extra <- if (any(positions > 9L)) "0" else as.character(positions)
  paste0(unique(c("c", standard, extra)), collapse = "")
}

read_selector <- function(file, specs, preserve_all_dimensions) {
  if (preserve_all_dimensions) return("*")
  header <- lidR::readLASheader(file)
  build_read_selector(inventory_read_dimensions(specs), extra_byte_names(header))
}

project_inventory_dimensions <- function(las, specs) {
  available <- names(las@data)
  keep <- unique(c(
    "X", "Y", "Z", intersect("Classification", available),
    intersect(inventory_read_dimensions(specs), available)
  ))
  dropped <- setdiff(available, keep)
  las@data <- las@data[, ..keep]

  extra_bytes <- las@header@VLR$Extra_Bytes[["Extra Bytes Description"]]
  if (!is.null(extra_bytes)) {
    for (dimension in setdiff(names(extra_bytes), keep)) {
      las@header@VLR$Extra_Bytes[["Extra Bytes Description"]][[dimension]] <- NULL
    }
  }
  attr(las, "dropped_dimensions") <- dropped
  las
}

read_point_cloud <- function(files, specs, preserve_all_dimensions = FALSE) {
  selectors <- vapply(
    files,
    read_selector,
    character(1),
    specs = specs,
    preserve_all_dimensions = preserve_all_dimensions
  )
  clouds <- Map(function(file, select) {
    cloud <- suppressWarnings(lidR::readLAS(file, select = select))
    if (is.null(cloud) || lidR::is.empty(cloud)) return(cloud)
    if (!preserve_all_dimensions) cloud <- project_inventory_dimensions(cloud, specs)
    cloud
  }, files, selectors)
  if (any(vapply(clouds, function(cloud) is.null(cloud) || lidR::is.empty(cloud), logical(1)))) {
    abort("At least one input point-cloud file is empty")
  }
  dropped_dimensions <- unique(unlist(lapply(clouds, attr, which = "dropped_dimensions")))
  if (is.null(dropped_dimensions)) dropped_dimensions <- character()
  if (length(files) == 1L) {
    las <- clouds[[1L]]
  } else {
    las <- CspStandSegmentation::las_merge(clouds, fill = TRUE)
  }
  if (is.null(las) || lidR::is.empty(las)) abort("The input point cloud is empty")
  attr(las, "read_selectors") <- stats::setNames(unname(selectors), basename(files))
  attr(las, "dropped_dimensions") <- dropped_dimensions
  las
}

has_valid_ground <- function(las) {
  if (!"Classification" %in% names(las@data)) return(FALSE)
  ground <- las@data[Classification == 2 & is.finite(X) & is.finite(Y)]
  nrow(ground) >= 3L && length(unique(ground$X)) >= 2L && length(unique(ground$Y)) >= 2L
}

make_dtm <- function(las, resolution) {
  if (has_valid_ground(las)) {
    grounded <- las
    method <- "classification_2"
  } else {
    grounded <- lidR::classify_ground(las, lidR::csf(), last_returns = FALSE)
    if (!has_valid_ground(grounded)) abort("Ground classification did not produce enough ground points for a DTM")
    method <- "csf"
  }
  dtm <- lidR::rasterize_terrain(grounded, res = resolution, algorithm = lidR::tin())
  if (is.null(dtm) || all(is.na(terra::values(dtm)))) abort("DTM generation produced no valid cells")
  list(las = grounded, dtm = dtm, method = method)
}

read_aoi <- function(path, las_crs) {
  if (is.null(path)) return(NULL)
  if (!file.exists(path)) abort(sprintf("AOI JSON does not exist: %s", path))
  aoi <- suppressWarnings(sf::st_read(path, quiet = TRUE))
  if (!nrow(aoi)) abort("AOI JSON contains no features")
  types <- unique(as.character(sf::st_geometry_type(aoi)))
  if (!all(types %in% c("POLYGON", "MULTIPOLYGON"))) {
    abort("AOI JSON must contain only Polygon or MultiPolygon geometries")
  }
  # AOI coordinates are defined by the CLI contract to already be in the
  # point-cloud CRS. GeoJSON readers otherwise label them as WGS84 even when
  # they contain projected native coordinates, so assign rather than transform.
  sf::st_crs(aoi) <- las_crs
  aoi <- sf::st_make_valid(aoi)
  union <- sf::st_union(sf::st_geometry(aoi))
  if (length(union) != 1L || sf::st_is_empty(union)) abort("AOI geometry is empty after normalization")
  union
}

polygon_area <- function(x, y) {
  if (length(x) < 3L) return(NA_real_)
  indices <- chull(x, y)
  xx <- x[indices]
  yy <- y[indices]
  0.5 * abs(sum(xx * c(yy[-1L], yy[1L]) - yy * c(xx[-1L], xx[1L])))
}

determine_area <- function(las, aoi) {
  if (!is.null(aoi)) {
    area <- as.numeric(sf::st_area(aoi))
    return(list(area_m2 = area, area_source = "aoi"))
  }
  list(
    area_m2 = polygon_area(las@data$X, las@data$Y),
    area_source = "point_cloud_convex_hull"
  )
}

validate_specs <- function(las, specs) {
  available <- names(las@data)
  for (spec in specs) {
    required <- c(spec$instance, spec$species, spec$species_prob)
    missing <- setdiff(required[!vapply(required, is.null, logical(1))], available)
    if (length(missing)) {
      abort(sprintf("Requested dimensions are missing for %s: %s", spec$instance, paste(missing, collapse = ", ")))
    }
  }
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
    inventory_quality_code = data.table::fcoalesce(as.integer(quality_flag), 9L)
  )]
  result[, inventory_quality_label := quality_label(inventory_quality_code, dbh_m)]

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
    "wood_share", "leaf_share", "inventory_quality_code", "inventory_quality_label"
  )
  result <- result[, intersect(keep, names(result)), with = FALSE]
  result <- filter_inventory_to_aoi(result, aoi, crs)
  attr(result, "segmentation_source") <- spec$source
  attr(result, "instance_dimension") <- instance
  result
}

safe_stat <- function(values, function_) {
  values <- values[is.finite(values)]
  if (length(values)) function_(values) else NA_real_
}

safe_quantile <- function(values, probability) {
  safe_stat(values, function(x) as.numeric(stats::quantile(x, probability, names = FALSE)))
}

stand_summary_for_inventory <- function(inventory, area) {
  area_ha <- area$area_m2 / 10000
  source <- if (nrow(inventory)) unique(inventory$segmentation_source) else attr(inventory, "segmentation_source")
  dimension <- if (nrow(inventory)) unique(inventory$instance_dimension) else attr(inventory, "instance_dimension")
  data.table::data.table(
    segmentation_source = source,
    instance_dimension = dimension,
    area_source = area$area_source,
    area_m2 = area$area_m2,
    tree_count = nrow(inventory),
    trees_per_ha = if (is.finite(area_ha) && area_ha > 0) nrow(inventory) / area_ha else NA_real_,
    basal_area_m2_per_ha = if (is.finite(area_ha) && area_ha > 0 && any(is.finite(inventory$dbh_m))) sum(pi * (inventory$dbh_m / 2)^2, na.rm = TRUE) / area_ha else NA_real_,
    mean_dbh_m = safe_stat(inventory$dbh_m, mean),
    median_dbh_m = safe_stat(inventory$dbh_m, stats::median),
    dbh_p05_m = safe_quantile(inventory$dbh_m, 0.05),
    dbh_p25_m = safe_quantile(inventory$dbh_m, 0.25),
    dbh_p75_m = safe_quantile(inventory$dbh_m, 0.75),
    dbh_p95_m = safe_quantile(inventory$dbh_m, 0.95),
    mean_height_m = safe_stat(inventory$height_m, mean),
    median_height_m = safe_stat(inventory$height_m, stats::median),
    height_p05_m = safe_quantile(inventory$height_m, 0.05),
    height_p25_m = safe_quantile(inventory$height_m, 0.25),
    height_p75_m = safe_quantile(inventory$height_m, 0.75),
    height_p95_m = safe_quantile(inventory$height_m, 0.95)
  )
}

species_composition <- function(inventories) {
  eligible <- Filter(function(inventory) "species_id" %in% names(inventory), inventories)
  if (!length(eligible)) return(NULL)
  combined <- data.table::rbindlist(eligible, fill = TRUE)
  combined <- combined[!is.na(species_id)]
  if (!nrow(combined)) {
    return(data.table::data.table(
      segmentation_source = character(),
      instance_dimension = character(),
      species_id = numeric(),
      tree_count = integer(),
      tree_proportion = numeric()
    ))
  }
  result <- combined[, .(tree_count = .N), by = .(segmentation_source, instance_dimension, species_id)]
  result[, tree_proportion := tree_count / sum(tree_count), by = .(segmentation_source, instance_dimension)]
  result[]
}

read_supplied_seeds <- function(path) {
  if (!file.exists(path)) abort(sprintf("Seed file does not exist: %s", path))
  seeds <- data.table::fread(path)
  required <- c("X", "Y", "Z", "TreeID")
  missing <- setdiff(required, names(seeds))
  if (length(missing)) abort(sprintf("Seed file is missing columns: %s", paste(missing, collapse = ", ")))
  seeds <- as.data.frame(seeds[, ..required])
  if (!nrow(seeds) || anyNA(seeds) || !all(vapply(seeds, is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(seeds)))) {
    abort("Seed columns must contain finite, non-missing numeric values")
  }
  if (any(seeds$TreeID <= 0) || any(seeds$TreeID != as.integer(seeds$TreeID)) || anyDuplicated(seeds$TreeID)) {
    abort("Seed TreeID values must be unique positive integers")
  }
  seeds
}

run_csp <- function(las, config, stage_dir) {
  seeds <- if (config$seed_mode == "supplied") {
    read_supplied_seeds(config$seed_file)
  } else {
    CspStandSegmentation::find_base_coordinates_raster(
      las,
      res = config$seed_resolution,
      zmin = config$seed_z_min,
      zmax = config$seed_z_max,
      q = config$seed_quantile,
      eps = config$seed_eps
    )
  }
  if (!nrow(seeds)) abort("CSP seed detection produced no seeds")
  if (nrow(seeds) < 2L) abort("CSP requires at least two effective seeds with upstream version 0.2.0")
  data.table::fwrite(seeds, file.path(stage_dir, "effective_seeds.tsv"), sep = "\t", na = "")

  weighted <- any(c(config$verticality_weight, config$linearity_weight, config$sphericity_weight) > 0)
  working <- las
  if (weighted) {
    working <- CspStandSegmentation::add_geometry(
      working,
      k = config$geometry_k,
      n_cores = config$geometry_threads
    )
  }
  working <- lidR::add_lasattribute(working, seq_len(nrow(working@data)), "CSPPointIndex", "Temporary point index")
  segmented <- CspStandSegmentation::csp_cost_segmentation(
    working,
    seeds,
    Voxel_size = config$voxel_resolution,
    V_w = config$verticality_weight,
    L_w = config$linearity_weight,
    S_w = config$sphericity_weight,
    N_cores = config$routing_workers,
    N_trees = 1
  )
  if (!all(c("TreeID", "CSPPointIndex") %in% names(segmented@data))) abort("CSP did not return point assignments")
  point_indices <- as.integer(round(segmented@data$CSPPointIndex))
  valid <- !is.na(point_indices) & point_indices >= 1L & point_indices <= nrow(las@data)
  if (!any(valid)) abort("CSP returned no assignments that could be mapped to input points")
  ids <- integer(nrow(las@data))
  ids[point_indices[valid]] <- as.integer(segmented@data$TreeID[valid])
  list(
    seeds = seeds,
    ids = ids,
    assigned_point_count = sum(ids > 0L),
    unassigned_point_count = sum(ids == 0L)
  )
}

write_csp_cloud <- function(original_las, csp_ids, path) {
  if ("PredInstance_CSP" %in% names(original_las@data)) abort("Input already contains PredInstance_CSP")
  output <- lidR::add_lasattribute(original_las, csp_ids, "PredInstance_CSP", "CSP tree instance ID")
  lidR::writeLAS(output, path)
  invisible(path)
}

process_peak_rss_kb <- function() {
  status <- tryCatch(readLines("/proc/self/status", warn = FALSE), error = function(error) character())
  line <- grep("^VmHWM:", status, value = TRUE)
  if (!length(line)) return(NA_real_)
  as.numeric(sub("^VmHWM:\\s+([0-9]+).*", "\\1", line[[1L]]))
}

write_outputs <- function(inventories, area, stage_dir) {
  inventory_dir <- file.path(stage_dir, "inventories")
  dir.create(inventory_dir, recursive = TRUE)
  for (name in names(inventories)) {
    data.table::fwrite(inventories[[name]], file.path(inventory_dir, paste0(name, ".tsv")), sep = "\t", na = "")
  }
  combined <- data.table::rbindlist(inventories, use.names = TRUE, fill = TRUE)
  has_species <- any(vapply(inventories, function(x) "species_id" %in% names(x), logical(1)))
  if (!has_species) {
    species_columns <- intersect(c("species_id", "species_prob", "species_conflict_fraction"), names(combined))
    if (length(species_columns)) combined[, (species_columns) := NULL]
  }
  data.table::fwrite(combined, file.path(stage_dir, "inventory_combined.tsv"), sep = "\t", na = "")
  summaries <- data.table::rbindlist(lapply(inventories, stand_summary_for_inventory, area = area), fill = TRUE)
  data.table::fwrite(summaries, file.path(stage_dir, "stand_summary.tsv"), sep = "\t", na = "")
  composition <- species_composition(inventories)
  if (!is.null(composition)) {
    data.table::fwrite(composition, file.path(stage_dir, "species_composition.tsv"), sep = "\t", na = "")
  }
}

publish_stage <- function(stage_dir, output_dir) {
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    abort(sprintf("Output directory is not empty: %s", output_dir))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  entries <- list.files(stage_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  if (!all(file.rename(entries, file.path(output_dir, basename(entries))))) abort("Failed to publish staged outputs")
  unlink(stage_dir, recursive = TRUE, force = TRUE)
}

point_cloud_bounds <- function(files) {
  headers <- lapply(files, lidR::readLASheader)
  list(
    xmin = min(vapply(headers, function(x) x@PHB[["Min X"]], numeric(1))),
    xmax = max(vapply(headers, function(x) x@PHB[["Max X"]], numeric(1))),
    ymin = min(vapply(headers, function(x) x@PHB[["Min Y"]], numeric(1))),
    ymax = max(vapply(headers, function(x) x@PHB[["Max Y"]], numeric(1))),
    headers = headers
  )
}

spatial_chunks <- function(bounds, size) {
  x_breaks <- seq(floor(bounds$xmin / size) * size, ceiling(bounds$xmax / size) * size, by = size)
  y_breaks <- seq(floor(bounds$ymin / size) * size, ceiling(bounds$ymax / size) * size, by = size)
  if (length(x_breaks) < 2L) x_breaks <- c(x_breaks, x_breaks + size)
  if (length(y_breaks) < 2L) y_breaks <- c(y_breaks, y_breaks + size)
  # Core reads are half-open on their upper edges to avoid double counting.
  # Add a final interval when the point-cloud maximum lies exactly on a grid
  # line so points on that maximum are still assigned to one chunk.
  if (tail(x_breaks, 1L) <= bounds$xmax) x_breaks <- c(x_breaks, tail(x_breaks, 1L) + size)
  if (tail(y_breaks, 1L) <= bounds$ymax) y_breaks <- c(y_breaks, tail(y_breaks, 1L) + size)
  grid <- expand.grid(x = seq_len(length(x_breaks) - 1L), y = seq_len(length(y_breaks) - 1L))
  lapply(seq_len(nrow(grid)), function(index) {
    x <- grid$x[[index]]
    y <- grid$y[[index]]
    c(xmin = x_breaks[[x]], xmax = x_breaks[[x + 1L]], ymin = y_breaks[[y]], ymax = y_breaks[[y + 1L]])
  })
}

read_spatial_chunk <- function(files, selectors, extent, specs = NULL) {
  filter <- sprintf(
    "-keep_xy %.10f %.10f %.10f %.10f",
    extent[["xmin"]], extent[["ymin"]], extent[["xmax"]], extent[["ymax"]]
  )
  clouds <- Map(function(file, select) {
    cloud <- suppressWarnings(lidR::readLAS(file, select = select, filter = filter))
    if (!is.null(cloud) && !lidR::is.empty(cloud) && !is.null(specs)) {
      cloud <- project_inventory_dimensions(cloud, specs)
    }
    cloud
  }, files, selectors)
  clouds <- Filter(function(x) !is.null(x) && !lidR::is.empty(x), clouds)
  if (!length(clouds)) return(NULL)
  if (length(clouds) == 1L) clouds[[1L]] else CspStandSegmentation::las_merge(clouds, fill = TRUE)
}

chunked_dtm <- function(files, chunks, buffer, resolution, work_dir, stage_dir) {
  partial_dir <- file.path(work_dir, "dtm")
  dir.create(partial_dir, recursive = TRUE)
  partials <- character()
  methods <- character()
  selectors <- rep("c", length(files))
  for (index in seq_along(chunks)) {
    core <- chunks[[index]]
    buffered <- core + c(xmin = -buffer, xmax = buffer, ymin = -buffer, ymax = buffer)
    las <- read_spatial_chunk(files, selectors, buffered)
    if (is.null(las)) next
    result <- tryCatch(make_dtm(las, resolution), error = function(error) NULL)
    if (is.null(result)) next
    tile <- tryCatch(
      terra::crop(result$dtm, terra::ext(core[["xmin"]], core[["xmax"]], core[["ymin"]], core[["ymax"]])),
      error = function(error) NULL
    )
    if (is.null(tile) || terra::ncell(tile) == 0L || all(is.na(terra::values(tile)))) next
    path <- file.path(partial_dir, sprintf("dtm_%05d.tif", index))
    terra::writeRaster(tile, path, overwrite = TRUE)
    partials <- c(partials, path)
    methods <- c(methods, result$method)
    rm(las, result, tile)
    gc(FALSE)
  }
  if (!length(partials)) abort("Chunked DTM generation produced no valid cells")
  dtm <- terra::vrt(partials)
  output <- file.path(stage_dir, "dtm_full.tif")
  terra::writeRaster(dtm, output, overwrite = TRUE)
  list(
    dtm = terra::rast(output),
    method = if (all(methods == "classification_2")) "chunked_classification_2" else "chunked_csf",
    chunk_count = length(chunks)
  )
}

update_streaming_hull <- function(hull, x, y) {
  candidate <- data.table::data.table(x = x, y = y)
  if (nrow(candidate) > 3L) candidate <- candidate[chull(x, y)]
  if (!is.null(hull)) candidate <- data.table::rbindlist(list(hull, candidate))
  if (nrow(candidate) > 3L) candidate <- candidate[chull(x, y)]
  candidate
}

append_partition <- function(path, values) {
  connection <- file(path, open = "ab")
  on.exit(close(connection))
  writeBin(as.double(t(as.matrix(values))), connection, size = 8L)
}

stream_inventory_partitions <- function(files, chunks, specs, config, work_dir) {
  selectors <- vapply(files, read_selector, character(1), specs = specs, preserve_all_dimensions = FALSE)
  partition_root <- file.path(work_dir, "partitions")
  dir.create(partition_root, recursive = TRUE)
  layouts <- list()
  for (spec in specs) {
    slug <- gsub("[^A-Za-z0-9_]+", "_", spec$source)
    columns <- unique(c(
      "X", "Y", "Z", spec$instance, spec$species, spec$species_prob,
      if (identical(spec$source, "FM")) c("PredScore_FM", "PredSemantic_FM")
    ))
    columns <- columns[!is.na(columns) & nzchar(columns)]
    directory <- file.path(partition_root, slug)
    dir.create(directory)
    layouts[[slug]] <- list(spec = spec, columns = columns, directory = directory)
  }

  point_count <- 0
  hull <- NULL
  loaded_dimensions <- character()
  dropped_dimensions <- character()
  for (core in chunks) {
    las <- read_spatial_chunk(files, selectors, core, specs)
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

    for (layout in layouts) {
      spec <- layout$spec
      instance <- spec$instance
      available_columns <- intersect(layout$columns, names(points))
      selected <- points[
        !is.na(get(instance)) & is.finite(get(instance)) &
          get(instance) >= 0 & !(get(instance) %in% config$non_tree_ids),
        ..available_columns
      ]
      if (!nrow(selected)) next
      missing <- setdiff(layout$columns, names(selected))
      for (dimension in missing) data.table::set(selected, j = dimension, value = NA_real_)
      layout_columns <- layout$columns
      selected <- selected[, ..layout_columns]
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
    rm(las, points)
    gc(FALSE)
  }
  list(
    layouts = layouts,
    point_count = point_count,
    hull = hull,
    loaded_dimensions = loaded_dimensions,
    dropped_dimensions = dropped_dimensions,
    selectors = selectors
  )
}

read_binary_partition <- function(path, columns) {
  count <- file.info(path)$size / 8
  values <- readBin(path, numeric(), n = count, size = 8L)
  if (length(values) %% length(columns) != 0L) abort(sprintf("Invalid inventory partition: %s", path))
  data.table::as.data.table(matrix(values, ncol = length(columns), byrow = TRUE, dimnames = list(NULL, columns)))
}

inventory_from_partitions <- function(stream, dtm, config, aoi, crs) {
  inventories <- list()
  for (slug in names(stream$layouts)) {
    layout <- stream$layouts[[slug]]
    paths <- sort(list.files(layout$directory, pattern = "\\.bin$", full.names = TRUE))
    partials <- list()
    for (path in paths) {
      points <- read_binary_partition(path, layout$columns)
      terrain <- terra::extract(dtm, data.frame(X = points$X, Y = points$Y), ID = FALSE)[[1L]]
      data.table::set(points, j = "Zref", value = points$Z)
      data.table::set(points, j = "Z", value = points$Z - terrain)
      points <- points[is.finite(Z)]
      if (!nrow(points)) next
      las <- lidR::LAS(as.data.frame(points))
      partials[[length(partials) + 1L]] <- inventory_for_spec(las, layout$spec, config, aoi, crs)
      rm(points, las)
      gc(FALSE)
    }
    if (!length(partials)) abort(sprintf("No valid inventory partitions for %s", layout$spec$instance))
    inventory <- data.table::rbindlist(partials, use.names = TRUE, fill = TRUE)
    data.table::setorder(inventory, instance_id)
    attr(inventory, "segmentation_source") <- layout$spec$source
    attr(inventory, "instance_dimension") <- layout$spec$instance
    inventories[[slug]] <- inventory
  }
  inventories
}

run_chunked_inventory <- function(config, files, stage_dir, output_parent, started) {
  work_dir <- tempfile("csp-chunks-", tmpdir = output_parent)
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
  bounds <- point_cloud_bounds(files)
  chunks <- spatial_chunks(bounds, config$read_chunk_size)
  available <- unique(c(
    "X", "Y", "Z", "Classification",
    unlist(lapply(bounds$headers, extra_byte_names), use.names = FALSE)
  ))
  for (spec in config$segmentation_specs) {
    required <- c(spec$instance, spec$species, spec$species_prob)
    missing <- setdiff(required[!is.na(required)], available)
    if (length(missing)) abort(sprintf("Requested dimensions are missing for %s: %s", spec$instance, paste(missing, collapse = ", ")))
  }
  las_crs <- sf::st_crs(bounds$headers[[1L]])
  aoi <- read_aoi(config$aoi_json, las_crs)
  dtm_result <- chunked_dtm(files, chunks, config$chunk_buffer, config$dtm_resolution, work_dir, stage_dir)
  if (!is.null(aoi)) {
    aoi_dtm <- terra::mask(terra::crop(dtm_result$dtm, terra::vect(aoi)), terra::vect(aoi))
    terra::writeRaster(aoi_dtm, file.path(stage_dir, "dtm_aoi.tif"), overwrite = TRUE)
  }
  stream <- stream_inventory_partitions(files, chunks, config$segmentation_specs, config, work_dir)
  area <- if (!is.null(aoi)) {
    list(area_m2 = as.numeric(sf::st_area(aoi)), area_source = "aoi")
  } else {
    list(area_m2 = polygon_area(stream$hull$x, stream$hull$y), area_source = "point_cloud_convex_hull")
  }
  inventories <- inventory_from_partitions(stream, dtm_result$dtm, config, aoi, las_crs)
  write_outputs(inventories, area, stage_dir)

  metadata <- list(
    tool = "3Dtrees: CSP StandSegmentation",
    package_version = as.character(utils::packageVersion("CspStandSegmentation")),
    input_files = unname(files), input_bytes = unname(sum(file.info(files)$size)),
    input_points = stream$point_count, input_dimensions = stream$loaded_dimensions,
    input_read_selectors = unname(stream$selectors),
    input_dimensions_dropped = stream$dropped_dimensions,
    segmentation_dimensions = vapply(config$segmentation_specs, `[[`, character(1), "instance"),
    species_dimensions_supplied = any(vapply(config$segmentation_specs, function(spec) !is.null(spec$species), logical(1))),
    csp_enabled = FALSE, csp_seed_count = 0L, csp_unassigned_point_count = 0L,
    read_mode = "spatial_chunks_and_instance_partitions",
    read_chunk_size_m = config$read_chunk_size,
    chunk_buffer_m = config$chunk_buffer,
    spatial_chunk_count = length(chunks),
    inventory_partitions = config$inventory_partitions,
    ground_method = dtm_result$method, dtm_resolution_m = config$dtm_resolution,
    area_source = area$area_source, area_m2 = area$area_m2,
    random_seed = config$random_seed, routing_workers = config$routing_workers,
    geometry_threads = config$geometry_threads, voxel_resolution_m = config$voxel_resolution
  )
  jsonlite::write_json(metadata, file.path(stage_dir, "run_metadata.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
  completed <- Sys.time()
  resources <- list(
    elapsed_seconds = as.numeric(difftime(completed, started, units = "secs")),
    process_peak_rss_kb = process_peak_rss_kb(), input_bytes = unname(sum(file.info(files)$size)),
    input_points = stream$point_count, routing_workers = config$routing_workers,
    geometry_threads = config$geometry_threads,
    note = "Process-local VmHWM only; container/cgroup telemetry is recorded by benchmark runs."
  )
  jsonlite::write_json(resources, file.path(stage_dir, "resource_summary.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
  resources
}

run_tool <- function(config) {
  started <- Sys.time()
  set.seed(config$random_seed)
  lidR::set_lidr_threads(config$geometry_threads)
  files <- resolve_input_files(config$input)
  output_parent <- dirname(config$output_dir)
  dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)
  stage_dir <- tempfile("csp-stage-", tmpdir = output_parent)
  dir.create(stage_dir)
  published <- FALSE
  on.exit(if (!published) unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)

  if (!config$enable_csp) {
    resources <- run_chunked_inventory(config, files, stage_dir, output_parent, started)
    publish_stage(stage_dir, config$output_dir)
    published <- TRUE
    message(sprintf("Completed CSP stand inventory in %.1f seconds", resources$elapsed_seconds))
    return(invisible(resources))
  }

  original_las <- read_point_cloud(
    files,
    config$segmentation_specs,
    preserve_all_dimensions = config$enable_csp
  )
  input_point_count <- nrow(original_las@data)
  input_dimensions <- names(original_las@data)
  input_read_selectors <- attr(original_las, "read_selectors")
  input_dimensions_dropped <- attr(original_las, "dropped_dimensions")
  original_data <- if (config$enable_csp) data.table::copy(original_las@data) else NULL
  validate_specs(original_las, config$segmentation_specs)
  las_crs <- sf::st_crs(original_las)
  aoi <- read_aoi(config$aoi_json, las_crs)
  area <- determine_area(original_las, aoi)

  dtm_result <- make_dtm(original_las, config$dtm_resolution)
  terra::writeRaster(dtm_result$dtm, file.path(stage_dir, "dtm_full.tif"), overwrite = TRUE)
  if (!is.null(aoi)) {
    aoi_dtm <- terra::mask(terra::crop(dtm_result$dtm, terra::vect(aoi)), terra::vect(aoi))
    terra::writeRaster(aoi_dtm, file.path(stage_dir, "dtm_aoi.tif"), overwrite = TRUE)
  }
  normalized <- lidR::normalize_height(dtm_result$las, lidR::tin(), dtm = dtm_result$dtm)

  specs <- config$segmentation_specs
  inventory_las <- normalized
  csp_seed_count <- 0L
  csp_unassigned_point_count <- 0L
  if (config$enable_csp) {
    csp <- run_csp(normalized, config, stage_dir)
    csp_seed_count <- nrow(csp$seeds)
    csp_unassigned_point_count <- csp$unassigned_point_count
    original_las@data <- original_data
    write_csp_cloud(original_las, csp$ids, file.path(stage_dir, "segmented_csp.laz"))
    inventory_las <- normalized
    inventory_las <- lidR::add_lasattribute(inventory_las, csp$ids, "PredInstance_CSP", "CSP tree instance ID")
    specs <- c(specs, list(list(instance = "PredInstance_CSP", source = "CSP", species = NULL, species_prob = NULL)))
  }

  validate_specs(inventory_las, specs)
  inventories <- list()
  for (spec in specs) {
    result <- inventory_for_spec(inventory_las, spec, config, aoi, las_crs)
    slug <- gsub("[^A-Za-z0-9_]+", "_", spec$source)
    inventories[[slug]] <- result
  }
  write_outputs(inventories, area, stage_dir)

  metadata <- list(
    tool = "3Dtrees: CSP StandSegmentation",
    package_version = as.character(utils::packageVersion("CspStandSegmentation")),
    input_files = unname(files),
    input_bytes = unname(sum(file.info(files)$size)),
    input_points = input_point_count,
    input_dimensions = input_dimensions,
    input_read_selectors = unname(input_read_selectors),
    input_dimensions_dropped = unname(input_dimensions_dropped),
    segmentation_dimensions = vapply(specs, `[[`, character(1), "instance"),
    species_dimensions_supplied = any(vapply(specs, function(spec) !is.null(spec$species), logical(1))),
    csp_enabled = config$enable_csp,
    csp_seed_count = csp_seed_count,
    csp_unassigned_point_count = csp_unassigned_point_count,
    ground_method = dtm_result$method,
    dtm_resolution_m = config$dtm_resolution,
    area_source = area$area_source,
    area_m2 = area$area_m2,
    random_seed = config$random_seed,
    routing_workers = config$routing_workers,
    geometry_threads = config$geometry_threads,
    voxel_resolution_m = config$voxel_resolution
  )
  jsonlite::write_json(metadata, file.path(stage_dir, "run_metadata.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
  completed <- Sys.time()
  resources <- list(
    elapsed_seconds = as.numeric(difftime(completed, started, units = "secs")),
    process_peak_rss_kb = process_peak_rss_kb(),
    input_bytes = unname(sum(file.info(files)$size)),
    input_points = input_point_count,
    routing_workers = config$routing_workers,
    geometry_threads = config$geometry_threads,
    note = "Process-local VmHWM only; container/cgroup telemetry is recorded by benchmark runs."
  )
  jsonlite::write_json(resources, file.path(stage_dir, "resource_summary.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")

  publish_stage(stage_dir, config$output_dir)
  published <- TRUE
  message(sprintf("Completed CSP stand inventory in %.1f seconds", resources$elapsed_seconds))
  invisible(config$output_dir)
}
