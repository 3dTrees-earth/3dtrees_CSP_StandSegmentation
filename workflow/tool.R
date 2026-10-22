.csp_workflow_dir <- dirname(normalizePath(sys.frame(1)$ofile, mustWork = TRUE))
for (.csp_module in c("inventory_policy.R", "instance_inventory.R", "partition_inventory.R", "owned_inventory.R",
                      "inventory_readers.R", "inventory_run.R")) {
  source(file.path(.csp_workflow_dir, .csp_module), local = TRUE)
}
rm(.csp_workflow_dir, .csp_module)

inventory_standard_attributes <- jsonlite::fromJSON(
  file.path(dirname(normalizePath(sys.frame(1)$ofile)), "inventory_standard_attributes.json"))

inventory_canonical_names <- function(names) {
  native <- vapply(inventory_standard_attributes, `[[`, character(1), "python")
  matched <- match(names, native)
  names[!is.na(matched)] <- names(inventory_standard_attributes)[matched[!is.na(matched)]]
  names
}

resolve_input_files <- function(path) {
  if (file.exists(path) && !dir.exists(path)) return(normalizePath(path, mustWork = TRUE))
  if (!dir.exists(path)) abort(sprintf("Input does not exist: %s", path))
  files <- list.files(path, pattern = "\\.(las|laz)$", full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) abort(sprintf("Input directory contains no LAS/LAZ files: %s", path))
  sort(normalizePath(files, mustWork = TRUE))
}

inventory_spec_dimensions <- function(spec) {
  unique(c(spec$instance, spec$species, spec$species_prob, spec$attributes))
}

inventory_read_dimensions <- function(specs) {
  unique(unlist(lapply(specs, inventory_spec_dimensions), use.names = FALSE))
}

extra_byte_names <- function(header) {
  extra_bytes <- header@VLR$Extra_Bytes[["Extra Bytes Description"]]
  if (is.null(extra_bytes)) character() else names(extra_bytes)
}

point_cloud_crs <- function(header) {
  parse_wkt <- function(records) {
    wkt <- records[["WKT OGC CS"]][["WKT OGC COORDINATE SYSTEM"]]
    if (!is.character(wkt) || length(wkt) != 1L || !nzchar(wkt)) return(NULL)
    tryCatch(suppressWarnings(sf::st_crs(wkt)), error = function(error) NULL)
  }
  candidates <- list(
    tryCatch(suppressWarnings(sf::st_crs(header)), error = function(error) NULL),
    parse_wkt(header@EVLR),
    parse_wkt(header@VLR)
  )
  valid <- Filter(function(crs) !is.null(crs) && !is.na(crs), candidates)
  if (!length(valid)) abort("Point cloud has no parseable CRS in its header or WKT records")
  if (any(vapply(valid[-1L], function(crs) crs != valid[[1L]], logical(1)))) {
    abort("Point cloud header contains conflicting valid CRS records")
  }
  valid[[1L]]
}

build_read_selector <- function(dimensions, extra_dimensions) {
  standard_codes <- vapply(inventory_standard_attributes, `[[`, character(1), "selector")
  standard <- unname(standard_codes[intersect(dimensions, names(standard_codes))])
  positions <- sort(match(intersect(dimensions, extra_dimensions), extra_dimensions))

  # LASlib can address only the first nine extra-byte records individually.
  # If a requested dimension is later, load all extra bytes rather than risk
  # silently omitting a requested inventory field.
  extra <- if (any(positions > 9L)) "0" else as.character(positions)
  paste0(unique(c("xyz", standard, extra)), collapse = "")
}

read_selector <- function(file, specs) {
  header <- lidR::readLASheader(file)
  build_read_selector(inventory_read_dimensions(specs), extra_byte_names(header))
}

project_inventory_dimensions <- function(las, specs) {
  validate_specs(las, specs)
  available <- names(las@data)
  keep <- unique(c(
    "X", "Y", "Z",
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

read_point_cloud <- function(files, specs) {
  selectors <- vapply(
    files,
    read_selector,
    character(1),
    specs = specs
  )
  clouds <- Map(function(file, select) {
    cloud <- suppressWarnings(lidR::readLAS(file, select = select))
    if (is.null(cloud) || lidR::is.empty(cloud)) return(cloud)
    cloud <- project_inventory_dimensions(cloud, specs)
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

load_supplied_dtm <- function(path, source_crs, bounds, stage_dir) {
  if (!file.exists(path) || dir.exists(path)) abort(sprintf("DTM does not exist: %s", path))
  dtm <- terra::rast(path)
  if (terra::nlyr(dtm) != 1L) abort("Supplied DTM must have exactly one raster band")
  raster_crs <- sf::st_crs(terra::crs(dtm))
  if (!same_horizontal_crs(source_crs, raster_crs, bounds)) {
    abort("Supplied DTM and point cloud must have the same valid CRS or equivalent horizontal CRS")
  }
  if (terra::xmin(dtm) > bounds$xmin || terra::xmax(dtm) < bounds$xmax ||
      terra::ymin(dtm) > bounds$ymin || terra::ymax(dtm) < bounds$ymax) {
    abort("Supplied DTM must cover the full point-cloud extent")
  }
  output <- file.path(stage_dir, "dtm_full.tif")
  if (!file.copy(path, output, overwrite = TRUE)) abort("Could not stage the supplied DTM")
  list(
    las = NULL, dtm = terra::rast(output), method = "supplied", strategy = "supplied",
    workers = 0L, valid_chunk_count = 0L, candidate_metadata = NULL,
    worker_peak_rss_sum_kb = 0, worker_peak_rss_max_kb = 0
  )
}

same_horizontal_crs <- function(source_crs, target_crs, bounds, tolerance_m = 0.001) {
  if (is.na(source_crs) || is.na(target_crs) ||
      isTRUE(sf::st_is_longlat(source_crs)) || isTRUE(sf::st_is_longlat(target_crs))) return(FALSE)
  if (source_crs == target_crs) return(TRUE)
  xy <- rbind(
    c(bounds$xmin, bounds$ymin), c(bounds$xmax, bounds$ymax),
    c((bounds$xmin + bounds$xmax) / 2, (bounds$ymin + bounds$ymax) / 2)
  )
  points <- sf::st_sfc(lapply(seq_len(nrow(xy)), function(index) sf::st_point(xy[index, ])),
                       crs = source_crs)
  transformed <- tryCatch(
    sf::st_coordinates(suppressWarnings(sf::st_transform(points, target_crs)))[, 1:2, drop = FALSE],
    error = function(error) NULL
  )
  !is.null(transformed) && all(is.finite(transformed)) &&
    max(abs(transformed - xy)) <= tolerance_m
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
    required <- inventory_spec_dimensions(spec)
    missing <- setdiff(required[!vapply(required, is.null, logical(1))], available)
    if (length(missing)) {
      abort(sprintf("Requested dimensions are missing: %s", paste(missing, collapse = ", ")))
    }
  }
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
  measured_dbh <- is.finite(inventory$dbh_m)
  source <- if (nrow(inventory)) unique(inventory$segmentation_source) else attr(inventory, "segmentation_source")
  dimension <- if (nrow(inventory)) unique(inventory$instance_dimension) else attr(inventory, "instance_dimension")
  data.table::data.table(
    segmentation_source = source,
    instance_dimension = dimension,
    area_source = area$area_source,
    area_m2 = area$area_m2,
    tree_count = nrow(inventory),
    dbh_measured_tree_count = sum(measured_dbh),
    dbh_coverage_fraction = if (nrow(inventory)) mean(measured_dbh) else NA_real_,
    trees_per_ha = if (is.finite(area_ha) && area_ha > 0) nrow(inventory) / area_ha else NA_real_,
    basal_area_m2_per_ha = if (is.finite(area_ha) && area_ha > 0 &&
                               nrow(inventory) > 0L && all(measured_dbh)) {
      sum(pi * (inventory$dbh_m / 2)^2) / area_ha
    } else NA_real_,
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
  crs <- lapply(headers, point_cloud_crs)
  if (any(vapply(crs[-1L], function(value) value != crs[[1L]], logical(1)))) {
    abort("Input point-cloud files have different CRS")
  }
  list(
    xmin = min(vapply(headers, function(x) x@PHB[["Min X"]], numeric(1))),
    xmax = max(vapply(headers, function(x) x@PHB[["Max X"]], numeric(1))),
    ymin = min(vapply(headers, function(x) x@PHB[["Min Y"]], numeric(1))),
    ymax = max(vapply(headers, function(x) x@PHB[["Max Y"]], numeric(1))),
    headers = headers,
    crs = crs[[1L]]
  )
}

point_cloud_bounds_python <- function(file, work_dir) {
  root <- if (exists("repository_root", inherits = TRUE)) {
    get("repository_root", inherits = TRUE)
  } else normalizePath(getwd(), mustWork = TRUE)
  script <- file.path(root, "workflow", "copc_header.py")
  output <- file.path(work_dir, "copc_header.json")
  log <- file.path(work_dir, "copc_header.log")
  status <- system2("python3", c(shQuote(script), "--input", shQuote(file),
                                "--output", shQuote(output)),
                    stdout = log, stderr = log)
  if (!identical(status, 0L) || !file.exists(output)) {
    detail <- if (file.exists(log)) paste(tail(readLines(log, warn = FALSE), 20L), collapse = "\n") else ""
    abort(sprintf("Cannot read 64-bit COPC header%s",
                  if (nzchar(detail)) paste0(":\n", detail) else ""))
  }
  header <- jsonlite::fromJSON(output, simplifyVector = TRUE)
  list(xmin = header$xmin, xmax = header$xmax,
       ymin = header$ymin, ymax = header$ymax,
       crs = sf::st_crs(header$crs_wkt),
       point_count = as.numeric(header$point_count),
       available_dimensions = header$available_dimensions,
       headers = NULL)
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

read_spatial_chunk <- function(files, selectors, extent, specs = NULL, crs = NULL) {
  filter <- sprintf(
    "-keep_xy %.10f %.10f %.10f %.10f",
    extent[["xmin"]], extent[["ymin"]], extent[["xmax"]], extent[["ymax"]]
  )
  clouds <- Map(function(file, select) {
    cloud <- suppressWarnings(lidR::readLAS(file, select = select, filter = filter))
    if (!is.null(cloud) && !lidR::is.empty(cloud) && !is.null(crs) &&
        is.na(suppressWarnings(sf::st_crs(cloud)))) sf::st_crs(cloud) <- crs
    if (!is.null(cloud) && !lidR::is.empty(cloud) && !is.null(specs)) {
      cloud <- project_inventory_dimensions(cloud, specs)
    }
    cloud
  }, files, selectors)
  clouds <- Filter(function(x) !is.null(x) && !lidR::is.empty(x), clouds)
  if (!length(clouds)) return(NULL)
  if (length(clouds) == 1L) clouds[[1L]] else CspStandSegmentation::las_merge(clouds, fill = TRUE)
}

parallel_chunk_map <- function(values, workers, fun) {
  workers <- min(as.integer(workers), length(values))
  if (workers <= 1L) return(lapply(values, fun))
  results <- parallel::mclapply(
    values,
    fun,
    mc.cores = workers,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
  failed <- vapply(results, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    abort(sprintf("Parallel DTM worker failed: %s", as.character(results[[which(failed)[[1L]]]])))
  }
  results
}

chunked_dtm <- function(files, chunks, buffer, resolution, work_dir, stage_dir, workers = 1L) {
  partial_dir <- file.path(work_dir, "dtm")
  dir.create(partial_dir, recursive = TRUE)
  bounds <- point_cloud_bounds(files)
  selectors <- rep("c", length(files))
  process_chunk <- function(index) {
    lidR::set_lidr_threads(1L)
    core <- chunks[[index]]
    buffered <- core + c(xmin = -buffer, xmax = buffer, ymin = -buffer, ymax = buffer)
    las <- read_spatial_chunk(files, selectors, buffered, crs = bounds$crs)
    if (is.null(las)) return(NULL)
    if (nrow(las@data) < 3L || length(unique(las@data$X)) < 2L || length(unique(las@data$Y)) < 2L) {
      return(NULL)
    }
    result <- make_dtm(las, resolution)
    tile <- terra::crop(
      result$dtm,
      terra::ext(core[["xmin"]], core[["xmax"]], core[["ymin"]], core[["ymax"]])
    )
    if (terra::ncell(tile) == 0L || all(is.na(terra::values(tile)))) return(NULL)
    path <- file.path(partial_dir, sprintf("dtm_%05d.tif", index))
    terra::writeRaster(tile, path, overwrite = TRUE)
    list(path = path, method = result$method, peak_rss_kb = process_peak_rss_kb())
  }
  results <- parallel_chunk_map(as.list(seq_along(chunks)), workers, process_chunk)
  results <- Filter(Negate(is.null), results)
  partials <- vapply(results, `[[`, character(1), "path")
  methods <- vapply(results, `[[`, character(1), "method")
  if (!length(partials)) abort("Chunked DTM generation produced no valid cells")
  dtm <- terra::vrt(partials)
  target <- terra::ext(
    floor(bounds$xmin / resolution) * resolution,
    ceiling(bounds$xmax / resolution) * resolution,
    floor(bounds$ymin / resolution) * resolution,
    ceiling(bounds$ymax / resolution) * resolution
  )
  dtm <- terra::extend(dtm, target)
  output <- file.path(stage_dir, "dtm_full.tif")
  terra::writeRaster(dtm, output, overwrite = TRUE)
  list(
    dtm = terra::rast(output),
    method = if (all(methods == "classification_2")) "chunked_classification_2" else "chunked_csf",
    chunk_count = length(chunks),
    valid_chunk_count = length(partials),
    workers = min(as.integer(workers), length(chunks)),
    worker_peak_rss_sum_kb = sum(vapply(results, `[[`, numeric(1), "peak_rss_kb")),
    worker_peak_rss_max_kb = max(vapply(results, `[[`, numeric(1), "peak_rss_kb"))
  )
}

header_point_count <- function(header) {
  value <- header@PHB[["Number of point records"]]
  if (is.null(value)) abort("LAS header does not declare a point count")
  as.numeric(value)
}

streaming_dtm_candidates <- function(files, config, work_dir) {
  candidate_dir <- file.path(work_dir, "dtm_candidates")
  root <- if (exists("repository_root", inherits = TRUE)) {
    get("repository_root", inherits = TRUE)
  } else {
    normalizePath(getwd(), mustWork = TRUE)
  }
  script <- file.path(root, "workflow", "stream_dtm_candidates.py")
  if (!file.exists(script)) abort(sprintf("Streaming DTM helper is missing: %s", script))
  log_path <- file.path(work_dir, "dtm_candidates.log")
  arguments <- c(
    unlist(lapply(files, function(path) c("--input", path)), use.names = FALSE),
    "--output-dir", candidate_dir,
    "--candidate-resolution", as.character(config$dtm_candidate_resolution),
    "--output-resolution", as.character(config$dtm_resolution),
    "--workers", as.character(config$dtm_workers)
  )
  status <- system2("python3", c(script, arguments), stdout = log_path, stderr = log_path)
  if (!identical(status, 0L)) {
    detail <- if (file.exists(log_path)) paste(tail(readLines(log_path, warn = FALSE), 20L), collapse = "\n") else ""
    abort(sprintf("Streaming DTM candidate generation failed%s", if (nzchar(detail)) paste0(":\n", detail) else ""))
  }
  metadata_path <- file.path(candidate_dir, "candidate_metadata.json")
  candidate_path <- file.path(candidate_dir, "minimum_candidates.laz")
  if (!file.exists(metadata_path) || !file.exists(candidate_path)) {
    abort("Streaming DTM candidate generation did not produce its declared outputs")
  }
  list(
    path = candidate_path,
    metadata = jsonlite::fromJSON(metadata_path, simplifyVector = TRUE)
  )
}

make_inventory_dtm <- function(files, bounds, config, work_dir, stage_dir) {
  if (!is.null(config$dtm)) {
    return(load_supplied_dtm(config$dtm, bounds$crs, bounds, stage_dir))
  }
  input_points <- if (!is.null(bounds$point_count)) bounds$point_count else
    sum(vapply(bounds$headers, header_point_count, numeric(1)))
  strategy <- config$dtm_strategy
  if (identical(strategy, "auto")) {
    strategy <- if (input_points >= config$dtm_streaming_threshold) "streaming" else "spatial"
  }
  if (identical(strategy, "spatial")) {
    result <- chunked_dtm(
      files, spatial_chunks(bounds, config$read_chunk_size), config$chunk_buffer, config$dtm_resolution,
      work_dir, stage_dir, config$dtm_workers
    )
    result$strategy <- "spatial"
    result$candidate_metadata <- NULL
    return(result)
  }

  candidates <- streaming_dtm_candidates(files, config, work_dir)
  candidate_bounds <- point_cloud_bounds(candidates$path)
  candidate_chunks <- spatial_chunks(candidate_bounds, config$read_chunk_size)
  result <- chunked_dtm(
    candidates$path, candidate_chunks, config$chunk_buffer,
    config$dtm_resolution, work_dir, stage_dir, config$dtm_workers
  )
  result$method <- paste0("streaming_candidates_", result$method)
  result$strategy <- "streaming"
  result$candidate_metadata <- candidates$metadata
  result
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

  source_crs <- point_cloud_bounds(files)$crs
  validate_inventory_crs(source_crs)

  read_specs <- c(config$segmentation_specs, list(list(attributes = unique(c(
    config$inventory_attributes, if (is.null(config$dtm)) "Classification"
  )))))
  original_las <- read_point_cloud(files, read_specs)
  if (is.na(sf::st_crs(original_las))) sf::st_crs(original_las) <- source_crs
  input_point_count <- nrow(original_las@data)
  input_dimensions <- names(original_las@data)
  input_read_selectors <- attr(original_las, "read_selectors")
  input_dimensions_dropped <- attr(original_las, "dropped_dimensions")
  original_data <- if (config$enable_csp) data.table::copy(original_las@data) else NULL
  validate_specs(original_las, config$segmentation_specs)
  las_crs <- source_crs
  aoi <- read_aoi(config$aoi_json, las_crs)
  area <- determine_area(original_las, aoi)

  dtm_result <- if (is.null(config$dtm)) {
    generated <- make_dtm(original_las, config$dtm_resolution)
    terra::writeRaster(generated$dtm, file.path(stage_dir, "dtm_full.tif"), overwrite = TRUE)
    generated
  } else {
    load_supplied_dtm(
      config$dtm, las_crs,
      list(
        xmin = min(original_las@data$X), xmax = max(original_las@data$X),
        ymin = min(original_las@data$Y), ymax = max(original_las@data$Y)
      ),
      stage_dir
    )
  }
  if (!is.null(aoi)) {
    aoi_dtm <- terra::mask(terra::crop(dtm_result$dtm, terra::vect(aoi)), terra::vect(aoi))
    terra::writeRaster(aoi_dtm, file.path(stage_dir, "dtm_aoi.tif"), overwrite = TRUE)
  }
  normalized <- lidR::normalize_height(
    if (is.null(dtm_result$las)) original_las else dtm_result$las,
    dtm_result$dtm
  )
  if (any(!is.finite(normalized@data$Z))) {
    abort("DTM has missing or non-finite cells under point-cloud points")
  }

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
  exports <- write_outputs(inventories, area, stage_dir, config, las_crs)

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
    output_formats = config$output_formats,
    ground_method = dtm_result$method,
    dtm_resolution_m = as.numeric(terra::res(dtm_result$dtm)),
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
  write_inventory_manifest(exports, config, las_crs, dtm_result, stage_dir)

  publish_stage(stage_dir, config$output_dir)
  published <- TRUE
  message(sprintf("Completed CSP stand inventory in %.1f seconds", resources$elapsed_seconds))
  invisible(config$output_dir)
}
