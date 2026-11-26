inventory_schema_version <- "1.0.0"

# LAS WKT and GeoTIFF may carry WGS 84 geocentric without the top-level EPSG ID.
is_inventory_standin <- function(crs) {
  isTRUE(crs$epsg == 4978L) ||
    (!is.na(crs) && identical(crs$proj4string, sf::st_crs(4978)$proj4string))
}

local_inventory_crs <- function(crs) is.na(crs) || is_inventory_standin(crs)

same_local_inventory_crs <- function(source_crs, raster_crs) {
  (is.na(source_crs) && is.na(raster_crs)) ||
    (is_inventory_standin(source_crs) && is_inventory_standin(raster_crs))
}

near_geocentric_origin <- function(bounds) {
  xyz <- unlist(bounds[c("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")], use.names = FALSE)
  if (length(xyz) != 6L || any(!is.finite(xyz))) return(FALSE)
  farthest <- apply(matrix(abs(xyz), nrow = 2L), 2L, max)
  sum(farthest^2) < 6000000^2
}

validate_inventory_crs <- function(crs, bounds) {
  if (is.na(crs)) {
    xy_spans <- c(bounds$xmax - bounds$xmin, bounds$ymax - bounds$ymin)
    if (max(xy_spans) < 1) {
      abort("Point cloud has no CRS and spans less than one unit in both X and Y; supply a verified projected metre CRS")
    }
    if (max(abs(c(bounds$zmin, bounds$zmax))) > 10000 || bounds$zmax - bounds$zmin > 10000) {
      abort("Point cloud has no CRS and elevations or elevation range exceed 10,000 units; supply a verified projected metre CRS")
    }
    centre <- c(mean(c(bounds$xmin, bounds$xmax)), mean(c(bounds$ymin, bounds$ymax)),
                mean(c(bounds$zmin, bounds$zmax)))
    radius <- sqrt(sum(centre^2))
    # A UTM northing can sit inside the ECEF radius band. The elevation
    # check above catches ECEF away from the equator. Near the equator,
    # prefer a plausible UTM frame when the coordinates are ambiguous.
    utm_like <- centre[[1L]] >= 0 && centre[[1L]] <= 1000000 && centre[[2L]] >= 0
    if (!utm_like && radius >= 6350000 && radius <= 6400000) {
      abort("Point cloud has no CRS and its coordinates resemble geocentric (ECEF) coordinates; supply a verified projected metre CRS")
    }
    return(invisible(crs))
  }
  if (is_inventory_standin(crs) && near_geocentric_origin(bounds)) return(invisible(crs))
  if (is_inventory_standin(crs) || is_geocentric_wkt(crs$wkt)) {
    abort("Point cloud has a geocentric CRS; only local EPSG:4978 stand-ins are supported")
  }
  if (isTRUE(sf::st_is_longlat(crs)) ||
      !tolower(crs$units_gdal %||% "") %in% c("metre", "meter", "metres", "meters")) {
    abort("Point cloud CRS must be projected with metre units")
  }
  invisible(crs)
}

inventory_crs_label <- function(crs) {
  if (is.na(crs)) return("local (assumed metres)")
  if (is_inventory_standin(crs)) return("local (EPSG:4978 stand-in)")
  if (is.na(crs$epsg)) crs$wkt else sprintf("EPSG:%d", crs$epsg)
}

`%||%` <- function(left, right) if (is.null(left) || is.na(left)) right else left

with_geographic_coordinates <- function(inventory, crs) {
  result <- data.table::copy(inventory)
  if (!all(c("x", "y", "z") %in% names(result))) abort("Inventory lacks source x, y, or z coordinates")
  if (any(!is.finite(result$x) | !is.finite(result$y) | !is.finite(result$z))) {
    abort("Inventory contains non-finite source coordinates")
  }
  if (local_inventory_crs(crs)) {
    result[, `:=`(longitude = NA_real_, latitude = NA_real_)]
    return(result)
  }
  geo <- sf::st_as_sf(as.data.frame(result), coords = c("x", "y"), crs = crs)
  transformed <- sf::st_transform(geo, 4326)
  coordinates <- sf::st_coordinates(transformed)
  result[, `:=`(longitude = coordinates[, 1L], latitude = coordinates[, 2L])]
  if (any(!is.finite(result$longitude) | !is.finite(result$latitude) |
          abs(result$longitude) > 180 | abs(result$latitude) > 90)) {
    abort("Inventory WGS84 coordinate transformation produced invalid positions")
  }
  result
}

read_species_lookup <- function(path) {
  if (is.null(path)) return(NULL)
  if (!file.exists(path)) abort(sprintf("Species lookup does not exist: %s", path))
  lookup <- data.table::fread(path)
  required <- c("segmentation_source", "species_id", "scientific_name", "taxonomy_version", "model_version")
  missing <- setdiff(required, names(lookup))
  if (length(missing)) abort(sprintf("Species lookup is missing columns: %s", paste(missing, collapse = ", ")))
  if (!"common_name" %in% names(lookup)) lookup[, common_name := NA_character_]
  if (anyDuplicated(lookup[, .(segmentation_source, species_id)])) {
    abort("Species lookup contains duplicate segmentation_source/species_id pairs")
  }
  for (column in c("scientific_name", "taxonomy_version", "model_version")) {
    if (anyNA(lookup[[column]]) || any(!nzchar(trimws(lookup[[column]])))) {
      abort(sprintf("Species lookup %s must be populated", column))
    }
  }
  lookup[, c(required, "common_name"), with = FALSE]
}

attach_species_metadata <- function(inventory, lookup) {
  if (is.null(lookup) || !"species_id" %in% names(inventory)) return(inventory)
  merge(inventory, lookup, by = c("segmentation_source", "species_id"), all.x = TRUE, sort = FALSE)
}

output_field_units <- function(fields) {
  units <- stats::setNames(rep("unitless", length(fields)), fields)
  units[intersect(fields, c("x", "y", "z", "dbh_m", "height_m", "dbh_spread_m",
                            "mean_dbh_m", "median_dbh_m",
                            "dbh_p05_m", "dbh_p25_m", "dbh_p75_m", "dbh_p95_m", "mean_height_m",
                            "median_height_m", "height_p05_m", "height_p25_m", "height_p75_m",
                            "height_p95_m"))] <- "m"
  units[intersect(fields, c("area_m2", "convex_hull_area_m2"))] <- "m2"
  units[intersect(fields, "basal_area_m2_per_ha")] <- "m2/ha"
  units[intersect(fields, c("longitude", "latitude"))] <- "degree"
  units[intersect(fields, c("trees_per_ha"))] <- "trees/ha"
  as.list(units)
}

evidence_columns <- c("dbh_fits", "dbh_slices", "crown_hull")

# inventory.csv and inventory.parquet hold the same table: one row per tree and
# segmentation, measurements first, then the DBH evidence as
# JSON (slices, runs and the crown hull as a GeoJSON geometry), which the viewer
# reads only for an opened tree. stand_summary.csv holds one row per
# segmentation and the CRS of every coordinate.
write_outputs <- function(inventories, area, stage_dir, config, crs) {
  lookup <- read_species_lookup(config$species_lookup)
  inventories <- lapply(inventories, function(inventory) {
    attach_species_metadata(with_geographic_coordinates(inventory, crs), lookup)
  })
  inventory <- data.table::rbindlist(inventories, use.names = TRUE, fill = TRUE)
  inventory[, instance_id := as.numeric(instance_id)]
  evidence <- tree_evidence(config, inventory)
  if (!is.null(evidence)) inventory <- merge(inventory, evidence, by = instance_sort_keys, all.x = TRUE, sort = FALSE)
  for (column in setdiff(evidence_columns, names(inventory))) data.table::set(inventory, j = column, value = NA_character_)
  leading <- c("segmentation_source", instance_sort_keys)
  data.table::setcolorder(inventory, c(leading, setdiff(names(inventory), c(leading, evidence_columns)),
                                       evidence_columns))
  seeded <- intersect(c("seed_x", "seed_y", "seed_height_m"), names(inventory))
  if (length(seeded)) data.table::setcolorder(inventory, seeded, after = "z")
  sorted_by <- sort_by_instance(inventory)
  summary <- data.table::rbindlist(lapply(inventories, stand_summary_for_inventory, area = area), fill = TRUE)
  sort_by_instance(summary)
  summary[, crs := inventory_crs_label(crs)]
  write_csv_table(inventory, file.path(stage_dir, "inventory.csv"))
  write_parquet_table(inventory, file.path(stage_dir, "inventory.parquet"))
  write_csv_table(summary, file.path(stage_dir, "stand_summary.csv"))
  artifacts <- list(
    list(key = "inventory.csv", format = "csv", table = "inventory", sorted_by = as.list(sorted_by)),
    list(key = "inventory.parquet", format = "parquet", table = "inventory", sorted_by = as.list(sorted_by),
         row_group_rows = parquet_row_group_rows),
    list(key = "stand_summary.csv", format = "csv", table = "stand_summary")
  )
  list(inventories = inventories, tables = list(inventory = inventory, stand_summary = summary),
       artifacts = artifacts)
}

# Detailed settings, timings and audits require --debug-outputs true.
write_debug_json <- function(config, value, path) {
  if (isTRUE(config$debug_outputs)) {
    jsonlite::write_json(value, path, auto_unbox = TRUE, pretty = TRUE, na = "null")
  }
}

write_inventory_run_metadata <- function(config, metadata, dtm, stage_dir) {
  debug <- isTRUE(config$debug_outputs)
  if (!debug && dtm$filled_cells == 0L) return(invisible(NULL))
  terrain <- inventory_terrain_metadata(config, dtm)
  prefix <- names(terrain) != "terrain_lookup"
  names(terrain)[prefix] <- paste0("dtm_", names(terrain)[prefix])
  fields <- if (debug) c(metadata, terrain) else terrain
  jsonlite::write_json(fields, file.path(stage_dir, "run_metadata.json"), auto_unbox = TRUE,
                       pretty = TRUE, digits = NA, na = "null")
}

manifest_scalar <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
}

validate_inventory_manifest <- function(manifest, stage_dir, expected_dimensions = NULL) {
  if (!identical(manifest$schema_version, inventory_schema_version)) {
    abort("Inventory manifest schema version is invalid")
  }
  coordinates <- manifest$coordinate_reference
  local <- identical(coordinates$source_frame, "local")
  source_crs <- sf::st_crs(coordinates$source_wkt %||% NA_character_)
  if (local != local_inventory_crs(source_crs) ||
      !identical(coordinates$source_xy_units, "m") ||
      !identical(coordinates$source_z_units, "m") ||
      !identical(coordinates$longitude_field, "longitude") ||
      !identical(coordinates$latitude_field, "latitude") ||
      !(if (local) is.null(coordinates$wgs84_epsg) else
          identical(as.integer(coordinates$wgs84_epsg), 4326L))) {
    abort("Inventory manifest has unsafe coordinate metadata")
  }
  if (!local) validate_inventory_crs(source_crs, NULL)
  if (!is.null(coordinates$source_epsg) &&
      (length(coordinates$source_epsg) != 1L ||
       (!is.na(coordinates$source_epsg) && source_crs != sf::st_crs(coordinates$source_epsg)))) {
    abort("Inventory manifest source EPSG and WKT disagree")
  }
  required_tables <- c("inventory", "stand_summary")
  if (!all(required_tables %in% names(manifest$tables)) || !length(manifest$segmentations)) {
    abort("Inventory manifest is missing required tables or segmentations")
  }
  if (any(vapply(manifest$tables[required_tables], function(table) !length(table$fields), logical(1)))) {
    abort("Inventory manifest has a table without fields")
  }
  if (any(vapply(manifest$segmentations, function(spec) {
    !manifest_scalar(spec$technical_name) || !manifest_scalar(spec$display_name) ||
      !identical(spec$technical_name, spec$potree_instance_dimension) ||
      !all(c("x", "y", "z", "longitude", "latitude", "dbh_quality_code", "dbh_quality_label") %in%
           names(spec$fields)) ||
      !identical(spec$fields$x, "m") || !identical(spec$fields$y, "m") ||
      !identical(spec$fields$z, "m") || !identical(spec$fields$longitude, "degree") ||
      !identical(spec$fields$latitude, "degree")
  }, logical(1)))) abort("Inventory manifest has an incomplete segmentation")
  dimensions <- vapply(manifest$segmentations, `[[`, character(1), "technical_name")
  if (anyDuplicated(dimensions) ||
      (!is.null(expected_dimensions) && !setequal(dimensions, expected_dimensions))) {
    abort("Inventory manifest segmentation dimensions do not match the outputs")
  }
  raster <- manifest$dtm
  dtm_raster <- terra::rast(file.path(stage_dir, "dtm_full.tif"), guessCRS = FALSE)
  dtm_bounds <- list(xmin = terra::xmin(dtm_raster), xmax = terra::xmax(dtm_raster),
                     ymin = terra::ymin(dtm_raster), ymax = terra::ymax(dtm_raster))
  raster_crs <- sf::st_crs(if (nzchar(raster$crs_wkt %||% "")) raster$crs_wkt else NA_character_)
  matching_crs <- if (local) same_local_inventory_crs(source_crs, raster_crs) else
    same_horizontal_crs(source_crs, raster_crs, dtm_bounds)
  if (!identical(raster$key, "dtm_full.tif") || !manifest_scalar(raster$source) ||
      !matching_crs ||
      length(raster$resolution_m) != 2L || any(!is.finite(raster$resolution_m)) ||
      any(raster$resolution_m <= 0)) {
    abort("Inventory manifest DTM metadata is invalid")
  }
  keys <- vapply(manifest$artifacts, `[[`, character(1), "key")
  actual <- setdiff(list.files(stage_dir, recursive = TRUE), "manifest.json")
  if (anyDuplicated(keys) || !setequal(keys, actual) ||
      any(!file.exists(file.path(stage_dir, keys))) || !"dtm_full.tif" %in% keys) {
    abort("Inventory manifest contains duplicate or missing artifact keys")
  }
  if (!all(c("inventory.csv", "inventory.parquet", "stand_summary.csv") %in% keys)) {
    abort("Inventory manifest is missing an inventory or stand summary file")
  }
  invisible(manifest)
}

write_inventory_manifest <- function(exports, config, crs, dtm, stage_dir) {
  specs <- lapply(exports$inventories, function(inventory) {
    list(
      technical_name = unique(inventory$instance_dimension),
      display_name = unique(inventory$segmentation_source),
      potree_instance_dimension = unique(inventory$instance_dimension),
      potree_species_dimension = {
        matching <- Filter(function(spec) identical(spec$instance, unique(inventory$instance_dimension)),
                           config$segmentation_specs)
        if (length(matching)) matching[[1L]]$species else NULL
      },
      fields = output_field_units(names(inventory))
    )
  })
  raster <- terra::rast(file.path(stage_dir, "dtm_full.tif"), guessCRS = FALSE)
  terrain <- inventory_terrain_metadata(config, dtm)
  manifest <- list(
    schema_version = inventory_schema_version,
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    coordinate_reference = list(
      source_epsg = crs$epsg, source_wkt = crs$wkt,
      source_xy_units = "m", source_z_units = "m",
      longitude_field = "longitude", latitude_field = "latitude",
      wgs84_epsg = if (!local_inventory_crs(crs)) 4326L
    ),
    segmentations = unname(specs),
    tables = lapply(exports$tables, function(table) list(fields = output_field_units(names(table)))),
    dtm = c(list(key = "dtm_full.tif", source = dtm$method,
                 crs_wkt = terra::crs(raster), resolution_m = as.numeric(terra::res(raster)),
                 aoi_key = if (file.exists(file.path(stage_dir, "dtm_aoi.tif"))) "dtm_aoi.tif" else NULL),
            terrain),
    artifacts = c(exports$artifacts, list(list(key = "dtm_full.tif", format = "geotiff")),
                  lapply(intersect(c("dtm_aoi.tif", "effective_seeds.tsv", "segmented_csp.laz", "instance_index.json",
                                     "run_metadata.json", "resource_summary.json"), list.files(stage_dir)),
                         function(key) list(key = key, format = tools::file_ext(key))))
  )
  if (local_inventory_crs(crs)) manifest$coordinate_reference$source_frame <- "local"
  expected_dimensions <- vapply(exports$inventories, function(inventory) unique(inventory$instance_dimension), character(1))
  validate_inventory_manifest(manifest, stage_dir, expected_dimensions)
  if (isTRUE(config$debug_outputs)) {
    jsonlite::write_json(manifest, file.path(stage_dir, "manifest.json"), auto_unbox = TRUE,
                         pretty = TRUE, digits = NA, na = "null", null = "null")
  }
  invisible(manifest)
}
