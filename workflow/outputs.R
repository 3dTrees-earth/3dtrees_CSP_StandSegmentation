inventory_schema_version <- "1.0.0"

validate_inventory_crs <- function(crs) {
  if (is.na(crs)) abort("Point cloud CRS is missing; metric inventory output needs a valid CRS")
  if (isTRUE(sf::st_is_longlat(crs)) ||
      !tolower(crs$units_gdal %||% "") %in% c("metre", "meter", "metres", "meters")) {
    abort("Point cloud CRS must be projected with metre units")
  }
  invisible(crs)
}

`%||%` <- function(left, right) if (is.null(left) || is.na(left)) right else left

with_geographic_coordinates <- function(inventory, crs) {
  result <- data.table::copy(inventory)
  if (!all(c("x", "y", "z") %in% names(result))) abort("Inventory lacks source x, y, or z coordinates")
  if (any(!is.finite(result$x) | !is.finite(result$y) | !is.finite(result$z))) {
    abort("Inventory contains non-finite source coordinates")
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
  units[intersect(fields, c("x", "y", "z", "dbh_m", "height_m", "dbh_spread_m", "dbh_band_voxel_m",
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

# One line per tree on what to know before using its values.
tree_diagnostics <- function(inventory) {
  method <- if ("inventory_method" %in% names(inventory)) inventory$inventory_method else rep(NA_character_, nrow(inventory))
  large <- !is.na(method) & method == "streaming_large_instance"
  voxel <- if ("dbh_band_voxel_m" %in% names(inventory)) inventory$dbh_band_voxel_m else rep(NA_real_, nrow(inventory))
  thinned <- is.finite(voxel) & voxel > 0
  notes <- cbind(
    ifelse(large, sprintf(
      "streaming_large_instance: %.0f points; species probability and FM score use a deterministic sample of at most %d points%s",
      inventory$point_count, inventory_large_instance_sample_limit,
      ifelse(is.na(inventory$dbh_m), "; x/y from that sample because no DBH could be fitted", "")), NA),
    ifelse(thinned, sprintf(
      "dbh_band_thinned: stem band kept %.0f of %.0f points, one per %.0f mm voxel",
      inventory$dbh_band_points, inventory$dbh_band_source_points, voxel * 1000), NA))
  apply(notes, 1L, function(row) if (all(is.na(row))) NA_character_ else paste(row[!is.na(row)], collapse = " | "))
}

evidence_columns <- c("dbh_fits", "dbh_slices", "crown_hull")

# inventory.csv and inventory.parquet hold the same table: one row per tree and
# segmentation, measurements first, then `diagnostics` and the evidence columns
# (JSON arrays and WKT) that the viewer reads only for an opened tree.
# stand_summary.csv holds one row per segmentation.
write_outputs <- function(inventories, area, stage_dir, config, crs) {
  validate_inventory_crs(crs)
  lookup <- read_species_lookup(config$species_lookup)
  inventories <- lapply(inventories, function(inventory) {
    attach_species_metadata(with_geographic_coordinates(inventory, crs), lookup)
  })
  inventory <- data.table::rbindlist(inventories, use.names = TRUE, fill = TRUE)
  inventory[, instance_id := as.numeric(instance_id)]
  inventory[, diagnostics := tree_diagnostics(inventory)]
  evidence <- tree_evidence(config, inventory)
  if (!is.null(evidence)) inventory <- merge(inventory, evidence, by = instance_sort_keys, all.x = TRUE, sort = FALSE)
  for (column in setdiff(evidence_columns, names(inventory))) data.table::set(inventory, j = column, value = NA_character_)
  leading <- c("segmentation_source", instance_sort_keys)
  data.table::setcolorder(inventory, c(leading, setdiff(names(inventory), c(leading, "diagnostics", evidence_columns)),
                                       "diagnostics", evidence_columns))
  sorted_by <- sort_by_instance(inventory)
  summary <- data.table::rbindlist(lapply(inventories, stand_summary_for_inventory, area = area), fill = TRUE)
  sort_by_instance(summary)
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

manifest_scalar <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
}

validate_inventory_manifest <- function(manifest, stage_dir, expected_dimensions = NULL) {
  if (!identical(manifest$schema_version, inventory_schema_version)) {
    abort("Inventory manifest schema version is invalid")
  }
  coordinates <- manifest$coordinate_reference
  if (!manifest_scalar(coordinates$source_wkt) ||
      !identical(coordinates$source_xy_units, "m") ||
      !identical(coordinates$source_z_units, "m") ||
      !identical(coordinates$longitude_field, "longitude") ||
      !identical(coordinates$latitude_field, "latitude") ||
      !identical(as.integer(coordinates$wgs84_epsg), 4326L)) {
    abort("Inventory manifest has unsafe coordinate metadata")
  }
  source_crs <- sf::st_crs(coordinates$source_wkt)
  validate_inventory_crs(source_crs)
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
  dtm_raster <- terra::rast(file.path(stage_dir, "dtm_full.tif"))
  dtm_bounds <- list(xmin = terra::xmin(dtm_raster), xmax = terra::xmax(dtm_raster),
                     ymin = terra::ymin(dtm_raster), ymax = terra::ymax(dtm_raster))
  if (!identical(raster$key, "dtm_full.tif") || !manifest_scalar(raster$source) ||
      !manifest_scalar(raster$crs_wkt) ||
      !same_horizontal_crs(source_crs, sf::st_crs(raster$crs_wkt), dtm_bounds) ||
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
  validate_inventory_crs(crs)
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
  raster <- terra::rast(file.path(stage_dir, "dtm_full.tif"))
  manifest <- list(
    schema_version = inventory_schema_version,
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    coordinate_reference = list(
      source_epsg = crs$epsg, source_wkt = crs$wkt,
      source_xy_units = "m", source_z_units = "m",
      longitude_field = "longitude", latitude_field = "latitude", wgs84_epsg = 4326L
    ),
    segmentations = unname(specs),
    tables = lapply(exports$tables, function(table) list(fields = output_field_units(names(table)))),
    dtm = list(key = "dtm_full.tif", source = dtm$method,
               crs_wkt = terra::crs(raster), resolution_m = as.numeric(terra::res(raster)),
               aoi_key = if (file.exists(file.path(stage_dir, "dtm_aoi.tif"))) "dtm_aoi.tif" else NULL),
    artifacts = c(exports$artifacts, list(
      list(key = "dtm_full.tif", format = "geotiff"),
      list(key = "run_metadata.json", format = "json"),
      list(key = "resource_summary.json", format = "json")
    ), lapply(intersect(c("dtm_aoi.tif", "effective_seeds.tsv", "segmented_csp.laz", "instance_index.json"),
                        list.files(stage_dir)), function(key) list(key = key, format = tools::file_ext(key))))
  )
  expected_dimensions <- vapply(exports$inventories, function(inventory) unique(inventory$instance_dimension), character(1))
  validate_inventory_manifest(manifest, stage_dir, expected_dimensions)
  jsonlite::write_json(manifest, file.path(stage_dir, "manifest.json"), auto_unbox = TRUE,
                       pretty = TRUE, na = "null", null = "null")
  invisible(manifest)
}
