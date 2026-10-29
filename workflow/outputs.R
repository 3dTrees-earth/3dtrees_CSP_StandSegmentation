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

write_gpkg_attribute <- function(connection, name, table) {
  DBI::dbWriteTable(connection, name, as.data.frame(table), overwrite = TRUE)
  DBI::dbExecute(
    connection,
    "INSERT OR REPLACE INTO gpkg_contents (table_name, data_type, identifier, srs_id) VALUES (?, 'attributes', ?, 0)",
    params = list(name, name)
  )
}

write_inventory_gpkg <- function(path, inventories, tables, crs) {
  layers <- list()
  for (name in names(inventories)) {
    inventory <- inventories[[name]]
    layer <- paste0("inventory_", name)
    points <- sf::st_sfc(lapply(seq_len(nrow(inventory)), function(index) {
      sf::st_point(as.numeric(unlist(inventory[index, .(x, y, z)], use.names = FALSE)), dim = "XYZ")
    }), crs = crs)
    spatial <- sf::st_sf(as.data.frame(inventory), geometry = points)
    sf::st_write(spatial, path, layer = layer, append = file.exists(path), quiet = TRUE)
    layers[[name]] <- layer
  }
  connection <- DBI::dbConnect(RSQLite::SQLite(), path)
  on.exit(DBI::dbDisconnect(connection), add = TRUE)
  for (name in names(tables)) write_gpkg_attribute(connection, name, tables[[name]])
  layers
}

output_field_units <- function(fields) {
  units <- stats::setNames(rep("unitless", length(fields)), fields)
  units[intersect(fields, c("x", "y", "z", "dbh_m", "height_m", "mean_dbh_m", "median_dbh_m",
                            "dbh_p05_m", "dbh_p25_m", "dbh_p75_m", "dbh_p95_m", "mean_height_m",
                            "median_height_m", "height_p05_m", "height_p25_m", "height_p75_m",
                            "height_p95_m"))] <- "m"
  units[intersect(fields, c("area_m2", "convex_hull_area_m2"))] <- "m2"
  units[intersect(fields, "basal_area_m2_per_ha")] <- "m2/ha"
  units[intersect(fields, c("longitude", "latitude"))] <- "degree"
  units[intersect(fields, c("trees_per_ha"))] <- "trees/ha"
  as.list(units)
}

write_outputs <- function(inventories, area, stage_dir, config, crs) {
  validate_inventory_crs(crs)
  lookup <- read_species_lookup(config$species_lookup)
  inventories <- lapply(inventories, function(inventory) {
    attach_species_metadata(with_geographic_coordinates(inventory, crs), lookup)
  })
  inventory_dir <- file.path(stage_dir, "inventories")
  dir.create(inventory_dir, recursive = TRUE)
  combined <- data.table::rbindlist(inventories, use.names = TRUE, fill = TRUE)
  summaries <- data.table::rbindlist(lapply(inventories, stand_summary_for_inventory, area = area), fill = TRUE)
  composition <- species_composition(inventories)
  if (is.null(composition)) {
    composition <- data.table::data.table(
      segmentation_source = character(), instance_dimension = character(),
      species_id = numeric(), tree_count = integer(), tree_proportion = numeric()
    )
  }
  if (!is.null(lookup) && nrow(composition)) {
    composition <- merge(composition, lookup, by = c("segmentation_source", "species_id"), all.x = TRUE, sort = FALSE)
  }
  validation <- data.table::rbindlist(lapply(inventories, function(inventory) {
    data.table::data.table(
      segmentation_source = unique(inventory$segmentation_source),
      instance_dimension = unique(inventory$instance_dimension),
      tree_count = nrow(inventory),
      dbh_unavailable_count = sum(is.na(inventory$dbh_m))
    )
  }))
  issues <- data.table::data.table(
    severity = character(), code = character(), segmentation_source = character(),
    instance_dimension = character(), instance_id = numeric(), message = character()
  )
  large <- if ("inventory_method" %in% names(combined)) {
    combined[inventory_method == "streaming_large_instance",
             .(segmentation_source, instance_dimension, instance_id, point_count,
               inventory_method, dbh_m)]
  } else {
    combined[0, .(segmentation_source, instance_dimension, instance_id, point_count,
                  inventory_method = character(), dbh_m = numeric())]
  }
  if (nrow(large)) {
    issues <- data.table::rbindlist(list(issues, large[, .(
      severity = "info", code = inventory_method,
      segmentation_source, instance_dimension, instance_id,
      message = sprintf(
        "%.0f points; species probability and FM score use a deterministic sample of at most %d points%s",
        point_count, inventory_large_instance_sample_limit,
        ifelse(is.na(dbh_m), "; x/y from that sample because no DBH could be fitted", "")
      )
    )]))
  }
  thinned <- if ("dbh_band_voxel_m" %in% names(combined)) {
    combined[is.finite(dbh_band_voxel_m) & dbh_band_voxel_m > 0]
  } else combined[0]
  if (nrow(thinned)) {
    issues <- data.table::rbindlist(list(issues, thinned[, .(
      severity = "info", code = "dbh_band_thinned",
      segmentation_source, instance_dimension, instance_id,
      message = sprintf(
        "Stem band kept %.0f of %.0f points, one per %.0f mm voxel, to fit the per-tree DBH budget",
        dbh_band_points, dbh_band_source_points, dbh_band_voxel_m * 1000
      )
    )]))
  }
  if (nrow(composition) && (is.null(lookup) || anyNA(composition$scientific_name))) {
    issues <- data.table::rbindlist(list(issues, data.table::data.table(
      severity = "warning", code = "species_lookup_missing", segmentation_source = "",
      instance_dimension = "", instance_id = NA_real_,
      message = "Supply --species-lookup to publish scientific names and taxonomy/model versions"
    )))
  }
  tables <- list(
    inventory_combined = combined, species = composition,
    stand_summary = summaries, validation = validation, issues = issues
  )
  # One row order for every format: by instance (species rows by species_id).
  inventory_order <- lapply(inventories, sort_by_instance)
  table_order <- lapply(names(tables), function(name) {
    sort_by_instance(tables[[name]], if (name == "species") "species_id")
  })
  names(table_order) <- names(tables)
  artifacts <- list()
  for (format in selected_table_formats(config)) {
    for (name in names(inventories)) {
      path <- file.path("inventories", paste0(name, ".", format))
      write_table_format(inventories[[name]], file.path(stage_dir, path), format)
      artifacts[[length(artifacts) + 1L]] <- table_artifact(path, format, "inventory", inventory_order[[name]],
                                                            segmentation = name)
    }
    for (name in names(tables)) {
      path <- paste0(name, ".", format)
      write_table_format(tables[[name]], file.path(stage_dir, path), format)
      artifacts[[length(artifacts) + 1L]] <- table_artifact(path, format, name, table_order[[name]])
    }
  }
  layers <- list()
  if ("gpkg" %in% config$output_formats) {
    layers <- write_inventory_gpkg(file.path(stage_dir, "inventory.gpkg"), inventories, tables, crs)
    artifacts[[length(artifacts) + 1L]] <- list(key = "inventory.gpkg", format = "gpkg", layers = layers,
                                                 tables = names(tables))
  }
  list(inventories = inventories, tables = tables, artifacts = artifacts, layers = layers)
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
  required_tables <- c("inventory_combined", "species", "stand_summary", "validation", "issues")
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
  segmentation_slugs <- unique(unlist(lapply(manifest$artifacts, `[[`, "segmentation")))
  parquet <- c(paste0(required_tables, ".parquet"),
               file.path("inventories", paste0(segmentation_slugs, ".parquet")))
  if (!length(segmentation_slugs) || !all(parquet %in% keys)) {
    abort("Inventory manifest is missing a Parquet table")
  }
  gpkg <- Filter(function(artifact) identical(artifact$format, "gpkg"), manifest$artifacts)
  if (length(gpkg)) {
    layers <- vapply(manifest$segmentations, function(spec) {
      if (!manifest_scalar(spec$layer)) abort("Inventory manifest is missing a GeoPackage layer")
      spec$layer
    }, character(1))
    if (anyDuplicated(layers) ||
        !all(c(layers, required_tables) %in% sf::st_layers(file.path(stage_dir, "inventory.gpkg"))$name)) {
      abort("Inventory manifest GeoPackage layers do not match the output")
    }
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
      layer = exports$layers[[unique(inventory$segmentation_source)]],
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
