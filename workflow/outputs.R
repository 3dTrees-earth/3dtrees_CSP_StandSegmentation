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

write_table_format <- function(table, path, format) {
  if (format == "csv") data.table::fwrite(table, path, sep = ",", na = "")
  else if (format == "tsv") data.table::fwrite(table, path, sep = "\t", na = "")
  else if (format == "json") {
    jsonlite::write_json(as.data.frame(table), path, dataframe = "rows", auto_unbox = TRUE, na = "null", pretty = TRUE)
  } else abort(sprintf("Unsupported table format: %s", format))
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
  artifacts <- list()
  for (format in setdiff(config$output_formats, "gpkg")) {
    for (name in names(inventories)) {
      path <- file.path("inventories", paste0(name, ".", format))
      write_table_format(inventories[[name]], file.path(stage_dir, path), format)
      artifacts[[length(artifacts) + 1L]] <- list(key = path, format = format, table = "inventory", segmentation = name)
    }
    for (name in names(tables)) {
      path <- paste0(name, ".", format)
      write_table_format(tables[[name]], file.path(stage_dir, path), format)
      artifacts[[length(artifacts) + 1L]] <- list(key = path, format = format, table = name)
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

validate_inventory_manifest <- function(manifest, stage_dir) {
  if (!identical(manifest$schema_version, inventory_schema_version)) {
    abort("Inventory manifest schema version is invalid")
  }
  required_tables <- c("inventory_combined", "species", "stand_summary", "validation", "issues")
  if (!all(required_tables %in% names(manifest$tables)) || !length(manifest$segmentations)) {
    abort("Inventory manifest is missing required tables or segmentations")
  }
  if (any(vapply(manifest$segmentations, function(spec) {
    !length(spec$technical_name) || !length(spec$display_name) || !length(spec$fields)
  }, logical(1)))) abort("Inventory manifest has an incomplete segmentation")
  keys <- vapply(manifest$artifacts, `[[`, character(1), "key")
  if (anyDuplicated(keys) || any(!file.exists(file.path(stage_dir, keys))) ||
      !"dtm_full.tif" %in% keys) {
    abort("Inventory manifest contains duplicate or missing artifact keys")
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
    ), lapply(intersect(c("dtm_aoi.tif", "effective_seeds.tsv", "segmented_csp.laz"),
                        list.files(stage_dir)), function(key) list(key = key, format = tools::file_ext(key))))
  )
  validate_inventory_manifest(manifest, stage_dir)
  jsonlite::write_json(manifest, file.path(stage_dir, "manifest.json"), auto_unbox = TRUE,
                       pretty = TRUE, na = "null", null = "null")
  invisible(manifest)
}
