# Table file formats. Parquet is the canonical format of every table and is
# always written; CSV, TSV and JSON are optional copies of the same rows
# (--output-formats), and the 3Dtrees platform makes CSV/GeoJSON/GeoPackage at
# download time from the Parquet files.
#
# Every table is sorted by (instance_dimension, instance_id) where it has those
# columns before it is written, and Parquet row groups hold at most
# `parquet_row_group_rows` rows, so a reader can fetch one tree's rows by range
# using the row-group min/max statistics.

parquet_row_group_rows <- 16384L
table_text_formats <- c("csv", "tsv", "json")
instance_sort_keys <- c("instance_dimension", "instance_id")

# Sort in place by the instance keys the table has, then `extra` keys (such as
# slice_id or run); returns the keys used. data.table sorts strings in C locale
# and is stable, so equal keys keep their measured order.
sort_by_instance <- function(table, extra = character()) {
  keys <- unique(c(intersect(instance_sort_keys, names(table)), intersect(extra, names(table))))
  if (length(keys) && nrow(table)) data.table::setorderv(table, keys)
  keys
}

geoparquet_geometry_types <- c(POINT = "Point", LINESTRING = "LineString", POLYGON = "Polygon",
                               MULTIPOINT = "MultiPoint", MULTILINESTRING = "MultiLineString",
                               MULTIPOLYGON = "MultiPolygon")

# Parquet, or GeoParquet 1.1 when `crs` is given: the `wkt` column becomes a WKB
# `geometry` column in the same position, with the PROJJSON CRS, geometry types
# and bbox in the "geo" schema metadata. GDAL in the image has no Parquet
# driver, so the metadata is written directly. arrow's R-specific "r" schema
# metadata is dropped so the file reads the same in every Parquet reader.
write_parquet_table <- function(table, path, crs = NULL) {
  frame <- as.data.frame(table)
  geo <- NULL
  if (!is.null(crs)) {
    geometry <- sf::st_as_sfc(frame$wkt, crs = crs)
    names(frame)[[match("wkt", names(frame))]] <- "geometry"
    frame$geometry <- arrow::Array$create(unclass(sf::st_as_binary(geometry)), type = arrow::binary())
    bbox <- as.numeric(sf::st_bbox(geometry))
    types <- unname(geoparquet_geometry_types[unique(as.character(sf::st_geometry_type(geometry)))])
    geo <- list(version = "1.1.0", primary_column = "geometry", columns = list(geometry = list(
      encoding = "WKB", geometry_types = as.list(types),
      crs = jsonlite::fromJSON(sf::st_crs(crs)$ProjJson, simplifyVector = FALSE),
      bbox = if (length(geometry) && all(is.finite(bbox))) bbox else NULL)))
  }
  arrow_table <- arrow::Table$create(frame)
  arrow_table$metadata$r <- NULL
  if (!is.null(geo)) {
    arrow_table$metadata$geo <- as.character(jsonlite::toJSON(geo, auto_unbox = TRUE, digits = NA, null = "null"))
  }
  arrow::write_parquet(arrow_table, path, chunk_size = parquet_row_group_rows, compression = "snappy")
}

# CSV/TSV are written by data.table::fwrite (doubles to 15 significant digits,
# NA as an empty field, empty strings as ""); JSON keeps full double precision.
write_table_format <- function(table, path, format, crs = NULL) {
  if (format == "parquet") write_parquet_table(table, path, crs)
  else if (format == "csv") data.table::fwrite(table, path, sep = ",", na = "")
  else if (format == "tsv") data.table::fwrite(table, path, sep = "\t", na = "")
  else if (format == "json") {
    jsonlite::write_json(as.data.frame(table), path, dataframe = "rows", auto_unbox = TRUE,
                         na = "null", digits = NA, pretty = TRUE)
  } else abort(sprintf("Unsupported table format: %s", format))
}

# The table formats a run writes: always Parquet, then the selected copies.
selected_table_formats <- function(config) {
  c("parquet", intersect(table_text_formats, config$output_formats))
}

# Manifest entry; Parquet entries also state the row order and row-group size.
table_artifact <- function(key, format, table, sorted_by, ...) {
  artifact <- list(key = key, format = format, table = table, ...)
  if (format == "parquet") {
    artifact$sorted_by <- as.list(sorted_by)
    artifact$row_group_rows <- parquet_row_group_rows
  }
  artifact
}
