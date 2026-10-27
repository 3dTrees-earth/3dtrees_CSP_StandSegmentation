# DBH evidence (--dbh-evidence true): the per-slice circles of the reported run,
# every run's DBH, and the crown hull polygon the crown area was measured from.
# Measurement code only records what it already computed; nothing here feeds
# back into an inventory value, so the inventories stay byte-identical.
#
# Instances are measured in forked workers, so each process appends its rows to
# its own file under the evidence directory; write_dbh_evidence() merges them.

dbh_slice_statuses <- c("fitted", "outlier_2sd", "growth_1_3x", "radius_below_0_04m",
                        "too_few_points", "no_fit")

evidence_enabled <- function(config) !is.null(config$dbh_evidence_dir)

record_evidence <- function(config, table_name, rows) {
  if (!evidence_enabled(config) || is.null(rows) || !nrow(rows)) return(invisible())
  path <- file.path(config$dbh_evidence_dir, sprintf("%s-%d.csv", table_name, Sys.getpid()))
  data.table::fwrite(rows, path, append = file.exists(path), col.names = !file.exists(path))
  invisible()
}

# The crown hull as one closed ring in WKT. `priority` lets the streaming summary
# of a large instance replace the hull its stem-band remeasurement wrote.
record_crown_hull <- function(config, instance, id, x, y, priority) {
  if (!evidence_enabled(config) || length(x) < 3L) return(invisible())
  ring <- chull(x, y)
  ring <- c(ring, ring[[1L]])
  record_evidence(config, "crown_hulls", data.table::data.table(
    instance_dimension = instance, instance_id = as.numeric(id), priority = priority,
    vertex_count = length(ring) - 1L,
    wkt = sprintf("POLYGON ((%s))", paste(sprintf("%.3f %.3f", x[ring], y[ring]), collapse = ", "))
  ))
}

# Called once per measured tree with every run's upstream result (`chosen` 0:
# every run failed, so only the hull is recorded).
record_tree_evidence <- function(config, instance, id, seed, fits, slices, chosen, points) {
  if (!evidence_enabled(config)) return(invisible())
  if (nrow(fits)) record_evidence(config, "dbh_runs", data.table::data.table(
    instance_dimension = instance, instance_id = as.numeric(id), run = seq_len(nrow(fits)),
    seed = (seed + seq_len(nrow(fits)) - 1L) %% .Machine$integer.max,
    dbh_m = fits$DBH, quality_code = as.integer(fits$quality_flag),
    reported = seq_len(nrow(fits)) == chosen
  ))
  reported <- if (chosen > 0L) slices[[chosen]] else NULL
  if (!is.null(reported) && nrow(reported)) {
    record_evidence(config, "dbh_slices", cbind(
      data.table::data.table(instance_dimension = instance, instance_id = as.numeric(id), run = chosen),
      reported[, !"tree_id"]))
  }
  # Upstream measures the crown above 0.5 m; the workflow falls back to all points.
  crown <- is.finite(points$Z) & points$Z > 0.5
  if (sum(crown) >= 3L) record_crown_hull(config, instance, id, points$X[crown], points$Y[crown], 1L)
  else record_crown_hull(config, instance, id, points$X, points$Y, 1L)
}

read_evidence <- function(directory, table_name) {
  files <- sort(list.files(directory, pattern = sprintf("^%s-[0-9]+[.]csv$", table_name), full.names = TRUE))
  if (!length(files)) return(NULL)
  data.table::rbindlist(lapply(files, data.table::fread, na.strings = ""), use.names = TRUE, fill = TRUE)
}

# GeoParquet 1.1 (WKB geometry, PROJJSON CRS) for the hulls; GDAL here has no
# Parquet driver, so the "geo" schema metadata is written directly.
write_parquet_table <- function(table, path, crs = NULL) {
  if (is.null(crs)) return(arrow::write_parquet(table, path))
  geometry <- sf::st_as_sfc(table$wkt, crs = crs)
  frame <- as.data.frame(table[, !"wkt"])
  frame$geometry <- arrow::Array$create(unclass(sf::st_as_binary(geometry)), type = arrow::binary())
  arrow_table <- arrow::Table$create(frame)
  bbox <- as.numeric(sf::st_bbox(geometry))
  geo <- list(version = "1.1.0", primary_column = "geometry", columns = list(geometry = list(
    encoding = "WKB", geometry_types = list("Polygon"),
    crs = jsonlite::fromJSON(sf::st_crs(crs)$ProjJson, simplifyVector = FALSE),
    bbox = if (all(is.finite(bbox))) bbox else NULL)))
  arrow_table$metadata$geo <- as.character(jsonlite::toJSON(geo, auto_unbox = TRUE, digits = NA, null = "null"))
  arrow::write_parquet(arrow_table, path)
}

# Keep rows of inventoried trees only (the AOI filter runs after measurement),
# one row per key, in a stable order; write CSV for people and Parquet for tools.
write_dbh_evidence <- function(stage_dir, inventories, config, crs) {
  if (!evidence_enabled(config)) return(list())
  keys <- unique(data.table::rbindlist(lapply(inventories, function(inventory)
    inventory[, .(instance_dimension, instance_id = as.numeric(instance_id))])))
  out_dir <- file.path(stage_dir, "dbh")
  dir.create(out_dir, showWarnings = FALSE)
  layouts <- list(
    dbh_slices = list(order = c("instance_dimension", "instance_id", "slice_id"), unique = c("instance_dimension", "instance_id", "slice_id")),
    dbh_runs = list(order = c("instance_dimension", "instance_id", "run"), unique = c("instance_dimension", "instance_id", "run")),
    crown_hulls = list(order = c("instance_dimension", "instance_id"), unique = c("instance_dimension", "instance_id"))
  )
  artifacts <- list()
  for (name in names(layouts)) {
    rows <- read_evidence(config$dbh_evidence_dir, name)
    if (is.null(rows)) next
    rows[, instance_id := as.numeric(instance_id)]
    rows <- rows[keys, on = c("instance_dimension", "instance_id"), nomatch = NULL]
    if (name == "crown_hulls") {
      data.table::setorderv(rows, c("instance_dimension", "instance_id", "priority"), c(1L, 1L, -1L))
      rows <- rows[, .SD[1L], by = c("instance_dimension", "instance_id")][, !"priority"]
    }
    # A tree measured twice (large-instance remeasurement) keeps its last rows.
    rows <- unique(rows, by = layouts[[name]]$unique, fromLast = TRUE)
    data.table::setorderv(rows, layouts[[name]]$order)
    if (name == "dbh_slices" && length(setdiff(unique(rows$status), dbh_slice_statuses))) {
      abort("Unknown DBH slice status in evidence")
    }
    csv <- file.path("dbh", paste0(name, ".csv"))
    parquet <- file.path("dbh", paste0(name, ".parquet"))
    data.table::fwrite(rows, file.path(stage_dir, csv), na = "")
    write_parquet_table(rows, file.path(stage_dir, parquet), if (name == "crown_hulls") crs)
    artifacts <- c(artifacts, list(
      list(key = csv, format = "csv", table = name),
      list(key = parquet, format = "parquet", table = name)))
  }
  unlink(config$dbh_evidence_dir, recursive = TRUE)
  artifacts
}
