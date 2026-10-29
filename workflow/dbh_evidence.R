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

# Column types are fixed so a table whose column happens to be all-NA in one run
# (no spline, no fitted run) keeps the same Parquet schema as every other run.
dbh_evidence_column_classes <- list(
  character = c("instance_dimension", "status", "wkt"),
  numeric = c("instance_id", "z_lo_m", "z_hi_m", "ground_z", "circle_x", "circle_y", "radius_m",
              "spline_x", "spline_y", "spline_r", "dbh_m"),
  integer = c("run", "slice_id", "n_points", "n_fit_points", "seed", "quality_code", "priority",
              "vertex_count"),
  logical = "reported"
)

read_evidence <- function(directory, table_name) {
  files <- sort(list.files(directory, pattern = sprintf("^%s-[0-9]+[.]csv$", table_name), full.names = TRUE))
  if (!length(files)) return(NULL)
  data.table::rbindlist(lapply(files, function(path) {
    header <- names(data.table::fread(path, nrows = 0L))
    classes <- lapply(dbh_evidence_column_classes, intersect, header)
    data.table::fread(path, na.strings = "", colClasses = classes[lengths(classes) > 0L])
  }), use.names = TRUE, fill = TRUE)
}

# Keep rows of inventoried trees only (the AOI filter runs after measurement),
# one row per key, sorted by (instance_dimension, instance_id, ...). Parquet
# (GeoParquet for the hulls) is always written; CSV/TSV/JSON copies follow
# --output-formats.
write_dbh_evidence <- function(stage_dir, inventories, config, crs) {
  if (!evidence_enabled(config)) return(list())
  keys <- unique(data.table::rbindlist(lapply(inventories, function(inventory)
    inventory[, .(instance_dimension, instance_id = as.numeric(instance_id))])))
  out_dir <- file.path(stage_dir, "dbh")
  dir.create(out_dir, showWarnings = FALSE)
  layouts <- list(
    dbh_slices = c("instance_dimension", "instance_id", "slice_id"),
    dbh_runs = c("instance_dimension", "instance_id", "run"),
    crown_hulls = c("instance_dimension", "instance_id")
  )
  artifacts <- list()
  for (name in names(layouts)) {
    rows <- read_evidence(config$dbh_evidence_dir, name)
    if (is.null(rows)) next
    rows <- rows[keys, on = c("instance_dimension", "instance_id"), nomatch = NULL]
    if (name == "crown_hulls") {
      data.table::setorderv(rows, c("instance_dimension", "instance_id", "priority"), c(1L, 1L, -1L))
      rows <- rows[, .SD[1L], by = c("instance_dimension", "instance_id")][, !"priority"]
    }
    # A tree measured twice (large-instance remeasurement) keeps its last rows.
    rows <- unique(rows, by = layouts[[name]], fromLast = TRUE)
    sorted_by <- sort_by_instance(rows, layouts[[name]])
    if (name == "dbh_slices" && length(setdiff(unique(rows$status), dbh_slice_statuses))) {
      abort("Unknown DBH slice status in evidence")
    }
    for (format in selected_table_formats(config)) {
      key <- file.path("dbh", paste0(name, ".", format))
      write_table_format(rows, file.path(stage_dir, key), format,
                         if (name == "crown_hulls" && format == "parquet") crs)
      artifacts[[length(artifacts) + 1L]] <- table_artifact(key, format, name, sorted_by)
    }
  }
  unlink(config$dbh_evidence_dir, recursive = TRUE)
  artifacts
}
