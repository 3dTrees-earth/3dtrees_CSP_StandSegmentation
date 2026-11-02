# DBH evidence (--dbh-evidence true): the per-slice circles of the reported run,
# every run's DBH, and the crown hull polygon the crown area was measured from.
# Measurement code only records what it already computed; nothing here feeds
# back into an inventory value.
#
# Instances are measured in forked workers, so each process appends its rows to
# its own file under the evidence directory; tree_evidence() merges them into
# one row per tree for the inventory table.

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

# Column types are fixed so a column that happens to be all-NA in one run (no
# spline, no fitted run) still writes numbers, not strings, into the JSON.
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

# One JSON object per row, built column by column; NA and non-finite values
# become null.
json_objects <- function(rows, columns) {
  fields <- lapply(columns, function(name) {
    value <- rows[[name]]
    text <- if (is.character(value)) {
      ifelse(is.na(value), "null", encodeString(value, quote = '"'))
    } else if (is.logical(value)) {
      ifelse(is.na(value), "null", tolower(as.character(value)))
    } else {
      ifelse(is.finite(value), as.character(value), "null")
    }
    paste0('"', name, '":', text)
  })
  paste0("{", do.call(paste, c(fields, sep = ",")), "}")
}

# A JSON array per tree from its rows, in `order`.
json_per_tree <- function(rows, order, column) {
  data.table::setorderv(rows, c(instance_sort_keys, order))
  columns <- setdiff(names(rows), instance_sort_keys)
  rows[, object := json_objects(.SD, columns), .SDcols = columns]
  rows[, stats::setNames(list(paste0("[", paste(object, collapse = ","), "]")), column),
       by = instance_sort_keys]
}

# The evidence of every inventoried tree as inventory columns: `dbh_fits` and
# `dbh_slices` (JSON arrays) and `crown_hull` (WKT). A tree measured twice
# (large-instance remeasurement) keeps its last rows; the streaming hull of a
# large instance wins over its stem-band hull.
tree_evidence <- function(config, keys) {
  if (!evidence_enabled(config)) return(NULL)
  keys <- unique(keys[, .(instance_dimension, instance_id = as.numeric(instance_id))])
  read <- function(name, unique_by) {
    rows <- read_evidence(config$dbh_evidence_dir, name)
    if (is.null(rows)) return(NULL)
    rows <- rows[keys, on = instance_sort_keys, nomatch = NULL]
    unique(rows, by = c(instance_sort_keys, unique_by), fromLast = TRUE)
  }
  evidence <- keys
  runs <- read("dbh_runs", "run")
  if (!is.null(runs)) evidence <- merge(evidence, json_per_tree(runs, "run", "dbh_fits"), all.x = TRUE)
  slices <- read("dbh_slices", "slice_id")
  if (!is.null(slices)) {
    if (length(setdiff(unique(slices$status), dbh_slice_statuses))) abort("Unknown DBH slice status in evidence")
    evidence <- merge(evidence, json_per_tree(slices[, !"run"], "slice_id", "dbh_slices"), all.x = TRUE)
  }
  hulls <- read_evidence(config$dbh_evidence_dir, "crown_hulls")
  if (!is.null(hulls)) {
    hulls <- hulls[keys, on = instance_sort_keys, nomatch = NULL]
    data.table::setorderv(hulls, c(instance_sort_keys, "priority"), c(1L, 1L, -1L))
    hulls <- unique(hulls, by = instance_sort_keys)
    evidence <- merge(evidence, hulls[, .(instance_dimension, instance_id, crown_hull = wkt)], all.x = TRUE)
  }
  unlink(config$dbh_evidence_dir, recursive = TRUE)
  evidence
}
