#!/usr/bin/env Rscript

# Parquet is the canonical table format; CSV/TSV are copies of the same rows.
# For every table this test asserts that the Parquet file and the CSV hold the
# same values (rewriting the Parquet rows with the CSV writer reproduces the CSV
# byte for byte), that columns keep their order and types (double/int32/bool/
# string, NA as null), that rows are sorted by (instance_dimension, instance_id)
# in row groups of at most 16,384 rows, and that crown hulls are GeoParquet 1.1.

arguments <- commandArgs(trailingOnly = FALSE)
test_file <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[[1L]])
root <- normalizePath(file.path(dirname(test_file), "..", ".."), mustWork = TRUE)
source(file.path(root, "workflow", "cli.R"))
source(file.path(root, "workflow", "tool.R"))
source(file.path(root, "workflow", "outputs.R"))

work <- tempfile("csp-parquet-")
dir.create(work)
crs <- sf::st_crs(32632)
set.seed(3)

# Two segmentations, 35,000 trees in shuffled order: the combined inventory and
# the evidence tables span several row groups.
make_inventory <- function(source, count, offset) {
  ids <- sample(seq_len(count) * 3 + offset)
  method <- rep("standard", count)
  method[seq(1L, count, by = 5000L)] <- "streaming_large_instance"
  data.table::data.table(
    segmentation_source = source, instance_dimension = paste0("PredInstance_", source),
    instance_id = as.numeric(ids),
    x = 500000 + stats::runif(count, 0, 100), y = 5300000 + stats::runif(count, 0, 100),
    z = 400 + stats::runif(count, 0, 10),
    dbh_m = ifelse(stats::runif(count) < 0.1, NA_real_, stats::runif(count, 0.05, 0.9)),
    height_m = stats::runif(count, 2, 40), convex_hull_area_m2 = stats::runif(count, 1, 80),
    point_count = sample.int(100000L, count, replace = TRUE),
    species_id = ifelse(stats::runif(count) < 0.05, NA_real_, sample(c(1, 2, 7), count, replace = TRUE)),
    species_prob = stats::runif(count),
    dbh_quality_code = sample(c(4L, 9L, NA_integer_), count, replace = TRUE),
    dbh_quality_label = sample(c("spline_fit", "dbh_unavailable", ""), count, replace = TRUE),
    dbh_band_points = sample.int(1000000L, count, replace = TRUE),
    dbh_band_source_points = 2000000L,
    dbh_band_voxel_m = ifelse(stats::runif(count) < 0.001, 0.002, 0),
    inventory_method = method,
    pred_score_mixed = sample(c(TRUE, FALSE, NA), count, replace = TRUE),
    pred_score_fm = NA_real_
  )
}
inventories <- list(SAT = make_inventory("SAT", 20000L, 0), FM = make_inventory("FM", 15000L, 1))

# DBH evidence as the forked workers record it: one CSV per table and process.
evidence_dir <- file.path(work, "stage", ".dbh-evidence")
dir.create(evidence_dir, recursive = TRUE)
config <- list(output_formats = c("parquet", "csv", "tsv"), species_lookup = NULL,
               dbh_evidence_dir = evidence_dir,
               segmentation_specs = list(list(instance = "PredInstance_SAT", species = "species_id_SAT"),
                                         list(instance = "PredInstance_FM", species = "species_id_FM")))
trees <- data.table::rbindlist(lapply(inventories, function(inventory)
  inventory[, .(instance_dimension, instance_id, x, y)]))
trees <- trees[sample(nrow(trees))]
slices <- trees[rep(seq_len(nrow(trees)), each = 2L)]
slices[, `:=`(run = 3L, slice_id = rep(2:1, nrow(trees)), z_lo_m = 1.2, z_hi_m = 1.3,
              ground_z = 400 + stats::runif(.N), n_points = 50L, n_fit_points = 45L,
              circle_x = x, circle_y = y, radius_m = stats::runif(.N, 0.03, 0.5),
              status = sample(dbh_slice_statuses, .N, replace = TRUE),
              spline_x = NA_real_, spline_y = NA_real_, spline_r = NA_real_)]
slices[, c("x", "y") := NULL]
record_evidence(config, "dbh_slices", slices)
runs <- trees[rep(seq_len(nrow(trees)), each = 2L), .(instance_dimension, instance_id)]
runs[, `:=`(run = rep(2:1, nrow(trees)), seed = 42L, dbh_m = NA_real_,
            quality_code = 9L, reported = rep(c(FALSE, TRUE), nrow(trees)))]
record_evidence(config, "dbh_runs", runs)
record_evidence(config, "crown_hulls", trees[, .(
  instance_dimension, instance_id, priority = 1L, vertex_count = 4L,
  wkt = sprintf("POLYGON ((%1$.3f %2$.3f, %3$.3f %2$.3f, %3$.3f %4$.3f, %1$.3f %4$.3f, %1$.3f %2$.3f))",
                x - 1, y - 1, x + 1, y + 1)
)])

stage <- file.path(work, "stage")
exports <- write_outputs(inventories, list(area_m2 = 10000, area_source = "aoi"), stage, config, crs)
evidence <- write_dbh_evidence(stage, exports$inventories, config, crs)
stopifnot(!dir.exists(evidence_dir), length(evidence) == 9L)

tables <- c("inventory_combined", "stand_summary", "species", "issues", "validation",
            file.path("inventories", c("SAT", "FM")), file.path("dbh", c("dbh_slices", "dbh_runs", "crown_hulls")))
arrow_types <- c(character = "string", numeric = "double", integer = "int32", logical = "bool")

hull_wkt <- function(geometry) {
  vapply(sf::st_as_sfc(structure(as.list(geometry), class = "WKB"), EWKB = FALSE), function(polygon) {
    ring <- polygon[[1L]]
    sprintf("POLYGON ((%s))", paste(sprintf("%.3f %.3f", ring[, 1L], ring[, 2L]), collapse = ", "))
  }, character(1))
}

for (table in tables) {
  parquet_path <- file.path(stage, paste0(table, ".parquet"))
  csv_path <- file.path(stage, paste0(table, ".csv"))
  arrow_table <- arrow::read_parquet(parquet_path, as_data_frame = FALSE)
  rows <- as.data.frame(arrow_table)
  is_hulls <- table == "dbh/crown_hulls"
  if (is_hulls) {
    position <- match("geometry", names(rows))
    rows$geometry <- hull_wkt(rows$geometry)
    names(rows)[[position]] <- "wkt"
  }
  # Same column order as the CSV.
  csv_header <- names(data.table::fread(csv_path, nrows = 0L))
  stopifnot(identical(names(rows), csv_header))
  # Types: numbers as double/int32, booleans as bool, text as string; R reads
  # each back as the class it was written from.
  schema <- arrow_table$schema
  for (column in setdiff(names(rows), "wkt")) {
    stopifnot(identical(schema[[column]]$type$ToString(), unname(arrow_types[[class(rows[[column]])[[1L]]]])))
  }
  # Values: the CSV writer applied to the Parquet rows reproduces the CSV and TSV.
  for (format in c("csv", "tsv")) {
    rewritten <- file.path(work, paste0("rewritten.", format))
    write_table_format(rows, rewritten, format)
    stopifnot(identical(unname(tools::md5sum(rewritten)),
                        unname(tools::md5sum(file.path(stage, paste0(table, ".", format))))))
  }
  # NA is null in Parquet and an empty field in the CSV; read with the Parquet
  # classes, the CSV holds the same values (doubles to fwrite's 15 digits).
  classes <- vapply(rows, function(column) class(column)[[1L]], character(1))
  csv <- as.data.frame(data.table::fread(csv_path, na.strings = "", colClasses = unname(classes)))
  stopifnot(identical(lapply(csv, is.na), lapply(rows, is.na)))
  for (column in names(rows)) {
    if (is.double(rows[[column]])) {
      stopifnot(isTRUE(all.equal(csv[[column]], rows[[column]], tolerance = 1e-14)))
    } else stopifnot(identical(csv[[column]], rows[[column]]))
  }
  # Sorted by instance, in row groups of at most 16,384 rows that never share a key range.
  keys <- intersect(c("instance_dimension", "instance_id"), names(rows))
  if (length(keys) && nrow(rows) > 1L) {
    ordered <- do.call(order, c(unname(as.list(rows[keys])), na.last = FALSE, method = "radix"))
    stopifnot(identical(ordered, seq_len(nrow(rows))))
  }
  reader <- arrow::ParquetFileReader$create(parquet_path)
  stopifnot(reader$num_row_groups == max(1L, ceiling(nrow(rows) / parquet_row_group_rows)))
  if (reader$num_row_groups > 1L) {
    last <- NULL
    for (group in seq_len(reader$num_row_groups) - 1L) {
      chunk <- as.data.frame(reader$ReadRowGroup(group))
      stopifnot(nrow(chunk) <= parquet_row_group_rows)
      first_key <- paste(chunk$instance_dimension[[1L]], sprintf("%012.0f", chunk$instance_id[[1L]]))
      if (!is.null(last)) stopifnot(first_key >= last)
      last <- paste(chunk$instance_dimension[[nrow(chunk)]], sprintf("%012.0f", chunk$instance_id[[nrow(chunk)]]))
    }
  }
  # No R-only schema metadata; the hulls carry GeoParquet metadata.
  stopifnot(is.null(arrow_table$metadata$r), identical(!is.null(arrow_table$metadata$geo), is_hulls))
}

# Row groups and sorting are real on the large tables.
stopifnot(arrow::ParquetFileReader$create(file.path(stage, "inventory_combined.parquet"))$num_row_groups == 3L,
          arrow::ParquetFileReader$create(file.path(stage, "dbh", "dbh_slices.parquet"))$num_row_groups == 5L)
combined <- arrow::read_parquet(file.path(stage, "inventory_combined.parquet"))
stopifnot(identical(combined$instance_dimension[[1L]], "PredInstance_FM"),
          identical(combined$instance_id[1:2], c(4, 7)))
# A single tree is a key range: one row in the inventory, its slices adjacent.
slices_parquet <- arrow::read_parquet(file.path(stage, "dbh", "dbh_slices.parquet"))
tree <- which(slices_parquet$instance_dimension == "PredInstance_SAT" & slices_parquet$instance_id == 300)
stopifnot(length(tree) == 2L, diff(tree) == 1L, identical(slices_parquet$slice_id[tree], 1:2))
# Types survive all-NA columns: no spline and no fitted DBH stay double.
slices_schema <- arrow::read_parquet(file.path(stage, "dbh", "dbh_slices.parquet"), as_data_frame = FALSE)$schema
runs_schema <- arrow::read_parquet(file.path(stage, "dbh", "dbh_runs.parquet"), as_data_frame = FALSE)$schema
stopifnot(identical(slices_schema$spline_r$type$ToString(), "double"),
          identical(runs_schema$dbh_m$type$ToString(), "double"),
          identical(runs_schema$reported$type$ToString(), "bool"),
          identical(arrow::read_parquet(file.path(stage, "inventory_combined.parquet"),
                                        as_data_frame = FALSE)$schema$pred_score_fm$type$ToString(), "double"))

# GeoParquet 1.1: WKB geometry named "geometry", PROJJSON CRS, Polygon type, bbox.
hulls <- arrow::read_parquet(file.path(stage, "dbh", "crown_hulls.parquet"), as_data_frame = FALSE)
geo <- jsonlite::fromJSON(hulls$metadata$geo, simplifyVector = FALSE)
column <- geo$columns$geometry
geometry <- sf::st_as_sfc(structure(as.list(as.data.frame(hulls)$geometry), class = "WKB"), EWKB = FALSE)
stopifnot(identical(geo$version, "1.1.0"), identical(geo$primary_column, "geometry"),
          identical(column$encoding, "WKB"), identical(unlist(column$geometry_types), "Polygon"),
          sf::st_crs(as.character(jsonlite::toJSON(column$crs, auto_unbox = TRUE, digits = NA))) == crs,
          isTRUE(all.equal(unlist(column$bbox), unname(as.numeric(sf::st_bbox(geometry))))),
          identical(hulls$schema$geometry$type$ToString(), "binary"),
          identical(names(hulls)[[length(names(hulls))]], "geometry"))

# The manifest names every Parquet artifact with its order and row-group size.
parquet_artifacts <- Filter(function(artifact) identical(artifact$format, "parquet"),
                            c(exports$artifacts, evidence))
sorted_by <- lapply(parquet_artifacts, function(artifact) unlist(artifact$sorted_by))
names(sorted_by) <- vapply(parquet_artifacts, `[[`, character(1), "key")
stopifnot(length(parquet_artifacts) == length(tables),
          all(vapply(parquet_artifacts, function(artifact)
            identical(artifact$row_group_rows, parquet_row_group_rows), logical(1))),
          identical(sorted_by[["inventory_combined.parquet"]], c("instance_dimension", "instance_id")),
          identical(sorted_by[["inventories/SAT.parquet"]], c("instance_dimension", "instance_id")),
          identical(sorted_by[["species.parquet"]], c("instance_dimension", "species_id")),
          identical(sorted_by[["stand_summary.parquet"]], "instance_dimension"),
          identical(sorted_by[["dbh/dbh_slices.parquet"]], c("instance_dimension", "instance_id", "slice_id")))

# Default: Parquet only, no CSV copies.
default_stage <- file.path(work, "default")
dir.create(default_stage)
default_evidence <- file.path(default_stage, ".dbh-evidence")
dir.create(default_evidence)
default_config <- modifyList(config, list(output_formats = "parquet", dbh_evidence_dir = default_evidence))
record_evidence(default_config, "dbh_runs", runs)
written <- c(write_outputs(inventories, list(area_m2 = 10000, area_source = "aoi"), default_stage,
                           default_config, crs)$artifacts,
             write_dbh_evidence(default_stage, inventories, default_config, crs))
stopifnot(all(vapply(written, `[[`, character(1), "format") == "parquet"),
          setequal(list.files(default_stage, recursive = TRUE),
                   c(paste0(c("inventory_combined", "stand_summary", "species", "issues", "validation",
                              "inventories/SAT", "inventories/FM", "dbh/dbh_runs"), ".parquet"))))

unlink(work, recursive = TRUE, force = TRUE)
cat("parquet output tests passed\n")
