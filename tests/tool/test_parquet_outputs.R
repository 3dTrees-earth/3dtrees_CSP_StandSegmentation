#!/usr/bin/env Rscript

# The inventory is one table in two files: inventory.parquet (viewer) and
# inventory.csv (people). This test asserts that both hold the same rows and
# columns (rewriting the Parquet rows with the CSV writer reproduces the CSV
# byte for byte), that rows are sorted by (instance_dimension, instance_id) in
# row groups of at most 2,048 rows, that the evidence columns (JSON and WKT)
# come last and parse per tree, and that the stand summary is one CSV.

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
    pred_ok = sample(c(TRUE, FALSE, NA), count, replace = TRUE),
    all_na_value = NA_real_
  )
}
inventories <- list(SAT = make_inventory("SAT", 20000L, 0), FM = make_inventory("FM", 15000L, 1))

# DBH evidence as the forked workers record it: one CSV per table and process.
evidence_dir <- file.path(work, "stage", ".dbh-evidence")
dir.create(evidence_dir, recursive = TRUE)
config <- list(species_lookup = NULL,
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
stopifnot(setequal(list.files(stage, recursive = TRUE, all.files = TRUE),
                   c("inventory.csv", "inventory.parquet", "stand_summary.csv")),
          !dir.exists(evidence_dir))

# One table, two files: the Parquet rows written as CSV are the CSV.
parquet <- arrow::read_parquet(file.path(stage, "inventory.parquet"))
rewritten <- file.path(work, "rewritten.csv")
write_csv_table(data.table::as.data.table(parquet), rewritten)
stopifnot(identical(unname(tools::md5sum(rewritten)), unname(tools::md5sum(file.path(stage, "inventory.csv")))),
          nrow(parquet) == 35000L)
arrow_table <- arrow::read_parquet(file.path(stage, "inventory.parquet"), as_data_frame = FALSE)
stopifnot(is.null(arrow_table$metadata$r),
          identical(tail(names(parquet), 4L), c("diagnostics", "dbh_fits", "dbh_slices", "crown_hull")),
          identical(arrow_table$schema$all_na_value$type$ToString(), "double"),
          identical(arrow_table$schema$pred_ok$type$ToString(), "bool"),
          identical(arrow_table$schema$dbh_slices$type$ToString(), "string"))

# Sorted rows in 2,048-row groups: one tree's evidence is one small range read.
stopifnot(arrow::ParquetFileReader$create(file.path(stage, "inventory.parquet"))$num_row_groups == 18L,
          identical(parquet$instance_dimension[[1L]], "PredInstance_FM"),
          identical(parquet$instance_id[1:2], c(4, 7)),
          !is.unsorted(paste(parquet$instance_dimension, sprintf("%09.0f", parquet$instance_id))))
tree <- parquet[parquet$instance_dimension == "PredInstance_SAT" & parquet$instance_id == 300, ]
fits <- jsonlite::fromJSON(tree$dbh_fits)
tree_slices <- jsonlite::fromJSON(tree$dbh_slices)
stopifnot(nrow(tree) == 1L, identical(fits$run, 1:2), identical(fits$reported, c(TRUE, FALSE)),
          identical(tree_slices$slice_id, 1:2), !"run" %in% names(tree_slices),
          all(tree_slices$status %in% dbh_slice_statuses), all(is.na(tree_slices$spline_r)),
          grepl("^POLYGON \\(\\(", tree$crown_hull),
          isTRUE(all.equal(tree_slices$circle_x, rep(tree$x, 2L))))

# Diagnostics name streamed and thinned trees only.
streamed <- parquet$inventory_method == "streaming_large_instance"
stopifnot(all(grepl("^streaming_large_instance: ", parquet$diagnostics[streamed])),
          all(is.na(parquet$diagnostics[!streamed & parquet$dbh_band_voxel_m == 0])),
          all(grepl("dbh_band_thinned: ", parquet$diagnostics[parquet$dbh_band_voxel_m > 0])))

# The stand summary: one CSV row per segmentation.
summary <- utils::read.csv(file.path(stage, "stand_summary.csv"))
stopifnot(identical(summary$instance_dimension, c("PredInstance_FM", "PredInstance_SAT")),
          identical(summary$tree_count, c(15000L, 20000L)))
artifacts <- vapply(exports$artifacts, `[[`, character(1), "key")
stopifnot(identical(artifacts, c("inventory.csv", "inventory.parquet", "stand_summary.csv")),
          identical(exports$artifacts[[2L]]$row_group_rows, 2048L),
          identical(unlist(exports$artifacts[[2L]]$sorted_by), c("instance_dimension", "instance_id")))

# --dbh-evidence false: the same columns, empty.
plain_stage <- file.path(work, "plain")
dir.create(plain_stage)
plain <- write_outputs(inventories, list(area_m2 = 10000, area_source = "aoi"), plain_stage,
                       modifyList(config, list(dbh_evidence_dir = NULL)), crs)
stopifnot(identical(names(plain$tables$inventory), names(parquet)),
          all(is.na(plain$tables$inventory$dbh_fits)), all(is.na(plain$tables$inventory$crown_hull)))

unlink(work, recursive = TRUE, force = TRUE)
cat("inventory table tests passed\n")
