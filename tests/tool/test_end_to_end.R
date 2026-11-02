#!/usr/bin/env Rscript

# Local, deterministic COPC integration test. No production service is used.
arguments <- commandArgs(trailingOnly = FALSE)
test_file <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[[1L]])
root <- normalizePath(file.path(dirname(test_file), "..", ".."), mustWork = TRUE)
source(file.path(root, "workflow", "cli.R"))
input <- file.path(root, "tests", "tool", "fixtures", "multi-instance.copc.laz")
runner <- file.path(root, "exec", "run.R")
work <- tempfile("csp-end-to-end-")
dir.create(work)
on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)

header <- lidR::readLASheader(input)
extra <- names(header@VLR$Extra_Bytes[["Extra Bytes Description"]])
stopifnot(match("PredInstance_ULS", extra) > 9L)
cloud <- suppressWarnings(lidR::readLAS(input, select = "*"))
stopifnot(nrow(cloud@data) == 29011L, sf::st_crs(cloud)$epsg == 32632L)

bounds <- header@PHB
dtm_path <- file.path(work, "provided.tif")
dtm <- terra::rast(
  xmin = floor(bounds[["Min X"]]) - 2, xmax = ceiling(bounds[["Max X"]]) + 2,
  ymin = floor(bounds[["Min Y"]]) - 2, ymax = ceiling(bounds[["Max Y"]]) + 2,
  resolution = 1, crs = sf::st_crs(cloud)$wkt
)
terra::values(dtm) <- 2
terra::writeRaster(dtm, dtm_path)

lookup_path <- file.path(work, "species.tsv")
data.table::fwrite(data.table::data.table(
  segmentation_source = "FM", species_id = 1:3,
  scientific_name = rep("Fagus sylvatica", 3L),
  common_name = rep("European beech", 3L),
  taxonomy_version = rep("fixture-taxonomy-1", 3L),
  model_version = rep("fixture-model-1", 3L)
), lookup_path, sep = "\t")

specifications <- c(
  "--instance-attribute", "PredInstance_SAT",
  "--instance-attribute", "PredInstance_FM",
  "--species-attribute", "species_id_FM", "--species-probability-attribute", "species_prob_FM",
  "--instance-attribute", "PredInstance_ULS"
)
dimensions <- c(SAT = "PredInstance_SAT", FM = "PredInstance_FM", ULS = "PredInstance_ULS")

# The same points as a plain (non-COPC) LAZ, as most platform inputs arrive.
plain_input <- file.path(work, "multi-instance.laz")
lidR::writeLAS(cloud, plain_input)
stopifnot(!input_is_copc(plain_input), input_is_copc(input))

# NULL leaves the reader at its default (auto).
run_case <- function(label, mode, csp = FALSE, chunk_size = 300,
                     inventory_read_mode = NULL, source = input) {
  output <- file.path(work, label)
  log <- file.path(work, paste0(label, ".log"))
  args <- c(runner, "--input", source, "--output-dir", output,
            specifications,
            "--species-lookup", lookup_path,
            "--inventory-attributes", "PredScore_FM,PredSemantic_FM",
            "--read-chunk-size", as.character(chunk_size),
            if (!is.null(inventory_read_mode)) c("--inventory-read-mode", inventory_read_mode),
            "--dtm-workers", "1")
  if (mode == "supplied") args <- c(args, "--dtm", dtm_path)
  if (csp) args <- c(args, "--enable-csp", "true")
  status <- suppressWarnings(system2("Rscript", args, stdout = log, stderr = log))
  if (!identical(status, 0L)) {
    stop(sprintf("%s failed:\n%s", label, paste(readLines(log, warn = FALSE), collapse = "\n")))
  }
  output
}

check_case <- function(output, mode, csp = FALSE, read_mode = NULL) {
  manifest <- jsonlite::fromJSON(file.path(output, "manifest.json"), simplifyVector = FALSE)
  metadata <- jsonlite::fromJSON(file.path(output, "run_metadata.json"))
  segmentations <- if (csp) c(dimensions, CSP = "PredInstance_CSP") else dimensions
  stopifnot(
    identical(manifest$schema_version, "1.0.0"),
    metadata$input_points == nrow(cloud@data),
    identical(unname(vapply(manifest$segmentations, `[[`, character(1), "technical_name")),
              unname(segmentations)),
    manifest$coordinate_reference$source_epsg == 32632L,
    manifest$coordinate_reference$wgs84_epsg == 4326L,
    identical(manifest$segmentations[[2L]]$potree_species_dimension, "species_id_FM"),
    is.null(read_mode) || identical(metadata$read_mode, read_mode)
  )
  keys <- vapply(manifest$artifacts, `[[`, character(1), "key")
  actual <- setdiff(list.files(output, recursive = TRUE), "manifest.json")
  stopifnot(!anyDuplicated(keys), setequal(keys, actual),
            all(c("inventory.csv", "inventory.parquet", "stand_summary.csv", "dtm_full.tif") %in% actual))
  if (mode == "supplied") {
    stopifnot(identical(manifest$dtm$source, "supplied"),
              identical(unname(tools::md5sum(file.path(output, "dtm_full.tif"))),
                        unname(tools::md5sum(dtm_path))))
  } else {
    stopifnot(!identical(manifest$dtm$source, "supplied"),
              any(is.finite(terra::values(terra::rast(file.path(output, "dtm_full.tif"))))))
  }
  # One table in two files.
  inventory <- data.table::as.data.table(utils::read.csv(file.path(output, "inventory.csv"), na.strings = ""))
  parquet <- data.table::as.data.table(arrow::read_parquet(file.path(output, "inventory.parquet")))
  columns <- c("instance_dimension", "instance_id", "point_count", "x", "y", "z", "longitude",
               "latitude", "dbh_m", "dbh_fits", "dbh_slices", "crown_hull")
  same <- function(column) {
    if (is.character(parquet[[column]])) return(identical(inventory[[column]], parquet[[column]]))
    isTRUE(all.equal(as.numeric(inventory[[column]]), as.numeric(parquet[[column]]), tolerance = 1e-6))
  }
  stopifnot(identical(names(inventory), names(parquet)), all(vapply(columns, same, logical(1))))
  summary <- utils::read.csv(file.path(output, "stand_summary.csv"))
  stopifnot(setequal(summary$instance_dimension, segmentations),
            identical(summary$tree_count[match(segmentations, summary$instance_dimension)],
                      vapply(segmentations, function(d) sum(inventory$instance_dimension == d), integer(1),
                             USE.NAMES = FALSE)))
  output_cloud <- if (csp) suppressWarnings(lidR::readLAS(
    file.path(output, "segmented_csp.laz"), select = "*"
  )) else NULL
  for (slug in names(segmentations)) {
    rows <- inventory[segmentation_source == slug]
    source_data <- if (slug == "CSP") output_cloud@data else cloud@data
    expected <- as.data.frame(table(source_data[[segmentations[[slug]]]]))
    expected <- expected[as.numeric(as.character(expected$Var1)) > 0, ]
    stopifnot(
      identical(as.numeric(rows$instance_id), as.numeric(as.character(expected$Var1))),
      identical(as.numeric(rows$point_count), as.numeric(expected$Freq)),
      all(is.finite(rows$longitude)), all(is.finite(rows$latitude)),
      all(grepl("^POLYGON", rows$crown_hull)),
      all(vapply(rows$dbh_fits, function(text) is.data.frame(jsonlite::fromJSON(text)), logical(1))),
      !any(c("inventory_quality_code", "inventory_quality_label") %in% names(rows))
    )
    if (slug == "FM") {
      stopifnot(all(rows$scientific_name == "Fagus sylvatica"),
                all(rows$taxonomy_version == "fixture-taxonomy-1"))
    }
  }
  stopifnot(sum(inventory$segmentation_source == "SAT") == 2L)
  if (csp) {
    stopifnot(nrow(output_cloud@data) == nrow(cloud@data),
              all(c(unname(dimensions), "PredInstance_CSP") %in% names(output_cloud@data)))
    for (dimension in dimensions) {
      stopifnot(identical(output_cloud@data[[dimension]], cloud@data[[dimension]]))
    }
  } else {
    if (identical(metadata$read_mode, "spatial_chunks_and_instance_partitions")) {
      stopifnot(identical(metadata$input_read_selectors, "xyz0"))
    }
    stopifnot(all(paste0("unused_", 1:4) %in% metadata$input_dimensions_dropped))
    stopifnot(!file.exists(file.path(output, "segmented_csp.laz")))
  }
  invisible(TRUE)
}

spatial <- run_case("supplied", "supplied", chunk_size = 10, inventory_read_mode = "spatial")
check_case(spatial, "supplied", read_mode = "spatial_chunks_and_instance_partitions")
cat("passed supplied (spatial reader)\n")

# Both readers must feed the same point order into the serial inventory so
# the bounded scan cannot change stochastic DBH measurements.
sequential <- run_case("supplied-laz", "supplied", chunk_size = 10, source = plain_input)
check_case(sequential, "supplied", read_mode = "sequential_scan_and_instance_partitions")
sequential_metadata <- jsonlite::fromJSON(file.path(sequential, "run_metadata.json"))
stopifnot(sequential_metadata$read_chunk_size_m == 10,
          sequential_metadata$spatial_bucket_count > 0)
stopifnot(identical(unname(tools::md5sum(file.path(spatial, "inventory.parquet"))),
                    unname(tools::md5sum(file.path(sequential, "inventory.parquet")))))
cat("passed supplied-laz (auto: sequential reader, same inventory as spatial)\n")

output <- run_case("generated", "generated")
check_case(output, "generated", read_mode = "copc_owned_windows")
cat("passed generated (auto: owned COPC reader)\n")

output <- run_case("generated-laz", "generated", source = plain_input)
check_case(output, "generated", read_mode = "sequential_scan_and_instance_partitions")
cat("passed generated-laz (auto: sequential reader)\n")

for (mode in c("supplied", "generated")) {
  label <- paste0(mode, "-csp")
  output <- run_case(label, mode, csp = TRUE)
  check_case(output, mode, csp = TRUE)
  cat(sprintf("passed %s\n", label))
}

cat("multi-instance COPC end-to-end tests passed\n")
