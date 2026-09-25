#!/usr/bin/env Rscript

# Local, deterministic COPC integration test. No production service is used.
arguments <- commandArgs(trailingOnly = FALSE)
test_file <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[[1L]])
root <- normalizePath(file.path(dirname(test_file), "..", ".."), mustWork = TRUE)
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
  "PredInstance_SAT",
  "PredInstance_FM,species_id_FM,species_prob_FM",
  "PredInstance_ULS"
)
dimensions <- c(SAT = "PredInstance_SAT", FM = "PredInstance_FM", ULS = "PredInstance_ULS")
formats <- c("csv", "tsv", "json", "gpkg")
tables <- c("inventory_combined", "species", "stand_summary", "validation", "issues")

run_case <- function(label, mode, selected, csp = FALSE, chunk_size = 300,
                     inventory_read_mode = "spatial") {
  output <- file.path(work, label)
  log <- file.path(work, paste0(label, ".log"))
  args <- c(runner, "--input", input, "--output-dir", output,
            as.vector(rbind(rep("--segmentation-spec", length(specifications)), specifications)),
            "--species-lookup", lookup_path,
            "--output-formats", paste(selected, collapse = ","),
            "--read-chunk-size", as.character(chunk_size),
            "--inventory-read-mode", inventory_read_mode,
            "--dtm-workers", "1")
  if (mode == "supplied") args <- c(args, "--dtm", dtm_path)
  if (csp) args <- c(args, "--enable-csp", "true")
  status <- suppressWarnings(system2("Rscript", args, stdout = log, stderr = log))
  if (!identical(status, 0L)) {
    stop(sprintf("%s failed:\n%s", label, paste(readLines(log, warn = FALSE), collapse = "\n")))
  }
  output
}

read_inventory <- function(output, slug, format) {
  if (format == "csv") return(data.table::fread(file.path(output, "inventories", paste0(slug, ".csv"))))
  if (format == "tsv") return(data.table::fread(file.path(output, "inventories", paste0(slug, ".tsv"))))
  if (format == "json") return(data.table::as.data.table(jsonlite::fromJSON(
    file.path(output, "inventories", paste0(slug, ".json"))
  )))
  data.table::as.data.table(sf::st_drop_geometry(sf::st_read(
    file.path(output, "inventory.gpkg"), layer = paste0("inventory_", slug), quiet = TRUE
  )))
}

read_table <- function(output, table, format) {
  if (format == "csv") return(data.table::fread(file.path(output, paste0(table, ".csv"))))
  if (format == "tsv") return(data.table::fread(file.path(output, paste0(table, ".tsv"))))
  if (format == "json") return(data.table::as.data.table(jsonlite::fromJSON(
    file.path(output, paste0(table, ".json"))
  )))
  connection <- DBI::dbConnect(RSQLite::SQLite(), file.path(output, "inventory.gpkg"))
  on.exit(DBI::dbDisconnect(connection), add = TRUE)
  data.table::as.data.table(DBI::dbReadTable(connection, table))
}

check_case <- function(output, mode, selected, csp = FALSE) {
  manifest <- jsonlite::fromJSON(file.path(output, "manifest.json"), simplifyVector = FALSE)
  metadata <- jsonlite::fromJSON(file.path(output, "run_metadata.json"))
  segmentations <- if (csp) c(dimensions, CSP = "PredInstance_CSP") else dimensions
  stopifnot(
    identical(manifest$schema_version, "1.0.0"),
    identical(unname(metadata$output_formats), selected),
    metadata$input_points == nrow(cloud@data),
    identical(unname(vapply(manifest$segmentations, `[[`, character(1), "technical_name")),
              unname(segmentations)),
    manifest$coordinate_reference$source_epsg == 32632L,
    manifest$coordinate_reference$wgs84_epsg == 4326L,
    identical(manifest$segmentations[[2L]]$potree_species_dimension, "species_id_FM"),
    file.exists(file.path(output, "dtm_full.tif"))
  )
  keys <- vapply(manifest$artifacts, `[[`, character(1), "key")
  actual <- setdiff(list.files(output, recursive = TRUE), "manifest.json")
  stopifnot(!anyDuplicated(keys), setequal(keys, actual))
  if (mode == "supplied") {
    stopifnot(identical(manifest$dtm$source, "supplied"),
              identical(unname(tools::md5sum(file.path(output, "dtm_full.tif"))),
                        unname(tools::md5sum(dtm_path))))
  } else {
    stopifnot(!identical(manifest$dtm$source, "supplied"),
              any(is.finite(terra::values(terra::rast(file.path(output, "dtm_full.tif"))))))
  }
  for (format in formats) {
    for (slug in names(segmentations)) {
      present <- file.exists(file.path(output, "inventories", paste0(slug, ".", format)))
      if (format != "gpkg") stopifnot(identical(present, format %in% selected))
    }
    if (format == "gpkg") {
      stopifnot(identical(file.exists(file.path(output, "inventory.gpkg")), "gpkg" %in% selected))
    } else {
      stopifnot(all(vapply(tables, function(table) {
        identical(file.exists(file.path(output, paste0(table, ".", format))), format %in% selected)
      }, logical(1))))
    }
  }
  if ("gpkg" %in% selected) {
    stopifnot(setequal(sf::st_layers(file.path(output, "inventory.gpkg"))$name,
                       c(paste0("inventory_", names(segmentations)), tables)))
  }
  for (table in tables) {
    per_format <- lapply(selected, function(format) read_table(output, table, format))
    expected_rows <- switch(table, inventory_combined = if (csp) 7L +
                              length(unique(read_inventory(output, "CSP", selected[[1L]])$instance_id)) else 7L,
                            species = 3L, stand_summary = length(segmentations),
                            validation = length(segmentations), issues = 0L)
    stopifnot(all(vapply(per_format, nrow, integer(1)) == expected_rows))
    if (table == "species") {
      stopifnot(all(vapply(per_format, function(rows) {
        all(rows$scientific_name == "Fagus sylvatica") &&
          all(rows$taxonomy_version == "fixture-taxonomy-1")
      }, logical(1))))
    }
    if (length(per_format) > 1L && table %in% c("inventory_combined", "species", "validation")) {
      columns <- switch(table,
                        inventory_combined = c("segmentation_source", "instance_id", "point_count", "x", "y", "longitude", "latitude"),
                        species = c("segmentation_source", "species_id", "tree_count", "tree_proportion"),
                        validation = c("segmentation_source", "tree_count", "dbh_unavailable_count"))
      selected_rows <- lapply(per_format, function(rows) {
        data.table::setorderv(rows, columns[1:2])
        rows[, ..columns]
      })
      for (rows in selected_rows[-1L]) {
        stopifnot(isTRUE(all.equal(selected_rows[[1L]], rows, tolerance = 1e-6,
                                   check.attributes = FALSE)))
      }
    }
  }
  output_cloud <- if (csp) suppressWarnings(lidR::readLAS(
    file.path(output, "segmented_csp.laz"), select = "*"
  )) else NULL
  canonical <- NULL
  for (slug in names(segmentations)) {
    per_format <- lapply(selected, function(format) read_inventory(output, slug, format))
    names(per_format) <- selected
    for (inventory in per_format) {
      source_data <- if (slug == "CSP") output_cloud@data else cloud@data
      expected <- as.data.frame(table(source_data[[segmentations[[slug]]]]))
      expected <- expected[as.numeric(as.character(expected$Var1)) > 0, ]
      data.table::setorder(inventory, instance_id)
      stopifnot(
        identical(as.numeric(inventory$instance_id), as.numeric(as.character(expected$Var1))),
        identical(as.numeric(inventory$point_count), as.numeric(expected$Freq)),
        all(c("x", "y", "z", "longitude", "latitude", "dbh_quality_code",
              "dbh_quality_label") %in% names(inventory)),
        all(is.finite(inventory$longitude)), all(is.finite(inventory$latitude)),
        !any(c("inventory_quality_code", "inventory_quality_label") %in% names(inventory))
      )
      if (slug == "FM") {
        stopifnot(all(inventory$scientific_name == "Fagus sylvatica"),
                  all(inventory$taxonomy_version == "fixture-taxonomy-1"))
      }
    }
    if (length(per_format) > 1L) {
      columns <- c("instance_id", "point_count", "x", "y", "z", "longitude", "latitude")
      reference <- per_format[[1L]][, ..columns]
      for (inventory in per_format[-1L]) {
        stopifnot(isTRUE(all.equal(reference, inventory[, ..columns], tolerance = 1e-6,
                                   check.attributes = FALSE)))
      }
    }
    if (slug == "SAT") canonical <- per_format[[1L]]
  }
  stopifnot(nrow(canonical) == 2L)
  if (csp) {
    stopifnot(nrow(output_cloud@data) == nrow(cloud@data),
              all(c(unname(dimensions), "PredInstance_CSP") %in% names(output_cloud@data)))
    for (dimension in dimensions) {
      stopifnot(identical(output_cloud@data[[dimension]], cloud@data[[dimension]]))
    }
  } else {
    stopifnot(identical(metadata$input_read_selectors, "c0"),
              all(paste0("unused_", 1:4) %in% metadata$input_dimensions_dropped))
    stopifnot(!file.exists(file.path(output, "segmented_csp.laz")))
  }
  invisible(TRUE)
}

for (mode in c("supplied", "generated")) {
  selections <- c(list(formats), as.list(formats))
  for (selected in selections) {
    label <- paste(mode, paste(selected, collapse = "-"), sep = "-")
    output <- run_case(label, mode, selected, chunk_size = if (length(selected) == 4L) 10 else 300)
    check_case(output, mode, selected)
    cat(sprintf("passed %s\n", label))
  }
  label <- paste0(mode, "-csp")
  output <- run_case(label, mode, formats, csp = TRUE)
  check_case(output, mode, formats, csp = TRUE)
  cat(sprintf("passed %s\n", label))
}

# Both readers must feed the same point order into the serial inventory so
# requesting the faster COPC scan cannot change stochastic DBH measurements.
spatial <- run_case("ordered-spatial", "supplied", formats, chunk_size = 10)
sequential <- run_case("ordered-sequential", "supplied", formats, chunk_size = 10,
                       inventory_read_mode = "sequential")
sequential_metadata <- jsonlite::fromJSON(file.path(sequential, "run_metadata.json"))
stopifnot(sequential_metadata$read_chunk_size_m == 10,
          sequential_metadata$spatial_bucket_count > 0)
data_files <- c(
  file.path("inventories", paste0(names(dimensions), ".csv")),
  paste0(tables, ".csv")
)
stopifnot(all(vapply(data_files, function(name) {
  identical(unname(tools::md5sum(file.path(spatial, name))),
            unname(tools::md5sum(file.path(sequential, name))))
}, logical(1))))

cat("multi-instance COPC end-to-end tests passed\n")
