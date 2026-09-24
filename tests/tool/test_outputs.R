#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = FALSE)
test_file <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[[1L]])
root <- normalizePath(file.path(dirname(test_file), "..", ".."), mustWork = TRUE)
source(file.path(root, "workflow", "cli.R"))
source(file.path(root, "workflow", "tool.R"))
source(file.path(root, "workflow", "outputs.R"))

expect_error <- function(expression, pattern) {
  message <- tryCatch({ force(expression); NULL }, error = conditionMessage)
  stopifnot(!is.null(message), grepl(pattern, message))
}

work <- tempfile("csp-outputs-")
dir.create(work)
crs <- sf::st_crs(32632)
bounds <- list(xmin = 500000, xmax = 500100, ymin = 5300000, ymax = 5300100)
compound <- sf::st_crs("EPSG:32632+5773")
stopifnot(same_horizontal_crs(compound, crs, bounds),
          !same_horizontal_crs(compound, sf::st_crs(32633), bounds))
provided <- file.path(work, "provided.tif")
raster <- terra::rast(xmin = 499900, xmax = 500200, ymin = 5299900, ymax = 5300200,
                      resolution = 10, crs = crs$wkt)
terra::values(raster) <- 100
terra::writeRaster(raster, provided)
stage <- file.path(work, "stage")
dir.create(stage)
dtm <- load_supplied_dtm(provided, crs,
                         bounds, stage)
stopifnot(identical(dtm$method, "supplied"), file.exists(file.path(stage, "dtm_full.tif")))
expect_error(load_supplied_dtm(provided, sf::st_crs(4326),
                               list(xmin = 500000, xmax = 500100, ymin = 5300000, ymax = 5300100),
                               stage), "same valid CRS")
expect_error(load_supplied_dtm(provided, crs,
                               list(xmin = 500000, xmax = 500300, ymin = 5300000, ymax = 5300100),
                               stage), "cover the full")

lookup_path <- file.path(work, "species.tsv")
data.table::fwrite(data.table::data.table(
  segmentation_source = "SAT", species_id = 10,
  scientific_name = "Fagus sylvatica", common_name = "European beech",
  taxonomy_version = "taxonomy-v1", model_version = "model-v1"
), lookup_path, sep = "\t")
inventory <- data.table::data.table(
  segmentation_source = "SAT", instance_dimension = "PredInstance_SAT",
  instance_id = c(1, 2), x = c(500000, 500010), y = c(5300000, 5300010),
  z = c(100, 101), dbh_m = c(0.3, NA_real_), height_m = c(20, 18),
  convex_hull_area_m2 = c(4, 5), point_count = c(20L, 30L),
  species_id = c(10, 10), species_prob = c(0.8, 0.9),
  dbh_quality_code = c(4L, 9L), dbh_quality_label = c("spline_fit", "dbh_unavailable")
)
config <- list(
  output_formats = c("csv", "tsv", "json", "gpkg"), species_lookup = lookup_path,
  segmentation_specs = list(list(instance = "PredInstance_SAT", source = "SAT",
                                 species = "species_id_SAT", species_prob = "species_prob_SAT"))
)
exports <- write_outputs(list(SAT = inventory), list(area_m2 = 10000, area_source = "aoi"),
                         stage, config, crs)
stopifnot(
  file.exists(file.path(stage, "inventory_combined.csv")),
  file.exists(file.path(stage, "inventory_combined.tsv")),
  file.exists(file.path(stage, "inventory_combined.json")),
  file.exists(file.path(stage, "inventory.gpkg")),
  identical(exports$inventories$SAT$scientific_name, rep("Fagus sylvatica", 2L)),
  all(is.finite(exports$inventories$SAT$longitude)),
  all(is.finite(exports$inventories$SAT$latitude))
)
gpkg <- file.path(stage, "inventory.gpkg")
stopifnot(all(c("inventory_SAT", "inventory_combined", "species", "stand_summary",
                "validation", "issues") %in% sf::st_layers(gpkg)$name))
spatial <- sf::st_read(gpkg, layer = "inventory_SAT", quiet = TRUE)
stopifnot(nrow(spatial) == 2L, ncol(sf::st_coordinates(spatial)) == 3L)
json <- jsonlite::fromJSON(file.path(stage, "inventory_combined.json"))
stopifnot(nrow(json) == 2L, all(c("x", "y", "z", "longitude", "latitude",
                                "dbh_quality_code", "dbh_quality_label") %in% names(json)))
writeLines("{}", file.path(stage, "run_metadata.json"))
writeLines("{}", file.path(stage, "resource_summary.json"))
written <- write_inventory_manifest(exports, config, crs, dtm, stage)
manifest <- jsonlite::fromJSON(file.path(stage, "manifest.json"), simplifyVector = FALSE)
stopifnot(identical(manifest$schema_version, "1.0.0"),
          identical(manifest$dtm$source, "supplied"),
          identical(manifest$segmentations[[1L]]$potree_instance_dimension, "PredInstance_SAT"))
invalid <- written
invalid$schema_version <- "0.0.0"
expect_error(validate_inventory_manifest(invalid, stage), "schema version")
invalid <- written
invalid$segmentations[[1L]]$potree_instance_dimension <- "missing_dimension"
expect_error(validate_inventory_manifest(invalid, stage), "incomplete segmentation")
invalid <- written
invalid$segmentations[[1L]]$technical_name <- "missing_dimension"
invalid$segmentations[[1L]]$potree_instance_dimension <- "missing_dimension"
expect_error(validate_inventory_manifest(invalid, stage, "PredInstance_SAT"), "do not match")
invalid <- written
invalid$coordinate_reference$wgs84_epsg <- 3857L
expect_error(validate_inventory_manifest(invalid, stage), "unsafe coordinate metadata")
invalid <- written
invalid$coordinate_reference$source_epsg <- 32633L
expect_error(validate_inventory_manifest(invalid, stage), "source EPSG")
invalid <- written
invalid$artifacts <- invalid$artifacts[-1L]
expect_error(validate_inventory_manifest(invalid, stage), "artifact keys")

unlink(work, recursive = TRUE, force = TRUE)
cat("output contract tests passed\n")
