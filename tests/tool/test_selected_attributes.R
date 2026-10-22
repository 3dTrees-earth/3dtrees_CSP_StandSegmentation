#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub('^--file=', '', grep('^--file=', args, value = TRUE)[[1L]])), '..', '..'))
repository_root <- root
source(file.path(root, 'workflow', 'cli.R'))
source(file.path(root, 'workflow', 'tool.R'))
source(file.path(root, 'workflow', 'outputs.R'))
work <- tempfile('csp-attributes-')
dir.create(work)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
input <- file.path(root, 'tests/tool/fixtures/multi-instance.copc.laz')
cloud <- suppressWarnings(lidR::readLAS(input, select = '*'))
dtm <- terra::rast(xmin = floor(min(cloud$X))-1, xmax = ceiling(max(cloud$X))+1,
                   ymin = floor(min(cloud$Y))-1, ymax = ceiling(max(cloud$Y))+1,
                   resolution = 1, crs = sf::st_crs(cloud)$wkt)
terra::values(dtm) <- 2
dtm_path <- file.path(work, 'dtm.tif')
terra::writeRaster(dtm, dtm_path)
base <- c('--input', input, '--segmentation-spec', 'PredInstance_FM',
          '--dtm', dtm_path, '--output-formats', 'csv', '--dtm-workers', '1')
expected <- cloud@data[PredInstance_FM > 0, .(point_count = .N), by = PredInstance_FM]
data.table::setorder(expected, PredInstance_FM)
config <- parse_cli_args(c(base, '--output-dir', file.path(work, 'extra'),
                           '--inventory-read-mode', 'copc-owned',
                           '--inventory-attributes', 'PredScore_FM,PredSemantic_FM,unused_1,Intensity,ReturnNumber,Synthetic_flag'))
selected <- read_point_cloud(input, config$segmentation_specs)
stopifnot(setequal(names(selected@data), c('X', 'Y', 'Z', 'PredInstance_FM', config$inventory_attributes)),
          identical(selected@data$unused_1, cloud@data$unused_1))
run_tool(config)
metadata <- jsonlite::fromJSON(file.path(config$output_dir, 'run_metadata.json'))
stopifnot(setequal(metadata$input_dimensions, names(selected@data)),
          !'Intensity' %in% metadata$input_dimensions_dropped)
rows <- data.table::fread(file.path(config$output_dir, 'inventories/FM.csv'))
stopifnot(any(is.finite(rows$pred_score_fm)), any(is.finite(rows$wood_point_count)))
# Missing explicit attributes must fail; do not silently substitute NA.
missing <- config
missing$output_dir <- file.path(work, 'missing')
missing$segmentation_specs[[1L]]$attributes <- 'does_not_exist'
error <- tryCatch({run_tool(missing); NULL}, error = conditionMessage)
stopifnot(!is.null(error), grepl('does_not_exist', error), !dir.exists(missing$output_dir))
# Full CSP output preserves requested attributes only, with its new instance ID.
for (generated in c(FALSE, TRUE)) {
  csp <- config
  csp$inventory_read_mode <- 'spatial'
  csp$enable_csp <- TRUE
  if (generated) csp$dtm <- NULL
  csp$output_dir <- file.path(work, paste0('csp-', generated))
  run_tool(csp)
  output <- suppressWarnings(lidR::readLAS(file.path(csp$output_dir, 'segmented_csp.laz'), select = '*'))
  stopifnot(nrow(output@data) == nrow(cloud@data),
            'PredInstance_CSP' %in% names(output@data),
            !any(c('unused_2', 'unused_3', 'unused_4') %in% names(output@data)),
            identical(output@data$unused_1, cloud@data$unused_1))
}
for (mode in c('copc-owned', 'sequential', 'spatial')) {
  config <- parse_cli_args(c(base, '--output-dir', file.path(work, mode), '--inventory-read-mode', mode))
  run_tool(config)
  metadata <- jsonlite::fromJSON(file.path(config$output_dir, 'run_metadata.json'))
  stopifnot(setequal(metadata$input_dimensions, c('X', 'Y', 'Z', 'PredInstance_FM')))
  rows <- data.table::fread(file.path(config$output_dir, 'inventories/FM.csv'))
  stopifnot(identical(rows$instance_id, expected$PredInstance_FM),
            identical(rows$point_count, expected$point_count),
            all(is.na(rows$pred_score_fm)), all(is.na(rows$wood_point_count)))
}
cat('Explicit attribute selection integration passed\n')
