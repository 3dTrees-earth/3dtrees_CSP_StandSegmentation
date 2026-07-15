abort <- function(message) {
  stop(message, call. = FALSE)
}

parse_bool <- function(value, name) {
  normalized <- tolower(value)
  if (normalized %in% c("true", "1", "yes")) return(TRUE)
  if (normalized %in% c("false", "0", "no")) return(FALSE)
  abort(sprintf("%s must be true or false", name))
}

parse_number <- function(value, name, minimum = -Inf, maximum = Inf, integer = FALSE) {
  parsed <- suppressWarnings(as.numeric(value))
  if (length(parsed) != 1L || is.na(parsed) || parsed < minimum || parsed > maximum) {
    abort(sprintf("%s must be between %s and %s", name, minimum, maximum))
  }
  if (integer && parsed != as.integer(parsed)) abort(sprintf("%s must be an integer", name))
  if (integer) as.integer(parsed) else parsed
}

parse_non_tree_ids <- function(value) {
  values <- strsplit(value, ",", fixed = TRUE)[[1L]]
  parsed <- suppressWarnings(as.numeric(trimws(values)))
  if (anyNA(parsed)) abort("--non-tree-ids must be a comma-separated list of numbers")
  unique(parsed)
}

validate_dimension <- function(value, label) {
  if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", value)) {
    abort(sprintf("%s is not a valid LAS dimension name: %s", label, value))
  }
  value
}

parse_segmentation_spec <- function(value) {
  parts <- strsplit(value, ",", fixed = TRUE)[[1L]]
  parts <- trimws(parts)
  if (length(parts) > 3L || !length(parts) || !nzchar(parts[[1L]])) {
    abort("--segmentation-spec must be INSTANCE[,SPECIES,SPECIES_PROB]")
  }
  parts <- c(parts, rep("", 3L - length(parts)))
  has_species <- nzchar(parts[[2L]])
  has_probability <- nzchar(parts[[3L]])
  if (xor(has_species, has_probability)) {
    abort(sprintf("Species and species-probability dimensions must be supplied together for %s", parts[[1L]]))
  }
  instance <- validate_dimension(parts[[1L]], "Instance dimension")
  source <- sub("^PredInstance_", "", instance)
  if (identical(source, instance)) source <- instance
  list(
    instance = instance,
    source = source,
    species = if (has_species) validate_dimension(parts[[2L]], "Species dimension") else NULL,
    species_prob = if (has_probability) validate_dimension(parts[[3L]], "Species probability dimension") else NULL
  )
}

usage <- function() {
  paste(
    "Usage: run.R --input FILE_OR_DIRECTORY --output-dir DIR [options]",
    "",
    "Common options:",
    "  --segmentation-spec INSTANCE[,SPECIES,SPECIES_PROB]  Repeat for each existing segmentation",
    "  --non-tree-ids IDS                  Comma-separated IDs to exclude (default: 0)",
    "  --enable-csp true|false             Create PredInstance_CSP (default: false)",
    "  --seed-mode automatic|supplied     CSP seed source (default: automatic)",
    "  --seed-file TSV                    X/Y/Z/TreeID table for supplied seeds",
    "  --aoi-json GEOJSON                 Optional AOI in the point-cloud CRS",
    "  --dtm-resolution METRES            DTM cell size (default: 0.2)",
    "  --read-chunk-size METRES            Inventory-only spatial chunk size (default: 25)",
    "  --chunk-buffer METRES               DTM chunk buffer (default: 5)",
    "  --inventory-partitions COUNT        Disk-backed instance partitions (default: 64)",
    "  --random-seed INTEGER              RANSAC seed (default: 42)",
    "",
    "Fine-tuning options:",
    "  --routing-workers INTEGER          CSP routing workers (default: 1)",
    "  --geometry-threads INTEGER         Geometry threads (default: 1)",
    "  --voxel-resolution METRES          CSP routing voxel (default: 0.3)",
    "  --geometry-k INTEGER               Geometry neighbours (default: 10)",
    "  --verticality-weight NUMBER        CSP V_w (default: 0)",
    "  --linearity-weight NUMBER          CSP L_w (default: 0)",
    "  --sphericity-weight NUMBER         CSP S_w (default: 0)",
    "  --seed-resolution METRES           Raster seed cell size (default: 0.1)",
    "  --seed-z-min METRES                Seed slice minimum (default: 0.5)",
    "  --seed-z-max METRES                Seed slice maximum (default: 2)",
    "  --seed-quantile NUMBER             Seed density quantile (default: 0.975)",
    "  --seed-eps METRES                  Seed clustering radius (default: 0.2)",
    "  --slice-min METRES                 Inventory slice minimum (default: 0.3)",
    "  --slice-max METRES                 Inventory slice maximum (default: 4)",
    "  --slice-increment METRES           Inventory slice step (default: 0.2)",
    "  --slice-width METRES               Inventory slice width (default: 0.1)",
    "  --max-dbh METRES                   Maximum accepted DBH (default: 1)",
    sep = "\n"
  )
}

parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  config <- list(
    input = NULL,
    output_dir = NULL,
    segmentation_specs = list(),
    non_tree_ids = 0,
    enable_csp = FALSE,
    seed_mode = "automatic",
    seed_file = NULL,
    aoi_json = NULL,
    dtm_resolution = 0.2,
    read_chunk_size = 25,
    chunk_buffer = 5,
    inventory_partitions = 64L,
    random_seed = 42L,
    routing_workers = 1L,
    geometry_threads = 1L,
    voxel_resolution = 0.3,
    geometry_k = 10L,
    verticality_weight = 0,
    linearity_weight = 0,
    sphericity_weight = 0,
    seed_resolution = 0.1,
    seed_z_min = 0.5,
    seed_z_max = 2,
    seed_quantile = 0.975,
    seed_eps = 0.2,
    slice_min = 0.3,
    slice_max = 4,
    slice_increment = 0.2,
    slice_width = 0.1,
    max_dbh = 1
  )

  if (any(args %in% c("--help", "-h"))) {
    cat(usage(), "\n")
    quit(status = 0L)
  }

  repeated <- character()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) abort(sprintf("Unexpected argument: %s", key))
    if (i == length(args)) abort(sprintf("Missing value for %s", key))
    value <- args[[i + 1L]]
    i <- i + 2L

    if (key == "--input") config$input <- value
    else if (key == "--output-dir") config$output_dir <- value
    else if (key == "--segmentation-spec") repeated <- c(repeated, value)
    else if (key == "--non-tree-ids") config$non_tree_ids <- parse_non_tree_ids(value)
    else if (key == "--enable-csp") config$enable_csp <- parse_bool(value, key)
    else if (key == "--seed-mode") config$seed_mode <- value
    else if (key == "--seed-file") config$seed_file <- value
    else if (key == "--aoi-json") config$aoi_json <- value
    else if (key == "--dtm-resolution") config$dtm_resolution <- parse_number(value, key, 0.001)
    else if (key == "--read-chunk-size") config$read_chunk_size <- parse_number(value, key, 1)
    else if (key == "--chunk-buffer") config$chunk_buffer <- parse_number(value, key, 0)
    else if (key == "--inventory-partitions") config$inventory_partitions <- parse_number(value, key, 1, 4096, TRUE)
    else if (key == "--random-seed") config$random_seed <- parse_number(value, key, 0, .Machine$integer.max, TRUE)
    else if (key == "--routing-workers") config$routing_workers <- parse_number(value, key, 1, 256, TRUE)
    else if (key == "--geometry-threads") config$geometry_threads <- parse_number(value, key, 1, 256, TRUE)
    else if (key == "--voxel-resolution") config$voxel_resolution <- parse_number(value, key, 0.001)
    else if (key == "--geometry-k") config$geometry_k <- parse_number(value, key, 1, 1000, TRUE)
    else if (key == "--verticality-weight") config$verticality_weight <- parse_number(value, key, 0, 1)
    else if (key == "--linearity-weight") config$linearity_weight <- parse_number(value, key, 0, 1)
    else if (key == "--sphericity-weight") config$sphericity_weight <- parse_number(value, key, 0, 1)
    else if (key == "--seed-resolution") config$seed_resolution <- parse_number(value, key, 0.001)
    else if (key == "--seed-z-min") config$seed_z_min <- parse_number(value, key)
    else if (key == "--seed-z-max") config$seed_z_max <- parse_number(value, key)
    else if (key == "--seed-quantile") config$seed_quantile <- parse_number(value, key, 0, 1)
    else if (key == "--seed-eps") config$seed_eps <- parse_number(value, key, 0.001)
    else if (key == "--slice-min") config$slice_min <- parse_number(value, key, 0)
    else if (key == "--slice-max") config$slice_max <- parse_number(value, key, 0)
    else if (key == "--slice-increment") config$slice_increment <- parse_number(value, key, 0.001)
    else if (key == "--slice-width") config$slice_width <- parse_number(value, key, 0.001)
    else if (key == "--max-dbh") config$max_dbh <- parse_number(value, key, 0.001)
    else abort(sprintf("Unknown option: %s", key))
  }

  if (is.null(config$input)) abort("--input is required")
  if (is.null(config$output_dir)) abort("--output-dir is required")
  if (!config$seed_mode %in% c("automatic", "supplied")) abort("--seed-mode must be automatic or supplied")
  if (config$enable_csp && config$seed_mode == "supplied" && is.null(config$seed_file)) {
    abort("--seed-file is required when --seed-mode supplied")
  }
  if (config$seed_z_min >= config$seed_z_max) abort("--seed-z-min must be lower than --seed-z-max")
  if (config$slice_min >= config$slice_max) abort("--slice-min must be lower than --slice-max")
  if (config$slice_width > config$slice_increment / 2) {
    abort("--slice-width may not exceed half of --slice-increment")
  }
  config$segmentation_specs <- lapply(repeated, parse_segmentation_spec)
  dimensions <- vapply(config$segmentation_specs, `[[`, character(1), "instance")
  if (anyDuplicated(dimensions)) abort("Each instance dimension may only be supplied once")
  if (config$enable_csp && "PredInstance_CSP" %in% dimensions) {
    abort("Do not supply PredInstance_CSP when CSP is enabled; the tool creates it")
  }
  if (!config$enable_csp && !length(config$segmentation_specs)) {
    abort("Provide at least one --segmentation-spec or enable CSP")
  }
  config
}
