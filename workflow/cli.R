source(file.path(dirname(normalizePath(sys.frame(1)$ofile)), "inventory_policy.R"), local = TRUE)

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

# One segmentation: its instance attribute and, optionally, the species and
# species-probability attributes that go with it (both or neither).
segmentation_spec <- function(instance, species = NULL, species_prob = NULL) {
  if (xor(is.null(species), is.null(species_prob))) {
    abort(sprintf("--species-attribute and --species-probability-attribute must be supplied together for %s", instance))
  }
  instance <- validate_dimension(instance, "Instance attribute")
  list(
    instance = instance,
    source = sub("^PredInstance_", "", instance),
    species = if (!is.null(species)) validate_dimension(species, "Species attribute"),
    species_prob = if (!is.null(species_prob)) validate_dimension(species_prob, "Species probability attribute")
  )
}

validate_owned_config <- function(config) {
  if (config$inventory_retained_points < inventory_large_instance_sample_limit + config$dbh_band_points) {
    abort("--inventory-retained-points must hold one large tree: sample limit + --dbh-band-points")
  }
  if (config$enable_csp) abort("Owned COPC inventory cannot be combined with --enable-csp")
  if (config$inventory_workers + config$inventory_decode_threads > 20L) {
    abort("Owned inventory workers + shared decode threads must be at most 20 CPUs")
  }
  invisible(config)
}

# COPC requires its info VLR first, right after the 375-byte LAS 1.4 header.
input_is_copc <- function(path) {
  bytes <- readBin(path, raw(), 393L)
  length(bytes) == 393L && identical(rawToChar(bytes[1:4]), "LASF") &&
    identical(rawToChar(bytes[378:381]), "copc")
}

input_point_count <- function(files) {
  sum(vapply(files, function(path) header_point_count(lidR::readLASheader(path)), numeric(1)))
}

# Partitions slow down superlinearly at about 7 M points each (925.7 M-point
# benchmark, workflow/resource_request.py), so aim for about 3 M.
default_inventory_partitions <- function(points) {
  wanted <- max(64, ceiling(points / 3e6))
  as.integer(min(4096, 2^ceiling(log2(wanted))))
}

# auto picks a reader whose memory does not grow with the cloud: the owned
# COPC windows for one COPC file, one sequential scan for one LAS/LAZ file.
# Several files and CSP segmentation keep the spatial reader.
resolve_inventory_read_mode <- function(config, files) {
  if (identical(config$inventory_read_mode, "auto")) {
    config$inventory_read_mode <- if (config$enable_csp || length(files) != 1L) "spatial" else
      if (input_is_copc(files[[1L]])) "copc-owned" else "sequential"
    if (config$inventory_read_mode == "copc-owned") validate_owned_config(config)
  }
  if (is.na(config$inventory_partitions)) {
    config$inventory_partitions <- default_inventory_partitions(input_point_count(files))
  }
  if (is.na(config$inventory_index_workers)) {
    config$inventory_index_workers <- if (config$inventory_read_mode == "copc-owned") 8L else 0L
  }
  if (config$inventory_index_workers > 0L && config$inventory_read_mode != "copc-owned") {
    abort("--inventory-index-workers requires a COPC input read with copc-owned")
  }
  config
}

usage <- function() {
  paste(
    "Usage: run.R --input MERGED.(copc.)laz --output-dir DIR [options]",
    "",
    "Common options:",
    "  --instance-attribute NAME          Instance attribute of an existing segmentation; repeat",
    "                                     for each segmentation (e.g. PredInstance_FM)",
    "  --species-attribute NAME           Species attribute of the preceding instance attribute",
    "                                     (e.g. species_id_FM; optional)",
    "  --species-probability-attribute NAME  Species probability attribute of the preceding",
    "                                     instance attribute (e.g. species_prob_FM; with species)",
    "  --inventory-attributes LIST        Additional numeric LAS attributes (comma-separated; default: none)",
    "  --non-tree-ids IDS                  Comma-separated IDs to exclude (default: 0)",
    "  --enable-csp true|false             Create PredInstance_CSP (default: false)",
    "  --seed-mode automatic|supplied     CSP seed source (default: automatic)",
    "  --seed-file TSV                    X/Y/Z/TreeID table for supplied seeds",
    "  --aoi-json GEOJSON                 Optional AOI in the point-cloud CRS",
    "  --dtm GEOTIFF                      Use this DTM; otherwise generate one",
    "  --species-lookup TSV               Species names and taxonomy/model versions",
    "  --dtm-resolution METRES            DTM cell size (default: 0.2)",
    "  --read-chunk-size METRES            Spatial DTM/inventory tile size (default: 300)",
    "  --chunk-buffer METRES               DTM chunk buffer (default: 5)",
    "  --dtm-workers INTEGER                Parallel DTM workers (default: 10)",
    "  --dtm-strategy auto|spatial|streaming DTM read strategy (default: auto)",
    "  --dtm-candidate-resolution METRES   Streaming low-surface grid (default: 0.1)",
    "  --dtm-streaming-threshold POINTS    Auto threshold (default: 50000000)",
    "  --inventory-partitions COUNT        Disk-backed instance partitions (default: about 3 M",
    "                                     points each, a power of two from 64 to 4096)",
    "  --inventory-read-mode auto|spatial|sequential|copc-owned  Inventory reader (default: auto:",
    "                                     copc-owned for one COPC file, sequential for one LAS/LAZ",
    "                                     file, spatial for several files or with --enable-csp)",
    "  --inventory-index-workers COUNT    COPC index workers (0: legacy; 1..8; default: 8 with",
    "                                     copc-owned, else 0)",
    "  --inventory-index-queue-mb MIB      Compressed index queue (32..1024; default: 256)",
    "  --inventory-core-size METRES        Owned core width (default: 50)",
    "  --inventory-buffer METRES           Owned window initial buffer (default: 20)",
    "  --inventory-batch-points COUNT      Max COPC node batch (default: 1000000)",
    "  --inventory-prefetch-mb MIB         Shared in-memory queue (default: 256; max: 16384)",
    "  --inventory-retained-points COUNT   Window retention budget (default: 40000000)",
    "  --inventory-decode-threads COUNT    COPC decode threads per worker (default: 4)",
    "  --inventory-workers INTEGER         Parallel instance partitions (default: 1)",
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
    "  --dbh-band-points N                Stem-band points kept per tree for DBH; above it",
    "                                     points are merged from 1 mm voxels up (default: 1000000)",
    "  --dbh-runs N                       Upstream DBH fits per tree; the median is reported (default: 5)",
    "  --dbh-evidence true|false          Fill the inventory's dbh_fits, dbh_slices and crown_hull",
    "                                     columns (default: true)",
    "  --scratch-dir DIR                  Intermediate files (default: $TMPDIR, else next to output)",
    sep = "\n"
  )
}

parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  config <- list(
    input = NULL,
    output_dir = NULL,
    segmentation_specs = list(),
    inventory_attributes = character(),
    non_tree_ids = 0,
    enable_csp = FALSE,
    seed_mode = "automatic",
    seed_file = NULL,
    aoi_json = NULL,
    dtm = NULL,
    species_lookup = NULL,
    dtm_resolution = 0.2,
    read_chunk_size = 300,
    chunk_buffer = 5,
    dtm_workers = 10L,
    dtm_strategy = "auto",
    dtm_candidate_resolution = 0.1,
    dtm_streaming_threshold = 50000000L,
    inventory_partitions = NA_integer_,
    inventory_read_mode = "auto",
    inventory_workers = 1L,
    inventory_index_workers = NA_integer_,
    inventory_index_queue_mb = 256L,
    inventory_core_size = 50,
    inventory_buffer = 20,
    inventory_batch_points = 2000000L,
    inventory_retained_points = 40000000L,
    inventory_decode_threads = 4L,
    inventory_prefetch_mb = 256L,
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
    max_dbh = 1,
    dbh_band_points = inventory_dbh_band_point_budget,
    dbh_runs = 5L,
    dbh_evidence = TRUE,
    scratch_dir = NULL
  )

  if (any(args %in% c("--help", "-h"))) {
    cat(usage(), "\n")
    quit(status = 0L)
  }

  segmentations <- list()
  # Species options belong to the most recent --instance-attribute.
  set_species <- function(field, key, value) {
    last <- length(segmentations)
    if (!last) abort(sprintf("%s must follow an --instance-attribute", key))
    if (!is.null(segmentations[[last]][[field]])) abort(sprintf("%s given twice for %s", key, segmentations[[last]]$instance))
    segmentations[[last]][[field]] <<- value
  }
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) abort(sprintf("Unexpected argument: %s", key))
    if (i == length(args)) abort(sprintf("Missing value for %s", key))
    value <- args[[i + 1L]]
    i <- i + 2L

    if (key == "--input") config$input <- value
    else if (key == "--output-dir") config$output_dir <- value
    else if (key == "--instance-attribute") segmentations[[length(segmentations) + 1L]] <- list(instance = value)
    else if (key == "--species-attribute") set_species("species", key, value)
    else if (key == "--species-probability-attribute") set_species("species_prob", key, value)
    else if (key == "--inventory-attributes") {
      values <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
      if (!length(values) || any(!nzchar(values)) || endsWith(value, ",")) abort("--inventory-attributes requires dimension names")
      values <- vapply(values, validate_dimension, character(1), label = "Inventory attribute", USE.NAMES = FALSE)
      if (any(values %in% c("Zref", "source_order__"))) abort("Reserved inventory attribute")
      config$inventory_attributes <- unique(c(config$inventory_attributes, values))
    }
    else if (key == "--non-tree-ids") config$non_tree_ids <- parse_non_tree_ids(value)
    else if (key == "--enable-csp") config$enable_csp <- parse_bool(value, key)
    else if (key == "--seed-mode") config$seed_mode <- value
    else if (key == "--seed-file") config$seed_file <- value
    else if (key == "--aoi-json") config$aoi_json <- value
    else if (key == "--dtm") config$dtm <- value
    else if (key == "--species-lookup") config$species_lookup <- value
    else if (key == "--dtm-resolution") config$dtm_resolution <- parse_number(value, key, 0.001)
    else if (key == "--read-chunk-size") config$read_chunk_size <- parse_number(value, key, 1)
    else if (key == "--chunk-buffer") config$chunk_buffer <- parse_number(value, key, 0)
    else if (key == "--dtm-workers") config$dtm_workers <- parse_number(value, key, 1, 256, TRUE)
    else if (key == "--dtm-strategy") config$dtm_strategy <- value
    else if (key == "--dtm-candidate-resolution") config$dtm_candidate_resolution <- parse_number(value, key, 0.001)
    else if (key == "--dtm-streaming-threshold") config$dtm_streaming_threshold <- parse_number(value, key, 1, Inf, TRUE)
    else if (key == "--inventory-partitions") config$inventory_partitions <- parse_number(value, key, 1, 4096, TRUE)
    else if (key == "--inventory-read-mode") config$inventory_read_mode <- value
    else if (key == "--inventory-workers") config$inventory_workers <- parse_number(value, key, 1, 20, TRUE)
    else if (key == "--inventory-index-workers") config$inventory_index_workers <- parse_number(value, key, 0, 8, TRUE)
    else if (key == "--inventory-index-queue-mb") config$inventory_index_queue_mb <- parse_number(value, key, 32, 1024, TRUE)
    else if (key == "--inventory-core-size") config$inventory_core_size <- parse_number(value, key, 0.001)
    else if (key == "--inventory-buffer") config$inventory_buffer <- parse_number(value, key, 0)
    else if (key == "--inventory-batch-points") config$inventory_batch_points <- parse_number(value, key, 1, 10000000, TRUE)
    else if (key == "--inventory-retained-points") config$inventory_retained_points <- parse_number(value, key, inventory_instance_point_limit, 100000000, TRUE)
    else if (key == "--inventory-prefetch-mb") config$inventory_prefetch_mb <- parse_number(value, key, 1, 16384, TRUE)
    else if (key == "--inventory-decode-threads") config$inventory_decode_threads <- parse_number(value, key, 1, 19, TRUE)
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
    else if (key == "--dbh-band-points") config$dbh_band_points <- parse_number(value, key, 1000, 100000000, TRUE)
    else if (key == "--dbh-runs") config$dbh_runs <- parse_number(value, key, 1, 99, TRUE)
    else if (key == "--dbh-evidence") config$dbh_evidence <- parse_bool(value, key)
    else if (key == "--scratch-dir") config$scratch_dir <- value
    else abort(sprintf("Unknown option: %s", key))
  }

  if (is.null(config$input)) abort("--input is required")
  if (is.null(config$output_dir)) abort("--output-dir is required")
  if (!config$dtm_strategy %in% c("auto", "spatial", "streaming")) {
    abort("--dtm-strategy must be auto, spatial, or streaming")
  }
  if (!config$inventory_read_mode %in% c("auto", "spatial", "sequential", "copc-owned")) {
    abort("--inventory-read-mode must be auto, spatial, sequential, or copc-owned")
  }
  if (isTRUE(config$inventory_index_workers > 0L) &&
      !config$inventory_read_mode %in% c("auto", "copc-owned")) {
    abort("--inventory-index-workers requires --inventory-read-mode copc-owned")
  }
  if (config$inventory_read_mode == "copc-owned") validate_owned_config(config)
  if (config$dtm_candidate_resolution > config$dtm_resolution) {
    abort("--dtm-candidate-resolution may not exceed --dtm-resolution")
  }
  if (!config$seed_mode %in% c("automatic", "supplied")) abort("--seed-mode must be automatic or supplied")
  if (config$enable_csp && config$seed_mode == "supplied" && is.null(config$seed_file)) {
    abort("--seed-file is required when --seed-mode supplied")
  }
  if (config$seed_z_min >= config$seed_z_max) abort("--seed-z-min must be lower than --seed-z-max")
  if (config$slice_min >= config$slice_max) abort("--slice-min must be lower than --slice-max")
  if (config$slice_width > config$slice_increment / 2) {
    abort("--slice-width may not exceed half of --slice-increment")
  }
  config$segmentation_specs <- lapply(segmentations, function(value) {
    spec <- segmentation_spec(value$instance, value$species, value$species_prob)
    spec$attributes <- config$inventory_attributes
    spec
  })
  dimensions <- vapply(config$segmentation_specs, `[[`, character(1), "instance")
  if (anyDuplicated(dimensions)) abort("Each instance dimension may only be supplied once")
  sources <- vapply(config$segmentation_specs, `[[`, character(1), "source")
  if (any(!nzchar(sources)) || anyDuplicated(sources)) {
    abort("Each segmentation must have a distinct, non-empty output source name")
  }
  if (config$enable_csp && "PredInstance_CSP" %in% dimensions) {
    abort("Do not supply PredInstance_CSP when CSP is enabled; the tool creates it")
  }
  if (!config$enable_csp && !length(config$segmentation_specs)) {
    abort("Provide at least one --instance-attribute or enable CSP")
  }
  config
}
