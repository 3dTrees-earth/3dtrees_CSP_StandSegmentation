#!/usr/bin/env Rscript

all_arguments <- commandArgs(trailingOnly = FALSE)
test_file <- sub("^--file=", "", grep("^--file=", all_arguments, value = TRUE)[[1L]])
root <- normalizePath(file.path(dirname(test_file), "..", ".."), mustWork = TRUE)
source(file.path(root, "workflow", "cli.R"))
source(file.path(root, "workflow", "tool.R"))

expect_error <- function(expression, pattern) {
  message <- tryCatch({ force(expression); NULL }, error = conditionMessage)
  stopifnot(!is.null(message), grepl(pattern, message))
}

spec <- segmentation_spec("PredInstance_FM", "species_id_FM", "species_prob_FM")
stopifnot(
  identical(spec$source, "FM"),
  identical(spec$species, "species_id_FM"),
  identical(spec$species_prob, "species_prob_FM")
)
expect_error(segmentation_spec("PredInstance_FM", "species_id_FM"), "supplied together")
expect_error(segmentation_spec("bad-name"), "valid LAS dimension")
# Species options belong to the --instance-attribute before them.
paired <- parse_cli_args(c("--input", "input.laz", "--output-dir", "output",
                           "--instance-attribute", "PredInstance_RCT",
                           "--species-attribute", "species_id_FM",
                           "--species-probability-attribute", "species_prob_FM",
                           "--instance-attribute", "PredInstance_SAT"))$segmentation_specs
stopifnot(identical(paired[[1L]]$species, "species_id_FM"), identical(paired[[1L]]$species_prob, "species_prob_FM"),
          is.null(paired[[2L]]$species), is.null(paired[[2L]]$species_prob))
expect_error(parse_cli_args(c("--input", "input.laz", "--output-dir", "output",
                              "--species-attribute", "species_id_FM")), "must follow an --instance-attribute")
expect_error(parse_cli_args(c("--input", "input.laz", "--output-dir", "output",
                              "--instance-attribute", "PredInstance_FM",
                              "--species-attribute", "species_id_FM")), "supplied together")

config <- parse_cli_args(c(
  "--input", "input.laz",
  "--output-dir", "output",
  "--instance-attribute", "PredInstance_SAT",
  "--instance-attribute", "PredInstance_FM", "--species-attribute", "species_id_FM", "--species-probability-attribute", "species_prob_FM",
  "--inventory-attributes", "PredScore_FM,PredSemantic_FM",
  "--dtm-resolution", "0.2",
  "--routing-workers", "1"
))
stopifnot(
  length(config$segmentation_specs) == 2L,
  identical(config$dtm_resolution, 0.2),
  identical(config$routing_workers, 1L),
  identical(config$read_chunk_size, 300),
  identical(config$chunk_buffer, 5),
  is.na(config$inventory_partitions),
  identical(config$inventory_read_mode, "auto"),
  is.na(config$inventory_index_workers),
  identical(config$inventory_workers, 1L),
  identical(config$dtm_workers, 10L),
  identical(config$dtm_strategy, "auto"),
  identical(config$dtm_candidate_resolution, 0.1),
  identical(config$dtm_streaming_threshold, 50000000L),
  isTRUE(config$dbh_evidence),
  !config$enable_csp
)
stopifnot(!parse_cli_args(c("--input", "input.laz", "--output-dir", "output",
                             "--instance-attribute", "PredInstance_SAT",
                             "--dbh-evidence", "false"))$dbh_evidence)

dtm_config <- parse_cli_args(c(
  "--input", "merged.copc.laz", "--output-dir", "output",
  "--instance-attribute", "PredInstance_SAT", "--dtm", "terrain.tif"
))
stopifnot(identical(dtm_config$dtm, "terrain.tif"))
# One inventory table in CSV and Parquet: there is no format choice any more.
expect_error(parse_cli_args(c("--input", "input.laz", "--output-dir", "output",
                              "--instance-attribute", "PredInstance_SAT", "--output-formats", "csv")),
             "Unknown option")

worker_pids <- unlist(parallel_chunk_map(
  as.list(seq_len(4L)),
  workers = 2L,
  fun = function(value) {
    Sys.sleep(0.05)
    Sys.getpid()
  }
))
stopifnot(length(unique(worker_pids)) == 2L)

edge_chunks <- spatial_chunks(list(xmin = 0, xmax = 25, ymin = 0, ymax = 25), 25)
stopifnot(
  length(edge_chunks) == 4L,
  any(vapply(edge_chunks, function(chunk) chunk[["xmin"]] == 25, logical(1))),
  any(vapply(edge_chunks, function(chunk) chunk[["ymin"]] == 25, logical(1)))
)

dimensions <- inventory_read_dimensions(config$segmentation_specs)
stopifnot(
  setequal(
    dimensions,
    c(
      "PredInstance_SAT", "PredInstance_FM", "species_id_FM",
      "species_prob_FM", "PredScore_FM", "PredSemantic_FM"
    )
  ),
  identical(
    build_read_selector(
      dimensions,
      c(
        "unused", "PredInstance_FM", "species_id_FM", "species_prob_FM",
        "PredScore_FM", "PredSemantic_FM", "PredInstance_SAT"
      )
    ),
    "xyz234567"
  ),
  identical(build_read_selector("late_dimension", c(rep("unused", 9L), "late_dimension")), "xyz0")
)

expect_error(
  parse_cli_args(c("--input", "input.laz", "--output-dir", "output")),
  "instance-attribute"
)
expect_error(
  parse_cli_args(c("--input", "input.laz", "--output-dir", "output",
                   "--instance-attribute", "PredInstance_SAT",
                   "--instance-attribute", "SAT")),
  "distinct"
)
expect_error(
  parse_cli_args(c(
    "--input", "input.laz", "--output-dir", "output",
    "--instance-attribute", "PredInstance_SAT", "--dtm-strategy", "invalid"
  )),
  "dtm-strategy"
)
stopifnot(identical(parse_cli_args(c(
  "--input", "input.laz", "--output-dir", "output",
  "--instance-attribute", "PredInstance_SAT", "--inventory-workers", "4"
))$inventory_workers, 4L))
expect_error(parse_cli_args(c(
  "--input", "input.laz", "--output-dir", "output",
  "--instance-attribute", "PredInstance_SAT", "--inventory-workers", "21"
)), "inventory-workers")
expect_error(parse_cli_args(c(
  "--input", "input.laz", "--output-dir", "output",
  "--instance-attribute", "PredInstance_SAT", "--inventory-read-mode", "invalid"
)), "inventory-read-mode")
expect_error(
  parse_cli_args(c(
    "--input", "input.laz", "--output-dir", "output",
    "--instance-attribute", "PredInstance_SAT", "--dtm-resolution", "0.2",
    "--dtm-candidate-resolution", "0.3"
  )),
  "candidate-resolution"
)
expect_error(
  parse_cli_args(c(
    "--input", "input.laz", "--output-dir", "output",
    "--enable-csp", "true", "--seed-mode", "supplied"
  )),
  "seed-file"
)

points <- data.table::data.table(
  PredInstance_FM = c(1, 1, 1, 2),
  species_id_FM = c(4, 4, 7, -1),
  species_prob_FM = c(0.8, 0.6, 0.9, -1),
  PredScore_FM = c(0.9, 0.9, 0.2, -1),
  PredSemantic_FM = c(1, 2, 2, 0)
)
species <- species_by_instance(points, spec)
tree_one <- species[PredInstance_FM == 1]
stopifnot(
  tree_one$species_id == 4,
  tree_one$species_prob == 0.7,
  abs(tree_one$species_conflict_fraction - 1 / 3) < 1e-12
)
fm <- fm_by_instance(points, "PredInstance_FM")
fm_one <- fm[PredInstance_FM == 1]
stopifnot(
  identical(fm_one$pred_score_fm_values, '[{"score":0.9,"point_count":2},{"score":0.2,"point_count":1}]'),
  is.na(fm[PredInstance_FM == 2]$pred_score_fm_values),
  fm_one$wood_point_count == 1,
  fm_one$leaf_point_count == 2,
  abs(fm_one$wood_share - 1 / 3) < 1e-12,
  abs(fm_one$leaf_share - 2 / 3) < 1e-12
)

header <- lidR::readLASheader(file.path(root, "tests", "tool", "fixtures", "multi-instance.copc.laz"))
valid_wkt <- header@VLR[["WKT OGC CS"]][["WKT OGC COORDINATE SYSTEM"]]
header@EVLR[["WKT OGC CS"]] <- header@VLR[["WKT OGC CS"]]
header@VLR[["WKT OGC CS"]][["WKT OGC COORDINATE SYSTEM"]] <- "COMPD_CS[\"broken\""
stopifnot(point_cloud_crs(header)$epsg == 32632L)
header@EVLR[["WKT OGC CS"]][["WKT OGC COORDINATE SYSTEM"]] <- sf::st_crs(32633)$wkt
header@VLR[["WKT OGC CS"]][["WKT OGC COORDINATE SYSTEM"]] <- valid_wkt
expect_error(point_cloud_crs(header), "conflicting valid CRS")

cat("tool helper tests passed\n")

owned_config <- parse_cli_args(c("--input", "fixture.copc.laz", "--output-dir", "out", "--instance-attribute", "PredInstance_FM",
                                  "--inventory-read-mode", "copc-owned", "--inventory-workers", "16",
                                  "--inventory-decode-threads", "4"))
stopifnot(owned_config$inventory_workers == 16L)
expect_error(parse_cli_args(c("--input", "fixture.copc.laz", "--output-dir", "out", "--instance-attribute", "PredInstance_FM",
                             "--inventory-read-mode", "copc-owned", "--inventory-workers", "17",
                             "--inventory-decode-threads", "4")), "at most 20")

expect_error(parse_cli_args(c("--input", "input.laz", "--output-dir", "output",
                               "--enable-csp", "true", "--inventory-read-mode", "copc-owned")),
             "cannot be combined")

minimal <- parse_cli_args(c("--input", "input.laz", "--output-dir", "out",
                           "--instance-attribute", "PredInstance_FM"))
stopifnot(identical(inventory_read_dimensions(minimal$segmentation_specs), "PredInstance_FM"),
          identical(inventory_layouts(minimal$segmentation_specs)$FM$columns,
                    c("X", "Y", "Z", "PredInstance_FM")))
expect_error(parse_cli_args(c("--input", "input.laz", "--output-dir", "out",
                             "--instance-attribute", "PredInstance_FM",
                             "--inventory-attributes", "Zref")), "Reserved")
expect_error(parse_cli_args(c("--input", "input.laz", "--output-dir", "out",
                             "--instance-attribute", "PredInstance_FM",
                             "--inventory-attributes", "Intensity,")), "requires dimension")

index_args <- c("--input", "input.copc.laz", "--output-dir", "out",
                "--instance-attribute", "PredInstance_FM")
index_config <- parse_cli_args(c(index_args, "--inventory-read-mode", "copc-owned",
                                 "--inventory-index-workers", "8", "--inventory-index-queue-mb", "256"))
stopifnot(index_config$inventory_index_workers == 8L, index_config$inventory_index_queue_mb == 256L,
          is.na(parse_cli_args(index_args)$inventory_index_workers))
expect_error(parse_cli_args(c(index_args, "--inventory-read-mode", "sequential",
                              "--inventory-index-workers", "8")), "requires")
expect_error(parse_cli_args(c(index_args, "--inventory-index-workers", "9")), "inventory-index-workers")
expect_error(parse_cli_args(c(index_args, "--inventory-index-queue-mb", "31")), "inventory-index-queue-mb")

# auto resolves to a reader whose memory does not grow with the cloud.
copc <- file.path(root, "tests", "tool", "fixtures", "multi-instance.copc.laz")
plain_dir <- tempfile("csp-plain-laz-")
dir.create(plain_dir)
plain <- file.path(plain_dir, "plain.laz")
lidR::writeLAS(suppressWarnings(lidR::readLAS(copc, select = "xyz")), plain)
stopifnot(input_is_copc(copc), !input_is_copc(plain))
auto <- parse_cli_args(c("--input", copc, "--output-dir", "out", "--instance-attribute", "PredInstance_FM"))
owned <- resolve_inventory_read_mode(auto, copc)
stopifnot(identical(owned$inventory_read_mode, "copc-owned"), owned$inventory_index_workers == 8L,
          owned$inventory_partitions == 64L)
stopifnot(default_inventory_partitions(29011) == 64L, default_inventory_partitions(192e6) == 64L,
          default_inventory_partitions(668e6) == 256L, default_inventory_partitions(4.26e9) == 2048L,
          default_inventory_partitions(1e11) == 4096L)
explicit_partitions <- parse_cli_args(c("--input", copc, "--output-dir", "out", "--instance-attribute", "PredInstance_FM",
                                        "--inventory-partitions", "16"))
stopifnot(resolve_inventory_read_mode(explicit_partitions, copc)$inventory_partitions == 16L)
sequential <- resolve_inventory_read_mode(auto, plain)
stopifnot(identical(sequential$inventory_read_mode, "sequential"), sequential$inventory_index_workers == 0L)
stopifnot(identical(resolve_inventory_read_mode(auto, c(plain, plain))$inventory_read_mode, "spatial"))
csp_auto <- parse_cli_args(c("--input", copc, "--output-dir", "out", "--enable-csp", "true"))
stopifnot(identical(resolve_inventory_read_mode(csp_auto, copc)$inventory_read_mode, "spatial"))
explicit <- parse_cli_args(c("--input", copc, "--output-dir", "out", "--instance-attribute", "PredInstance_FM",
                             "--inventory-read-mode", "spatial"))
stopifnot(identical(resolve_inventory_read_mode(explicit, copc)$inventory_read_mode, "spatial"))
indexed <- parse_cli_args(c("--input", plain, "--output-dir", "out", "--instance-attribute", "PredInstance_FM",
                            "--inventory-index-workers", "4"))
expect_error(resolve_inventory_read_mode(indexed, plain), "requires a COPC input")
# A DTM tile whose points all lie in its buffer is skipped, not cropped.
grid <- expand.grid(X = seq(0, 9.95, by = 0.05), Y = seq(0, 9.95, by = 0.05))
ground <- lidR::LAS(data.table::data.table(grid, Z = 0, Classification = 2L))
sf::st_crs(ground) <- 32632
ground_file <- file.path(plain_dir, "ground.laz")
lidR::writeLAS(ground, ground_file)
dtm_work <- file.path(plain_dir, "dtm-work")
dtm_stage <- file.path(plain_dir, "dtm-stage")
dir.create(dtm_work)
dir.create(dtm_stage)
tiles <- list(c(xmin = 0, xmax = 10, ymin = 0, ymax = 10), c(xmin = 10, xmax = 20, ymin = 0, ymax = 10))
buffered <- chunked_dtm(ground_file, tiles, 2, 0.5, dtm_work, dtm_stage)
stopifnot(buffered$chunk_count == 2L, buffered$valid_chunk_count == 1L)
unlink(plain_dir, recursive = TRUE)
cat("read mode resolution tests passed\n")

# lidR's spatial index snapped 4021833.4 to 4021833.4000000004 at 0.1 m and
# rejected the point that set its own bounds (patches/lidR-grid-bounds.patch).
set.seed(1)
edge <- data.table::data.table(X = 4021833.4 + c(0, runif(499, 0, 2)), Y = 5300000 + runif(500, 0, 2),
                               Z = runif(500, 0, 3))
edge <- suppressMessages(lidR::LAS(edge))
stopifnot(sprintf("%.6f", min(edge@data$X)) == "4021833.400000",
          nrow(CspStandSegmentation::add_geometry(edge, n_cores = 1)@data) == 500L)
cat("spatial index bounds test passed\n")
