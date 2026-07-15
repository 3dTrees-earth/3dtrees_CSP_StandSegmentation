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

spec <- parse_segmentation_spec("PredInstance_FM,species_id_FM,species_prob_FM")
stopifnot(
  identical(spec$source, "FM"),
  identical(spec$species, "species_id_FM"),
  identical(spec$species_prob, "species_prob_FM")
)
expect_error(parse_segmentation_spec("PredInstance_FM,species_id_FM"), "supplied together")
expect_error(parse_segmentation_spec("bad-name"), "valid LAS dimension")

config <- parse_cli_args(c(
  "--input", "input.laz",
  "--output-dir", "output",
  "--segmentation-spec", "PredInstance_SAT",
  "--segmentation-spec", "PredInstance_FM,species_id_FM,species_prob_FM",
  "--dtm-resolution", "0.2",
  "--routing-workers", "1"
))
stopifnot(
  length(config$segmentation_specs) == 2L,
  identical(config$dtm_resolution, 0.2),
  identical(config$routing_workers, 1L),
  identical(config$read_chunk_size, 300),
  identical(config$chunk_buffer, 5),
  identical(config$inventory_partitions, 64L),
  identical(config$dtm_workers, 10L),
  identical(config$dtm_strategy, "auto"),
  identical(config$dtm_candidate_resolution, 0.1),
  identical(config$dtm_streaming_threshold, 50000000L),
  !config$enable_csp
)

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
  identical(
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
    "c234567"
  ),
  identical(build_read_selector("late_dimension", c(rep("unused", 9L), "late_dimension")), "c0")
)

expect_error(
  parse_cli_args(c("--input", "input.laz", "--output-dir", "output")),
  "segmentation-spec"
)
expect_error(
  parse_cli_args(c(
    "--input", "input.laz", "--output-dir", "output",
    "--segmentation-spec", "PredInstance_SAT", "--dtm-strategy", "invalid"
  )),
  "dtm-strategy"
)
expect_error(
  parse_cli_args(c(
    "--input", "input.laz", "--output-dir", "output",
    "--segmentation-spec", "PredInstance_SAT", "--dtm-resolution", "0.2",
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
  fm_one$pred_score_fm == 0.9,
  isTRUE(fm_one$pred_score_mixed),
  fm_one$wood_point_count == 1,
  fm_one$leaf_point_count == 2,
  abs(fm_one$wood_share - 1 / 3) < 1e-12,
  abs(fm_one$leaf_share - 2 / 3) < 1e-12
)

cat("tool helper tests passed\n")
