#!/usr/bin/env Rscript

# Local benchmark helper: consume completed FM partitions and a supplied DTM.
args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(3L, 4L))) {
  stop("Usage: run_partition_inventory.R PARTITIONS DTM OUTPUT_CSV [PARTITION_BASENAME]")
}
arguments <- commandArgs(trailingOnly = FALSE)
test_file <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[[1L]])
root <- normalizePath(file.path(dirname(test_file), "..", ".."), mustWork = TRUE)
source(file.path(root, "workflow", "cli.R"))
tool_file <- if (identical(Sys.getenv("CSP_BASELINE_TOOL"), "1")) {
  "/opt/CspStandSegmentation/workflow/tool.R"
} else {
  file.path(root, "workflow", "tool.R")
}
source(tool_file)
lidR::set_lidr_threads(1L)
set.seed(42L)

partition_root <- normalizePath(args[[1L]], mustWork = TRUE)
layout_dir <- file.path(partition_root, "FM")
if (length(args) == 4L) {
  selected <- args[[4L]]
  if (basename(selected) != selected || !grepl("^[0-9]{5}\\.bin$", selected)) {
    stop("PARTITION_BASENAME must be a five-digit .bin name")
  }
  source_partition <- file.path(layout_dir, selected)
  if (!file.exists(source_partition)) stop("Selected partition does not exist")
  layout_dir <- tempfile("csp-selected-partition-")
  dir.create(layout_dir)
  if (!file.symlink(source_partition, file.path(layout_dir, selected))) {
    stop("Could not select partition for benchmark")
  }
}
dtm <- terra::rast(args[[2L]])
summary <- jsonlite::fromJSON(file.path(partition_root, "partition_summary.json"),
                               simplifyVector = TRUE)
spec <- list(source = "FM", instance = "PredInstance_FM",
             species = "species_id_FM", species_prob = "species_prob_FM")
columns <- c("X", "Y", "Z", "PredInstance_FM", "species_id_FM",
             "species_prob_FM", "PredScore_FM", "PredSemantic_FM")
stream <- list(layouts = list(FM = list(spec = spec, columns = columns,
                                       directory = layout_dir)),
               largest_instances = summary$largest_instances,
               oversized_instances = summary$oversized_instances)
workers <- as.integer(Sys.getenv("CSP_INVENTORY_WORKERS", "1"))
if (is.na(workers) || workers < 1L || workers > 20L) stop("Invalid CSP_INVENTORY_WORKERS")
config <- list(inventory_workers = workers, random_seed = 42L, non_tree_ids = 0,
               slice_min = 0.3, slice_max = 4, slice_increment = 0.2,
               slice_width = 0.1, max_dbh = 1, geometry_threads = 1L)
started <- Sys.time()
inventories <- inventory_from_partitions(stream, dtm, config, NULL,
                                         sf::st_crs(terra::crs(dtm)))
data.table::fwrite(inventories$FM, args[[3L]])
cat(sprintf("inventory_elapsed_seconds=%.3f trees=%d points=%.0f\n",
            as.numeric(difftime(Sys.time(), started, units = "secs")),
            nrow(inventories$FM), sum(inventories$FM$point_count)))
if (file.exists("/sys/fs/cgroup/memory.peak")) {
  cat(sprintf("cgroup_memory_peak_bytes=%s\n",
              readLines("/sys/fs/cgroup/memory.peak", n = 1L)))
}
if (file.exists("/sys/fs/cgroup/cpu.stat")) {
  cpu <- readLines("/sys/fs/cgroup/cpu.stat")
  cat(cpu[grepl("^usage_usec ", cpu)], "\n", sep = "")
}
