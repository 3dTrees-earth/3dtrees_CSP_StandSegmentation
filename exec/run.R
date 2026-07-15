#!/usr/bin/env Rscript

script_path <- function() {
  arguments <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", arguments, value = TRUE)
  if (!length(file_argument)) return(getwd())
  dirname(normalizePath(sub("^--file=", "", file_argument[[1L]]), mustWork = TRUE))
}

repository_root <- normalizePath(file.path(script_path(), ".."), mustWork = TRUE)
source(file.path(repository_root, "workflow", "cli.R"))
source(file.path(repository_root, "workflow", "tool.R"))

if ("--version" %in% commandArgs(trailingOnly = TRUE)) {
  description <- read.dcf(file.path(repository_root, "DESCRIPTION"))
  cat(sprintf("CspStandSegmentation %s\n", description[1L, "Version"]))
  quit(status = 0L)
}

required_packages <- c("CspStandSegmentation", "data.table", "jsonlite", "lidR", "sf", "terra")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  abort(sprintf("Missing required R packages: %s", paste(missing_packages, collapse = ", ")))
}

tryCatch(
  run_tool(parse_cli_args()),
  error = function(error) {
    message(sprintf("ERROR: %s", conditionMessage(error)))
    quit(status = 1L)
  }
)
