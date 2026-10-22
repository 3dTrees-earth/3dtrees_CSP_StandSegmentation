#!/usr/bin/env Rscript

point_count <- 20L
las <- lidR::LAS(data.frame(
  X = as.numeric(seq_len(point_count)),
  Y = as.numeric(seq_len(point_count) %% 3L),
  Z = as.numeric(seq_len(point_count)) / 10,
  Classification = rep(2L, point_count)
))

expected <- paste0("extra_", seq_len(12L))
for (index in seq_along(expected)) {
  las <- lidR::add_lasattribute(
    las,
    rep(as.numeric(index), point_count),
    expected[[index]],
    expected[[index]]
  )
}

path <- tempfile(fileext = ".las")
on.exit(unlink(path), add = TRUE)
lidR::writeLAS(las, path)
loaded <- lidR::readLAS(path, select = "c0")

stopifnot(all(expected %in% names(loaded@data)))
cat(sprintf("patched rlas loaded all %d extra-byte dimensions\n", length(expected)))
