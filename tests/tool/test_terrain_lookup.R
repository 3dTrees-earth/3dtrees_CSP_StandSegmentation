#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = FALSE)
test_file <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[[1L]])
root <- normalizePath(file.path(dirname(test_file), "..", ".."))
source(file.path(root, "workflow", "cli.R"))
source(file.path(root, "workflow", "terrain.R"))

dtm <- terra::rast(xmin = 0, xmax = 2, ymin = 0, ymax = 2, nrows = 2, ncols = 2, crs = "")
terra::values(dtm) <- c(0, 4, 8, 12)
points <- data.table::data.table(X = c(0, 0.5, 0.75, 1, 2),
                                Y = c(2, 1.5, 1.25, 1, 0), Z = 20)
nearest <- normalize_inventory_points(data.table::copy(points), dtm)
bilinear <- normalize_inventory_points(data.table::copy(points), dtm, "bilinear")
stopifnot(identical(nearest$Z, c(20, 20, 20, 8, 8)),
          identical(bilinear$Z, c(20, 20, 17, 14, 8)), all(bilinear$Zref == 20))
terra::values(dtm) <- c(0, NA, 8, 12)
stopifnot(identical(bilinear_terrain(c(0.5, 0.75, 1.25, -1), c(1.5, 1.25, 1.25, 1), dtm,
                                   terrain_values(dtm)), c(0, 0, NA_real_, NA_real_)))
fill <- new_terrain_fill(dtm, list())
filled <- bilinear_terrain(c(0.5, 0.75, 1.25, -1), c(1.5, 1.25, 1.25, 1), dtm,
                           terrain_values(dtm), fill)
stopifnot(identical(filled, c(0, 0, 6, NA_real_)),
          identical(terrain_fill_distances(fill), c(`2` = 1)))
cat("Nearest and bilinear terrain lookup checks passed\n")
