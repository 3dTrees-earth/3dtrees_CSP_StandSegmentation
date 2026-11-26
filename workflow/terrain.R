# DTM cell values, read once per file-backed DTM. terra::extract needs about
# 1 GB per million points; the cell index gives the same values for a fraction.
terrain_cache <- new.env()
terrain_values <- function(dtm) {
  key <- paste(terra::sources(dtm), collapse = "|")
  if (!nzchar(key)) return(terra::values(dtm, mat = FALSE))
  if (!identical(terrain_cache$key, key)) {
    terrain_cache$values <- terra::values(dtm, mat = FALSE)
    terrain_cache$key <- key
  }
  terrain_cache$values
}

terrain_fill_cap <- function(dtm, config) {
  min(config$dtm_fill_max_distance %||% 5,
      if (is.null(config$dtm_fill_radius)) Inf else config$dtm_fill_radius * min(terra::res(dtm)))
}

new_terrain_fill <- function(dtm, config) {
  cap <- terrain_fill_cap(dtm, config)
  if (cap == 0) return(NULL)
  list(cap = cap, estimates = new.env(parent = emptyenv()))
}

# Find the closest original finite centres in growing circles. Scan one row
# at a time to bound temporary memory even for a finely spaced raster.
estimate_terrain <- function(dtm, values, cell, cap) {
  resolution <- terra::res(dtm)
  cols <- ncol(dtm)
  row <- (cell - 1) %/% cols + 1
  col <- (cell - 1) %% cols + 1
  radius <- min(cap, min(resolution))
  repeat {
    half_width <- ceiling(radius / resolution[[1L]])
    half_height <- ceiling(radius / resolution[[2L]])
    columns <- seq.int(max(1, col - half_width), min(cols, col + half_width))
    dx2 <- ((columns - col) * resolution[[1L]])^2
    # Include cells exactly at the cap when resolution arithmetic rounds up.
    best2 <- (radius * (1 + 8 * .Machine$double.eps))^2
    elevations <- numeric()
    for (r in seq.int(max(1, row - half_height), min(nrow(dtm), row + half_height))) {
      dy2 <- ((r - row) * resolution[[2L]])^2
      if (dy2 > best2) next
      offsets2 <- dx2 + dy2
      candidates <- values[(r - 1) * cols + columns]
      available <- which(offsets2 <= best2 & is.finite(candidates))
      if (!length(available)) next
      closest2 <- min(offsets2[available])
      nearest <- candidates[available[offsets2[available] == closest2]]
      elevations <- if (closest2 < best2) nearest else c(elevations, nearest)
      best2 <- closest2
    }
    if (length(elevations)) return(list(value = mean(elevations), distance = sqrt(best2)))
    if (radius == cap) {
      xy <- terra::xyFromCell(dtm, cell)[1L, ]
      abort(sprintf(paste("DTM NoData hole at cell %.0f (x=%.6f, y=%.6f) extends more than %g m",
                          "to the nearest original finite cell centre (--dtm-fill-max-distance or --dtm-fill-radius)"),
                    cell, xy[[1L]], xy[[2L]], cap))
    }
    radius <- min(cap, 2 * radius)
  }
}

terrain_under_cells <- function(dtm, values, cells, fill) {
  terrain <- values[cells]
  if (is.null(fill)) return(terrain)
  missing <- which(!is.finite(terrain))
  missing <- missing[!is.na(cells[missing])]
  if (!length(missing)) return(terrain)
  needed <- unique(cells[missing])
  estimated <- vapply(needed, function(cell) {
    key <- sprintf("%.0f", cell)
    if (!exists(key, envir = fill$estimates, inherits = FALSE)) {
      fill$estimates[[key]] <- estimate_terrain(dtm, values, cell, fill$cap)
    }
    fill$estimates[[key]]$value
  }, numeric(1))
  terrain[missing] <- estimated[match(cells[missing], needed)]
  terrain
}

terrain_fill_distances <- function(fill) {
  vapply(as.list(fill$estimates), `[[`, numeric(1), "distance")
}

merge_terrain_fills <- function(fills) {
  distances <- unlist(unname(fills))
  distances[!duplicated(names(distances))]
}

bilinear_terrain <- function(x, y, dtm, values, fill = NULL) {
  cols <- ncol(dtm)
  rows <- nrow(dtm)
  resolution <- terra::res(dtm)
  xmin <- terra::xmin(dtm)
  ymax <- terra::ymax(dtm)
  xmax <- terra::xmax(dtm)
  ymin <- terra::ymin(dtm)
  terrain <- rep(NA_real_, length(x))
  if (!length(x)) return(terrain)
  # Bound temporary vectors even when rebuilding a very large tree's grid.
  batch <- 250000L
  for (start in seq.int(1L, length(x), by = batch)) {
    i <- seq.int(start, min(start + batch - 1L, length(x)))
    inside <- is.finite(x[i]) & is.finite(y[i]) &
      x[i] >= xmin & x[i] <= xmax & y[i] >= ymin & y[i] <= ymax
    i <- i[inside]
    # Replicate the outermost centre between that centre and the raster edge.
    col <- pmin(cols, pmax(1, (x[i] - xmin) / resolution[[1L]] + 0.5))
    row <- pmin(rows, pmax(1, (ymax - y[i]) / resolution[[2L]] + 0.5))
    left <- floor(col)
    top <- floor(row)
    right <- ceiling(col)
    bottom <- ceiling(row)
    wx <- col - left
    wy <- row - top
    upper <- values[(top - 1) * cols + left] * (1 - wx) + values[(top - 1) * cols + right] * wx
    lower <- values[(bottom - 1) * cols + left] * (1 - wx) + values[(bottom - 1) * cols + right] * wx
    terrain[i] <- upper * (1 - wy) + lower * wy
    # Preserve the containing-cell fallback at NoData boundaries. Only
    # estimate terrain if the containing cell itself has no elevation.
    missing <- i[!is.finite(terrain[i])]
    if (length(missing)) {
      cells <- terra::cellFromXY(dtm, cbind(x[missing], y[missing]))
      terrain[missing] <- terrain_under_cells(dtm, values, cells, fill)
    }
  }
  terrain
}

normalize_inventory_points <- function(points, dtm, lookup = "nearest", fill = NULL) {
  values <- terrain_values(dtm)
  terrain <- if (identical(lookup, "bilinear")) {
    bilinear_terrain(points$X, points$Y, dtm, values, fill)
  } else {
    cells <- terra::cellFromXY(dtm, cbind(points$X, points$Y))
    terrain_under_cells(dtm, values, cells, fill)
  }
  if (any(!is.finite(terrain))) abort("DTM has missing or non-finite cells under inventoried points")
  data.table::set(points, j = "Zref", value = points$Z)
  data.table::set(points, j = "Z", value = points$Z - terrain)
  points
}
