#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = FALSE)
test_file <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[[1L]])
repository_root <- normalizePath(file.path(dirname(test_file), "..", ".."))
source(file.path(repository_root, "workflow", "cli.R"))
source(file.path(repository_root, "workflow", "tool.R"))
source(file.path(repository_root, "workflow", "outputs.R"))

work <- tempfile("inventory-local-crs-")
dir.create(work)
script <- file.path(work, "fixture.py")
wkt_path <- file.path(work, "standin.wkt")
bound_ecef_wkt <- sf::st_crs("+proj=geocent +ellps=GRS80 +towgs84=1,2,3,0,0,0,0 +units=m")$wkt
bound_ecef_path <- file.path(work, "bound-ecef.wkt")
writeLines(bound_ecef_wkt, bound_ecef_path)
writeLines(sub(',\\s*ID\\["EPSG",4978\\]\\]$', "]", sf::st_crs(4978)$wkt), wkt_path)
writeLines(c(
  "import copy, struct, sys",
  "from pathlib import Path",
  "import laspy",
  "import numpy as np",
  "from itertools import product",
  "root = Path(sys.argv[1])",
  "header = laspy.LasHeader(point_format=6, version='1.4')",
  "header.scales = [0.0001] * 3",
  "header.add_extra_dim(laspy.ExtraBytesParams(name='PredInstance_FM', type='int32'))",
  "ground = list(product(range(13), range(13), [0.]))",
  "tree = list(product([4.8, 5.0, 5.2], [4.8, 5.0, 5.2], np.arange(0.5, 12., 0.05)))",
  "xyz = np.array(ground + tree)",
  "for name in ('empty', 'standin', 'standin-wkt', 'standin-wrapped', 'user-projected',",
  "             'ecef', 'bound-ecef', 'geographic', 'degrees', 'high', 'range', 'unlabelled-ecef', 'utm57'):",
  "    out = laspy.LasData(copy.deepcopy(header))",
  "    x, y, z = xyz.T.copy()",
  "    if name == 'degrees':",
  "        x, y = 7.8 + x / 1000, 48.0 + y / 1000",
  "    if name == 'high':",
  "        z += 12000",
  "    if name == 'range':",
  "        z = (z - 6) * 1000",
  "    if name == 'unlabelled-ecef':",
  "        out.header.offsets = [6378137., 0., 0.]",
  "        x += 6378137",
  "    if name == 'utm57':",
  "        out.header.offsets = [400000., 6379994., 0.]",
  "        x, y, z = x + 400000, y + 6379994, z * (60 / z.max())",
  "    if name in ('ecef', 'bound-ecef'):",
  "        out.header.offsets = [0., 0., 6356752.]",
  "        z += 6356752",
  "    if name == 'standin-wrapped':",
  "        wrap_y = 14 * 2**32 * header.scales[1]",
  "        out.header.offsets = [0., wrap_y, 0.]",
  "        y += wrap_y",
  "    out.x, out.y, out.z = x, y, z",
  "    out.classification = [2] * len(ground) + [1] * len(tree)",
  "    out.PredInstance_FM = [0] * len(ground) + [1] * len(tree)",
  "    if name in ('standin-wkt', 'bound-ecef'):",
  "        wkt = Path(sys.argv[2 if name == 'standin-wkt' else 3]).read_bytes()",
  "        out.vlrs.append(laspy.VLR('LASF_Projection', 2112, record_data=wkt + b'\\x00'))",
  "    elif name in ('empty', 'degrees', 'high', 'range', 'unlabelled-ecef', 'utm57'):",
  "        out.vlrs.append(laspy.VLR('LASF_Projection', 2112, record_data=b''))",
  "    elif name == 'user-projected':",
  "        keys = struct.pack('<16H', 1, 1, 0, 3, 1024, 0, 1, 1, 2048, 0, 1, 4269, 3072, 0, 1, 32767)",
  "        out.vlrs.append(laspy.VLR('LASF_Projection', 34735, record_data=keys))",
  "    else:",
  "        model, epsg = (2, 4326) if name == 'geographic' else (3, 4978)",
  "        keys = struct.pack('<12H', 1, 1, 0, 2, 1024, 0, 1, model, 2048, 0, 1, epsg)",
  "        out.vlrs.append(laspy.VLR('LASF_Projection', 34735, record_data=keys))",
  "    out.write(root / (name + '.las'))"
), script)
stopifnot(system2("python3", shQuote(c(script, work, wkt_path, bound_ecef_path))) == 0L)

dtm <- terra::rast(xmin = -2, xmax = 20, ymin = -2, ymax = 20, resolution = 1, crs = "")
terra::values(dtm) <- 0
dtm_path <- file.path(work, "supplied.tif")
terra::writeRaster(dtm, dtm_path)

local_config <- function(case, output, ...) {
  input <- file.path(work, paste0(case, ".las"))
  parse_cli_args(c("--input", input, "--output-dir", output,
                   "--instance-attribute", "PredInstance_FM",
                   "--inventory-workers", "1", "--dtm-workers", "1",
                   "--debug-outputs", "true", ...))
}

check_local <- function(case, supplied, label, epsg) {
  output <- file.path(work, paste0(case, if (supplied) "-supplied" else "-generated"))
  config <- local_config(case, output, if (supplied) c("--dtm", dtm_path))
  run_tool(config)
  rows <- utils::read.csv(file.path(output, "inventory.csv"))
  summary <- utils::read.csv(file.path(output, "stand_summary.csv"))
  raster <- terra::rast(file.path(output, "dtm_full.tif"), guessCRS = FALSE)
  manifest <- jsonlite::fromJSON(file.path(output, "manifest.json"))
  stopifnot(nrow(rows) > 0L, all(is.na(rows$longitude)), all(is.na(rows$latitude)),
            identical(summary$crs, label),
            identical(manifest$coordinate_reference$source_frame, "local"),
            is.null(manifest$coordinate_reference$wgs84_epsg))
  stopifnot(if (is.na(epsg)) terra::crs(raster) == "" else is_inventory_standin(sf::st_crs(terra::crs(raster))))
}
for (supplied in c(TRUE, FALSE)) {
  check_local("empty", supplied, "local (assumed metres)", NA_integer_)
  check_local("standin", supplied, "local (EPSG:4978 stand-in)", 4978L)
}

# This northing is in the Earth's radius band. Its easting and elevations
# are plausible UTM coordinates, so validate them through the CLI workflow.
check_local("utm57", FALSE, "local (assumed metres)", NA_integer_)
utm57 <- point_cloud_bounds_python(file.path(work, "utm57.las"), work)
stopifnot(mean(c(utm57$ymin, utm57$ymax)) == 6380000, utm57$zmin == 0, utm57$zmax == 60)

projected <- dtm
terra::crs(projected) <- "EPSG:25833"
projected_path <- file.path(work, "projected.tif")
terra::writeRaster(projected, projected_path)
for (case in c("empty", "standin")) {
  output <- file.path(work, paste0(case, "-projected"))
  config <- local_config(case, output, "--dtm", projected_path)
  refusal <- tryCatch({ run_tool(config); "" }, error = conditionMessage)
  stopifnot(grepl("--dtm-crs-local-ok true", refusal, fixed = TRUE), !dir.exists(output))
  config$dtm_crs_local_ok <- TRUE
  logs <- capture.output(run_tool(config), type = "message")
  stopifnot(any(startsWith(logs, "WARNING: --dtm-crs-local-ok dropped the projected DTM CRS")))
  manifest <- jsonlite::fromJSON(file.path(output, "manifest.json"))
  metadata <- jsonlite::fromJSON(file.path(output, "run_metadata.json"))
  stopifnot(sf::st_crs(manifest$dtm$dropped_crs_wkt) == sf::st_crs(25833),
            identical(metadata$dtm_dropped_crs_wkt, manifest$dtm$dropped_crs_wkt))
}

standin_dtm <- dtm
terra::crs(standin_dtm) <- paste(readLines(wkt_path), collapse = "\n")
standin_dtm_path <- file.path(work, "standin.tif")
terra::writeRaster(standin_dtm, standin_dtm_path)
aoi_path <- file.path(work, "aoi.geojson")
writeLines('{"type":"Polygon","coordinates":[[[1,1],[11,1],[11,11],[1,11],[1,1]]]}', aoi_path)
for (case in c("empty", "standin-wkt")) {
  aoi_output <- file.path(work, paste0(case, "-aoi"))
  withCallingHandlers(
    run_tool(local_config(case, aoi_output, "--dtm", standin_dtm_path, "--aoi-json", aoi_path)),
    warning = function(warning) {
      if (grepl("CRS do not match", conditionMessage(warning), fixed = TRUE)) stop(warning)
    }
  )
  cropped <- terra::rast(file.path(aoi_output, "dtm_aoi.tif"), guessCRS = FALSE)
  summary <- utils::read.csv(file.path(aoi_output, "stand_summary.csv"))
  stopifnot(nrow(cropped) == 10L, ncol(cropped) == 10L,
            identical(summary$area_source, "aoi"),
            identical(summary$crs, if (case == "empty") "local (assumed metres)" else "local (EPSG:4978 stand-in)"))
}

refusals <- c(ecef = "geocentric CRS", `bound-ecef` = "geocentric CRS",
              geographic = "projected with metre units",
              `standin-wrapped` = "geocentric CRS",
              `user-projected` = "without a supported projected EPSG code or WKT",
              degrees = "spans less than one unit", high = "exceed 10,000 units",
              range = "exceed 10,000 units", `unlabelled-ecef` = "resemble geocentric (ECEF)")
for (case in names(refusals)) {
  output <- file.path(work, paste0(case, "-refused"))
  refusal <- tryCatch({ run_tool(local_config(case, output)); "" }, error = conditionMessage)
  stopifnot(grepl(refusals[[case]], refusal, fixed = TRUE), !dir.exists(output))
}
wrapped_bounds <- point_cloud_bounds_python(file.path(work, "standin-wrapped.las"), work)
stopifnot(is_inventory_standin(wrapped_bounds$crs), wrapped_bounds$ymin > 6000000)

bound_3d <- sf::st_crs("+proj=utm +zone=32 +ellps=GRS80 +towgs84=1,2,3,0,0,0,0 +units=m +vunits=m")$wkt
bound_2d <- horizontal_wkt(bound_3d)
ecef_wkt <- sf::st_crs(4978)$wkt
stopifnot(startsWith(bound_3d, "BOUNDCRS["), grepl("CS[Cartesian,3]", bound_3d, fixed = TRUE),
          startsWith(bound_2d, "BOUNDCRS["), grepl("CS[Cartesian,2]", bound_2d, fixed = TRUE),
          !grepl("ORDER[3]", bound_2d, fixed = TRUE), !is.na(sf::st_crs(bound_2d)),
          identical(horizontal_wkt(ecef_wkt), ecef_wkt),
          startsWith(bound_ecef_wkt, "BOUNDCRS["),
          identical(horizontal_wkt(bound_ecef_wkt), bound_ecef_wkt))
cat("Local CRS supplied/generated terrain and genuine CRS refusal checks passed\n")
