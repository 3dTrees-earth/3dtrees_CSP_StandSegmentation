# Shared measurement limits; values live only in inventory_policy.json.
inventory_policy <- jsonlite::fromJSON(file.path(dirname(normalizePath(sys.frame(1)$ofile)), "inventory_policy.json"))
stopifnot(inventory_policy$schema_version == 2L,
          all(vapply(inventory_policy, function(x) length(x) == 1L && is.finite(x) && x > 0 && x == floor(x), logical(1))))
inventory_instance_point_limit <- inventory_policy$instance_point_limit
# Stem-band points kept per tree for DBH; above it, points are merged on a
# voxel grid starting at the initial size and doubling until the band fits.
inventory_dbh_band_point_budget <- inventory_policy$dbh_band_point_budget
inventory_dbh_band_initial_voxel_m <- inventory_policy$dbh_band_initial_voxel_um / 1e6
inventory_large_instance_sample_limit <- inventory_policy$large_instance_sample_limit
inventory_consumer_batch_points <- inventory_policy$consumer_batch_points
