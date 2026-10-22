# Shared measurement limits; values live only in inventory_policy.json.
inventory_policy <- jsonlite::fromJSON(file.path(dirname(normalizePath(sys.frame(1)$ofile)), "inventory_policy.json"))
stopifnot(inventory_policy$schema_version == 1L,
          all(vapply(inventory_policy, function(x) length(x) == 1L && is.finite(x) && x > 0 && x == floor(x), logical(1))))
inventory_instance_point_limit <- inventory_policy$instance_point_limit
inventory_dbh_geometry_limit <- inventory_policy$dbh_geometry_limit
inventory_large_instance_sample_limit <- inventory_policy$large_instance_sample_limit
inventory_consumer_batch_points <- inventory_policy$consumer_batch_points
