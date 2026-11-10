# Output tables: inventory.csv and inventory.parquet with the same rows and
# columns, and stand_summary.csv. Rows are sorted by
# (instance_dimension, instance_id) and Parquet row groups hold at most
# `parquet_row_group_rows` rows, so the viewer fetches one tree's row group by
# range using the row-group min/max statistics.

parquet_row_group_rows <- 2048L
instance_sort_keys <- c("instance_dimension", "instance_id")

# Sort in place by the instance keys the table has, then `extra` keys; returns
# the keys used. data.table sorts strings in C locale and is stable, so equal
# keys keep their measured order.
sort_by_instance <- function(table, extra = character()) {
  keys <- unique(c(intersect(instance_sort_keys, names(table)), intersect(extra, names(table))))
  if (length(keys) && nrow(table)) data.table::setorderv(table, keys)
  keys
}

# Snappy, not zstd: the viewer's Parquet reader has no zstd decoder. arrow's
# R-specific "r" schema metadata is dropped so every reader sees the same file.
write_parquet_table <- function(table, path) {
  arrow_table <- arrow::Table$create(as.data.frame(table))
  arrow_table$metadata$r <- NULL
  arrow::write_parquet(arrow_table, path, chunk_size = parquet_row_group_rows, compression = "snappy")
}

# data.table::fwrite writes doubles to 15 significant digits and NA as "".
write_csv_table <- function(table, path) data.table::fwrite(table, path, sep = ",", na = "")
