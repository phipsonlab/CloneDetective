#' Calculate clone barcode with at least x number of cells
#'
#' This function calculates the number of clones with at least a specified number of cells.
#' Clones with fewer cells than the threshold are excluded from the count.
#'
#' @param count_data A data.table where each row corresponds to a clone barcode.
#' @param thresholds A numeric vector listing the minimum cell counts required for a clone to be included in the count.
#' The calculation will be applied for each value in the vector.
#' @param count_column A character indicating the column storing the count of the clone barcodes.
#' @param grouping_col A character string, default NA. If provided, the function repeats the calculation for each subsets of `count_data` defined by this column, effectively a group-by operation.
#'
#' @return A data.table containing the counts of clones that meet or exceed the specified threshold.
#' @import data.table
#'
#' @examples
#' library(data.table)
#'
#' toy_clone_counts <- data.table(
#' sample_name = c(rep("test1", 2), rep("test2", 3)),
#' clone_barcode = c("A", "B", "C", "D", "E"),
#' read_count = c(1, 5, 10, 15, 20)
#' )
#'
#' count_retained_clones(
#'   count_data = toy_clone_counts,
#'   thresholds = c(2, 10),
#'   count_column = "read_count",
#'   grouping_col = "sample_name"
#' )
#'
#' @export
count_retained_clones <- function(count_data, thresholds, count_column,
                                  grouping_col = NA) {
  count_data_dt <- as.data.table(count_data)

  clone_counts <- lapply(thresholds, function(thres) {
    count_data_subset <- count_data_dt[get(count_column) >= thres, ]

    if (!is.na(grouping_col)) {
      n_clones <- count_data_subset[, .(cnt = .N), by = c(grouping_col)]

      # Be explicit when the count is zero
      samples <- data.table(fake_samp = unique(count_data_dt[[grouping_col]]))
      setnames(samples, "fake_samp", grouping_col)

      n_clones <- merge.data.table(samples, n_clones, all.x = TRUE)
      n_clones[is.na(n_clones)] <- 0L

    } else {
      n_clones <- count_data_subset[, .(cnt = .N), ]

      # Be explicit when the count is zero
      if (nrow(n_clones) == 0) {
        n_clones < data.table(cnt = c(0L))
      }

      # have to do this so the merge later can find a common column
      n_clones[, col_for_merging := rep("dummy", nrow(n_clones))]
    }

    setnames(n_clones, "cnt", paste0("at_least_", thres, "_cells"))


    return(n_clones)
  })
  clone_counts <- Reduce(function(...) merge(..., all = TRUE), clone_counts)

  # if there no such column, it will just do nothing
  if (is.na(grouping_col)) {
    clone_counts[, col_for_merging := NULL]
  }

  return(clone_counts)
}


#' Remove clones with counts below a specified threshold
#'
#' Removes clones with counts strictly less than the specified threshold.
#'
#' @param count_data A data.table containing clone barcodes and their counts.
#' One row per clone barcode.
#' @param threshold A numeric value indicating the minimum count required for a
#' clone barcode to be retained.
#' @param count_column A character string specifying the column name for
#' clone barcode counts.
#'
#' @return A data.table with clones above the threshold.
#' @import data.table
#'
#' @examples
#' library(data.table)
#'
#' toy_clone_counts <- data.table(
#' sample_name = c(rep("test1", 2), rep("test2", 3)),
#' read_count = c(1, 5, 10, 15, 20),
#' clone_barcodes = c("ACGT", "CATG", "CATG", "ACGT", "ATGC")
#' )
#'
#' remove_clones_below_threshold(
#'   count_data = toy_clone_counts,
#'   threshold = 2,
#'   count_column = "read_count"
#' )
#'
#' @export
#'
remove_clones_below_threshold <- function(count_data, threshold, count_column) {
  count_data_dt <- as.data.table(count_data)
  res <- count_data_dt[get(count_column) >= threshold, ]
  return(res)
}


#' Get top barcodes and cumulative sum of the counts.
#'
#' Filter the count data to retain only the top x barcodes, where x is defined by top_threshold.
#'
#' @param count_data A data.table where each row corresponds to a clone barcode.
#' @param top_threshold A numeric value indicating how many of the top clone barcodes to retain.
#' @param count_column A character string indicating the column storing the count of the clone barcodes.
#' @param grouping_col A character string, default NA. If provided, the function repeats the calculation for each subset of `count_data` defined by this column, effectively performing a group-by operation.
#'
#' @return A filtered data.table with only the top clone barcodes.
#' @export
#'
#' @examples
#' library(data.table)
#'
#' toy_clone_counts <- data.table(
#' sample_name = c(rep("test1", 3), rep("test2", 3)),
#' read_count = c(1, 5, 20, 10, 15, 12),
#' clone_barcodes = c("ACGT", "CATG", "CATG", "ACGT", "ATGC", "TCGT")
#' )
#'
#' get_top_barcodes_and_cum_sum(
#'   count_data = toy_clone_counts,
#'   top_threshold = 2,
#'   count_column = "read_count",
#'   grouping_col = "sample_name"
#' )

get_top_barcodes_and_cum_sum <- function(count_data, top_threshold, count_column,
                             grouping_col = NA) {

  if (is.na(grouping_col)) {
    count_data_top <- count_data[order(-get(count_column)), head(.SD, top_threshold)]
    count_data_top[, cum_sum := cumsum(get(count_column))]
    count_data_top[, barcode_rank := seq_len(.N)]
  } else {
    count_data_top <- count_data[order(-get(count_column)), head(.SD, top_threshold), by=c(grouping_col)]
    count_data_top[, cum_sum := cumsum(get(count_column)), by=c(grouping_col)]
    count_data_top[, barcode_rank := seq_len(.N), by=c(grouping_col)]
  }

  setnames(count_data_top, "cum_sum", paste0("cum_sum_", count_column))

  return(count_data_top)
}






