#' Calculate clone projection
#'
#' Predicts the number of cells tagged with each clone barcode for a given projection amount.
#' This function is useful for planning single-cell experiments.
#' It allows you to estimate the number of cells associated with specific clone barcodes
#' and anticipate the abundance of clone barcodes (i.e., how many clone barcodes
#' with a substantial number of cells).
#'
#' How this projection is done: For each clone barcode, we compute their proportion,
#' and then multiply the proportion by the number of cells to project to
#' (specified by the `project_amnt` parameter).
#' If `grouping_col` is specified, the proportion will be calculated with respect
#' to the total number of cells in each group.
#'
#' @param count_data A data.table containing clone barcode count data.
#' @param grouping_col A character string, default NA. If provided, the projection
#' is repeated for each group defined by this column.
#' @param count_column A character string specifying the column storing clone barcode counts.
#' @param project_amnt A numeric value indicating the number of cells to project to.
#'
#' @return A modified data.table with projected cell counts.
#' @export
#'
#' @examples
#' toy_clone_counts <- data.table(
#' sample_name = c(rep("test1", 3), rep("test2", 3)),
#' read_count = c(1, 5, 20, 10, 15, 12),
#' clone_barcodes = c("ACGT", "CATG", "CATG", "ACGT", "ATGC", "TCGT")
#' )
#'
#' projecting_clones(
#'     count_data = toy_clone_counts,
#'     project_amnt = 100,
#'     count_column = "read_count",
#'     grouping_col = "sample_name"
#' )

projecting_clones <- function(count_data, count_column,
                              grouping_col = NA,
                              project_amnt = 10000) {

    count_data_proportion <- convert_count_to_proportion(
        count_data = count_data,
        grouping_col = grouping_col,
        count_column = count_column
    )
    count_data_proportion[, projection := read_proportion * project_amnt]
    setnames(count_data_proportion, "read_proportion", paste0(count_column, "_proportion"))
    setnames(count_data_proportion, "projection", paste0("projected_to_", project_amnt))

    return(count_data_proportion)
}

#' Convert count to proportion
#'
#' For each clone barcode, convert the absolute count to proportion.
#' If `grouping_col` is specified, the proportion will be calculated with respect
#' to the total number of cells in each group.
#'
#' @param count_data A data.table containing clone barcode count data.
#' @param grouping_col A character string, default NA. If provided, the projection
#' is repeated for each group defined by this column.
#' @param count_column A character string specifying the column storing clone barcode counts.
#'
#' @return A modified data.table with a new column representing the proportion.
#' @export
#'
#' @examples
#' toy_clone_counts <- data.table(
#' sample_name = c(rep("test1", 2), rep("test2", 3)),
#' read_count = c(1, 5, 10, 15, 20)
#' )
#'
#' convert_count_to_proportion(
#'     count_data = toy_clone_counts,
#'     grouping_col = "sample_name",
#'     count_column = "read_count"
#' )
convert_count_to_proportion <- function(count_data, count_column, grouping_col = NA) {
    res <- copy(count_data)

    if (is.na(grouping_col)) {
        res[, read_proportion := get(count_column) / sum(get(count_column))]
    } else {
        res[, read_proportion := get(count_column) / sum(get(count_column)), by = c(grouping_col)]
    }

    return(res)
}
