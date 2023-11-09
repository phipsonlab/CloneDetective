#' Calculate clone projection
#'
#' Predicts the number of cells tagged with each clone barcode for a given projection amount.
#' This function is useful for planning single-cell experiments.
#' It allows you to estimate the number of cells associated with specific clone barcodes
#' and anticipate the abundance of clone barcodes (i.e., how many clone barcodes
#' with a substantial number of cells).
#'
#' @details
#' How this projection is done: For each clone barcode, we compute their proportion,
#' and then multiply the proportion by the number of cells to project to
#' (specified by the `project_amnt` parameter).
#' If `grouping_col` is specified, the proportion will be calculated with respect
#' to the total number of cells in each group.
#'
#' The parameter `project_amnt` allows you to specify the desired number of cells
#' for projection.
#' You can provide a numeric vector containing one or more values, indicating the
#' number of cells you want to project your data onto.
#' For example, if you want to project to 10,000 and 20,000 cells,
#' you can create a vector `c(10000, 20000)` and pass it as the value for `project_amnt`.
#'
#'
#' @param count_data A data.table containing clone barcode count data.
#' @param grouping_col A character string, default NA. If provided, the projection
#' is repeated for each group defined by this column.
#' @param count_column A character string specifying the column storing clone barcode counts.
#' @param project_amnt A vector of numeric indicating the number of cells to project to.
#' @param confidence_threshold A numeric value indicating the confidence level for
#' clone barcode detection in cells.
#' It should be between 0 (exclusive) and 1 (inclusive).
#' This parameter allows you to specify how confident you are that you will be able to
#' detect clone barcodes for the cells in future scRNAseq experiment.
#' It ranges from 0 (indicating no confidence) to 1 (indicating complete confidence).
#' For example, if you are 70\% confident that all cells will have their clone barcodes
#' detected, you should set the `confidence_threshold` to 0.7.
#'
#'
#' @return A modified data.table with projected cell counts.
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
#' res <- project_clones(
#'     count_data = toy_clone_counts,
#'     project_amnt = c(100),
#'     count_column = "read_count",
#'     grouping_col = "sample_name"
#' )
#'
#' res
#'
#' @import stats
#'

project_clones <- function(count_data, count_column,
                           grouping_col = NA,
                           project_amnt = c(10000),
                           confidence_threshold = 1.0) {

    count_data_proportion <- convert_count_to_proportion(
        count_data = count_data,
        grouping_col = grouping_col,
        count_column = count_column
    )


    for (amnt in project_amnt) {
        count_data_proportion[, projection := read_proportion * amnt]

        # sample the count using binomial sampling with confidence of confidence_threshold
        # if it is less than 1.0.
        # this is to simulate cases where there are cells with no clones detected.

        if (confidence_threshold < 1.0) {
            # need to do rounding as rbinom is stupid and need integers
            count_data_proportion[, projection := round(projection)]
            projected_with_conf <- rbinom(
                n = nrow(count_data_proportion),
                size = count_data_proportion$projection,
                prob = confidence_threshold
            )
            count_data_proportion$projection <- projected_with_conf
        }

        setnames(
            x = count_data_proportion,
            old = "projection",
            new = paste0("projected_", amnt, "_confidence_", gsub("\\.", "_", as.character(confidence_threshold)))
        )
    }

    setnames(count_data_proportion, "read_proportion", paste0(count_column, "_proportion"))

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
#' library(data.table)
#'
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
#'
convert_count_to_proportion <- function(count_data, count_column, grouping_col = NA) {
    res <- copy(count_data)

    if (is.na(grouping_col)) {
        res[, read_proportion := get(count_column) / sum(get(count_column))]
    } else {
        res[, read_proportion := get(count_column) / sum(get(count_column)), by = c(grouping_col)]
    }

    return(res)
}

