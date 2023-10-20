#' Draw elbow plot
#'
#' The elbow plot shows the proportion of barcodes' frequency.
#' We name this elbow plot as the final plot resembles an elbow.
#'
#' @param count_data A data.table where each row corresponds to a clone barcode.
#' @param y_axis_column The column which data to display in y-axis.
#' Note, the numeric barcode ID will also be ordered based on this.
#' @param colour_column A character string, default NA. If provided, indicates the column to colour the elbow by.
#' @param facet_column A character string, default NA. If provided, indicates the column which value separates the plots into column.
#' @param facet_row A character string, default NA. If provided, indicates the column which value separates the plots into row
#'
#' @import data.table
#' @import ggplot2
#' @importFrom scales pretty_breaks
#' @return A ggplot object
#' @export
#'
#' @examples
#' toy_clone_counts <- data.table(
#'     treatment = c(rep("treat1", 4), rep("treat2", 4)),
#'     replicate = rep(c(rep("rep1", 2), rep("rep2", 2)), 2),
#'     read_count = c(20, 10, 15, 2, 3, 30, 15, 20)
#' )
#' draw_elbow_plot(
#'     count_data = toy_clone_counts,
#'     count_column = "read_count",
#'     elbow_grouping_col = "replicate",
#'     plot_grouping_col = "treatment",
#'     y_axis_column = 'read_count'
#' )
#'
draw_elbow_plot <- function(count_data, y_axis_column,
                            colour_column=NA,
                            facet_column=NA,
                            facet_row=NA) {

    # OMG this is ridiculous...
    # TODO fix me if possible
    if (is.na(facet_column) & is.na(facet_row)) {
        count_data_with_numeric_barcode_id <- assign_numeric_barcode_id(
            count_data = count_data,
            elbow_grouping_col = colour_column,
            plot_grouping_col = NA,
            count_column = y_axis_column
        )
    } else if (!is.na(facet_column) & is.na(facet_row)) {
        count_data_with_numeric_barcode_id <- assign_numeric_barcode_id(
            count_data = count_data,
            elbow_grouping_col = colour_column,
            plot_grouping_col = facet_column,
            count_column = y_axis_column
        )
    } else if (is.na(facet_column) & !is.na(facet_row)) {
        count_data_with_numeric_barcode_id <- assign_numeric_barcode_id(
            count_data = count_data,
            elbow_grouping_col = colour_column,
            plot_grouping_col = facet_row,
            count_column = y_axis_column
        )
    } else {
        # create one column concatenating the facet_row and facet_group
        count_data_copy <- copy(count_data)
        count_data[, concat_facet := paste(get(facet_column), get(facet_row), sep = "_")]
        count_data_with_numeric_barcode_id <- assign_numeric_barcode_id(
            count_data = count_data,
            elbow_grouping_col = colour_column,
            plot_grouping_col = "concat_facet",
            count_column = y_axis_column
        )
    }

    if (is.na(colour_column)) {
        plt <- ggplot(count_data_with_numeric_barcode_id,
                      aes(x = barcode_id, y = .data[[y_axis_column]]))
    } else {
        plt <- ggplot(count_data_with_numeric_barcode_id,
                      aes(x = barcode_id,
                          y = .data[[y_axis_column]],
                          colour = .data[[colour_column]]))
    }

    # Add the dots in
    plt <- plt + geom_point(size = 0.5)

    # Separate the plots into different facets if given
    # The else is essentially default where not facets are given/required.
    if (!is.na(facet_column) & !is.na(facet_row)) {
        plt <- plt + facet_grid(as.formula(paste(facet_row, " ~", facet_column)))
    } else if (is.na(facet_column) & !is.na(facet_row)) {
        plt <- plt + facet_grid(as.formula(paste(facet_row, " ~.")))
    } else if (!is.na(facet_column) & is.na(facet_row)) {
        plt <- plt + facet_grid(as.formula(paste(". ~", facet_column)))
    }

    # Makes the plot prettier.
    plt <- plt +
        scale_x_continuous(breaks = pretty_breaks(n = 10)) +
        scale_y_continuous(breaks = pretty_breaks(n = 10)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        scale_colour_viridis_d(option = "turbo") +
        labs(
            x = "Barcode ID",
            y = paste(y_axis_column, "assigned to each barcode"),
            title = paste(y_axis_column, "assigned to each barcode"),
            subtitle = paste("Barcode IDs are numerically assigned in order of", y_axis_column)
        ) +
        guides(colour = guide_legend(override.aes = list(size=6)))

    return(plt)
}




#' Assign numeric barcode id
#'
#' For each clone barcode assign a numeric id.
#' IDs are ordered based on the number of cells for each clone (clone with the largest count has small ID).
#'
#' @param count_data A data.table where each row corresponds to a clone barcode.
#' @param count_column A character indicating the column storing the count of the clone barcodes.
#' @param plot_grouping_col A character string, default NA, indicating the column which values separates each plot.
#' @param elbow_grouping_col A character string, default NA, indicating the column which values separates each elbow.
#'
#' @return A data.table with numerical ID for each clone barcode.
#'
assign_numeric_barcode_id <- function(count_data, count_column,
                                      elbow_grouping_col = NA,
                                      plot_grouping_col = NA) {

    if (is.na(elbow_grouping_col)) {
        # Means we only need 1 plot and 1 elbow.
        if (is.na(plot_grouping_col)) {
            res <- copy(count_data)
            res <- res[order(-get(count_column))]
            res$barcode_id <- seq(0, nrow(res)-1)
        } else {
            res <- lapply(split(count_data, by = c(plot_grouping_col)), function(chunk) {
                ord_chunk <- chunk[order(-get(count_column))]
                ord_chunk$barcode_id <- seq(0, nrow(ord_chunk)-1)
                return(ord_chunk)
            })
            res <- rbindlist(res)
        }

    } else {
        # Means we have multiple elbows
        if (is.na(plot_grouping_col)) {
            # We only have 1 plot
            # So we need to break per elbow, but for the next elbow we need to
            # continue the ID from the previous elbow.
            res <- assign_numeric_barcode_id_multiple_elbows(
                all_elbow_data = count_data,
                elbow_grouping_col = elbow_grouping_col,
                count_column = count_column
            )
        } else {
            # the worst nightmare, multiple elbows and multiple plots
            split_plot_data <- split(count_data, by = c(plot_grouping_col))

            res <- list()

            for (i in seq(length(split_plot_data))) {
                a_plot_data <- split_plot_data[[i]]
                res_all_elbows <- assign_numeric_barcode_id_multiple_elbows(
                    all_elbow_data = a_plot_data,
                    elbow_grouping_col = elbow_grouping_col,
                    count_column = count_column
                )
                res[[i]] <- res_all_elbows
            }
            res <- rbindlist(res)
        }
    }
    return(res)
}

#' Assign numeric barcode id for all elbows
#'
#' If a plot has multiple elbows, assign the numerical barcode id such that
#' for the first elbow, the barcode id starts from 0,
#' then the subsequent elbows, the barcode id starts from the end of the previous elbow.
#'
#' @param all_elbow_data A data.table storing the data for all elbows.
#' @param elbow_grouping_col A column in the data.table that separate each elbow.
#' @param count_column A character indicating the column storing the count of the clone barcodes.
#'
#' @return A new data.table
#'
assign_numeric_barcode_id_multiple_elbows <- function(all_elbow_data, elbow_grouping_col, count_column) {

    # We need to break per elbow, but for the next elbow we need to
    # continue the ID from the previous elbow.

    start_id <- 0

    split_elbow_data <- split(all_elbow_data, by = c(elbow_grouping_col))

    res <- list()

    for (i in seq_len(length(split_elbow_data))) {
        an_elbow_data <- split_elbow_data[[i]]

        ordered_elbow_data <- an_elbow_data[order(-get(count_column))]
        end_id <- start_id + nrow(ordered_elbow_data) - 1
        ordered_elbow_data[, barcode_id := seq(start_id, end_id)]
        start_id <- end_id + 1
        res[[i]] <- ordered_elbow_data

    }
    res <- rbindlist(res)
}



