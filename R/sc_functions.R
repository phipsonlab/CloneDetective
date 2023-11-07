#' Title
#'
#' @param cell_clone_bcode_dt
#' @param cell_bcode_col
#' @param clone_bcode_col
#' @param sample_col
#'
#' @return
#' @export
#'
#' @examples
get_clone_barcodes_exp <- function(cell_clone_bcode_dt,
                                   cell_bcode_col = "CellBarcode",
                                   clone_bcode_col = "CloneBarcode",
                                   sample_col = NA) {

    if (is.na(sample_col)) {
        grouping_col <- c(cell_bcode_col, clone_bcode_col)
    } else {
        grouping_col <- c(cell_bcode_col, clone_bcode_col, sample_col)
    }

    clone_exp <- cell_clone_bcode_dt[, .(n_reads = .N), by=grouping_col]

    return(clone_exp)
}

#' Title
#'
#' @param cell_clone_bcode_exp
#' @param valid_cells_bcodes
#' @param cell_bcode_col
#'
#' @return
#' @export
#'
#' @examples
#'
#' @import treemapify
draw_treemap <- function(cell_clone_bcode_exp, valid_cells_bcodes,
                         cell_bcode_col = "CellBarcode") {


    n_cells_with_clone_count <- count_clones_per_cells_for_treemap(
        valid_cells_bcodes = valid_cells_bcodes,
        cell_clone_bcode_exp = cell_clone_bcode_exp,
        cell_bcode_col = cell_bcode_col)

    plt <- ggplot(n_cells_with_clone_count,
                  aes(fill = n_clone_barcode, area = n_cells, label = perc_and_n_cells)) +
        geom_treemap() +
        geom_treemap_text(colour = "black",
                          place = "centre") +
        labs(
            title = "Distribution of Cells by No. of Detected Clone Barcodes",
            subtitle = paste("Of the",
                             prettyNum(length(valid_cells_bcodes), big.mark = ","),
                             "cells, how many have 0, 1, 2, etc. clone barcodes"),
            fill = "No. of clone barcodes found"
        ) +
        theme(legend.position = "bottom") +
        scale_fill_brewer(palette = "Set3")

    return(plt)
}


#' Title
#'
#' @param valid_cells_bcodes
#' @param cell_clone_bcode_exp
#' @param cell_bcode_col
#'
#' @return
#' @export
#'
#' @examples
count_clones_per_cells_for_treemap <- function(valid_cells_bcodes, cell_clone_bcode_exp,
                                       cell_bcode_col) {
    n_valid_cells <- length(valid_cells_bcodes)

    # get the clone barcodes for all cells that are deemed valid.
    # for most cases, these are barcodes for cells which cellranger
    # has deemed valid.


    clone_count_per_cell <- count_clones_per_cells(
        cell_clone_bcode_exp = cell_clone_bcode_exp,
        cell_bcode_col = cell_bcode_col,
        valid_cells_bcodes = valid_cells_bcodes
    )

    # Convert the n_clone_barcode that are > 10 to actually > 10
    clone_count_per_cell_inc_zeros <- copy(clone_count_per_cell)
    clone_count_per_cell_inc_zeros[, n_clone_barcode := ifelse(n_clone_barcode > 10, "> 10", as.character(n_clone_barcode))]

    # we will add in those valid cells which we can't find any clone barcodes for, and
    # assign n_clone_barcode as 0
    valid_cells_without_clones <- setdiff(
        x = valid_cells_bcodes,
        y = unique(cell_clone_bcode_exp[[cell_bcode_col]])
    )
    valid_cells_without_clones_dt <- data.table(
        CellBarcode = valid_cells_without_clones,
        n_clone_barcode = rep(0, length(valid_cells_without_clones))
    )
    setnames(valid_cells_without_clones_dt, "CellBarcode", cell_bcode_col)

    # combine everything
    clone_count_per_cell_inc_zeros <- rbind(
        clone_count_per_cell_inc_zeros,
        valid_cells_without_clones_dt
    )

    # set the level for the clone count per cell
    clone_count_per_cell_level <- c(as.character(seq(0, 10)), "> 10")
    clone_count_per_cell_inc_zeros[, n_clone_barcode := factor(n_clone_barcode, levels = clone_count_per_cell_level)]

    # have to actually count the number of cells with 0, 1, 2, .. clone barcodes
    # in order to get percentage that we draw the treemap for
    n_cells_with_clone_count <- clone_count_per_cell_inc_zeros[, .(n_cells = .N), by=n_clone_barcode]
    n_cells_with_clone_count <- n_cells_with_clone_count[order(n_clone_barcode)]

    n_cells_with_clone_count[, perc_and_n_cells := paste0(
        prettyNum(n_cells, big.mark = ","), " (",
        round(100 * n_cells / n_valid_cells, 1),
        "%)")]
}

count_clones_per_cells <- function(cell_clone_bcode_exp, cell_bcode_col, valid_cells_bcodes) {
    valid_cells_with_clones <- intersect(
        x = unique(cell_clone_bcode_exp[[cell_bcode_col]]),
        y = valid_cells_bcodes
    )

    valid_cell_clone_barcodes <- cell_clone_bcode_exp[get(cell_bcode_col) %in% valid_cells_with_clones]

    # count how many clones found per cell.
    clone_count_per_cell <- valid_cell_clone_barcodes[, .(n_clone_barcode = .N), by=cell_bcode_col]

    return(clone_count_per_cell)
}

#' Title
#'
#' @param cell_by_gene_mat
#' @param cell_clone_bcode_exp
#' @param cell_clone_reads_dt
#' @param cell_bcode_col
#' @param barcode_edit_dist_col
#' @param clone_bcode_col
#' @param export_as_hdf5
#'
#' @return
#' @export
#'
#' @examples
#'
#' @import HDF5Array
#'
embed_clone_ids_and_export <- function(cell_by_gene_mat,
                                       cell_clone_bcode_exp,
                                       cell_clone_reads_dt,
                                       export_as_hdf5 = TRUE,
                                       cell_bcode_col = "CellBarcode",
                                       barcode_edit_dist_col = "BarcodeEditDist",
                                       clone_bcode_col = "CloneBarcode") {

    valid_cells_bcodes <- colData(cell_by_gene_mat)$Barcode

    cell_clone_assignments <- get_cells_clones_assignment(
        cell_clone_bcode_exp = cell_clone_bcode_exp,
        valid_cells_bcodes = valid_cells_bcodes,
        cell_clone_reads_dt = cell_clone_reads_dt,
        cell_bcode_col = cell_bcode_col,
        barcode_edit_dist_col = barcode_edit_dist_col,
        clone_bcode_col = clone_bcode_col
    )

    # Attach clone bcode to SCE
    cell_barcodes <- as.data.table(colData(cell_by_gene_mat))
    # set the order of clone barcodes so they match the SCE object
    cell_clone_assignments <- cell_clone_assignments[order(match(get(cell_bcode_col), colData(cell_by_gene_mat)$Barcode))]

    colData(cell_by_gene_mat)$clone_barcode <- cell_clone_assignments[[clone_bcode_col]]

    if (export_as_hdf5) {
        saveHDF5SummarizedExperiment(
            x = cell_by_gene_mat,
            dir = "CloneDetective_output"
        )
    }
}

get_cells_clones_assignment <- function(valid_cells_bcodes,
                                        cell_clone_reads_dt,
                                        cell_bcode_col,
                                        barcode_edit_dist_col,
                                        clone_bcode_col,
                                        most_dominant_threshold = 0.5) {

    # convert the cell barcode reads data.table to cell barcode expression table.
    cell_clone_bcode_exp <- get_clone_barcodes_exp(
        cell_clone_bcode_dt = cell_clone_reads_dt,
        cell_bcode_col = cell_bcode_col,
        clone_bcode_col = clone_bcode_col
    )

    # first count how many clones were detected per cell
    clone_count_per_cell <- count_clones_per_cells(
        cell_clone_bcode_exp = cell_clone_bcode_exp,
        cell_bcode_col = cell_bcode_col,
        valid_cells_bcodes = valid_cells_bcodes
    )

    # isolate the easiest case where we found only 1 clone per cell
    single_clone_cells <- clone_count_per_cell[n_clone_barcode == 1,][[cell_bcode_col]]
    single_clone_cell_clone_ids <- cell_clone_bcode_exp[get(cell_bcode_col) %in% single_clone_cells]

    # get multiclone cells and rank their reads
    multiclone_cell_read_counts <- find_and_rank_clones_in_multiclone_cells(
        clone_count_per_cell = clone_count_per_cell,
        cell_bcode_col = cell_bcode_col,
        cell_clone_bcode_exp = cell_clone_bcode_exp
    )

    # find 2nd case, most abundant clone made up more than 50% of reads
    cells_assigned_to_majority_clone <- multiclone_cell_read_counts[ranking == 1 & prop_reads > most_dominant_threshold]

    # remove those cells which most abundant clones are dominant from cells_assigned_to_majority_clone,
    # then proceed to deal with the 3rd case, where we need to look at average edit distance
    # to determine the clone ID
    cells_with_determined_clones <- cells_assigned_to_majority_clone[[cell_bcode_col]]
    multiclone_cell_read_counts <- multiclone_cell_read_counts[! get(cell_bcode_col) %in% cells_with_determined_clones]
    # save memory
    rm(cells_with_determined_clones)
    cells_to_find_clones_for <- unique(multiclone_cell_read_counts[[cell_bcode_col]])

    clones_identified_via_avg_edit_dist <- assign_clone_by_edit_distance(
        cells_to_find = cells_to_find_clones_for,
        cell_bcode_col = cell_bcode_col,
        cell_clone_reads_dt = cell_clone_reads_dt,
        barcode_edit_dist_col = barcode_edit_dist_col,
        clone_bcode_col = clone_bcode_col
    )

    # Add in the cells which we can't find the clones for.
    cells_without_clones <- setdiff(
        x = valid_cells_bcodes,
        y = unique(cell_clone_bcode_exp[[cell_bcode_col]])
    )
    cells_without_clones <- data.table(cellbcode = cells_without_clones)
    setnames(cells_without_clones, "cellbcode", cell_bcode_col)

    # R is messy with NA. If a column is a character it shows <NA>.
    # If it is numeric, it is NA. Idiots..
    clone_barcodes <- rbindlist(
        list(
            single_clone_cell_clone_ids,
            cells_assigned_to_majority_clone,
            clones_identified_via_avg_edit_dist,
            cells_without_clones
        ), fill = TRUE
    )

    # No need to return the ranking
    clone_barcodes[, ranking := NULL]

    return(clone_barcodes)
}

#' Title
#'
#' @param clone_count_per_cell
#' @param cell_bcode_col
#' @param cell_clone_bcode_exp
#'
#' @return
#' @export
#'
#' @examples
find_and_rank_clones_in_multiclone_cells <- function(clone_count_per_cell, cell_bcode_col, cell_clone_bcode_exp) {
    multiclone_cells <- clone_count_per_cell[n_clone_barcode > 1,][[cell_bcode_col]]

    # The number of reads for each clone, for cells with multiple clones
    multiclone_cell_read_counts <- cell_clone_bcode_exp[get(cell_bcode_col) %in% multiclone_cells]

    # Order the table so for each cell barcode, we have clone with the most reads at the top
    multiclone_cell_read_counts <- multiclone_cell_read_counts[order(get(cell_bcode_col), -n_reads)]

    # Calculate proportion of reads for each clone
    multiclone_cell_read_counts[, prop_reads := n_reads / sum(n_reads), by=get(cell_bcode_col)]

    # Add ranking to the clones based on the number of reads
    multiclone_cell_read_counts[, ranking := frank(-n_reads, ties.method = "random"), by=c(cell_bcode_col)]

    return(multiclone_cell_read_counts)
}

#' Title
#'
#' @param cells_to_find
#' @param cell_clone_reads_dt
#' @param cell_bcode_col
#' @param barcode_edit_dist_col
#' @param clone_bcode_col
#'
#' @return
#' @export
#'
#' @examples
assign_clone_by_edit_distance <- function(cells_to_find,
                                          cell_clone_reads_dt,
                                          cell_bcode_col,
                                          barcode_edit_dist_col,
                                          clone_bcode_col) {

    edit_distances <- cell_clone_reads_dt[get(cell_bcode_col) %in% cells_to_find]

    # calculate average edit distance for each clone barcode and cell
    # can just convert the mean using mean function, but still need to calculate n reads, so no point.
    # just as complicated really.
    edit_distances_avg <- edit_distances[, .(sum_bcode_edit_dist = sum(get(barcode_edit_dist_col)), n_reads = .N), by = c(cell_bcode_col, clone_bcode_col)]
    edit_distances_avg[, mean_bcode_edit_dist := sum_bcode_edit_dist / n_reads]
    edit_distances_avg[, prop_reads := n_reads / sum(n_reads), by=c(cell_bcode_col)]

    # Rank by average barcode edit distance
    edit_distances_avg[, ranking := frank(list(mean_bcode_edit_dist, -prop_reads), ties.method = "random"), by=c(cell_bcode_col)]

    # order by cell barcode and then ranking
    # edit_distances_avg <- edit_distances_avg[order(get(cell_bcode_col), ranking)]

    res <-  edit_distances_avg[ranking == 1, ]
    return(res)
}

