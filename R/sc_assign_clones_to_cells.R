#' Assign Clones to Cells and Embed the Assignment to SingleCellExperiment Object
#'
#' This function assigns clone barcodes to cells and optionally embeds the
#' assignment into a SingleCellExperiment object. It calculates a cell-by-clone
#' matrix from a provided data.table where each row corresponds to a read
#' associated with a cell, UMI, and clone barcode.
#' See details for the assignment method.
#'
#' @param cell_by_gene_mat A SingleCellExperiment object containing the
#' cell-by-gene matrix. It must have a 'Barcode' column in 'colData' that
#' uniquely identifies each cell. The barcodes should match those in the
#' `cell_clone_reads_dt` data table.
#'
#' @param cell_clone_reads_dt A data.table object representing the reads. Each
#' row includes information about a cell, UMI, and clone barcode. For data
#' produced by NextClone, use `fread` from the data.table package and pass the
#' resulting object to this parameter.
#'
#' @param cell_bcode_col Name of the column in `cell_clone_reads_dt` that
#' indicates the cell barcode for each read.
#'
#' @param barcode_edit_dist_col Name of the column in `cell_clone_reads_dt` that
#' indicates the edit distance for the clone barcode of each read.
#'
#' @param clone_bcode_col Name of the column in `cell_clone_reads_dt` that
#' specifies the clone barcode for each read.
#'
#' @param umi_col Name of the column in `cell_clone_reads_dt` that specifies the
#' UMI barcode for each read.
#'
#' @param most_dominant_threshold A numeric value between 0 and 1 specifying the
#' proportion threshold of reads that must be associated with a clone barcode to
#' assign it to a cell in cases where multiple clone barcodes are detected.
#'
#' @param umi_clone_consensus_threshold A numeric value between 0 and 1 specifying
#' the proportion threshold of reads for collapsing UMIs when computing the
#' cell-by-clone matrix.
#' See details for more information.
#'
#' @param embed_to_mat A boolean indicating whether to embed the clone barcode
#' assignment directly into the `cell_by_gene_mat` object.
#'
#'
#' @return Depending on the `embed_to_mat` parameter, this function returns either
#' an updated SingleCellExperiment object with clone barcode assignments embedded
#' or a data.table with the cell-clone assignments.
#'
#' @details
#' Clone barcode assignment to cells follows a tiered approach:
#' Cells with a single detected clone barcode are straightforwardly assigned to that clone.
#' In cases where multiple clone barcodes are present, the most dominant clone barcode
#' (constituting over 50 percent of reads, adjustable via `most_dominant_threshold`) is assigned.
#' Remaining cells, where no clone barcode is sufficiently dominant, are assigned based
#' on the lowest average barcode edit distance. If edit distances are equal, the clone
#' barcode with the higher read count prevails.
#' To use the default 50 percent threshold, set `most_dominant_threshold` to 0.5.
#'
#' The cell-by-clone matrix construction first collapses reads with the same cell
#' and UMI barcodes
#' For a group of reads have the same cell barcode and UMI barcode,
#' if the reads are mapped to several clone barcodes,
#' by default, they are collapsed into one read and assigned to the clone barcode
#' comprising 70 percent or more of its group's reads.
#' This threshold modifiable by the `umi_clone_consensus_threshold` parameter.
#' To apply the default threshold of 70 percent, set this parameter to 0.7.
#'
#'
#' @examples
#' library(scater)
#' library(data.table)
#'
#' set.seed(42)
#'
#' sce <- mockSCE(ncells = 4, ngenes = 100)
#' colData(sce)$Barcode <- colnames(sce)
#'
#' cell_clone_reads_dt <- data.table(
#'     CellBarcode = c(
#'         rep("Cell_001", 6),
#'         rep("Cell_002", 7),
#'         rep("Cell_003", 4),
#'         rep("Cell_004", 5)
#'     ),
#'     CloneBarcode = c(
#'         rep("AA", 6),
#'         rep("BB1", 4), rep("BB2", 3),
#'         rep("CC1", 2), rep("CC2", 2),
#'         rep("CC1", 2), "DD1", "DD2", "DD3"
#'     ),
#'     BarcodeEditDist = c(
#'         rep(0L, 6),
#'         rep(0L, 7),
#'         c(0, 1, 2, 1),
#'         rep(0, 5)
#'
#'     )
#' )
#' cell_clone_reads_dt[, UMI := seq_len(nrow(cell_clone_reads_dt))]
#'
#' sce_with_clone <- assign_and_embed_clones(
#'     cell_by_gene_mat = sce,
#'     cell_clone_reads_dt = cell_clone_reads_dt,
#'     cell_bcode_col = "CellBarcode",
#'     barcode_edit_dist_col = "BarcodeEditDist",
#'     clone_bcode_col = "CloneBarcode",
#'     umi_col = "UMI"
#' )
#' colData(sce_with_clone)
#'
#' @import scater
#'
#' @export
assign_and_embed_clones <- function(cell_by_gene_mat,
                                    cell_clone_reads_dt,
                                    cell_bcode_col = "CellBarcode",
                                    barcode_edit_dist_col = "BarcodeEditDist",
                                    clone_bcode_col = "CloneBarcode",
                                    umi_col = "UMI",
                                    most_dominant_threshold = 0.5,
                                    umi_clone_consensus_threshold = 0.7,
                                    embed_to_mat = TRUE) {

    valid_cells_bcodes <- colData(cell_by_gene_mat)$Barcode

    cell_clone_assignments <- get_cells_clones_assignment(
        valid_cells_bcodes = valid_cells_bcodes,
        cell_clone_reads_dt = cell_clone_reads_dt,
        cell_bcode_col = cell_bcode_col,
        barcode_edit_dist_col = barcode_edit_dist_col,
        umi_col = umi_col,
        clone_bcode_col = clone_bcode_col
    )

    # Attach clone bcode to SCE
    if (embed_to_mat) {
        cell_barcodes <- as.data.table(colData(cell_by_gene_mat))
        # set the order of clone barcodes so they match the SCE object
        cell_clone_assignments <- cell_clone_assignments[order(match(get(cell_bcode_col), colData(cell_by_gene_mat)$Barcode))]

        colData(cell_by_gene_mat)$clone_barcode <- cell_clone_assignments[[clone_bcode_col]]
        colData(cell_by_gene_mat)$clone_barcode_criteria <- cell_clone_assignments$criteria

        return(cell_by_gene_mat)
    } else {
        cell_clone_assignments <- cell_clone_assignments[, c(cell_bcode_col, clone_bcode_col, "criteria"), with=FALSE]
        return(cell_clone_assignments)
    }


}

#' Get clone assignments.
#'
#' Internal function which apply the stratified protocol described in
#' `assign_and_embed_clones` function.
#'
#' @param valid_cells_bcodes A vector of cells barcodes to get the clone
#' assignments for.
#'
#' @param cell_clone_reads_dt A data.table object representing the reads. Each
#' row includes information about a cell, UMI, and clone barcode. For data
#' produced by NextClone, use `fread` from the data.table package and pass the
#' resulting object to this parameter.
#'
#' @param cell_bcode_col Name of the column in `cell_clone_reads_dt` that
#' indicates the cell barcode for each read.
#'
#' @param barcode_edit_dist_col Name of the column in `cell_clone_reads_dt` that
#' indicates the edit distance for the clone barcode of each read.
#'
#' @param umi_col Name of the column in `cell_clone_reads_dt` that specifies the
#' UMI barcode for each read.
#'
#' @param most_dominant_threshold A numeric value between 0 and 1 specifying the
#' proportion threshold of reads that must be associated with a clone barcode to
#' assign it to a cell in cases where multiple clone barcodes are detected.
#'
#' @param umi_clone_consensus_threshold A numeric value between 0 and 1 specifying
#' the proportion threshold of reads for collapsing UMIs when computing the
#' cell-by-clone matrix.
#' See details for more information.
#'
#' @param clone_bcode_col Name of the column in `cell_clone_reads_dt` that
#' specifies the clone barcode for each read.
#'
#' @return A data.table of cell and clone assignment.
#'
#' @keywords internal
#'
get_cells_clones_assignment <- function(valid_cells_bcodes,
                                        cell_clone_reads_dt,
                                        cell_bcode_col = "CellBarcode",
                                        barcode_edit_dist_col = "BarcodeEditDist",
                                        clone_bcode_col = "CloneBarcode",
                                        umi_col = "UMI",
                                        most_dominant_threshold = 0.5,
                                        umi_clone_consensus_threshold = 0.7) {

    # convert the cell barcode reads data.table to cell barcode expression table.
    cell_by_clone_matrix <- generate_cell_clone_barcode_matrix(
        cell_clone_bcode_dt = cell_clone_reads_dt,
        cell_bcode_col = cell_bcode_col,
        clone_bcode_col = clone_bcode_col,
        umi_col = umi_col,
        umi_clone_consensus_threshold = umi_clone_consensus_threshold
    )

    # first count how many clones were detected per cell
    clone_count_per_cell <- count_clones_per_cells(
        cell_by_clone_matrix = cell_by_clone_matrix,
        cell_bcode_col = cell_bcode_col,
        valid_cells_bcodes = valid_cells_bcodes
    )

    # isolate the easiest case where we found only 1 clone per cell
    single_clone_cells <- clone_count_per_cell[n_clone_barcode == 1,][[cell_bcode_col]]
    single_clone_cell_clone_ids <- cell_by_clone_matrix[get(cell_bcode_col) %in% single_clone_cells]
    single_clone_cell_clone_ids[, criteria := "single_clone"]

    # get multiclone cells and rank their reads
    multiclone_cell_read_counts <- find_and_rank_clones_in_multiclone_cells(
        clone_count_per_cell = clone_count_per_cell,
        cell_bcode_col = cell_bcode_col,
        cell_by_clone_matrix = cell_by_clone_matrix
    )

    # find 2nd case, most abundant clone made up more than 50% of reads
    cells_assigned_to_majority_clone <- multiclone_cell_read_counts[ranking == 1 & prop_reads > most_dominant_threshold]
    cells_assigned_to_majority_clone[, criteria := paste0("dominant_clone_moreThan_", gsub("\\.", "_", as.character(most_dominant_threshold))) ]

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
    clones_identified_via_avg_edit_dist[, criteria := "clone_from_edit_distance"]

    # Add in the cells which we can't find the clones for.
    cells_without_clones <- setdiff(
        x = valid_cells_bcodes,
        y = unique(cell_by_clone_matrix[[cell_bcode_col]])
    )
    cells_without_clones <- data.table(cellbcode = cells_without_clones)
    setnames(cells_without_clones, "cellbcode", cell_bcode_col)
    cells_without_clones[, criteria := "no_clones_found"]

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

    # Turn the criteria into factor
    clone_barcodes[, criteria := as.factor(criteria)]

    return(clone_barcodes)
}

#' Count number of clones per cells.
#'
#' Internal function to compute the number of clones detected per cell.
#'
#' @param valid_cells_bcodes A vector of cells barcodes to get the clone
#' assignments for.
#'
#' @param cell_by_clone_matrix A cell-by-clone matrix generated by
#' `generate_cell_clone_barcode_matrix` function.
#'
#' @param cell_bcode_col Name of the column in `cell_by_clone_matrix` that
#' indicates the cell barcode.
#'
#' @return A data.table denoting the number of clones detected per cell.
#'
#' @keywords internal
#'
count_clones_per_cells <- function(cell_by_clone_matrix, cell_bcode_col, valid_cells_bcodes) {
    valid_cells_with_clones <- intersect(
        x = unique(cell_by_clone_matrix[[cell_bcode_col]]),
        y = valid_cells_bcodes
    )

    valid_cell_clone_barcodes <- cell_by_clone_matrix[get(cell_bcode_col) %in% valid_cells_with_clones]

    # count how many clones found per cell.
    clone_count_per_cell <- valid_cell_clone_barcodes[, .(n_clone_barcode = .N), by=cell_bcode_col]

    return(clone_count_per_cell)
}

#' Count number of clones per cells.
#'
#' Internal function to find cells with multiple clones detected,
#' and rank the clones based on the proportion of reads detected for them.
#'
#' @param clone_count_per_cell A data.table denoting the number of clones detected
#' per cell. Generated by `count_clones_per_cells` function.
#'
#' @param cell_by_clone_matrix A cell-by-clone matrix generated by
#' `generate_cell_clone_barcode_matrix` function.
#'
#' @param cell_bcode_col Name of the column in `cell_by_clone_matrix` that
#' indicates the cell barcode.
#'
#' @return A data.table denoting the ranking of clones for multiclone cells.
#'
#' @keywords internal
#'
find_and_rank_clones_in_multiclone_cells <- function(clone_count_per_cell, cell_bcode_col, cell_by_clone_matrix) {
    multiclone_cells <- clone_count_per_cell[n_clone_barcode > 1,][[cell_bcode_col]]

    # The number of reads for each clone, for cells with multiple clones
    multiclone_cell_read_counts <- cell_by_clone_matrix[get(cell_bcode_col) %in% multiclone_cells]

    # Order the table so for each cell barcode, we have clone with the most reads at the top
    multiclone_cell_read_counts <- multiclone_cell_read_counts[order(get(cell_bcode_col), -n_reads)]

    # Calculate proportion of reads for each clone
    multiclone_cell_read_counts[, prop_reads := n_reads / sum(n_reads), by=get(cell_bcode_col)]

    # Add ranking to the clones based on the number of reads
    multiclone_cell_read_counts[, ranking := frank(-n_reads, ties.method = "random"), by=c(cell_bcode_col)]

    return(multiclone_cell_read_counts)
}

#' Assign clones by edit distance.
#'
#' Internal function to assign clones to multiclone cells based on the edit
#' distance of the clone barcode reads
#'
#' @param cells_to_find A vector of multiclone cells to do the assignment.
#'
#' @param cell_clone_reads_dt A data.table object representing the reads. Each
#' row includes information about a cell, UMI, and clone barcode. For data
#' produced by NextClone, use `fread` from the data.table package and pass the
#' resulting object to this parameter.
#'
#' @param cell_bcode_col Name of the column in `cell_clone_reads_dt` that
#' indicates the cell barcode for each read.
#'
#' @param barcode_edit_dist_col Name of the column in `cell_clone_reads_dt` that
#' indicates the edit distance for the clone barcode of each read.
#'
#' @param clone_bcode_col Name of the column in `cell_clone_reads_dt` that
#' specifies the clone barcode for each read.
#'
#' @return A data.table denoting the clones for multiclone cells.
#'
#' @keywords internal
#'
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

