#' Generate cell by clone matrix.
#'
#' Count how many reads are mapped to each clone/cell pair, and generate a table
#' that resemles cell-by-clone matrix.
#'
#' @param cell_clone_reads_dt A data.table object representing the reads. Each
#' row includes information about a cell, UMI, and clone barcode. For data
#' produced by NextClone, use `fread` from the data.table package and pass the
#' resulting object to this parameter.
#'
#' @param cell_bcode_col Name of the column in `cell_clone_reads_dt` that
#' indicates the cell barcode for each read.
#'
#' @param clone_bcode_col Name of the column in `cell_clone_reads_dt` that
#' specifies the clone barcode for each read.
#'
#' @param umi_col Name of the column in `cell_clone_reads_dt` that specifies the
#' UMI barcode for each read.
#'
#' @param umi_clone_consensus_threshold A numeric value between 0 and 1 specifying
#' the proportion threshold of reads for collapsing UMIs when computing the
#' cell-by-clone matrix.
#' See details for more information.
#'
#' @return A data.table resembling cell-by-clone matrix.
#'
#'
#' @details
#' The cell-by-clone matrix construction first collapses reads with the same cell
#' and UMI barcodes
#' For a group of reads have the same cell barcode and UMI barcode,
#' if the reads are mapped to several clone barcodes,
#' by default, they are collapsed into one read and assigned to the clone barcode
#' comprising 70% or more of its group's reads.
#' This threshold modifiable by the `umi_clone_consensus_threshold` parameter.
#' To apply the default threshold of 70%, set this parameter to 0.7.
#'
#' @examples
#' library(data.table)
#'
#' cell_clone_bcode_dt <- data.table(
#'        CellBarcode = c(
#'            rep("A", 2),
#'            rep("B", 4),
#'            rep("C", 5),
#'            rep("D", 7)
#'        ),
#'        CloneBarcode = c(
#'            rep("AA", 2),
#'            rep("AA", 3), "BB",
#'            rep("CC", 4), "DD",
#'            rep("XX", 2), rep("YY", 3), rep("ZZ", 2)
#'        ),
#'        UMI = c(
#'            rep("A1", 2),
#'            rep("B1", 2), "B2", "B3",
#'            rep("C1", 3), "C2", "C1",
#'            rep("D1", 5), rep("D2", 2)
#'        )
#'    )
#' res <- generate_cell_clone_barcode_matrix(cell_clone_bcode_dt)
#'
#' @importFrom purrr pmap
#'
#' @export
generate_cell_clone_barcode_matrix <- function(cell_clone_bcode_dt,
                                               cell_bcode_col = "CellBarcode",
                                               clone_bcode_col = "CloneBarcode",
                                               umi_col = "UMI",
                                               umi_clone_consensus_threshold = 0.7) {


    n_reads_per_cell_clone_umi <- cell_clone_bcode_dt[, .(n_reads = .N), by=c(cell_bcode_col, clone_bcode_col, umi_col)]

    # to deal with problematic cells, that is cells that have
    # multiple reads with the same UMI and but they all don't mapped to the same clone.
    n_clones_per_cell_umi <- cell_clone_bcode_dt[, .(n_clones_found = uniqueN(CloneBarcode)), by = c(cell_bcode_col, umi_col)]
    problematic_cells <- n_clones_per_cell_umi[n_clones_found > 1]

    # so we can see how many reads are mapped to which clone for all the reads
    # that have the same UMI but for the same cell.
    problematic_cells_info <- merge.data.table(
        x = problematic_cells,
        y = n_reads_per_cell_clone_umi
    )


    # let's retain only clones which have clear consensus.
    # i.e. 6 reads, 5 reads map to A, 1 read maps to B, clearly A is the winner
    # as 80% of the reads are mapped to A.
    # Set 80% as parameter.
    problematic_cells_info[, prop_reads := n_reads / sum(n_reads), by = c(cell_bcode_col, umi_col)]
    problematic_cells_to_retain <- problematic_cells_info[prop_reads > umi_clone_consensus_threshold,]
    problematic_cells_to_retain[, prop_reads := NULL]

    # easier ones where reads with the same umi all map to the same clone.
    non_problematic_cells <- n_clones_per_cell_umi[n_clones_found == 1]
    non_problematic_cells_info <- merge.data.table(
        x = non_problematic_cells,
        y = n_reads_per_cell_clone_umi
    )

    n_reads_per_cell_clone_umi_clean <- rbind(
        unique(non_problematic_cells_info, by = c(cell_bcode_col, umi_col)),
        problematic_cells_to_retain
    )

    clones_exp <- n_reads_per_cell_clone_umi_clean[, .(n_reads = .N), by=c(cell_bcode_col, clone_bcode_col)]

    return(clones_exp)
}
