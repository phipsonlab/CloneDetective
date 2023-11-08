test_that("get clone barcode expression works", {
    cell_clone_bcode_dt <- data.table(
        CellBarcode = c(
            rep("A", 2),
            rep("B", 4),
            rep("C", 5),
            rep("D", 7)
        ),
        CloneBarcode = c(
            rep("AA", 2),
            rep("AA", 3), "BB",
            rep("CC", 4), "DD",
            rep("XX", 2), rep("YY", 3), rep("ZZ", 2)
        ),
        UMI = c(
            rep("A1", 2),
            rep("B1", 2), "B2", "B3",
            rep("C1", 3), "C2", "C1",
            rep("D1", 5), rep("D2", 2)
        )
    )

    expected_res <- data.table(
        CellBarcode = c("A", "B", "B", "C", "D"),
        CloneBarcode = c("AA", "AA", "BB", "CC", "ZZ"),
        n_reads = as.integer(c(1,2,1,2,1))
    )

    res <- get_clone_barcodes_exp(cell_clone_bcode_dt)

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)

})

test_that("count clone barcodes per cell works", {

    cell_clone_bcode_exp <- data.table(
        CellBarcode = c("A", "B", "C", "C"),
        CloneBarcode = c("AA", "AA", "BB", "CC"),
        n_reads = as.integer(c(2, 2, 2, 1))
    )

    res <- count_clones_per_cells(
        valid_cells_bcodes = c("B", "C", "D"),
        cell_clone_bcode_exp = cell_clone_bcode_exp,
        cell_bcode_col = "CellBarcode"
    )

    expected_res <- data.table(
        CellBarcode = c("B", "C"),
        n_clone_barcode = c(1L,2L)
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)

})

test_that("count clone barcodes per cell for treemap works", {

    cell_clone_bcode_exp <- data.table(
        CellBarcode = c("A", "B", "C", "C"),
        CloneBarcode = c("AA", "AA", "BB", "CC"),
        n_reads = as.integer(c(2, 2, 2, 1))
    )

    res <- count_clones_per_cells_for_treemap(
        valid_cells_bcodes = c("B", "C", "D"),
        cell_clone_bcode_exp = cell_clone_bcode_exp,
        cell_bcode_col = "CellBarcode"
    )

    expected_res <- data.table(
        n_clone_barcode = factor(c(0, 1, 2)),
        n_cells = as.integer(rep(1, 3)),
        perc_and_n_cells = as.character(rep("1 (33.3%)", 3))
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)

})

test_that("find_and_rank_clones_in_multiclone_cells works", {

    cell_clone_bcode_exp <- data.table(
        CellBarcode = c("A", "B", "B", "C", "C", "C"),
        CloneBarcode = c("AA", "AA", "BB", "CC", "DD", "EE"),
        n_reads = as.integer(c(2, 2, 1, 1, 2, 3))
    )

    clone_cnts_per_cell <- count_clones_per_cells(
        valid_cells_bcodes = c("A", "B", "C"),
        cell_clone_bcode_exp = cell_clone_bcode_exp,
        cell_bcode_col = "CellBarcode"
    )

    res <- find_and_rank_clones_in_multiclone_cells(
        clone_count_per_cell = clone_cnts_per_cell,
        cell_bcode_col = "CellBarcode",
        cell_clone_bcode_exp = cell_clone_bcode_exp
    )

    expected_res <- data.table(
        CellBarcode = c("B", "B", "C", "C", "C"),
        CloneBarcode = c("AA", "BB", "EE", "DD", "CC"),
        n_reads = as.integer(c(2, 1, 3, 2, 1)),
        prop_reads = c(2/3, 1/3, 3/6, 2/6, 1/6),
        ranking = as.integer(c(1, 2, 1, 2, 3))
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)

})

test_that("determine_clone_by_edit_distance works", {
    # A has 3 clone barcodes, AA, BB, CC, all same read proportion (2/6), but AA has the smallest avg edit dist.
    # B has 2 clone barcodes, DD, EE, same avg edit dist of 0, but EE has more proportion of reads.
    cell_clone_reads_dt <- data.table(
        CellBarcode = c(
            rep("A", 6),
            rep("B", 5)
        ),
        CloneBarcode = c(
            "AA", "AA", "BB", "BB", "CC", "CC",
            "DD", "DD", "EE", "EE", "EE"
        ),
        BarcodeEditDist = c(
            0, 0, 1, 0, 1, 1,
            rep(0, 5)
        )
    )

    res <- assign_clone_by_edit_distance(
        cells_to_find = c("A", "B"),
        cell_clone_reads_dt = cell_clone_reads_dt,
        cell_bcode_col = "CellBarcode",
        barcode_edit_dist_col = "BarcodeEditDist",
        clone_bcode_col = "CloneBarcode"
    )

    expected_res <- data.table(
        CellBarcode = c("A", "B"),
        CloneBarcode = c("AA", "EE"),
        sum_bcode_edit_dist = rep(0, 2),
        n_reads = c(2L, 3L),
        mean_bcode_edit_dist = rep(0, 2),
        prop_reads = c(2/6, 3/5),
        ranking = c(1L, 1L)
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)
})

test_that("get_cells_clones_assignment works", {
    # A mix of everything, 1 clone, most abundant is dominant, edit distance,
    # tie edit distance then proportion of reads, and no clone.
    cell_clone_reads_dt <- data.table(
        CellBarcode = c(
            rep("A", 6),
            rep("B", 7),
            rep("C", 4),
            rep("D", 5)
        ),
        CloneBarcode = c(
            rep("AA", 6),
            rep("BB1", 4), rep("BB2", 3),
            rep("CC1", 2), rep("CC2", 2),
            rep("CC1", 2), "DD1", "DD2", "DD3"
        ),
        BarcodeEditDist = c(
            rep(0L, 6),
            rep(0L, 7),
            c(0, 1, 2, 1),
            rep(0, 5)

        )
    )
    cell_clone_reads_dt[, UMI := seq_len(nrow(cell_clone_reads_dt))]

    res <- get_cells_clones_assignment(
        valid_cells_bcodes = c("A", "B", "C", "D", "E"),
        cell_clone_reads_dt = cell_clone_reads_dt,
        cell_bcode_col = "CellBarcode",
        barcode_edit_dist_col = "BarcodeEditDist",
        clone_bcode_col = "CloneBarcode"
    )

    expected_res <- data.table(
        CellBarcode = c("A", "B", "C", "D", "E"),
        CloneBarcode = c("AA", "BB1", "CC1", "CC1", NA),
        n_reads = c(6L, 4L, 2L, 2L, NA),
        prop_reads = c(NA, 4/7, 2/4, 2/5, NA),
        sum_bcode_edit_dist = c(NA, NA, 1, 0, NA),
        mean_bcode_edit_dist = c(NA, NA, 1/2, 0, NA)
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)
})

test_that("get_cells_clones_assignment with more than 0.5 threshold works", {
    # A mix of everything, 1 clone, most abundant is dominant, edit distance,
    # tie edit distance then proportion of reads, and no clone.
    # Same as above basically, but now we will have a case where the definition of
    # dominant is altered, and depending on the threshold, cell E may or may not
    # be assigned clone based on edit distance.
    cell_clone_reads_dt <- data.table(
        CellBarcode = c(
            rep("A", 6),
            rep("B", 7),
            rep("C", 4),
            rep("D", 5),
            rep("E", 7)

        ),
        CloneBarcode = c(
            rep("AA", 6),
            rep("BB1", 6), "BB2",
            rep("CC1", 2), rep("CC2", 2),
            rep("CC1", 2), "DD1", "DD2", "DD3",
            rep("EE1", 4), rep("EE2", 3)
        ),
        BarcodeEditDist = c(
            rep(0L, 6),
            rep(0L, 7),
            c(0, 1, 2, 1),
            rep(0, 5),
            rep(1, 4), rep(0, 3)
        )
    )
    cell_clone_reads_dt[, UMI := seq_len(nrow(cell_clone_reads_dt))]

    res <- get_cells_clones_assignment(
        valid_cells_bcodes = c("A", "B", "C", "D", "E", "F"),
        cell_clone_reads_dt = cell_clone_reads_dt,
        cell_bcode_col = "CellBarcode",
        barcode_edit_dist_col = "BarcodeEditDist",
        clone_bcode_col = "CloneBarcode",
        most_dominant_threshold = 0.7
    )

    expected_res <- data.table(
        CellBarcode = c("A", "B", "C", "D", "E", "F"),
        CloneBarcode = c("AA", "BB1", "CC1", "CC1", "EE2", NA),
        n_reads = c(6L, 6L, 2L, 2L, 3L, NA),
        prop_reads = c(NA, 6/7, 2/4, 2/5, 3/7, NA),
        sum_bcode_edit_dist = c(NA, NA, 1, 0, 0, NA),
        mean_bcode_edit_dist = c(NA, NA, 1/2, 0, 0, NA)
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)
})



