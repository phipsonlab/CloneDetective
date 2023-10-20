test_that("one plot one elbow NO plot_grouping_col numerical barcode id assignment works", {

    # 1 plot, 1 "elbow"
    # Essentially treating everything in the data as 1 plot

    toy_clone_counts <- data.table(
        sample_name = rep("sample1", 6),
        read_count = c(1, 5, 20, 10, 15, 12),
        clone_barcodes = c("ACGT", "CATG", "CATG", "ACGT", "ATGC", "TCGT")
    )

    expected_res <- data.table(
        sample_name = rep("sample1", 6),
        read_count = c(20, 15, 12, 10, 5, 1),
        clone_barcodes = c("CATG", "ATGC", "TCGT", "ACGT", "CATG", "ACGT"),
        barcode_id = seq(0, 5)
    )

    res <- assign_numeric_barcode_id(
        count_data = toy_clone_counts,
        count_column = "read_count"
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)
})

test_that("one plot one elbow WITH plot_grouping_col numerical barcode id assignment works", {

    # 1 plot, 1 "elbow"
    # Essentially treating everything in the data as 1 plot

    toy_clone_counts <- data.table(
        sample_name = rep("sample1", 6),
        read_count = c(1, 5, 20, 10, 15, 12),
        clone_barcodes = c("ACGT", "CATG", "CATG", "ACGT", "ATGC", "TCGT")
    )

    expected_res <- data.table(
        sample_name = rep("sample1", 6),
        read_count = c(20, 15, 12, 10, 5, 1),
        clone_barcodes = c("CATG", "ATGC", "TCGT", "ACGT", "CATG", "ACGT"),
        barcode_id = seq(0, 5)
    )

    res <- assign_numeric_barcode_id(
        count_data = toy_clone_counts,
        count_column = "read_count",
        elbow_grouping_col = NA,
        plot_grouping_col = "sample_name"
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)
})

test_that("multiple plots one elbow numerical barcode id assignment works", {

    # 2 plots, 1 "elbow" each

    toy_clone_counts <- data.table(
        sample_name = c(rep("sample1", 3),rep("sample2", 3)),
        read_count = c(1, 5, 20, 10, 15, 12)
    )

    expected_res <- data.table(
        sample_name = c(rep("sample1", 3),rep("sample2", 3)),
        read_count = c(20, 5, 1, 15, 12, 10),
        barcode_id = rep(seq(0, 2), 2)
    )

    res <- assign_numeric_barcode_id(
        count_data = toy_clone_counts,
        count_column = "read_count",
        elbow_grouping_col = NA,
        plot_grouping_col = "sample_name"
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)
})

test_that("multiple plots, multiple elbows numerical barcode id assignment works", {

    # multiple plots, multiple "elbow"s each.
    # This is pretending as if in one ggplot plot (perhaps a treatment)
    # there are multiple elbows because for each treatment, we have multiple replicates
    # so 1 elbow per replicate

    toy_clone_counts <- data.table(
        treatment = c(rep("treat1", 4), rep("treat2", 4)),
        replicate = rep(c(rep("rep1", 2), rep("rep2", 2)), 2),
        read_count = c(20, 10, 15, 2, 3, 30, 15, 20)
    )

    expected_res <- data.table(
        treatment = c(rep("treat1", 4), rep("treat2", 4)),
        replicate = rep(c(rep("rep1", 2), rep("rep2", 2)), 2),
        read_count = c(20, 10, 15, 2, 30, 3, 20, 15),
        barcode_id = rep(seq(0, 3), 2)
    )

    res <- assign_numeric_barcode_id(
        count_data = toy_clone_counts,
        count_column = "read_count",
        elbow_grouping_col = "replicate",
        plot_grouping_col = "treatment"
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)
})

