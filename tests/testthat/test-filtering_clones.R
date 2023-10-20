library(data.table)

test_that("simple counting works correctly", {
    toy_clone_counts <- data.table(
        sample_name = rep("test", 5),
        read_count = c(1, 5, 10, 15, 20)
    )
    thresholds <- c(2, 10)

    expected_res <- data.table(
        at_least_2_cells = c(4L),
        at_least_10_cells = c(3L)
    )

    res <- count_retained_clones(
        count_data = toy_clone_counts,
        thresholds = thresholds,
        count_column = "read_count"
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)

})

test_that("simple counting no results", {
    toy_clone_counts <- data.table(
        sample_name = rep("test", 5),
        read_count = c(1, 5, 10, 15, 20)
    )
    thresholds <- c(30, 50)

    expected_res <- data.table(
        at_least_30_cells = c(0L),
        at_least_50_cells = c(0L)
    )

    res <- count_retained_clones(
        count_data = toy_clone_counts,
        thresholds = thresholds,
        count_column = "read_count"
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)

})

test_that("counting with group by works correctly", {
    toy_clone_counts <- data.table(
        sample_name = c(rep("test1", 2), rep("test2", 3)),
        read_count = c(1, 5, 10, 15, 20)
    )
    thresholds <- c(2, 10)

    expected_res <- data.table(
        sample_name = c("test1", "test2"),
        at_least_2_cells = c(1L, 3L),
        at_least_10_cells = c(0L, 3L)
    )

    res <- count_retained_clones(
        count_data = toy_clone_counts,
        thresholds = thresholds,
        count_column = "read_count",
        grouping_col = "sample_name"
    )
    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)

})

test_that("counting with group by a sample do not pass any thresholds", {
    toy_clone_counts <- data.table(
        sample_name = c(rep("test1", 2), rep("test2", 3)),
        read_count = c(1, 5, 10, 15, 20)
    )
    thresholds <- c(6, 10)

    expected_res <- data.table(
        sample_name = c("test1", "test2"),
        at_least_6_cells = c(0L, 3L),
        at_least_10_cells = c(0L, 3L)
    )

    res <- count_retained_clones(
        count_data = toy_clone_counts,
        thresholds = thresholds,
        count_column = "read_count",
        grouping_col = "sample_name"
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)

})

test_that("counting with group by no sample not pass any thresholds", {
    toy_clone_counts <- data.table(
        sample_name = c(rep("test1", 2), rep("test2", 3)),
        read_count = c(1, 5, 10, 15, 20)
    )
    thresholds <- c(100, 150)

    expected_res <- data.table(
        sample_name = c("test1", "test2"),
        at_least_100_cells = c(0L, 0L),
        at_least_150_cells = c(0L, 0L)
    )

    res <- count_retained_clones(
        count_data = toy_clone_counts,
        thresholds = thresholds,
        count_column = "read_count",
        grouping_col = "sample_name"
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)

})

test_that("filtering clone barcodes works correctly", {
    toy_clone_counts <- data.table(
        sample_name = c(rep("test1", 2), rep("test2", 3)),
        read_count = c(1, 5, 10, 15, 20),
        clone_barcodes = c("ACGT", "CATG", "CATG", "ACGT", "ATGC")
    )

    expected_res <- data.table(
        sample_name = c(rep("test1", 1), rep("test2", 3)),
        read_count = c(5, 10, 15, 20),
        clone_barcodes = c("CATG", "CATG", "ACGT", "ATGC")
    )

    res <- remove_clones_below_threshold(
        count_data = toy_clone_counts,
        threshold = 2,
        count_column = "read_count"
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)

})


test_that("get top clone barcodes works correctly", {
    toy_clone_counts <- data.table(
        sample_name = c(rep("test1", 3), rep("test2", 3)),
        read_count = c(1, 5, 20, 10, 15, 12),
        clone_barcodes = c("ACGT", "CATG", "CATG", "ACGT", "ATGC", "TCGT")
    )

    expected_res <- data.table(
        sample_name = c(rep("test1", 2), rep("test2", 2)),
        read_count = c(20, 5, 15, 12),
        clone_barcodes = c("CATG", "CATG", "ATGC", "TCGT")
    )

    res <- get_top_barcodes(
        count_data = toy_clone_counts,
        top_threshold = 2,
        count_column = "read_count",
        grouping_col = "sample_name"
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)

})

test_that("get top clone barcodes no group by works correctly", {
    toy_clone_counts <- data.table(
        sample_name = c(rep("test1", 3), rep("test2", 3)),
        read_count = c(1, 5, 20, 10, 15, 12),
        clone_barcodes = c("ACGT", "CATG", "CATG", "ACGT", "ATGC", "TCGT")
    )

    expected_res <- data.table(
        sample_name = c("test1", "test2","test2"),
        read_count = c(20, 15, 12),
        clone_barcodes = c("CATG", "ATGC", "TCGT")
    )

    res <- get_top_barcodes(
        count_data = toy_clone_counts,
        top_threshold = 3,
        count_column = "read_count"
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)

})



