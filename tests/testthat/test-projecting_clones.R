test_that("clone projection works", {
    toy_clone_counts <- data.table(
        sample_name = c(rep("test1", 3), rep("test2", 3)),
        read_count = c(1, 5, 20, 10, 15, 12),
        clone_barcodes = c("ACGT", "CATG", "CATG", "ACGT", "ATGC", "TCGT")
    )

    expected_res <- data.table(
        sample_name = c(rep("test1", 3), rep("test2", 3)),
        read_count = c(1, 5, 20, 10, 15, 12),
        clone_barcodes = c("ACGT", "CATG", "CATG", "ACGT", "ATGC", "TCGT"),
        read_count_proportion = c(1/26, 5/26, 20/26, 10/37, 15/37, 12/37)
    )
    expected_res[, projected_to_100 := 100 * read_count_proportion]
    expected_res[, projected_to_200 := 200 * read_count_proportion]

    res <- project_clones(
        count_data = toy_clone_counts,
        project_amnt = c(100, 200),
        count_column = "read_count",
        grouping_col = "sample_name"
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)
})

test_that("clone projection no grouping works", {

    read_counts <- c(1, 5, 20, 10, 15, 12)
    sum_read_counts <- sum(read_counts)

    toy_clone_counts <- data.table(
        sample_name = c(rep("test1", 3), rep("test2", 3)),
        read_count = read_counts,
        clone_barcodes = c("ACGT", "CATG", "CATG", "ACGT", "ATGC", "TCGT")
    )

    expected_res <- copy(toy_clone_counts)
    expected_res[, read_count_proportion := read_count / sum_read_counts]
    expected_res[, projected_to_100 := 100 * read_count_proportion]

    res <- project_clones(
        count_data = toy_clone_counts,
        project_amnt = c(100),
        count_column = "read_count"
    )

    expect_equal(names(res), names(expected_res))
    expect_identical(res, expected_res, ignore_attr = TRUE)
})
