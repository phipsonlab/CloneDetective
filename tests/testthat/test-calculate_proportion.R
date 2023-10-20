library(data.table)

test_that("proportion is calculated correctly", {
    toy_clone_counts <- data.table(
        sample_name = c(rep("test1", 2), rep("test2", 3)),
        read_count = c(1, 5, 10, 15, 20)
    )

    expected_res <- data.table(
        sample_name = c(rep("test1", 2), rep("test2", 3)),
        read_count = c(1, 5, 10, 15, 20),
        read_proportion = c(1/6, 5/6, 10/45, 15/45, 20/45)
    )

    res <- convert_count_to_proportion(
        count_data = toy_clone_counts,
        grouping_col = "sample_name",
        count_column = "read_count"
    )

    expect_identical(res, expected_res, ignore_attr = TRUE)


})

test_that("proportion no grouping is calculated correctly", {

    read_counts <- c(1, 5, 10, 15, 20)
    sum_read_counts <- sum(read_counts)

    toy_clone_counts <- data.table(
        sample_name = c(rep("test1", 2), rep("test2", 3)),
        read_count = read_counts
    )

    expected_res <- copy(toy_clone_counts)
    expected_res[, read_proportion := read_count / sum_read_counts]

    res <- convert_count_to_proportion(
        count_data = toy_clone_counts,
        count_column = "read_count"
    )

    expect_identical(res, expected_res, ignore_attr = TRUE)


})
