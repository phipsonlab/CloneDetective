% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_clones.R
\name{remove_clones_below_threshold}
\alias{remove_clones_below_threshold}
\title{Remove clones with counts below a specified threshold}
\usage{
remove_clones_below_threshold(count_data, threshold, count_column)
}
\arguments{
\item{count_data}{A data.table containing clone barcodes and their counts.
One row per clone barcode.}

\item{threshold}{A numeric value indicating the minimum count required for a
clone barcode to be retained.}

\item{count_column}{A character string specifying the column name for
clone barcode counts.}
}
\value{
A data.table with clones above the threshold.
}
\description{
Removes clones with counts strictly less than the specified threshold.
}
\examples{
library(data.table)

toy_clone_counts <- data.table(
sample_name = c(rep("test1", 2), rep("test2", 3)),
read_count = c(1, 5, 10, 15, 20),
clone_barcodes = c("ACGT", "CATG", "CATG", "ACGT", "ATGC")
)

remove_clones_below_threshold(
  count_data = toy_clone_counts,
  threshold = 2,
  count_column = "read_count"
)

}
