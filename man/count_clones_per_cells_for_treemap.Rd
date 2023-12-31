% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_draw_treemap.R
\name{count_clones_per_cells_for_treemap}
\alias{count_clones_per_cells_for_treemap}
\title{Count number of clones per cells.}
\usage{
count_clones_per_cells_for_treemap(
  valid_cells_bcodes,
  cell_by_clone_matrix,
  cell_bcode_col
)
}
\arguments{
\item{valid_cells_bcodes}{A vector of cells barcodes to get the clone
assignments for.}

\item{cell_by_clone_matrix}{A cell-by-clone matrix generated by
`generate_cell_clone_barcode_matrix` function.}

\item{cell_bcode_col}{Name of the column in `cell_by_clone_matrix` that
indicates the cell barcode.}
}
\value{
A data.table denoting the number of clones detected per cell.
}
\description{
Internal function to compute the number of clones detected per cell
for drawing treemap.
}
\keyword{internal}
