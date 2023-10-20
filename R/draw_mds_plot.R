#' Draw MDS plot using ggplot.
#'
#' The default plotMDS from edgeR uses R's `plot` function which is hard
#' to manipulate. This function will draw MDS plot as ggplot and arrange
#' it into a grid broken down by library sizes.
#'
#' @param dge_list A list of DGEList objects, 1 per library.
#' @param lib_sizes A character vector of library sizes.
#' 1 indexed element corresponds to 1 indexed element in `dge_list`.
#'
plot_mds_as_ggplot <- function(dge_list, lib_sizes) {
    mds_objs <- lapply(dge_list, function(y) {
        mds = plotMDS(y, plot=FALSE)
        return(mds)
    })

    # Convert the object to dataframe so we can use ggplot!
    mds_df <- pmap(
        list(lib_sizes, mds_objs, dge_list),
        function(lib_size, mds_obj, dge) {
            lib_label <- paste0(
                lib_size,
                " (dim1 var=", round(mds_obj$var.explained[1] * 100), "%", ",",
                " dim2 var=", round(mds_obj$var.explained[2] * 100), "%", ")"
            )
            df <- data.frame(
                x=mds_obj$x,
                y=mds_obj$y,
                treatment=dge$samples$group,
                replicate=str_split_i(rownames(dge$samples), "_", 2),
                library_size=lib_label
            )
            return(df)
        })
    mds_df <- bind_rows(mds_df)

    axis_lab <- mds_objs[[1]]$axislabel

    ggplot(mds_df, aes(x=x, y=y)) +
        geom_point(aes(colour=treatment, shape=replicate), size=3) +
        facet_wrap(~ library_size) +
        labs(
            x = paste(axis_lab, "1"),
            y = paste(axis_lab, "2"),
            color = "Treatment",
            shape = "Replicate",
            title = "MDS plot for each library for NGS data v2"
        ) +
        scale_y_continuous(breaks = pretty_breaks(n=10)) +
        scale_x_continuous(breaks = pretty_breaks(n=5)) +
        scale_color_viridis_d(option = 'turbo') +
        theme_bw() +
        theme(text=element_text(size=14))
}
