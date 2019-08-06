monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

#' Plot a dataset and trajectory in 3 dimensions
#'
#' @param cds cell_data_set to plot
#' @param dims numeric vector that indicates the dimensions used to create the
#'   3D plot, by default it is the first three dimensions.
#' @param reduction_method string indicating the reduction method to plot.
#' @param color_cells_by the cell attribute (e.g. the column of colData(cds))
#'   to map to each cell's color. Default is cluster.
#' @param genes a gene name or gene id to color the plot by.
#' @param show_trajectory_graph a logical used to indicate whether to graph the
#'   principal graph backbone. Default is TRUE.
#' @param trajectory_graph_color the color of graph backbone. Default is black.
#' @param trajectory_graph_segment_size numeric indicating the width of the
#'   graph backbone. Default is 5.
#' @param norm_method string indicating the method used to transform gene
#'   expression when gene markers are provided. Default is "log". "size_only"
#'   is also supported.
#' @param cell_size numeric indicating the size of the point to be plotted.
#'   Default is 25.
#' @param alpha numeric indicating the alpha value of the plotted cells.
#'   Default is 1.
#' @param min_expr numeric indicating the minimum marker gene value to be
#'   colored. Default is 0.1.
#' @return a plotly plot object
#' @export
#' @examples
#' \dontrun{
#' plot_cells_3d(cds, markers=c("Rbfox3, Neurod1", "Sox2"))
#' }
#'
#' @export
plot_cells_3d <- function(cds,
                          dims = c(1,2,3),
                          reduction_method = c("UMAP", "tSNE", "PCA", "LSI"),
                          color_cells_by="cluster",
                          #group_cells_by=c("cluster", "partition"), #
                          genes=NULL,
                          show_trajectory_graph=TRUE,
                          trajectory_graph_color="black",
                          trajectory_graph_segment_size=5,
                          norm_method = c("log", "size_only"),
                          #label_cell_groups = TRUE,#
                          #label_groups_by_cluster=TRUE,#
                          #group_label_size=2,#
                          #labels_per_group=1,#
                          #label_branch_points=TRUE,#
                          #label_roots=TRUE,#
                          #label_leaves=TRUE,#
                          #graph_label_size=2,#
                          cell_size=25,
                          alpha = 1,
                          min_expr=0.1) {

  reduction_method <- match.arg(reduction_method)
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste("No dimensionality reduction for",
                                      reduction_method, "calculated.",
                                      "Please run reduce_dimensions with",
                                      "reduction_method =", reduction_method,
                                      "before attempting to plot."))
  low_dim_coords <- reducedDims(cds)[[reduction_method]]
  if(!is.null(color_cells_by)) {
    assertthat::assert_that(color_cells_by %in% c("cluster", "partition",
                                                  "pseudotime") |
                              color_cells_by %in% names(colData(cds)),
                            msg = paste("color_cells_by must be a column in",
                                        "the colData table."))
  }

  assertthat::assert_that(!is.null(color_cells_by) || !is.null(markers),
                          msg = paste("Either color_cells_by or markers must",
                                      "be NULL, cannot color by both!"))
  norm_method = match.arg(norm_method)

  if (show_trajectory_graph &&
      is.null(principal_graph(cds)[[reduction_method]])) {
    message("No trajectory to plot. Has learn_graph() been called yet?")
    show_trajectory_graph = FALSE
  }

  gene_short_name <- NA
  sample_name <- NA

  x <- dims[[1]]
  y <- dims[[2]]
  z <- dims[[3]]

  S_matrix <- reducedDims(cds)[[reduction_method]]
  data_df <- data.frame(S_matrix[,c(dims)])

  colnames(data_df) <- c("data_dim_1", "data_dim_2", "data_dim_3")
  data_df$sample_name <- row.names(data_df)

  data_df <- as.data.frame(cbind(data_df, colData(cds)))

  if (color_cells_by == "cluster"){
    data_df$cell_color <- tryCatch({
      clusters(cds, reduction_method = reduction_method)[data_df$sample_name]},
      error = function(e) {NULL})
  } else if (color_cells_by == "partition") {
    data_df$cell_color <- tryCatch({
      partitions(cds,
                 reduction_method = reduction_method)[data_df$sample_name]},
      error = function(e) {NULL})
  } else if (color_cells_by == "pseudotime") {
    data_df$cell_color <- tryCatch({
      pseudotime(cds,
                 reduction_method = reduction_method)[data_df$sample_name]},
      error = function(e) {NULL})
  } else{
    data_df$cell_color <- colData(cds)[data_df$sample_name,color_cells_by]
  }

  ## Marker genes
  markers_exprs <- NULL
  if (!is.null(genes)) {
    if ((is.null(dim(genes)) == FALSE) && dim(genes) >= 2){
      markers <- unlist(genes[,1], use.names=FALSE)
    } else {
      markers <- genes
    }
    markers_rowData <-
      as.data.frame(subset(rowData(cds), gene_short_name %in% markers |
                             row.names(rowData(cds)) %in% markers))
    if (nrow(markers_rowData) >= 1) {
      cds_exprs <- SingleCellExperiment::counts(cds)[row.names(markers_rowData), ,drop=FALSE]
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds))

      if ((is.null(dim(genes)) == FALSE) && dim(genes) >= 2){
        genes <- as.data.frame(genes)
        row.names(genes) <- genes[,1]
        genes <- genes[row.names(cds_exprs),]
        agg_mat <-
          as.matrix(my.aggregate.Matrix(cds_exprs,
                                        as.factor(genes[,2]),
                                        fun="sum"))
        agg_mat <- t(scale(t(log10(agg_mat + 1))))
        agg_mat[agg_mat < -2] <- -2
        agg_mat[agg_mat > 2] <- 2
        markers_exprs <- agg_mat
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')

        markers_exprs$feature_label <- markers_exprs$feature_id
        #markers_linear <- TRUE
      } else {
        cds_exprs@x <- round(cds_exprs@x)
        markers_exprs <- matrix(cds_exprs, nrow=nrow(markers_rowData))
        colnames(markers_exprs) <- colnames(SingleCellExperiment::counts(cds))
        row.names(markers_exprs) <- row.names(markers_rowData)
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
        markers_exprs <- merge(markers_exprs, markers_rowData,
                               by.x = "feature_id", by.y="row.names")
        markers_exprs$feature_label <-
          as.character(markers_exprs$gene_short_name)
        markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <-
          markers_exprs$feature_id
        markers_exprs$feature_label <- factor(markers_exprs$feature_label,
                                              levels = markers)
      }
    }
  }

  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    data_df <- merge(data_df, markers_exprs, by.x="sample_name",
                     by.y="cell_id")
    data_df$expression <- with(data_df, ifelse(value >= min_expr, value, NA))
    sub1 <- data_df[!is.na(data_df$expression),]
    sub2 <- data_df[is.na(data_df$expression),]
    if(norm_method == "size_only"){

      p <- plotly::plot_ly(sub1, x = ~data_dim_1, y = ~data_dim_2,
                           z = ~data_dim_3, type = 'scatter3d',
                           color = ~expression, size=I(cell_size),
                           mode="markers", name = "Expression",
                           alpha = I(alpha)) %>%
        plotly::add_markers(x = sub2$data_dim_1, y = sub2$data_dim_2,
                            color = I("lightgrey"), z = sub2$data_dim_3,
                            marker=list(opacity = .4), showlegend=FALSE)
    } else {
      sub1$log10_expression <- log10(sub1$expression + min_expr)
      p <- plotly::plot_ly(sub1, x = ~data_dim_1, y = ~data_dim_2,
                           z = ~data_dim_3, type = 'scatter3d',
                           color = ~log10_expression, mode="markers",
                           size=I(cell_size), name = "Log10\nExpression",
                           alpha = I(alpha)) %>%
        plotly::add_markers(x = sub2$data_dim_1, y = sub2$data_dim_2,
                            z = sub2$data_dim_3, color = I("lightgrey"),
                            marker=list(opacity = .4), showlegend=FALSE)
    }
  } else {
    if(color_cells_by %in% c("cluster", "partition")){
      if (is.null(data_df$cell_color)){
        p <- plotly::plot_ly(data_df, x = ~data_dim_1, y = ~data_dim_2,
                             z = ~data_dim_3, type = 'scatter3d',
                             size=I(cell_size), color=I("gray"),
                             mode="markers", alpha = I(alpha))
        message(paste("cluster_cells() has not been called yet, can't color",
                      "cells by cluster or partition"))
      } else{
        p <- plotly::plot_ly(data_df, x = ~data_dim_1, y = ~data_dim_2,
                             z = ~data_dim_3, type = 'scatter3d',
                             size=I(cell_size), color=~cell_color,
                             mode="markers", alpha = I(alpha))
      }
    } else if(class(data_df$cell_color) == "numeric") {

      p <- plotly::plot_ly(data_df) %>%
        plotly::add_trace(x = ~data_dim_1, y = ~data_dim_2, z = ~data_dim_3,
                          type = 'scatter3d', size=I(cell_size), alpha = I(alpha),
                          mode="markers", marker=list(
                            colorbar = list(title = color_cells_by, len=0.5),
                            color=~cell_color,
                            line=list(width = 1,
                                      color = ~cell_color,
                                      colorscale="Viridis"),
                            colorscale="Viridis"))
    } else {
      p <- plotly::plot_ly(data_df, x = ~data_dim_1, y = ~data_dim_2,
                           z = ~data_dim_3, type = 'scatter3d',
                           size=I(cell_size), color=~cell_color,
                           mode="markers", alpha = I(alpha))
    }
  }
  p <- p %>%
    plotly::layout(scene = list(xaxis=list(title=paste("Component", x)),
                                yaxis=list(title=paste("Component", y)),
                                zaxis=list(title=paste("Component", z))))
  ## Graph info
  if (show_trajectory_graph) {

    ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
      as.data.frame() %>%
      dplyr::select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y,
                     prin_graph_dim_3 = z) %>%
      dplyr::mutate(sample_name = rownames(.),
                    sample_state = rownames(.))

    dp_mst <- cds@principal_graph[[reduction_method]]

    edge_df <- dp_mst %>%
      igraph::as_data_frame() %>%
      dplyr::select_(source = "from", target = "to") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select_(source="sample_name",
                                        source_prin_graph_dim_1="prin_graph_dim_1",
                                        source_prin_graph_dim_2="prin_graph_dim_2",
                                        source_prin_graph_dim_3="prin_graph_dim_3"),
                       by = "source") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select_(target="sample_name",
                                        target_prin_graph_dim_1="prin_graph_dim_1",
                                        target_prin_graph_dim_2="prin_graph_dim_2",
                                        target_prin_graph_dim_3="prin_graph_dim_3"),
                       by = "target")

    for (i in 1:nrow(edge_df)) {
      p <- p %>%
        plotly::add_trace(
          x = as.vector(t(edge_df[i, c("source_prin_graph_dim_1",
                                       "target_prin_graph_dim_1")])),
          y = as.vector(t(edge_df[i, c("source_prin_graph_dim_2",
                                       "target_prin_graph_dim_2")])),
          z = as.vector(t(edge_df[i, c("source_prin_graph_dim_3",
                                       "target_prin_graph_dim_3")])),
          color = trajectory_graph_color,
          line = list(color = I(trajectory_graph_color),
                      width = trajectory_graph_segment_size), mode = 'lines',
          type = 'scatter3d', showlegend = FALSE)
    }
  }
  p
}

#' Plots the cells along with their trajectories.
#'
#' @param cds cell_data_set for the experiment
#' @param x the column of reducedDims(cds) to plot on the horizontal axis
#' @param y the column of reducedDims(cds) to plot on the vertical axis
#' @param cell_size The size of the point for each cell
#' @param cell_stroke The stroke used for plotting each cell - default is 1/2
#'   of the cell_size
#' @param reduction_method The lower dimensional space in which to plot cells.
#'   Must be one of "UMAP", "tSNE", "PCA" and "LSI".
#' @param color_cells_by What to use for coloring the cells. Must be either the
#'   name of a column of colData(cds), or one of "clusters", "partitions", or
#'   "pseudotime".
#' @param group_cells_by How to group cells when labeling them. Must be either
#'   the name of a column of colData(cds), or one of "clusters" or "partitions".
#'   If a column in colData(cds), must be a categorical variable.
#' @param genes Facet the plot, showing the expression of each gene in a facet
#'   panel. Must be either a list of gene ids (or short names), or a dataframe
#'   with two columns that groups the genes into modules that will be
#'   aggregated prior to plotting. If the latter, the first column must be gene
#'   ids, and the second must the group for each gene.
#' @param show_trajectory_graph Whether to render the principal graph for the
#'   trajectory. Requires that learn_graph() has been called on cds.
#' @param trajectory_graph_color The color to be used for plotting the
#'   trajectory graph.
#' @param trajectory_graph_segment_size The size of the line segments used for
#'   plotting the trajectory graph.
#' @param norm_method How to normalize gene expression scores prior to plotting
#'   them. Must be one of "log" or "size_only".
#' @param label_cell_groups Whether to label cells in each group (as specified
#'   by group_cells_by) according to the most frequently occurring label(s) (as
#'   specified by color_cells_by) in the group. If false, plot_cells() simply
#'   adds a traditional color legend.
#' @param label_groups_by_cluster Instead of labeling each cluster of cells,
#'   place each label once, at the centroid of all cells carrying that label.
#' @param group_label_size Font size to be used for cell group labels.
#' @param labels_per_group How many labels to plot for each group of cells.
#'   Defaults to 1, which plots only the most frequent label per group.
#' @param label_branch_points Whether to plot a label for each branch point in
#'   the principal graph.
#' @param label_roots Whether to plot a label for each root in the principal
#'   graph.
#' @param label_leaves Whether to plot a label for each leaf node in the
#'   principal graph.
#' @param graph_label_size How large to make the branch, root, and leaf labels.
#' @param alpha Alpha for the cells. Useful for reducing overplotting.
#' @param min_expr Minimum expression threshold for plotting genes
#' @param rasterize Whether to plot cells as a rastered bitmap. Requires the
#'   ggrastr package.
#'
#' @return a ggplot2 plot object
#' @export
#' @examples
#' \dontrun{
#' lung <- load_A549()
#' plot_cells(lung)
#' plot_cells(lung, color_cells_by="log_dose")
#' plot_cells(lung, markers="GDF15")
#' }
plot_cells <- function(cds,
                       x=1,
                       y=2,
                       reduction_method = c("UMAP", "tSNE", "PCA", "LSI"),
                       color_cells_by="cluster",
                       group_cells_by=c("cluster", "partition"),
                       genes=NULL,
                       show_trajectory_graph=TRUE,
                       trajectory_graph_color="grey28",
                       trajectory_graph_segment_size=0.75,
                       norm_method = c("log", "size_only"),
                       label_cell_groups = TRUE,
                       label_groups_by_cluster=TRUE,
                       group_label_size=2,
                       labels_per_group=1,
                       label_branch_points=TRUE,
                       label_roots=TRUE,
                       label_leaves=TRUE,
                       graph_label_size=2,
                       cell_size=0.35,
                       cell_stroke= I(cell_size / 2),
                       alpha = 1,
                       min_expr=0.1,
                       rasterize=FALSE) {
  reduction_method <- match.arg(reduction_method)
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste("No dimensionality reduction for",
                                      reduction_method, "calculated.",
                                      "Please run reduce_dimensions with",
                                      "reduction_method =", reduction_method,
                                      "before attempting to plot."))
  low_dim_coords <- reducedDims(cds)[[reduction_method]]
  assertthat::assert_that(ncol(low_dim_coords) >=max(x,y),
                          msg = paste("x and/or y is too large. x and y must",
                                      "be dimensions in reduced dimension",
                                      "space."))
  if(!is.null(color_cells_by)) {
    assertthat::assert_that(color_cells_by %in% c("cluster", "partition",
                                                  "pseudotime") |
                              color_cells_by %in% names(colData(cds)),
                            msg = paste("color_cells_by must one of",
                                        "'cluster', 'partition', 'pseudotime,",
                                        "or a column in the colData table."))

    if(color_cells_by == "pseudotime") {
      tryCatch({pseudotime(cds, reduction_method = reduction_method)},
               error = function(x) {
                 stop(paste("No pseudotime for", reduction_method,
                            "calculated. Please run order_cells with",
                            "reduction_method =", reduction_method,
                            "before attempting to color by pseudotime."))})

    }
  }
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(markers),
                          msg = paste("Either color_cells_by or markers must",
                                      "be NULL, cannot color by both!"))

  norm_method = match.arg(norm_method)
  group_cells_by=match.arg(group_cells_by)
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(genes),
                          msg = paste("Either color_cells_by or genes must be",
                                      "NULL, cannot color by both!"))

  if (show_trajectory_graph &&
      is.null(principal_graph(cds)[[reduction_method]])) {
    message("No trajectory to plot. Has learn_graph() been called yet?")
    show_trajectory_graph = FALSE
  }

  gene_short_name <- NA
  sample_name <- NA
  #sample_state <- colData(cds)$State
  data_dim_1 <- NA
  data_dim_2 <- NA
  if (rasterize){
    plotting_func <- ggrastr::geom_point_rast
  }else{
    plotting_func <- ggplot2::geom_point
  }

  S_matrix <- reducedDims(cds)[[reduction_method]]
  data_df <- data.frame(S_matrix[,c(x,y)])

  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)

  data_df <- as.data.frame(cbind(data_df, colData(cds)))
  if (group_cells_by == "cluster"){
    data_df$cell_group <-
      tryCatch({clusters(cds,
                         reduction_method = reduction_method)[
                           data_df$sample_name]},
               error = function(e) {NULL})
  } else if (group_cells_by == "partition") {
    data_df$cell_group <-
      tryCatch({partitions(cds,
                           reduction_method = reduction_method)[
                             data_df$sample_name]},
               error = function(e) {NULL})
  } else{
    stop("Error: unrecognized way of grouping cells.")
  }

  if (color_cells_by == "cluster"){
    data_df$cell_color <-
      tryCatch({clusters(cds,
                         reduction_method = reduction_method)[
                           data_df$sample_name]},
               error = function(e) {NULL})
  } else if (color_cells_by == "partition") {
    data_df$cell_color <-
      tryCatch({partitions(cds,
                           reduction_method = reduction_method)[
                             data_df$sample_name]},
               error = function(e) {NULL})
  } else if (color_cells_by == "pseudotime") {
    data_df$cell_color <-
      tryCatch({pseudotime(cds,
                           reduction_method = reduction_method)[
                             data_df$sample_name]}, error = function(e) {NULL})
  } else{
    data_df$cell_color <- colData(cds)[data_df$sample_name,color_cells_by]
  }

  ## Graph info
  if (show_trajectory_graph) {

    ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
      as.data.frame() %>%
      dplyr::select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
      dplyr::mutate(sample_name = rownames(.),
                    sample_state = rownames(.))

    dp_mst <- cds@principal_graph[[reduction_method]]

    edge_df <- dp_mst %>%
      igraph::as_data_frame() %>%
      dplyr::select_(source = "from", target = "to") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select_(
                           source="sample_name",
                           source_prin_graph_dim_1="prin_graph_dim_1",
                           source_prin_graph_dim_2="prin_graph_dim_2"),
                       by = "source") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select_(
                           target="sample_name",
                           target_prin_graph_dim_1="prin_graph_dim_1",
                           target_prin_graph_dim_2="prin_graph_dim_2"),
                       by = "target")
  }

  ## Marker genes
  markers_exprs <- NULL
  expression_legend_label <- NULL
  if (!is.null(genes)) {
    if (!is.null(dim(genes)) && dim(genes) >= 2){
      markers = unlist(genes[,1], use.names=FALSE)
    } else {
      markers = genes
    }
    markers_rowData <- as.data.frame(subset(rowData(cds),
                                            gene_short_name %in% markers |
                                              row.names(rowData(cds)) %in%
                                              markers))
    if (nrow(markers_rowData) == 0) {
      stop("None of the provided genes were found in the cds")
    }
    if (nrow(markers_rowData) >= 1) {
      cds_exprs <- SingleCellExperiment::counts(cds)[row.names(markers_rowData), ,drop=FALSE]
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds))

      if (!is.null(dim(genes)) && dim(genes) >= 2){
        genes = as.data.frame(genes)
        row.names(genes) = genes[,1]
        genes = genes[row.names(cds_exprs),]

        agg_mat = as.matrix(my.aggregate.Matrix(cds_exprs, as.factor(genes[,2]),
                                                fun="sum"))

        agg_mat = t(scale(t(log10(agg_mat + 1))))
        agg_mat[agg_mat < -2] = -2
        agg_mat[agg_mat > 2] = 2
        markers_exprs = agg_mat
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
        if (is.factor(genes[,2]))
          markers_exprs$feature_id = factor(markers_exprs$feature_id,
                                            levels=levels(genes[,2]))

        markers_exprs$feature_label <- markers_exprs$feature_id
        norm_method = "size_only"
        expression_legend_label = "Expression score"
      } else {
        cds_exprs@x = round(cds_exprs@x)
        markers_exprs = matrix(cds_exprs, nrow=nrow(markers_rowData))
        colnames(markers_exprs) = colnames(SingleCellExperiment::counts(cds))
        row.names(markers_exprs) = row.names(markers_rowData)
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
        markers_exprs <- merge(markers_exprs, markers_rowData,
                               by.x = "feature_id", by.y="row.names")
        markers_exprs$feature_label <-
          as.character(markers_exprs$gene_short_name)

        markers_exprs$feature_label <- ifelse(is.na(markers_exprs$feature_label) | !as.character(markers_exprs$feature_label) %in% markers,
                                              as.character(markers_exprs$feature_id),
                                              as.character(markers_exprs$feature_label))

        markers_exprs$feature_label <- factor(markers_exprs$feature_label,
                                              levels = markers)
        if (norm_method == "size_only")
          expression_legend_label = "Expression"
        else
          expression_legend_label = "log10(Expression)"
      }
    }
  }

  if (label_cell_groups && is.null(color_cells_by) == FALSE){
    if (is.null(data_df$cell_color)){
      if (is.null(genes)){
        message(paste(color_cells_by, "not found in colData(cds), cells will",
                      "not be colored"))
      }
      text_df = NULL
      label_cell_groups = FALSE
    }else{
      if(is.character(data_df$cell_color) || is.factor(data_df$cell_color)) {

        if (label_groups_by_cluster && is.null(data_df$cell_group) == FALSE){
          text_df = data_df %>%
            dplyr::group_by(cell_group) %>%
            dplyr::mutate(cells_in_cluster= dplyr::n()) %>%
            dplyr::group_by(cell_color, add=TRUE) %>%
            dplyr::mutate(per=dplyr::n()/cells_in_cluster)
          median_coord_df = text_df %>%
            dplyr::summarize(fraction_of_group = dplyr::n(),
                             text_x = stats::median(x = data_dim_1),
                             text_y = stats::median(x = data_dim_2))
          text_df = suppressMessages(text_df %>% dplyr::select(per) %>%
                                       dplyr::distinct())
          text_df = suppressMessages(dplyr::inner_join(text_df,
                                                       median_coord_df))
          text_df = text_df %>% dplyr::group_by(cell_group) %>%
            dplyr::top_n(labels_per_group, per)
        } else {
          text_df = data_df %>% dplyr::group_by(cell_color) %>%
            dplyr::mutate(per=1)
          median_coord_df = text_df %>%
            dplyr::summarize(fraction_of_group = dplyr::n(),
                             text_x = stats::median(x = data_dim_1),
                             text_y = stats::median(x = data_dim_2))
          text_df = suppressMessages(text_df %>% dplyr::select(per) %>%
                                       dplyr::distinct())
          text_df = suppressMessages(dplyr::inner_join(text_df,
                                                       median_coord_df))
          text_df = text_df %>% dplyr::group_by(cell_color) %>%
            dplyr::top_n(labels_per_group, per)
        }

        text_df$label = as.character(text_df %>% dplyr::pull(cell_color))
        # I feel like there's probably a good reason for the bit below, but I
        # hate it and I'm killing it for now.
        # text_df$label <- paste0(1:nrow(text_df))
        # text_df$process_label <- paste0(1:nrow(text_df), '_',
        # as.character(as.matrix(text_df[, 1])))
        # process_label <- text_df$process_label
        # names(process_label) <- as.character(as.matrix(text_df[, 1]))
        # data_df[, group_by] <-
        #  process_label[as.character(data_df[, group_by])]
        # text_df$label = process_label
      } else {
        message(paste("Cells aren't colored in a way that allows them to",
                      "be grouped."))
        text_df = NULL
        label_cell_groups = FALSE
      }
    }
  }

  if (!is.null(markers_exprs) && nrow(markers_exprs) > 0){
    data_df <- merge(data_df, markers_exprs, by.x="sample_name",
                     by.y="cell_id")
    data_df$value <- with(data_df, ifelse(value >= min_expr, value, NA))
    na_sub <- data_df[is.na(data_df$value),]
    if(norm_method == "size_only"){
      g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) +
        plotting_func(aes(data_dim_1, data_dim_2), size=I(cell_size),
                      stroke = I(cell_stroke), color = "grey80",
                      data = na_sub) +
        plotting_func(aes(color=value), size=I(cell_size),
                      stroke = I(cell_stroke), na.rm = TRUE) +
        viridis::scale_color_viridis(option = "viridis",
                                     name = expression_legend_label,
                                     na.value = "grey80", end = 0.8) +
        guides(alpha = FALSE) + facet_wrap(~feature_label)
    } else {
      g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) +
        plotting_func(aes(data_dim_1, data_dim_2), size=I(cell_size),
                      stroke = I(cell_stroke), color = "grey80",
                      data = na_sub) +
        plotting_func(aes(color=log10(value+min_expr)),
                      size=I(cell_size), stroke = I(cell_stroke),
                      na.rm = TRUE) +
        viridis::scale_color_viridis(option = "viridis",
                                     name = expression_legend_label,
                                     na.value = "grey80", end = 0.8) +
        guides(alpha = FALSE) + facet_wrap(~feature_label)
    }
  } else {
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2))

    # We don't want to force users to call order_cells before even being able
    # to look at the trajectory, so check whether it's null and if so, just
    # don't color the cells
    if(color_cells_by %in% c("cluster", "partition")){
      if (is.null(data_df$cell_color)){
        g <- g + geom_point(color=I("gray"), size=I(cell_size),
                            stroke = I(cell_stroke), na.rm = TRUE,
                            alpha = I(alpha))
        message(paste("cluster_cells() has not been called yet, can't",
                      "color cells by cluster"))
      } else{
        g <- g + geom_point(aes(color = cell_color), size=I(cell_size),
                            stroke = I(cell_stroke), na.rm = TRUE,
                            alpha = alpha)
      }
      g <- g + guides(color = guide_legend(title = color_cells_by,
                                           override.aes = list(size = 4)))
    } else if (class(data_df$cell_color) == "numeric"){
      g <- g + geom_point(aes(color = cell_color), size=I(cell_size),
                          stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
      g <- g + viridis::scale_color_viridis(name = color_cells_by, option="C")
    } else {
      g <- g + geom_point(aes(color = cell_color), size=I(cell_size),
                          stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
      g <- g + guides(color = guide_legend(title = color_cells_by,
                                           override.aes = list(size = 4)))
    }

  }
  if (show_trajectory_graph){
    g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1",
                                     y="source_prin_graph_dim_2",
                                     xend="target_prin_graph_dim_1",
                                     yend="target_prin_graph_dim_2"),
                          size=trajectory_graph_segment_size,
                          color=I(trajectory_graph_color),
                          linetype="solid",
                          na.rm=TRUE,
                          data=edge_df)


    if (label_branch_points){
      mst_branch_nodes <- branch_nodes(cds)
      branch_point_df <- ica_space_df %>%
        dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
        dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))

      g <- g +
        geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                   shape = 21, stroke=I(trajectory_graph_segment_size),
                   color="white",
                   fill="black",
                   size=I(graph_label_size * 1.5),
                   na.rm=TRUE, branch_point_df) +
        geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                             label="branch_point_idx"),
                  size=I(graph_label_size), color="white", na.rm=TRUE,
                  branch_point_df)
    }

    if (label_leaves){
      mst_leaf_nodes <- leaf_nodes(cds)
      leaf_df <- ica_space_df %>%
        dplyr::slice(match(names(mst_leaf_nodes), sample_name)) %>%
        dplyr::mutate(leaf_idx = seq_len(dplyr::n()))

      g <- g +
        geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                   shape = 21, stroke=I(trajectory_graph_segment_size),
                   color="black",
                   fill="lightgray",
                   size=I(graph_label_size * 1.5),
                   na.rm=TRUE,
                   leaf_df) +
        geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                             label="leaf_idx"),
                  size=I(graph_label_size), color="black", na.rm=TRUE, leaf_df)
    }

    if (label_roots){
      mst_root_nodes <- root_nodes(cds)
      root_df <- ica_space_df %>%
        dplyr::slice(match(names(mst_root_nodes), sample_name)) %>%
        dplyr::mutate(root_idx = seq_len(dplyr::n()))

      g <- g +
        geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                   shape = 21, stroke=I(trajectory_graph_segment_size),
                   color="black",
                   fill="white",
                   size=I(graph_label_size * 1.5),
                   na.rm=TRUE,
                   root_df) +
        geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                             label="root_idx"),
                  size=I(graph_label_size), color="black", na.rm=TRUE, root_df)
    }
  }

  if(label_cell_groups) {
    g <- g + ggrepel::geom_text_repel(data = text_df,
                                      mapping = aes_string(x = "text_x",
                                                           y = "text_y",
                                                           label = "label"),
                                      size=I(group_label_size))
    # If we're coloring by gene expression, don't hide the legend
    if (is.null(markers_exprs))
      g <- g + theme(legend.position="none")
  }

  g <- g +
    #scale_color_brewer(palette="Set1") +
    monocle_theme_opts() +
    xlab(paste(reduction_method, x)) +
    ylab(paste(reduction_method, y)) +
    #guides(color = guide_legend(label.position = "top")) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white'))
  g
}



#' Plots expression for one or more genes as a function of pseudotime
#'
#' @param cds_subset subset cell_data_set including only the genes to be
#'   plotted.
#' @param min_expr the minimum (untransformed) expression level to plot.
#' @param cell_size the size (in points) of each cell used in the plot.
#' @param nrow the number of rows used when laying out the panels for each
#'   gene's expression.
#' @param ncol the number of columns used when laying out the panels for each
#'   gene's expression
#' @param panel_order vector of gene names indicating the order in which genes
#'   should be laid out (left-to-right, top-to-bottom). If
#'   \code{label_by_short_name = TRUE}, use gene_short_name values, otherwise
#'   use feature IDs.
#' @param color_cells_by the cell attribute (e.g. the column of colData(cds))
#'   to be used to color each cell.
#' @param trend_formula the model formula to be used for fitting the expression
#'   trend over pseudotime.
#' @param label_by_short_name label figure panels by gene_short_name (TRUE) or
#'   feature ID (FALSE).
#' @param vertical_jitter A value passed to ggplot to jitter the points in the
#'   vertical dimension. Prevents overplotting, and is particularly helpful for
#'   rounded transcript count data.
#' @param horizontal_jitter A value passed to ggplot to jitter the points in
#'   the horizontal dimension. Prevents overplotting, and is particularly
#'   helpful for rounded transcript count data.
#' @return a ggplot2 plot object
#' @export
plot_genes_in_pseudotime <-function(cds_subset,
                                    min_expr=NULL,
                                    cell_size=0.75,
                                    nrow=NULL,
                                    ncol=1,
                                    panel_order=NULL,
                                    color_cells_by="pseudotime",
                                    trend_formula="~ splines::ns(pseudotime, df=3)",
                                    label_by_short_name=TRUE,
                                    vertical_jitter=NULL,
                                    horizontal_jitter=NULL){
  assertthat::assert_that(methods::is(cds_subset, "cell_data_set"))
  tryCatch({pseudotime(cds_subset)}, error = function(x) {
    stop(paste("No pseudotime calculated. Must call order_cells first."))})
  colData(cds_subset)$pseudotime <- pseudotime(cds_subset)
  if(!is.null(min_expr)) {
    assertthat::assert_that(assertthat::is.number(min_expr))
  }
  assertthat::assert_that(assertthat::is.number(cell_size))

  if(!is.null(nrow)) {
    assertthat::assert_that(assertthat::is.count(nrow))
  }

  assertthat::assert_that(assertthat::is.count(ncol))
  assertthat::assert_that(is.logical(label_by_short_name))
  if (label_by_short_name) {
    assertthat::assert_that("gene_short_name" %in% names(rowData(cds_subset)),
                            msg = paste("When label_by_short_name = TRUE,",
                                        "rowData must have a column of gene",
                                        "names called gene_short_name."))
  }
  assertthat::assert_that(color_cells_by %in% c("cluster", "partition") |
                            color_cells_by %in% names(colData(cds_subset)),
                          msg = paste("color_cells_by must be a column in the",
                                      "colData table."))

  if(!is.null(panel_order)) {
    if (label_by_short_name) {
      assertthat::assert_that(all(panel_order %in%
                                    rowData(cds_subset)$gene_short_name))
    } else {
      assertthat::assert_that(all(panel_order %in%
                                    row.names(rowData(cds_subset))))
    }
  }
  assertthat::assert_that(nrow(rowData(cds_subset)) <= 100,
                          msg = paste("cds_subset has more than 100 genes -",
                                      "pass only the subset of the CDS to be",
                                      "plotted."))

  assertthat::assert_that(methods::is(cds_subset, "cell_data_set"))
  assertthat::assert_that("pseudotime" %in% names(colData(cds_subset)),
                          msg = paste("pseudotime must be a column in",
                                      "colData. Please run order_cells",
                                      "before running",
                                      "plot_genes_in_pseudotime."))
  if(!is.null(min_expr)) {
    assertthat::assert_that(assertthat::is.number(min_expr))
  }
  assertthat::assert_that(assertthat::is.number(cell_size))
  assertthat::assert_that(!is.null(size_factors(cds_subset)))
  if(!is.null(nrow)) {
    assertthat::assert_that(assertthat::is.count(nrow))
  }

  assertthat::assert_that(assertthat::is.count(ncol))
  assertthat::assert_that(is.logical(label_by_short_name))
  if (label_by_short_name) {
    assertthat::assert_that("gene_short_name" %in% names(rowData(cds_subset)),
                            msg = paste("When label_by_short_name = TRUE,",
                                        "rowData must have a column of gene",
                                        "names called gene_short_name."))
  }
  assertthat::assert_that(color_cells_by %in% c("cluster", "partition") |
                            color_cells_by %in% names(colData(cds_subset)),
                          msg = paste("color_cells_by must be a column in the",
                                      "colData table."))

  if(!is.null(panel_order)) {
    if (label_by_short_name) {
      assertthat::assert_that(all(panel_order %in%
                                    rowData(cds_subset)$gene_short_name))
    } else {
      assertthat::assert_that(all(panel_order %in%
                                    row.names(rowData(cds_subset))))
    }
  }
  assertthat::assert_that(nrow(rowData(cds_subset)) <= 100,
                          msg = paste("cds_subset has more than 100 genes -",
                                      "pass only the subset of the CDS to be",
                                      "plotted."))

  f_id <- NA
  Cell <- NA
  cds_subset = cds_subset[,is.finite(colData(cds_subset)$pseudotime)]

  cds_exprs <- SingleCellExperiment::counts(cds_subset)
  cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds_subset))
  cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))

  if (is.null(min_expr)) {
    min_expr <- 0
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_colData <- colData(cds_subset)
  cds_rowData <- rowData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_rowData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_colData, by.x = "Cell", by.y = "row.names")

  cds_exprs$adjusted_expression <- cds_exprs$expression

  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$f_id <- as.character(cds_exprs$f_id)
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)


  new_data <- data.frame(pseudotime = colData(cds_subset)$pseudotime)
  model_tbl = fit_models(cds_subset, model_formula_str = trend_formula)

  model_expectation <- model_predictions(model_tbl,
                                         new_data = colData(cds_subset))

  colnames(model_expectation) <- colnames(cds_subset)
  expectation <- plyr::ddply(cds_exprs, plyr::.(f_id, Cell),
                             function(x) {
                               data.frame(
                                 "expectation"=model_expectation[x$f_id,
                                                                 x$Cell])
                             })
  cds_exprs <- merge(cds_exprs, expectation)

  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
  if (!is.null(panel_order)) {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label,
                                      levels = panel_order)
  }
  q <- ggplot(aes(pseudotime, expression), data = cds_exprs)


  if (!is.null(color_cells_by)) {
    q <- q + geom_point(aes_string(color = color_cells_by),
                        size = I(cell_size),
                        position=position_jitter(horizontal_jitter,
                                                 vertical_jitter))
    if (class(colData(cds_subset)[,color_cells_by]) == "numeric"){
      q <- q + viridis::scale_color_viridis(option="C")
    }
  }
  else {
    q <- q + geom_point(size = I(cell_size),
                        position=position_jitter(horizontal_jitter,
                                                 vertical_jitter))
  }

  q <- q + geom_line(aes(x = pseudotime, y = expectation), data = cds_exprs)

  q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow,
                                        ncol = ncol, scales = "free_y")
  if (min_expr < 1) {
    q <- q + expand_limits(y = c(min_expr, 1))
  }

  q <- q + ylab("Expression")

  q <- q + xlab("pseudotime")
  q <- q + monocle_theme_opts()
  q
}

#' Plots the percentage of variance explained by the each component based on
#' PCA from the normalized expression data determined using preprocess_cds.
#'
#' @param cds cell_data_set of the experiment.
#' @return ggplot object.
#' @export
#' @examples
#' cds <- load_a549()
#' cds <- preprocess_cds(cds)
#' plot_pc_variance_explained(cds)
plot_pc_variance_explained <- function(cds) {
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(!is.null(reducedDims(cds)[["PCA"]]),
                          msg = paste("Data has not been preprocessed with",
                                      "PCA. Please run preprocess_cds with",
                                      "method = 'PCA' before running",
                                      "plot_pc_variance_explained."))
  prop_varex <- cds@preprocess_aux$prop_var_expl

  p <- qplot(1:length(prop_varex), prop_varex, alpha = I(0.5)) +
    monocle_theme_opts() +
    theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
    theme(panel.background = element_rect(fill='white')) +
    xlab('PCA components') +
    ylab('Variance explained \n by each component')
  return(p)
}

#' Plot expression for one or more genes as a violin plot
#'
#' @description Accepts a subset of a cell_data_set and an attribute to group
#' cells by, and produces a ggplot2 object that plots the level of expression
#' for each group of cells.
#'
#' @param cds_subset Subset cell_data_set to be plotted.
#' @param group_cells_by NULL of the cell attribute (e.g. the column of
#'   colData(cds)) to group cells by on the horizontal axis. If NULL, all cells
#'   are plotted together.
#' @param min_expr the minimum (untransformed) expression level to be plotted.
#'   Default is 0.
#' @param nrow the number of panels per row in the figure.
#' @param ncol the number of panels per column in the figure.
#' @param panel_order the order in which genes should be laid out
#'   (left-to-right, top-to-bottom). Should be gene_short_name if
#'   \code{label_by_short_name = TRUE} or feature ID if
#'   \code{label_by_short_name = FALSE}.
#' @param label_by_short_name label figure panels by gene_short_name (TRUE) or
#'   feature id (FALSE). Default is TRUE.
#' @param normalize Logical, whether or not to normalize expression by size
#'   factor. Default is TRUE.
#' @param log_scale Logical, whether or not to scale data logarithmically.
#'   Default is TRUE.
#' @return a ggplot2 plot object
#' @import ggplot2
#' @export
#' @examples
#' cds <- load_a549()
#' cds_subset <- cds[row.names(subset(rowData(cds),
#'                  gene_short_name %in% c("ACTA1", "ID1", "CCNB2"))),]
#' plot_genes_violin(cds_subset, group_cells_by="culture_plate", ncol=2,
#'                   min_expr=0.1)
#'
plot_genes_violin <- function (cds_subset,
                               group_cells_by = NULL,
                               min_expr = 0,
                               nrow = NULL,
                               ncol = 1,
                               panel_order = NULL,
                               label_by_short_name = TRUE,
                               normalize = TRUE,
                               log_scale = TRUE) {

  assertthat::assert_that(methods::is(cds_subset, "cell_data_set"))

  if(!is.null(group_cells_by)) {
    assertthat::assert_that(group_cells_by %in% names(colData(cds_subset)),
                            msg = paste("group_cells_by must be a column in",
                                        "the colData table"))
  }

  assertthat::assert_that(assertthat::is.number(min_expr))

  if(!is.null(nrow)) {
    assertthat::assert_that(assertthat::is.count(nrow))
  }

  assertthat::assert_that(assertthat::is.count(ncol))

  assertthat::assert_that(is.logical(label_by_short_name))
  if (label_by_short_name) {
    assertthat::assert_that("gene_short_name" %in% names(rowData(cds_subset)),
                            msg = paste("When label_by_short_name = TRUE,",
                                        "rowData must have a column of gene",
                                        "names called gene_short_name."))
  }
  if(!is.null(panel_order)) {
    if (label_by_short_name) {
      assertthat::assert_that(all(panel_order %in%
                                    rowData(cds_subset)$gene_short_name))
    } else {
      assertthat::assert_that(all(panel_order %in%
                                    row.names(rowData(cds_subset))))
    }
  }

  assertthat::assert_that(is.logical(normalize))
  assertthat::assert_that(is.logical(log_scale))

  assertthat::assert_that(nrow(rowData(cds_subset)) <= 100,
                          msg = paste("cds_subset has more than 100 genes -",
                                      "pass only the subset of the CDS to be",
                                      "plotted."))

  if (normalize) {
    cds_exprs <- SingleCellExperiment::counts(cds_subset)
    cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds_subset))
    cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
  } else {
    cds_exprs <- SingleCellExperiment::counts(cds_subset)
    cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
  }

  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr

  cds_exprs <- merge(cds_exprs, rowData(cds_subset), by.x = "f_id",
                     by.y = "row.names")
  cds_exprs <- merge(cds_exprs, colData(cds_subset), by.x = "Cell",
                     by.y = "row.names")

  if (label_by_short_name) {
    if (!is.null(cds_exprs$gene_short_name)) {
      cds_exprs$feature_label <- cds_exprs$gene_short_name
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    } else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  } else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }

  if (!is.null(panel_order)) {
    cds_exprs$feature_label = factor(cds_exprs$feature_label,
                                     levels = panel_order)
  }

  cds_exprs[,group_cells_by] <- as.factor(cds_exprs[,group_cells_by])

  q <- ggplot(aes_string(x = group_cells_by, y = "expression"),
              data = cds_exprs) +
    monocle_theme_opts()

  cds_exprs[,group_cells_by] <- as.factor(cds_exprs[,group_cells_by])
  q <- q + geom_violin(aes_string(fill = group_cells_by), scale="width") +
    guides(fill=FALSE)
  q <- q + stat_summary(fun.y=mean, geom="point", size=1, color="black")
  q <- q + facet_wrap(~feature_label, nrow = nrow,
                      ncol = ncol, scales = "free_y")
  if (min_expr < 1) {
    q <- q + expand_limits(y = c(min_expr, 1))
  }

  q <- q + ylab("Expression") + xlab(group_cells_by)

  if (log_scale){
    q <- q + scale_y_log10()
  }
  q
}


#' Plots the number of cells expressing one or more genes above a given value
#' as a barplot
#'
#'  @description Accepts a subset cell_data_set and the parameter
#'  \code{group_cells_by}, used for dividing cells into groups. Returns one or
#'  more bar graphs (one graph for each gene in the cell_data_set). Each graph
#'  shows the percentage (or number) of cells that express a gene in each
#'  sub-group in the cell_data_set.
#'
#' @param cds_subset Subset cell_data_set to be plotted.
#' @param group_cells_by the cell attribute (e.g. the column of colData(cds))
#'   to group cells by on the horizontal axis. If NULL, all cells plotted as
#'   one group.
#' @param min_expr the minimum (untransformed) expression level to consider the
#'   gene 'expressed'. Default is 0.
#' @param nrow the number of panels per row in the figure.
#' @param ncol the number of panels per column in the figure.
#' @param panel_order the order in which genes should be laid out
#'   (left-to-right, top-to-bottom). Should be gene_short_name if
#'   \code{label_by_short_name = TRUE} or feature ID if
#'   \code{label_by_short_name = FALSE}.
#' @param plot_as_count Logical, whether to plot as a count of cells rather
#'   than a percent. Default is FALSE.
#' @param label_by_short_name label figure panels by gene_short_name (TRUE) or
#'   feature id (FALSE). Default is TRUE.
#' @param normalize Logical, whether or not to normalize expression by size
#'   factor. Default is TRUE.
#' @param plot_limits A pair of number specifying the limits of the y axis. If
#'   \code{NULL}, scale to the range of the data. Example \code{c(0,100)}.
#' @param bootstrap_samples The number of bootstrap replicates to generate when
#'   plotting error bars. Default is 100.
#' @param conf_int_alpha The size of the confidence interval to use when plotting
#'   error bars. Default is 0.95.
#' @return a ggplot2 plot object
#' @import ggplot2
#' @export
#' @examples
#' cds <- load_a549()
#' cds_subset <- cds[row.names(subset(rowData(cds),
#'                                   gene_short_name %in% c("NDRG4", "HBG2"))),]
#' plot_percent_cells_positive(cds_subset, group_cells_by="culture_plate")
plot_percent_cells_positive <- function(cds_subset,
                                        group_cells_by = NULL,
                                        min_expr = 0,
                                        nrow = NULL,
                                        ncol = 1,
                                        panel_order = NULL,
                                        plot_as_count = FALSE,
                                        label_by_short_name=TRUE,
                                        normalize = TRUE,
                                        plot_limits = NULL,
                                        bootstrap_samples=100,
                                        conf_int_alpha = .95){

  assertthat::assert_that(methods::is(cds_subset, "cell_data_set"))

  if(!is.null(group_cells_by)) {
    assertthat::assert_that(group_cells_by %in% names(colData(cds_subset)),
                            msg = paste("group_cells_by must be a column in",
                                        "the colData table"))
  }
  assertthat::assert_that(assertthat::is.number(min_expr))

  if(!is.null(nrow)) {
    assertthat::assert_that(assertthat::is.count(nrow))
  }

  assertthat::assert_that(assertthat::is.count(ncol))
  assertthat::assert_that(is.logical(plot_as_count))
  assertthat::assert_that(is.logical(label_by_short_name))
  assertthat::assert_that(is.logical(normalize))

  assertthat::assert_that(nrow(rowData(cds_subset)) <= 100,
                          msg = paste("cds_subset has more than 100 genes -",
                                      "pass only the subset of the CDS to be",
                                      "plotted."))

  marker_exprs <- SingleCellExperiment::counts(cds_subset)

  if (normalize) {
    marker_exprs <- Matrix::t(Matrix::t(marker_exprs)/size_factors(cds_subset))
    marker_exprs_melted <- reshape2::melt(round(as.matrix(marker_exprs)))
  } else {
    marker_exprs_melted <- reshape2::melt(as.matrix(marker_exprs))
  }

  colnames(marker_exprs_melted) <- c("f_id", "Cell", "expression")

  marker_exprs_melted <- merge(marker_exprs_melted, colData(cds_subset),
                               by.x="Cell", by.y="row.names")
  marker_exprs_melted <- merge(marker_exprs_melted, rowData(cds_subset),
                               by.x="f_id", by.y="row.names")

  if (label_by_short_name) {
    if (!is.null(marker_exprs_melted$gene_short_name)){
      marker_exprs_melted$feature_label <- marker_exprs_melted$gene_short_name
      marker_exprs_melted$feature_label[
        is.na(marker_exprs_melted$feature_label)] <- marker_exprs_melted$f_id
    } else {
      marker_exprs_melted$feature_label <- marker_exprs_melted$f_id
    }
  } else {
    marker_exprs_melted$feature_label <- marker_exprs_melted$f_id
  }

  if (!is.null(panel_order)) {
    marker_exprs_melted$feature_label <-
      factor(marker_exprs_melted$feature_label, levels=panel_order)
  }

  if(is.null(group_cells_by)) {
    marker_exprs_melted$all_cell <- "All"
    group_cells_by <- "all_cell"
  }

  # marker_counts <-
  #   plyr::ddply(marker_exprs_melted,
  #               c("feature_label", group_cells_by),
  #               function(x) {
  #                 data.frame(target = sum(x$expression > min_expr),
  #                            target_fraction = sum(x$expression >
  #                                                    min_expr)/nrow(x))
  #               })

  marker_counts_bootstrap = rsample::bootstraps(marker_exprs_melted, times = bootstrap_samples)

  group_mean_bootstrap <- function(split) {
    rsample::analysis(split) %>%
      dplyr::group_by_("feature_label", group_cells_by) %>%
      dplyr::summarize(target = sum(expression > min_expr),
                       target_fraction = sum(expression > min_expr)/dplyr::n())
  }

  marker_counts <-
    marker_counts_bootstrap %>%
    dplyr::mutate(summary_stats = purrr::map(splits, group_mean_bootstrap)) %>%
    tidyr::unnest(summary_stats)
  marker_counts <- marker_counts %>% dplyr::ungroup() %>%
    dplyr::group_by_("feature_label", group_cells_by) %>%
    dplyr::summarize(target_mean = mean(target),
                     target_fraction_mean = mean(target_fraction),
                     target_low = quantile(target, conf_int_alpha / 2),
                     target_high = quantile(target, 1 - conf_int_alpha / 2),
                     target_fraction_low = quantile(target_fraction, (1 - conf_int_alpha) / 2),
                     target_fraction_high = quantile(target_fraction, 1 - (1 - conf_int_alpha) / 2))


  # marker_counts <-
  #   marker_exprs_melted %>% dplyr::group_by_("feature_label", group_cells_by) %>%
  #   dplyr::summarize(target = sum(expression > min_expr),
  #             target_fraction = sum(expression > min_expr)/dplyr::n())

  if (!plot_as_count){
    marker_counts$target_fraction_mean <- marker_counts$target_fraction_mean * 100
    marker_counts$target_fraction_low <- marker_counts$target_fraction_low * 100
    marker_counts$target_fraction_high <- marker_counts$target_fraction_high * 100
    qp <- ggplot(aes_string(x=group_cells_by, y="target_fraction_mean",
                            fill=group_cells_by),
                 data=marker_counts) +
      ylab("Cells (percent)")
  } else {
    qp <- ggplot(aes_string(x=group_cells_by, y="target_mean", fill=group_cells_by),
                 data=marker_counts) +
      ylab("Cells")
  }
  if (group_cells_by == "all_cell") {
    group_cells_by <- ""
    qp <- qp + theme(legend.title = element_blank())
  }
  qp <- qp + xlab(group_cells_by)

  if (is.null(plot_limits) == FALSE) {
    qp <- qp + scale_y_continuous(limits=plot_limits)
  }

  qp <- qp + facet_wrap(~feature_label, nrow=nrow, ncol=ncol, scales="free_y")
  qp <-  qp +
    geom_bar(stat="identity") +
    geom_linerange(aes(ymin=target_fraction_low, ymax=target_fraction_high))+
    monocle_theme_opts()

  return(qp)
}

#' Create a dot plot to visualize the mean gene expression and percentage of
#' expressed cells in each group of cells
#'
#' @param cds A cell_data_set for plotting.
#' @param markers A list of gene ids (or short names) to show in the plot
#' @param group_cells_by How to group cells when labeling them. Must be either
#'   the name of a column of colData(cds), or one of "clusters" or "partitions".
#'   If a column in colData(cds), must be a categorical variable.
#' @param reduction_method The dimensionality reduction method used for clusters
#'   and partitions.
#' @param norm_method Determines how to transform expression values prior to
#'   plotting. Options are "log" and "size_only". Default is "log".
#' @param lower_threshold The lowest gene expressed treated as expressed. By
#'   default, zero.
#' @param max.size The maximum size of the dot. By default, it is 10.
#' @param ordering_type How to order the genes / groups on the dot plot. Only
#'   accepts 'cluster_row_col' (use biclustering to cluster the rows and
#'   columns), 'maximal_on_diag' (position each column so that the maximal color
#'   shown on each column on the diagonal, if the current maximal is used in
#'   earlier columns, the next largest one is position), and 'none' (preserve
#'   the ordering from the input gene or alphabetical ordering of groups).
#'   Default is 'cluster_row_col'.
#' @param axis_order Whether to put groups on x-axis, genes on y-axis (option
#'   'group_marker') or the reverse order (option 'marker_group'). Default is
#'   "group_marker".
#' @param flip_percentage_mean Logical indicating whether to use color of the
#'   dot to represent the percentage (by setting flip_percentage_mean = FALSE,
#'   default) and size of the dot the mean expression, or the opposite (by
#'   setting flip_percentage_mean = TRUE).
#' @param pseudocount A pseudo-count added to the average gene expression.
#' @param scale_max The maximum value (in standard deviations) to show in the
#'   heatmap. Values larger than this are set to the max.
#' @param scale_min The minimum value (in standard deviations) to show in the
#'   heatmap. Values smaller than this are set to the min.
#'
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom reshape2 dcast
#' @importFrom viridis scale_color_viridis
#' @export
plot_genes_by_group <- function(cds,
                                markers,
                                group_cells_by="cluster",
                                reduction_method = "UMAP",
                                norm_method = c("log", "size_only"),
                                lower_threshold = 0,
                                max.size = 10,
                                # maybe be also do the maximum color on the
                                # diagonal; the axis change be switched too
                                ordering_type = c('cluster_row_col',
                                                  'maximal_on_diag',
                                                  'none'),
                                axis_order = c('group_marker', 'marker_group'),
                                flip_percentage_mean = FALSE,
                                pseudocount = 1,
                                scale_max = 3,
                                scale_min = -3) {

  assertthat::assert_that(methods::is(cds, "cell_data_set"))

  if(!is.null(group_cells_by)) {
    assertthat::assert_that(group_cells_by %in% c("cluster", "partition") |
                              group_cells_by %in% names(colData(cds)),
                            msg = paste("group_cells_by must be a column in",
                                        "the colData table."))
  }

  norm_method = match.arg(norm_method)
  #group_cells_by=match.arg(group_cells_by)

  gene_ids = as.data.frame(fData(cds)) %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(rowname %in% markers | gene_short_name %in% markers) %>%
    dplyr::pull(rowname)
  if(length(gene_ids) < 1)
    stop(paste('Please make sure markers are included in the gene_short_name",
               "column of the fData!'))

  if(flip_percentage_mean == FALSE){
    major_axis <- 1
    minor_axis <- 2
  } else if (flip_percentage_mean == TRUE){
    major_axis <- 2
    minor_axis <- 1
  }

  exprs_mat <- t(as.matrix(exprs(cds)[gene_ids, ]))
  exprs_mat <- reshape2::melt(exprs_mat)
  colnames(exprs_mat) <- c('Cell', 'Gene', 'Expression')
  exprs_mat$Gene <- as.character(exprs_mat$Gene)


  if (group_cells_by == "cluster"){
    cell_group <- tryCatch({clusters(cds,
                                     reduction_method = reduction_method)},
                           error = function(e) {NULL})
  } else if (group_cells_by == "partition") {
    cell_group <- tryCatch({partitions(cds,
                                       reduction_method = reduction_method)},
                           error = function(e) {NULL})
  } else{
    cell_group <- colData(cds)[,group_cells_by]
  }
  names(cell_group) = colnames(cds)

  exprs_mat$Group <- cell_group[exprs_mat$Cell]
  exprs_mat = exprs_mat %>% dplyr::filter(is.na(Group) == FALSE)
  ExpVal <- exprs_mat %>% dplyr::group_by(Group, Gene) %>%
    dplyr::summarize(mean = mean(log(Expression + pseudocount)),
                     percentage = sum(Expression > lower_threshold) /
                       length(Expression))
  ExpVal$mean <- ifelse(ExpVal$mean < scale_min, scale_min, ExpVal$mean)
  ExpVal$mean <- ifelse(ExpVal$mean > scale_max, scale_max, ExpVal$mean)

  ExpVal$Gene <- fData(cds)[ExpVal$Gene, 'gene_short_name']

  res <- reshape2::dcast(ExpVal[, 1:4], Group ~ Gene,
                         value.var = colnames(ExpVal)[2 + major_axis])
  group_id <- res[, 1]
  res <- res[, -1]
  row.names(res) <- group_id

  if(ordering_type == 'cluster_row_col') {
    row_dist <- stats::as.dist((1 - stats::cor(t(res[, -1])))/2)
    row_dist[is.na(row_dist)] <- 1

    col_dist <- stats::as.dist((1 - stats::cor(res[, -1]))/2)
    col_dist[is.na(col_dist)] <- 1

    ph <- pheatmap::pheatmap(res[, -1],
                             useRaster = T,
                             cluster_cols=TRUE,
                             cluster_rows=TRUE,
                             show_rownames=F,
                             show_colnames=F,
                             clustering_distance_cols=col_dist,
                             clustering_distance_rows=row_dist,
                             clustering_method = 'ward.D2',
                             silent=TRUE,
                             filename=NA)

    ExpVal$Gene <- factor(ExpVal$Gene,
                          levels = colnames(res)[-1][ph$tree_col$order])
    ExpVal$Group <- factor(ExpVal$Group,
                           levels = row.names(res)[ph$tree_row$order])

  } else if(ordering_type == 'maximal_on_diag'){

    order_mat <- t(apply(res, major_axis, order))
    max_ind_vec <- c()
    for(i in 1:nrow(order_mat)) {
      tmp <- max(which(!(order_mat[i, ] %in% max_ind_vec)))
      max_ind_vec <- c(max_ind_vec, order_mat[i, tmp])
    }
    max_ind_vec <- max_ind_vec[!is.na(max_ind_vec)]

    if(major_axis == 1){
      max_ind_vec <- c(max_ind_vec, setdiff(1:length(markers), max_ind_vec))
      ExpVal$Gene <- factor(ExpVal$Gene ,
                            levels = dimnames(res)[[2]][max_ind_vec])
    }
    else{
      max_ind_vec <- c(max_ind_vec, setdiff(1:length(unique(exprs_mat$Group)),
                                            max_ind_vec))
      ExpVal$Group <- factor(ExpVal$Group,
                             levels = dimnames(res)[[1]][max_ind_vec])
    }
  } else if(ordering_type == 'none'){
    ExpVal$Gene <- factor(ExpVal$Gene, levels = markers)
  }

  if(flip_percentage_mean){
    g <- ggplot(ExpVal, aes(y = Gene,  x = Group)) +
      geom_point(aes(colour = percentage,  size = mean)) +
      viridis::scale_color_viridis(name = 'percentage') +
      scale_size(name = 'log(mean + 0.1)', range = c(0, max.size))
  } else {
    g <- ggplot(ExpVal, aes(y = Gene,  x = Group)) +
      geom_point(aes(colour = mean,  size = percentage)) +
      viridis::scale_color_viridis(name = 'log(mean + 0.1)') +
      scale_size(name = 'percentage', range = c(0, max.size))
  }

  if (group_cells_by == "cluster"){
    g <- g + xlab("Cluster")
  } else if (group_cells_by == "partition") {
    g <- g + xlab("Partition")
  } else{
    g <- g + xlab(group_cells_by)
  }

  g <- g + ylab("Gene") + monocle_theme_opts() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  if(axis_order == 'marker_group') {
    g <- g + coord_flip()
  }

  g
}

