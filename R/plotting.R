monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

#' Plot a dataset and trajectory in 3 dimensions
#'
#' @param cds cell_data_set for the experiment
#' @param dim the dimensions used to create the 3D plot, by default it is the first three dimensions
#' @param color_by the cell attribute (e.g. the column of colData(cds)) to map to each cell's color
#' @param markers a gene name or gene id to use for setting the size of each cell in the plot
#' @param markers_linear a boolean used to indicate whether you want to scale the markers logarithimically or linearly
#' @param webGL_filename the name of a file to which you'd to write a webGL file containing the plot
#' @param image_filename the name of a file to which you'd to write an animated GIF containing the plot
#' @param obj_filename The filename that you'd like to save the 3d plot to
#' @param view_matrix A 4x4 matrix signifying your point of view on the 3d plot
#' @param scale_expr whether to tranform the log expression values to z scores
#' @param palette the color palette used for plotting
#' @param width the width of the plot in pixels
#' @param height the height of the plot in pixels
#' @param useNULL_GLdev if TRUE, don't show the plot on the screen (to be used with webGL or movie output)
#' @param show_group_labels Whether or not to show labels of groups
#' @param show_branch_points Whether to show icons for each branch point (only available after running assign_cell_states)
#' @param cell_size size of cells, default value is 5
#' @param backbone_segment_color color of backbone, value is a string containing a hex color value
#' @param backbone_vertex_color Color for the vertex on the principal graph backbone
#' @param cell_alpha Alpha value for the vertex in the 3D plot
#' @param ... Extra arguments passed to the function
#' @return a ggplot2 plot object
#' @export
#' @examples
#' \dontrun{
#' plot_cells_3d(cds, markers=c("Rbfox3, Neurod1", "Sox2"))
#' }
#'
plot_cells_3d <- function(cds,
                          reduction_method = "UMAP",
                          downsample_size = NULL,
                          dim=c(1, 2, 3),
                          color_by=NULL,
                          palette = NULL,
                          markers=NULL,
                          markers_linear=FALSE,
                          cell_size=5,
                          cell_alpha=0.5,
                          backbone_segment_color="#77B6EA",
                          backbone_segment_width=2,
                          backbone_vertex_color=NULL,
                          show_group_labels=TRUE,
                          show_branch_points=TRUE,
                          webGL_filename=NULL,
                          image_filename=NULL,
                          obj_filename=NULL,
                          scale_expr=TRUE,
                          width=800,
                          height=600,
                          view_matrix=NULL,
                          useNULL_GLdev = !interactive(),
                          text_cex = 1,
                          use_plotly = FALSE,
                          selfcontained=FALSE,
                          ...){
  gene_short_name <- NA
  sample_name <- NA

  #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
  if(is.null(downsample_size)) {
    cell_sampled <- 1:ncol(cds)
  } else {
    if(class(downsample_size) %in% c('integer', 'numeric')) {
      if(downsample_size < 2) {
        stop('Error: downsample_size must be at least larger 2!')
      }
      cell_sampled <- sample(1:ncol(cds), downsample_size)
    } else {
      stop('Error: downsample_size must be an integer!')
    }
  }

  lib_info_with_pseudo <- colData(cds)[cell_sampled, ]
  sample_state <- colData(cds)$State[cell_sampled]

  reduced_dim_coords <- cds@principal_graph_aux[[reduction_method]]$dp_mst

  if(length(dim) != 3)
    dim <- 1:3

  ica_space_df <- data.frame(reduced_dim_coords[,dim])
  colnames(ica_space_df) <- c("prin_graph_dim_1", "prin_graph_dim_2",  "prin_graph_dim_3")

  ica_space_df$sample_name <- row.names(ica_space_df)
  ica_space_df$sample_state <- row.names(ica_space_df)
  #ica_space_with_state_df <- merge(ica_space_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")
  #print(ica_space_with_state_df)
  dp_mst <- principal_graph(cds)[[reduction_method]]

  if (is.null(dp_mst)){
    stop("You must first call orderCells() before using this function")
  }

  edge_list <- as.data.frame(igraph::get.edgelist(dp_mst))
  colnames(edge_list) <- c("source", "target")

  edge_df <- merge(edge_list, ica_space_df, by.x="source", by.y="sample_name", all=F)
  edge_df <- plyr::rename(edge_df, c("prin_graph_dim_1"="source_prin_graph_dim_1", "prin_graph_dim_2"="source_prin_graph_dim_2", "prin_graph_dim_3"="source_prin_graph_dim_3"))
  edge_df <- merge(edge_df, ica_space_df[,c("sample_name", "prin_graph_dim_1", "prin_graph_dim_2", "prin_graph_dim_3")], by.x="target", by.y="sample_name", all=F)
  edge_df <- plyr::rename(edge_df, c("prin_graph_dim_1"="target_prin_graph_dim_1", "prin_graph_dim_2"="target_prin_graph_dim_2", "prin_graph_dim_3"="target_prin_graph_dim_3"))

  S_matrix <- reducedDims(cds)[[reduction_method]][cell_sampled,]
  data_df <- data.frame(S_matrix[,dim])
  #data_df <- cbind(data_df, sample_state)
  colnames(data_df) <- c("data_dim_1", "data_dim_2", "data_dim_3")
  data_df$sample_name <- row.names(data_df)
  data_df <- merge(data_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")

  markers_exprs = NULL
  if (is.null(markers) == FALSE){
    markers_rowData <- subset(rowData(cds), gene_short_name %in% markers)
    if (nrow(markers_rowData) >= 1){

      markers_expr_val <- counts(cds)[row.names(markers_rowData), cell_sampled]
      markers_expr_val <- Matrix::colSums(markers_expr_val)
      markers_expr_val <- markers_expr_val / colData(cds)$Size_Factor
      nz_points = markers_expr_val != 0
      if (scale_expr){
        markers_expr_val[nz_points] <- scale(log10(markers_expr_val[nz_points]))
        markers_expr_val[markers_expr_val < -3] = -3
        markers_expr_val[markers_expr_val > 3] = 3
      }
      markers_expr_val[markers_expr_val == 0] = NA

      markers_exprs <- data.frame(cell_id = colnames(cds)[cell_sampled])
      markers_exprs$value <- markers_expr_val
    }
  }

  point_colors_df <- data.frame(sample_name = data_df$sample_name,
                                point_colors = "darkgray")

  map2color<-function(x, limits=NULL){
    inf_vals = is.infinite(x)

    if (is.null(limits) == FALSE){
      x[x < limits[1]] = limits[1]
      x[x > limits[2]] = limits[2]
    }
    x[inf_vals] = stats::median(x)
    ii <- cut(x, breaks = seq(min(x, na.rm=T), max(x, na.rm=T), len = 100),
              include.lowest = TRUE)
    #colors <- colorRampPalette(c("white", "#E05263"))(99)[ii]
    colors = viridis::viridis(99)[ii]
    colors[inf_vals] = "darkgray"
    return(colors)
  }

  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    data_df <- merge(data_df, markers_exprs, by.x="sample_name", by.y="cell_id")

    if(markers_linear || scale_expr){
      point_colors_df$point_colors = map2color(data_df$value, c(-3, 3))
    } else {
      point_colors_df$point_colors = map2color(log10(data_df$value+0.1))
    }
  }else if (is.null(color_by) == FALSE && color_by %in% colnames(colData(cds))){
    point_colors_df = dplyr::data_frame(sample_name=data_df$sample_name,
                                        color_by=as.character(colData(cds)[data_df$sample_name,color_by]))

    if (class(colData(cds)[,color_by]) == "numeric"){
      point_colors_df$point_colors = map2color(colData(cds)[data_df$sample_name,color_by])
      point_colors_df$point_colors[is.na(point_colors_df$point_colors)] = "darkgray"
    }else{
      gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
      }

      if (is.null(palette)){

        if(is.factor(colData(cds)[,color_by])) {
          num_colors = length(levels(colData(cds)[,color_by]))
          point_colors = gg_color_hue(num_colors)
          names(point_colors) = levels(colData(cds)[,color_by])
        } else {
          num_colors = length(unique(colData(cds)[,color_by]))
          point_colors = gg_color_hue(num_colors)
          names(point_colors) = unique(colData(cds)[,color_by])
        }
      }else{
        point_colors = palette
      }

      point_colors_df$point_colors = point_colors[point_colors_df$color_by]
      #point_colors = colors[as.character(colData(cds)[data_df$sample_name,color_by])]
      point_colors_df$point_colors[is.na(point_colors_df$point_colors)] = "darkgray"
    }
  }

  if(!use_plotly) {
    rgl::open3d(windowRect=c(0,0,width,height),
                useNULL=useNULL_GLdev)
    if (is.null(view_matrix) == FALSE){
      rgl::view3d(userMatrix=view_matrix)
    }
    if (is.null(backbone_segment_color) == FALSE){
      rgl::segments3d(matrix(as.matrix(t(edge_df[,c(3,4,5, 7,8,9)])), ncol=3, byrow=T), lwd=backbone_segment_width,
                      col=backbone_segment_color,
                      line_antialias=TRUE)
      if (is.null(backbone_vertex_color) == FALSE)
        rgl::points3d(Matrix::t(reduced_dim_coords[dim,]), col=backbone_vertex_color)
    }
  }

  point_colors_df$point_alpha = cell_alpha
  point_colors_df$point_alpha[is.na(point_colors_df$point_colors)] = 0

  if(!use_plotly) {
    rgl::points3d(data_df[,c("data_dim_1", "data_dim_2", "data_dim_3")],
                  size = cell_size,
                  col=point_colors_df$point_colors,
                  alpha=point_colors_df$point_alpha,
                  shininess=75,
                  point_antialias=TRUE)
  }
  #bg3d(fogtype="linear")
  medoid_df <- NULL
  if (is.null(point_colors_df$color_by) == FALSE){
    point_colors_df = dplyr::inner_join(point_colors_df, data_df)

    medoid_df = point_colors_df %>% dplyr::group_by(color_by, point_colors) %>% dplyr::summarize(mean_d1 = stats::median(data_dim_1),
                                                                                                 mean_d2 = stats::median(data_dim_2),
                                                                                                 mean_d3 = stats::median(data_dim_3))
    if (show_group_labels && color_by %in% colnames(colData(cds)) && class(colData(cds)[,color_by]) != "numeric"){
      if(!use_plotly) {
        rgl::text3d(x=medoid_df$mean_d1, y=medoid_df$mean_d2, z=medoid_df$mean_d3, texts=as.character(medoid_df$color_by))
      }
    }
  }
  branch_point_df <- NULL
  if (show_branch_points){
    mst_branch_nodes <- cds@principal_graph_aux[[reduction_method]]$branch_points
    branch_point_df_source <- edge_df %>%
      dplyr::slice(match(mst_branch_nodes, source)) %>%
      dplyr::mutate(branch_point_idx = which(mst_branch_nodes == source))

    branch_point_df_target <- edge_df %>%
      dplyr::slice(match(mst_branch_nodes, target)) %>%
      dplyr::mutate(branch_point_idx = which(mst_branch_nodes == target))

    if(!use_plotly) {
      rgl::points3d(branch_point_df_source[,c("source_prin_graph_dim_1", "source_prin_graph_dim_2", "source_prin_graph_dim_3")],
                    size = 20,
                    col='black',
                    alpha=0.3,
                    shininess=75,
                    point_antialias=TRUE)

      rgl::text3d(x=branch_point_df_source$source_prin_graph_dim_1, y=branch_point_df_source$source_prin_graph_dim_2, z=branch_point_df_source$source_prin_graph_dim_3,
                  texts=as.character(branch_point_df_source$branch_point_idx), color = 'red', cex = text_cex)

      rgl::points3d(branch_point_df_target[,c("target_prin_graph_dim_1", "target_prin_graph_dim_2", "target_prin_graph_dim_3")],
                    size = 20,
                    col='black',
                    alpha=0.3,
                    shininess=75,
                    point_antialias=TRUE)

      if(show_group_labels){
        rgl::text3d(x=branch_point_df_target$target_prin_graph_dim_1, y=branch_point_df_target$target_prin_graph_dim_2, z=branch_point_df_target$target_prin_graph_dim_3,
                    texts=as.character(branch_point_df_target$branch_point_idx), color = 'red', cex = text_cex)
      }
    }
  }

  if(use_plotly){
    p <- plotly::plot_ly(data_df, x = ~data_dim_1, y = ~data_dim_2, z = ~data_dim_3)  %>%
      plotly::add_markers(
        type = "scatter", mode = "markers", text = ~sample_name,
        opacity = cell_alpha,
        color = I(point_colors_df$point_colors), marker  = list(size = cell_size), hoverinfo = 'text')  #%>%, showlegend=FALSE, showscale = FALSE
    # hide_guides

    if(nrow(edge_df) > 2000 | !is.null(backbone_segment_color)) {
      for (i in 1:nrow(edge_df)) { # add_segments(x = ~x, xend = ~x, y = ~ymin, yend = ~ymax)
        p <- p %>% plotly::add_trace(x = as.vector(t(edge_df[i, c(3, 7)])), y = as.vector(t(edge_df[i, c(4, 8)])), z = as.vector(t(edge_df[i, c(5, 9)])),
                                     line = list(color = I(backbone_segment_color), width = 4),
                                     marker = list(color = I(backbone_vertex_color), size = 5),
                                     mode = 'lines+markers', type = 'scatter3d')
      }
    }
    p <- plotly::hide_legend(p)

    axx <- list(
      title = paste0('component ', 1) #, showgrid = F
    )

    axy <- list(
      title = paste0('component ', 2) #, showgrid = F
    )

    axz <- list(
      title = paste0('component ', 3) #, showgrid = F
    )

    p <- p %>% plotly::layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
    if(is.null(point_colors_df$color_by) == FALSE) {
      p <- p %>% plotly::add_trace(x = medoid_df$mean_d1, y = medoid_df$mean_d2, z = medoid_df$mean_d3, type = "scatter3d", text = as.character(medoid_df$color_by), mode = "text")
    }

    if(nrow(branch_point_df_source) > 0) {
      p <- p %>% plotly::add_trace(x = branch_point_df_source$source_prin_graph_dim_1, y = branch_point_df_source$source_prin_graph_dim_2, z = branch_point_df_source$source_prin_graph_dim_3, type = "scatter3d", text = as.character(branch_point_df_source$branch_point_idx), mode = "text", color = I('black'))
      p <- p %>% plotly::add_trace(x = branch_point_df_source$source_prin_graph_dim_1, y = branch_point_df_source$source_prin_graph_dim_2, z = branch_point_df_source$source_prin_graph_dim_3,
                                   type = "scatter3d", mode = "markers", marker = list(color = I('red'), size = 10))
    }

    if(nrow(branch_point_df_target) > 0) {
      p <- p %>% plotly::add_trace(x = branch_point_df_target$target_prin_graph_dim_1, y = branch_point_df_target$target_prin_graph_dim_2, z = branch_point_df_target$target_prin_graph_dim_3, type = "scatter3d", text = as.character(branch_point_df_target$branch_point_idx), mode = "text", color = I('black'))
      p <- p %>% plotly::add_trace(x = branch_point_df_target$target_prin_graph_dim_1, y = branch_point_df_target$target_prin_graph_dim_2, z = branch_point_df_target$target_prin_graph_dim_3,
                                   type = "scatter3d", mode = "markers",  marker = list(color = I('red'), size = 10))
    }

    print(p)
  }

  if (is.null(obj_filename) == FALSE){
    if(!use_plotly) {
      writeOBJ(obj_filename)
    } else {
      warning("plotly doesn't support writting to a object file!")
    }
  }

  if (is.null(image_filename) == FALSE){
    if(!use_plotly) {
      rgl::rgl.snapshot(image_filename)
      grDevices::graphics.off()
    } else {
      plotly::orca(p, image_filename)
    }
  }

  widget = NULL
  if (is.null(webGL_filename) == FALSE){
    if(!use_plotly) {
      widget <- rgl::rglwidget(elementId = "example",
                               controllers = "player",
                               sizingPolicy = htmlwidgets::sizingPolicy(
                                 browser.fill = TRUE
                               ))
      htmlwidgets::saveWidget(widget, webGL_filename, selfcontained=selfcontained)
    } else {
      htmlwidgets::saveWidget(p, webGL_filename, selfcontained=selfcontained)
    }
  }

  return(widget)
}


#' Plots the cells along with their trajectories.
#'
#' @param cds cell_data_set for the experiment
#' @param x the column of reducedDimS(cds) to plot on the horizontal axis
#' @param y the column of reducedDimS(cds) to plot on the vertical axis
#' @param color_by the cell attribute (e.g. the column of colData(cds)) to map to each cell's color
#' @param show_backbone whether to show the diameter path of the MST used to order the cells
#' @param backbone_color the color used to render the backbone.
#' @param markers a gene name or gene id to use for setting the size of each cell in the plot
#' @param markers_linear a boolean used to indicate whether you want to scale the markers logarithimically or linearly
#' @param show_cell_names draw the name of each cell in the plot
#' @param show_state_number show state number
#' @param cell_size The size of the point for each cell
#' @param cell_link_size The size of the line segments connecting cells (when used with ICA) or the principal graph (when used with DDRTree)
#' @param cell_name_size the size of cell name labels
#' @param state_number_size the size of the state number
#' @param show_branch_points Whether to show icons for each branch point (only available after running assign_cell_states)
#' @param theta How many degrees you want to rotate the trajectory
#' @param alpha The alpha aesthetics for the original cell points, useful to highlight the learned principal graph
#' @param ... Additional arguments passed into scale_color_viridis function
#' @return a ggplot2 plot object
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung()
#' plot_cells(lung)
#' plot_cells(lung, color_by="Pseudotime", show_backbone=FALSE)
#' plot_cells(lung, markers="MYH3")
#' }
plot_cells <- function(cds,
                       x=1,
                       y=2,
                       reduction_method = c("UMAP", "tSNE"),
                       color_cells_by="cluster",
                       group_cells_by=c("cluster", "partition"),
                       genes=NULL,
                       show_trajectory_graph=TRUE,
                       trajectory_graph_color="black",
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
                       alpha = 1,
                       min_expr=0.1,
                       rasterize=FALSE) {
  reduction_method <- match.arg(reduction_method)
  assertthat::assert_that(is(cds, "cell_data_set"))
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
    assertthat::assert_that(color_cells_by %in% c("cluster", "partition") |
                              color_cells_by %in% names(colData(cds)),
                            msg = paste("color_cells_by must be a column in the",
                                        "colData table."))
  }

  assertthat::assert_that(!is.null(color_cells_by) || !is.null(markers),
                          msg = paste("Either color_cells_by or markers must be",
                                      "NULL, cannot color by both!"))

  #if (!is.null(color_cells_by) && color_cells_by == "cluster" && length(clusters(cds, reduction_method = reduction_method)) == 0){
  #  stop("Error: Clustering is not performed yet. Please call clusterCells() before calling this function.")
  #}
  norm_method = match.arg(norm_method)
  group_cells_by=match.arg(group_cells_by)
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(genes),
                          msg = paste("Either color_cells_by or genes must be",
                                      "NULL, cannot color by both!"))

  if (show_trajectory_graph && is.null(principal_graph(cds)[[reduction_method]])) {
    message("No trajectory to plot. Has learn_graph() been called yet?")
    show_trajectory_graph = FALSE
  }

  gene_short_name <- NA
  sample_name <- NA
  sample_state <- colData(cds)$State
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
    data_df$cell_group = tryCatch({clusters(cds, reduction_method = reduction_method)[data_df$sample_name]}, error = function(e) {NULL})
  } else if (group_cells_by == "partition") {
    data_df$cell_group = tryCatch({partitions(cds, reduction_method = reduction_method)[data_df$sample_name]}, error = function(e) {NULL})
  } else{
    stop("Error: unrecognized way of grouping cells.")
  }

  if (color_cells_by == "cluster"){
    data_df$cell_color = tryCatch({clusters(cds, reduction_method = reduction_method)[data_df$sample_name]}, error = function(e) {NULL})
  } else if (color_cells_by == "partition") {
    data_df$cell_color = tryCatch({partitions(cds, reduction_method = reduction_method)[data_df$sample_name]}, error = function(e) {NULL})
  } else{
    data_df$cell_color = colData(cds)[data_df$sample_name,color_cells_by]
  }

  ## Graph info
  if (show_trajectory_graph) {

    ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
      as.data.frame() %>%
      dplyr::select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
      dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))

    dp_mst <- cds@principal_graph[[reduction_method]]

    edge_df <- dp_mst %>%
      igraph::as_data_frame() %>%
      dplyr::select_(source = "from", target = "to") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select_(source="sample_name",
                                        source_prin_graph_dim_1="prin_graph_dim_1",
                                        source_prin_graph_dim_2="prin_graph_dim_2"),
                       by = "source") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select_(target="sample_name",
                                        target_prin_graph_dim_1="prin_graph_dim_1",
                                        target_prin_graph_dim_2="prin_graph_dim_2"),
                       by = "target")
  }

  ## Marker genes
  markers_exprs <- NULL
  if (!is.null(genes)) {
    if ((is.null(dim(genes)) == FALSE) && dim(genes) >= 2){
      markers = unlist(genes[,1], use.names=FALSE)
    } else {
      markers = genes
    }
    markers_rowData <- as.data.frame(subset(rowData(cds),
                                            gene_short_name %in% markers |
                                              rownames(rowData(cds)) %in% markers))
    if (nrow(markers_rowData) >= 1) {
      cds_exprs <- counts(cds)[row.names(markers_rowData), ,drop=FALSE]
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds))

      if ((is.null(dim(genes)) == FALSE) && dim(genes) >= 2){
        genes = as.data.frame(genes)
        row.names(genes) = genes[,1]
        genes = genes[row.names(cds_exprs),]
        agg_mat = as.matrix(Matrix.utils::aggregate.Matrix(cds_exprs, as.factor(genes[,2]), fun="sum"))
        agg_mat = t(scale(t(log10(agg_mat + 1))))
        agg_mat[agg_mat < -2] = -2
        agg_mat[agg_mat > 2] = 2
        markers_exprs = agg_mat
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
        #markers_exprs <- merge(markers_exprs, markers_rowData, by.x = "feature_id", by.y="row.names")
        #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
        markers_exprs$feature_label <- markers_exprs$feature_id
        markers_linear=TRUE
      } else {
        cds_exprs@x = round(cds_exprs@x)
        markers_exprs = matrix(cds_exprs, nrow=nrow(markers_rowData))
        colnames(markers_exprs) = colnames(counts(cds))
        row.names(markers_exprs) = row.names(markers_rowData)
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
        markers_exprs <- merge(markers_exprs, markers_rowData, by.x = "feature_id", by.y="row.names")
        #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
        markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
        markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$feature_id
        markers_exprs$feature_label <- factor(markers_exprs$feature_label,levels = markers)
        #markers_exprs = with(markers_exprs, ifelse(value > min_expr, value, NA))
      }
    }
  }

  if (label_cell_groups && is.null(color_cells_by) == FALSE){
    if (is.null(data_df$cell_color)){
      message(paste(color_cells_by, "not found in colData(cds), cells will not be colored"))
      text_df = NULL
      label_cell_groups = FALSE
    }else{
      if(is.character(data_df$cell_color) || is.factor(data_df$cell_color)) {

        if (label_groups_by_cluster && is.null(data_df$cell_group) == FALSE){
          #text_df = data_df %>% dplyr::color_cells_by_("cluster", color_cells_by)
          text_df = data_df %>%
            dplyr::group_by(cell_group) %>%
            dplyr::mutate(cells_in_cluster= dplyr::n()) %>%
            dplyr::group_by(cell_color, add=TRUE) %>%
            dplyr::mutate(per=dplyr::n()/cells_in_cluster)
          median_coord_df = text_df %>% dplyr::summarize(fraction_of_group = n(),
                                                         text_x = median(x = data_dim_1),
                                                         text_y = median(x = data_dim_2))
          text_df = text_df %>% dplyr::select(per) %>% dplyr::distinct()
          text_df = dplyr::inner_join(text_df, median_coord_df)
          text_df = text_df %>% dplyr::group_by(cell_group) %>% dplyr::top_n(labels_per_group, per)
        } else {
          text_df = data_df %>% dplyr::group_by(cell_color) %>% dplyr::mutate(per=1)
          median_coord_df = text_df %>% dplyr::summarize(fraction_of_group = n(),
                                                         text_x = median(x = data_dim_1),
                                                         text_y = median(x = data_dim_2))
          text_df = text_df %>% dplyr::select(per) %>% dplyr::distinct()
          text_df = dplyr::inner_join(text_df, median_coord_df)
          text_df = text_df %>% dplyr::group_by(cell_color) %>% dplyr::top_n(labels_per_group, per)
        }

        text_df$label = as.character(text_df %>% dplyr::pull(cell_color))
        # I feel like there's probably a good reason for the bit below, but I hate it and I'm killing it for now.
        # text_df$label <- paste0(1:nrow(text_df))
        # text_df$process_label <- paste0(1:nrow(text_df), '_', as.character(as.matrix(text_df[, 1])))
        # process_label <- text_df$process_label
        # names(process_label) <- as.character(as.matrix(text_df[, 1]))
        # data_df[, group_by] <- process_label[as.character(data_df[, group_by])]
        # text_df$label = process_label
      } else {
        message("Cells aren't colored in a way that allows them to be grouped.")
        text_df = NULL
        label_cell_groups = FALSE
      }
    }
  }

  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    data_df <- merge(data_df, markers_exprs, by.x="sample_name", by.y="cell_id")
    data_df$value <- with(data_df, ifelse(value >= min_expr, value, NA))
    if(norm_method == "size_only"){
      g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + plotting_func(aes(color=value), size=I(cell_size), stroke = I(cell_size / 2), na.rm = TRUE) +
        viridis::scale_color_viridis(option = "viridis", name = "Expression", na.value = "grey80", end = 0.8) +
        guides(alpha = FALSE) + facet_wrap(~feature_label)
      # g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + geom_point(aes(color= value), size=I(cell_size), na.rm = TRUE, alpha = alpha) +
      #     viridis::scale_color_viridis(name = paste0("value"), ...) + facet_wrap(~feature_label)
    } else {
      g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + plotting_func(aes(color=log10(value+min_expr), alpha = ifelse(!is.na(value), "2", "1")), size=I(cell_size), stroke = I(cell_size / 2), na.rm = TRUE) +
        viridis::scale_color_viridis(option = "viridis", name = "log10(Expression)", na.value = "grey80", end = 0.8) +
        guides(alpha = FALSE) + facet_wrap(~feature_label)
      # g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + geom_point(aes(color=log10(value + 0.1)), size=I(cell_size), na.rm = TRUE, alpha = alpha) +
      #     viridis::scale_color_viridis(name = paste0("log10(value + 0.1)"), ...) + facet_wrap(~feature_label)
    }
  } else {
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2))

    # We don't want to force users to call order_cells before even being able to look at the trajectory,
    # so check whether it's null and if so, just don't color the cells
    if (color_cells_by == "Pseudotime" & is.null(data_df$Pseudotime)){
      g <- g + geom_point(color=I("gray"), size=I(cell_size), na.rm = TRUE, alpha = I(alpha))
      message("order_cells() has not been called yet, can't color cells by Pseudotime")
    } else if(color_cells_by %in% c("cluster", "partition")){
      if (is.null(data_df$cell_color)){
        g <- g + geom_point(color=I("gray"), size=I(cell_size), na.rm = TRUE, alpha = I(alpha))
        message("cluster_cells() has not been called yet, can't color cells by cluster")
      } else{
        g <- g + geom_point(aes(color = cell_color), size=I(cell_size), na.rm = TRUE, alpha = alpha)
      }
    } else if (class(data_df$cell_color) == "numeric"){
      g <- g + geom_point(aes(color = cell_color), size=I(cell_size), na.rm = TRUE, alpha = alpha)
      g <- g + viridis::scale_color_viridis(option="C")
    } else {
      g <- g + geom_point(aes(color = cell_color), size=I(cell_size), na.rm = TRUE, alpha = alpha)
    }
  }
  if (show_trajectory_graph){
    g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1",
                                     y="source_prin_graph_dim_2",
                                     xend="target_prin_graph_dim_1",
                                     yend="target_prin_graph_dim_2"),
                          size=trajectory_graph_segment_size,
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
        geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2", label="branch_point_idx"),
                  size=I(graph_label_size), color="white", na.rm=TRUE, branch_point_df)
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
        geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2", label="leaf_idx"),
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
        geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2", label="root_idx"),
                  size=I(graph_label_size), color="black", na.rm=TRUE, root_df)
    }
  }

  if(label_cell_groups) {
    g <- g + ggrepel::geom_text_repel(data = text_df,
                                      mapping = aes_string(x = "text_x", y = "text_y", label = "label"),
                                      size=I(group_label_size)) + theme(legend.position="none")
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
#' @description Plots expression for one or more genes as a function of pseudotime.
#' Plotting allows you determine if the ordering produced by orderCells() is correct
#' and it does not need to be flipped using the "reverse" flag in orderCells
#'
#' @param cds_subset cell_data_set for the experiment
#' @param min_expr the minimum (untransformed) expression level to use in plotted the genes.
#' @param cell_size the size (in points) of each cell used in the plot
#' @param nrow the number of rows used when laying out the panels for each gene's expression
#' @param ncol the number of columns used when laying out the panels for each gene's expression
#' @param panel_order the order in which genes should be layed out (left-to-right, top-to-bottom)
#' @param color_by the cell attribute (e.g. the column of colData(cds)) to be used to color each cell
#' @param trend_formula the model formula to be used for fitting the expression trend over pseudotime
#' @param label_by_short_name label figure panels by gene_short_name (TRUE) or feature id (FALSE)
#' @param vertical_jitter A value passed to ggplot to jitter the points in the vertical dimension. Prevents overplotting, and is particularly helpful for rounded transcript count data.
#' @param horizontal_jitter A value passed to ggplot to jitter the points in the horizontal dimension. Prevents overplotting, and is particularly helpful for rounded transcript count data.
#' @return a ggplot2 plot object
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' my_genes <- row.names(subset(rowData(HSMM), gene_short_name %in% c("CDK1", "MEF2C", "MYH3")))
#' cds_subset <- HSMM[my_genes,]
#' plot_genes_in_pseudotime(cds_subset, color_by="Time")
#' }
plot_genes_in_pseudotime <-function(cds_subset,
                                    min_expr=NULL,
                                    cell_size=0.75,
                                    nrow=NULL,
                                    ncol=1,
                                    panel_order=NULL,
                                    color_by="Pseudotime",
                                    trend_formula="~ splines::ns(Pseudotime, df=3)",
                                    label_by_short_name=TRUE,
                                    vertical_jitter=NULL,
                                    horizontal_jitter=NULL){

  f_id <- NA
  Cell <- NA
  cds_subset = cds_subset[,is.finite(colData(cds_subset)$Pseudotime)]

  cds_exprs <- counts(cds_subset)
  if (is.null(size_factors(cds_subset))) {
    stop("Error: to call this function, you must call estimate_size_factors() first")
  }
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
  #cds_exprs$f_id <- as.character(cds_exprs$f_id)
  #cds_exprs$Cell <- as.character(cds_exprs$Cell)

  cds_exprs$adjusted_expression <- cds_exprs$expression

  # trend_formula <- paste("adjusted_expression", trend_formula,
  #     sep = "")
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


  new_data <- data.frame(Pseudotime = colData(cds_subset)$Pseudotime)
  model_tbl = fit_models(cds_subset, model_formula_str = trend_formula)

  model_expectation <- model_predictions(model_tbl,
                                         new_data = colData(cds_subset))

  colnames(model_expectation) <- colnames(cds_subset)
  expectation <- plyr::ddply(cds_exprs, plyr::.(f_id, Cell), function(x) data.frame("expectation"=model_expectation[x$f_id, x$Cell]))
  cds_exprs <- merge(cds_exprs, expectation)
  #cds_exprs$expectation <- expectation#apply(cds_exprs,1, function(x) model_expectation[x$f_id, x$Cell])

  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
  if (is.null(panel_order) == FALSE) {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label,
                                      levels = panel_order)
  }
  q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
  if (is.null(color_by) == FALSE) {
    q <- q + geom_point(aes_string(color = color_by), size = I(cell_size), position=position_jitter(horizontal_jitter, vertical_jitter))
    if (class(colData(cds_subset)[,color_by]) == "numeric"){
      q <- q + viridis::scale_color_viridis(option="C")
    }
  }
  else {
    q <- q + geom_point(size = I(cell_size), position=position_jitter(horizontal_jitter, vertical_jitter))
  }

  q <- q + geom_line(aes(x = Pseudotime, y = expectation), data = cds_exprs)

  q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow,
                                        ncol = ncol, scales = "free_y")
  if (min_expr < 1) {
    q <- q + expand_limits(y = c(min_expr, 1))
  }

  q <- q + ylab("Expression")

  q <- q + xlab("Pseudo-time")
  q <- q + monocle_theme_opts()
  q
}

#' Plots the percentage of variance explained by the each component based on PCA from the normalized expression
#' data using the same procedure used in reduceDimension function.
#'
#' @param cds CellDataSet for the experiment after running reduceDimension with reduction_method as tSNE
#' @param max_components Maximum number of components shown in the scree plot (variance explained by each component)
#' @param norm_method Determines how to transform expression values prior to reducing dimensionality
#' @param residual_model_formula_str A model formula specifying the effects to subtract from the data before clustering.
#' @param pseudo_count amount to increase expression values before dimensionality reduction
#' @param return_all A logical argument to determine whether or not the variance of each component is returned
#' @param use_existing_pc_variance Whether to plot existing results for variance explained by each PC
#' @param verbose Whether to emit verbose output during dimensionality reduction
#' @param ... additional arguments to pass to the dimensionality reduction function
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' plot_pc_variance_explained(HSMM)
#' }
plot_pc_variance_explained <- function(cds,
                                       max_components=100,
                                       norm_method = c("log", "vstExprs", "none"),
                                       residual_model_formula_str=NULL,
                                       pseudo_count=NULL,
                                       return_all = F,
                                       use_existing_pc_variance=FALSE,
                                       verbose=FALSE,
                                       ...){
  norm_method <- match.arg(norm_method)
  set.seed(2016)
  if(!is.null(int_metadata(cds)$tsne_variance_explained) & use_existing_pc_variance == T){
    prop_varex <-int_metadata(cds)$tsne_variance_explained
  } else if(is.null(assays(cds)$normalized_exprs)) {
    stop(paste("You must call preprocess(cds) before running this function"))
  } else {
    FM <- assays(cds)$normalized_exprs
    if (!is.null(rowData(cds)$use_for_ordering) &&
        nrow(subset(rowData(cds), use_for_ordering == TRUE)) > 0) {
      FM <- FM[rowData(cds)$use_for_ordering, ]
    }

    xm <- Matrix::rowMeans(FM)
    xsd <- sqrt(Matrix::rowMeans((FM - xm)^2))
    FM <- FM[xsd > 0,]

    if (is.null(residual_model_formula_str) == FALSE) {
      if (verbose)
        message("Removing batch effects")
      X.model_mat <- Matrix::sparse.model.matrix(stats::as.formula(residual_model_formula_str),
                                                 data = colData(cds), drop.unused.levels = TRUE)

      fit <- limma::lmFit(FM, X.model_mat, ...)
      beta <- fit$coefficients[, -1, drop = FALSE]
      beta[is.na(beta)] <- 0
      FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
    } else {
      X.model_mat <- NULL
    }

    if (nrow(FM) == 0) {
      stop("Error: all rows have standard deviation zero")
    }

    irlba_res <- sparse_prcomp_irlba(Matrix::t(FM), n = min(max_components, min(dim(FM)) - 1),
                                     center = TRUE, scale. = TRUE)
    prop_varex <- irlba_res$sdev^2 / sum(irlba_res$sdev^2)

  }

  p <- qplot(1:length(prop_varex), prop_varex, alpha = I(0.5)) +  monocle_theme_opts() +
    theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
    theme(panel.background = element_rect(fill='white')) + xlab('components') +
    ylab('Variance explained \n by each component')

  int_metadata(cds)$tsne_variance_explained <- prop_varex # update CDS slot for variance_explained

  if(return_all) {
    return(list(variance_explained = prop_varex, p = p))
  }
  else
    return(p)
}

#' @title Plot expression for one or more genes as a violin plot
#'
#' @description Accepts a subset of a cell_data_set and an attribute to group
#' cells by, and produces a ggplot2 object that plots the level of expression
#' for each group of cells.
#'
#' @param cds_subset Subset cell_data_set to be plotted.
#' @param grouping the cell attribute (e.g. the column of colData(cds)) to
#'   group cells by on the horizontal axis.
#' @param min_expr the minimum (untransformed) expression level to use when
#'   plotted the genes. If \code{NULL},
#'   zero is used.
#' @param nrow the number of panels per row in the figure.
#' @param ncol the number of panels per column in the figure.
#' @param panel_order the order in which genes should be layed out
#'   (left-to-right, top-to-bottom). Should be gene_short_name, if
#'   \code{label_by_short_name = TRUE} or gene ID if
#'   \code{label_by_short_name = FALSE}.
#' @param label_by_short_name label figure panels by gene_short_name (TRUE) or
#'   feature id (FALSE).
#' @param log_scale Logical, whether or not to scale data logarithmically.
#' @return a ggplot2 plot object
#' @import ggplot2
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' my_genes <- HSMM[row.names(subset(rowData(HSMM),
#'                  gene_short_name %in% c("ACTA1", "ID1", "CCNB2"))),]
#' plot_genes_violin(my_genes, grouping="Hours", ncol=2, min_expr=0.1)
#' }
plot_genes_violin <- function (cds_subset,
                               grouping = "State",
                               min_expr = NULL,
                               nrow = NULL,
                               ncol = 1,
                               panel_order = NULL,
                               label_by_short_name = TRUE,
                               normalize = TRUE,
                               log_scale = TRUE) {

  assertthat::assert_that(is(cds_subset, "cell_data_set"))

  assertthat::assert_that(grouping %in% names(colData(cds_subset)),
                          msg = paste("grouping must be a column in the",
                                      "colData table"))
  if(!is.null(min_expr)) {
    assertthat::assert_that(assertthat::is.number(min_expr))
  }

  if(!is.null(nrow)) {
    assertthat::assert_that(assertthat::is.count(nrow))
  }

  assertthat::assert_that(assertthat::is.count(ncol))

  assertthat::assert_that(is.logical(label_by_short_name))
  assertthat::assert_that(is.logical(normalize))
  assertthat::assert_that(is.logical(log_scale))

  if (normalize) {
    cds_exprs <- counts(cds_subset)
    cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds_subset))
    cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
  } else {
    cds_exprs <- counts(cds_subset)
    cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
  }
  if (is.null(min_expr)) {
    min_expr <- 0
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

  cds_exprs[,grouping] <- as.factor(cds_exprs[,grouping])

  q <- ggplot(aes_string(x = grouping, y = "expression"), data = cds_exprs) +
    monocle_theme_opts()

  cds_exprs[,grouping] <- as.factor(cds_exprs[,grouping])
  q <- q + geom_violin(aes_string(fill = grouping), scale="width") + guides(fill=FALSE)

  q <- q + facet_wrap(~feature_label, nrow = nrow,
                      ncol = ncol, scales = "free_y")
  if (min_expr < 1) {
    q <- q + expand_limits(y = c(min_expr, 1))
  }

  q <- q + ylab("Expression") + xlab(grouping)

  if (log_scale){
    q <- q + scale_y_log10()
  }
  q
}


#' Plots the number of cells expressing one or more genes above a given value
#' as a barplot
#'
#'  @description Accepts a cell_data_set and the parameter "grouping", used for
#'  dividing cells into groups. Returns one or more bar graphs (one graph for
#'  each gene in the cell_data_set). Each graph shows the percentage of cells
#'  that express a gene in each sub-group in the cell_data_set.
#'
#'  As an example, let's say the cell_data_set passed into the function as
#'  cds_subset included genes A, B, and C and the grouping parameter divided
#'  the cells into three groups called X, Y, and Z. Then three graphs would be
#'  produced called A, B, and C. In each graph there would be three bars one
#'  for X, one for Y, and one for Z. The X bar in the A graph would show the
#'  percentage of cells in the X group that express gene A.
#'
#' @param cds_subset Subset cell_data_set to be plotted.
#' @param grouping the cell attribute (e.g. the column of colData(cds)) to
#'   group cells by on the horizontal axis
#' @param min_expr the minimum (untransformed) expression level to consider the
#'   gene 'expressed'. If \code{NULL},
#'  zero is used.
#' @param nrow the number of panels per row in the figure.
#' @param ncol the number of panels per column in the figure.
#' @param panel_order the order in which genes should be layed out
#'   (left-to-right, top-to-bottom)
#' @param plot_as_count whether to plot as a count of cells rather than a
#'   percent.
#' @param label_by_short_name label figure panels by gene_short_name (TRUE) or
#'   feature id (FALSE).
#' @param plot_limits A pair of number specifying the limits of the y axis. If
#'   \code{NULL}, scale to the range of the data. Example \code{c(0,100)}.
#' @return a ggplot2 plot object
#' @import ggplot2
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' MYOG_ID1 <- HSMM[row.names(subset(rowData(HSMM),
#'                                   gene_short_name %in% c("MYOG", "ID1"))),]
#' plot_percent_cells_positive(MYOG_ID1, grouping="Media", ncol=2)
#' }
plot_percent_cells_positive <- function(cds_subset,
                                        grouping = "State",
                                        min_expr = NULL,
                                        nrow = NULL,
                                        ncol = 1,
                                        panel_order = NULL,
                                        plot_as_count = FALSE,
                                        label_by_short_name=TRUE,
                                        normalize = TRUE,
                                        plot_limits = NULL){

  assertthat::assert_that(is(cds_subset, "cell_data_set"))

  assertthat::assert_that(grouping %in% names(colData(cds_subset)),
                          msg = paste("grouping must be a column in the",
                                      "colData table"))
  if(!is.null(min_expr)) {
    assertthat::assert_that(assertthat::is.number(min_expr))
  }

  if(!is.null(nrow)) {
    assertthat::assert_that(assertthat::is.count(nrow))
  }

  assertthat::assert_that(assertthat::is.count(ncol))
  assertthat::assert_that(is.logical(plot_as_count))
  assertthat::assert_that(is.logical(label_by_short_name))
  assertthat::assert_that(is.logical(normalize))

  if (is.null(min_expr)) {
    min_expr <- 0
  }

  marker_exprs <- counts(cds_subset)

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

  marker_counts <- plyr::ddply(marker_exprs_melted, c("feature_label",
                                                      grouping),
                               function(x) {
                                 data.frame(target = sum(x$expression > min_expr),
                                            target_fraction = sum(x$expression > min_expr)/nrow(x)) })

  if (!plot_as_count){
    marker_counts$target_fraction <- marker_counts$target_fraction * 100
    qp <- ggplot(aes_string(x=grouping, y="target_fraction", fill=grouping),
                 data=marker_counts) +
      ylab("Cells (percent)")
  } else {
    qp <- ggplot(aes_string(x=grouping, y="target", fill=grouping),
                 data=marker_counts) +
      ylab("Cells")
  }

  if (is.null(plot_limits) == FALSE) {
    qp <- qp + scale_y_continuous(limits=plot_limits)
  }

  qp <- qp + facet_wrap(~feature_label, nrow=nrow, ncol=ncol, scales="free_y")
  qp <-  qp + geom_bar(stat="identity") + monocle_theme_opts()

  return(qp)
}


#' Create a dot plot to visualize the mean gene expression and percentage of expressed cells in each group of cells
#'
#' @param cds CellDataSet for the experiment
#' @param markers a gene name use for visualize the dot plot
#' @param group_by the cell attribute (e.g. the column of pData(cds)) to group cells
#' @param lower_threshold The lowest gene expressed treated as expressed. By default, it is cds@lowerDetectionLimit.
#' @param max.size The maximum size of the dot. By default, it is 10.
#' @param ordering_type How to order the genes / groups on the dot plot. Only accept 'cluster_row_col' (use biclustering to cluster the rows and columns),
#' 'maximal_on_diag' (position each column so that the maximal color shown on each column on the diagonal, if the current maximal is used in earlier columns, the next largest one is position),
#' 'none' (preserve the ordering from the input gene or alphabetical ordering of groups)
#' @param axis_order Wheter to put groups on x-axis, genes on y-axis (option 'group_marker') or the reverse order (option 'marker_group')
#' @param flip_percentage_mean Whether to use color of the dot to represent the percentage (by setting flip_percentage_mean = FALSE, default) and size of the dot the mean expression or the opposite (by setting flip_percentage_mean = T)
#' @param pseudocount A pseudo-count added to the average gene expression
#' @param scale_max The maximum value (in standard deviations) to show in the heatmap. Values larger than this are set to the max.
#' @param scale_min The minimum value (in standard deviations) to show in the heatmap. Values smaller than this are set to the min.
#' @param ... additional arguments passed into the function (not used for now)
#' @return a ggplot2 plot object
#' @import ggplot2
#' @import pheatmap
#' @importFrom reshape2 melt
#' @importFrom reshape2 dcast
#' @importFrom viridis scale_color_viridis
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' HSMM <- reduceDimension(HSMM, reduction_method = 'tSNE')
#' HSMM <- clusterCells(HSMM)
#' plot_gene_by_group(HSMM, get_classic_muscle_markers())
#' }
plot_genes_by_group <- function(cds,
                                markers,
                                group_cells_by="cluster",
                                reduction_method = "UMAP",
                                norm_method = c("log", "size_only"),
                                lower_threshold = 0,
                                max.size = 10,
                                ordering_type = c('cluster_row_col', 'maximal_on_diag', 'none'), # maybe be also do the maximum color on the diagonal; the axis change be switched too
                                axis_order = c('group_marker', 'marker_group'),
                                flip_percentage_mean = FALSE,
                                pseudocount = 1,
                                scale_max = 3,
                                scale_min = -3,
                                ...) {
  #if(!(group_by %in% colnames(pData(cds))))
  #  stop(paste0(group_by, ' is not in the pData, please first perform cell clustering (clusterCells) before running this function!'))

  assertthat::assert_that(is(cds, "cell_data_set"))

  if(!is.null(group_cells_by)) {
    assertthat::assert_that(group_cells_by %in% c("cluster", "partition") | group_cells_by %in% names(colData(cds)),
                            msg = paste("group_cells_by must be a column in the",
                                        "colData table."))
  }

  norm_method = match.arg(norm_method)
  #group_cells_by=match.arg(group_cells_by)

  gene_ids = as.data.frame(fData(cds)) %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(rowname %in% markers | gene_short_name %in% markers) %>% dplyr::pull(rowname)
  if(length(gene_ids) < 1)
    stop('Please make sure markers are included in the gene_short_name column of the fData!')

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
    cell_group <- tryCatch({clusters(cds, reduction_method = reduction_method)}, error = function(e) {NULL})
  } else if (group_cells_by == "partition") {
    cell_group <- tryCatch({partitions(cds, reduction_method = reduction_method)}, error = function(e) {NULL})
  } else{
    cell_group <- colData(cds)[,group_cells_by]
  }
  names(cell_group) = colnames(cds)

  exprs_mat$Group <- cell_group[exprs_mat$Cell]
  exprs_mat = exprs_mat %>% filter(is.na(Group) == FALSE)
  ExpVal <- exprs_mat %>% dplyr::group_by(Group, Gene) %>% dplyr::summarize(mean = mean(log(Expression + pseudocount)), percentage = sum(Expression > lower_threshold) / length(Expression))
  ExpVal$mean <- ifelse(ExpVal$mean < scale_min, scale_min, ExpVal$mean)
  ExpVal$mean <- ifelse(ExpVal$mean > scale_max, scale_max, ExpVal$mean)

  ExpVal$Gene <- fData(cds)[ExpVal$Gene, 'gene_short_name']

  res <- reshape2::dcast(ExpVal[, 1:4], Group ~ Gene, value.var = colnames(ExpVal)[2 + major_axis])
  group_id <- res[, 1]
  res <- res[, -1]
  row.names(res) <- group_id

  if(ordering_type == 'cluster_row_col') {
    row_dist <- as.dist((1 - cor(t(res[, -1])))/2)
    row_dist[is.na(row_dist)] <- 1

    col_dist <- as.dist((1 - cor(res[, -1]))/2)
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

    ExpVal$Gene <- factor(ExpVal$Gene, levels = colnames(res)[-1][ph$tree_col$order])
    ExpVal$Group <- factor(ExpVal$Group, levels = row.names(res)[ph$tree_row$order])

  } else if(ordering_type == 'maximal_on_diag'){

    order_mat <- t(apply(res, major_axis, order))
    max_ind_vec <- c()
    for(i in 1:nrow(order_mat)) {
      tmp <- max(which(!(order_mat[i, ] %in% max_ind_vec)))
      max_ind_vec <- c(max_ind_vec, order_mat[i, tmp])
    }
    max_ind_vec <- max_ind_vec[!is.na(max_ind_vec)]

    if(major_axis == 1){
      # ExpVal$Gene <- factor(ExpVal %>% dplyr::pull(colnames(ExpVal)[minor_axis]) , levels = dimnames(res)[[minor_axis]][max_ind_vec])
      max_ind_vec <- c(max_ind_vec, setdiff(1:length(markers), max_ind_vec))
      ExpVal$Gene <- factor(ExpVal$Gene , levels = dimnames(res)[[2]][max_ind_vec])
    }
    else{
      max_ind_vec <- c(max_ind_vec, setdiff(1:length(unique(exprs_mat$Group)), max_ind_vec))
      ExpVal$Group <- factor(ExpVal$Group , levels = dimnames(res)[[1]][max_ind_vec])
    }
  } else if(ordering_type == 'none'){
    ExpVal$Gene <- factor(ExpVal$Gene, levels = markers)
  }

  # + scale_color_gradient(low ="blue",   high = "red", limits=c(1, max(ExpVal$mean) ))
  if(flip_percentage_mean){
    g <- ggplot(ExpVal, aes(y = Gene,  x = Group)) + geom_point(aes(colour = percentage,  size = mean)) + viridis::scale_color_viridis(name = 'percentage') + scale_size(name = 'log(mean + 0.1)', range = c(0, max.size))
  } else {
    g <- ggplot(ExpVal, aes(y = Gene,  x = Group)) + geom_point(aes(colour = mean,  size = percentage)) + viridis::scale_color_viridis(name = 'log(mean + 0.1)') + scale_size(name = 'percentage', range = c(0, max.size))
  }

  g <- g + xlab("Cluster") + ylab("Gene") + monocle3:::monocle_theme_opts() + xlab("Cluster") + ylab("Gene") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  if(axis_order == 'marker_group') {
    g <- g + coord_flip()
  }

  g
}

