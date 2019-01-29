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

#' Plots the minimum spanning tree on cells.
#'
#' @param cds cell_data_set for the experiment
#' @param x the column of reducedDimS(cds) to plot on the horizontal axis
#' @param y the column of reducedDimS(cds) to plot on the vertical axis
#' @param color_by the cell attribute (e.g. the column of pData(cds)) to map to each cell's color
#' @param show_backbone whether to show the diameter path of the MST used to order the cells
#' @param backbone_color the color used to render the backbone.
#' @param markers a gene name or gene id to use for setting the size of each cell in the plot
#' @param use_color_gradient Whether or not to use color gradient instead of cell size to show marker expression level
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
#' @import ggplot2
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung()
#' plot_cell_trajectory(lung)
#' plot_cell_trajectory(lung, color_by="Pseudotime", show_backbone=FALSE)
#' plot_cell_trajectory(lung, markers="MYH3")
#' }
plot_cell_trajectory <- function(cds,
                                 x=1,
                                 y=2,
                                 color_by="Pseudotime",
                                 show_backbone=TRUE,
                                 backbone_color="black",
                                 markers=NULL,
                                 use_color_gradient = FALSE,
                                 markers_linear = FALSE,
                                 show_cell_names=FALSE,
                                 show_state_number = FALSE,
                                 cell_size=1.5,
                                 cell_link_size=0.75,
                                 cell_name_size=2,
                                 state_number_size = 2.9,
                                 show_branch_points=TRUE,
                                 theta = 0,
                                 alpha = 1,
                                 ...) {
  requireNamespace("igraph")
  gene_short_name <- NA
  sample_name <- NA
  sample_state <- pData(cds)$State
  data_dim_1 <- NA
  data_dim_2 <- NA

  #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
  lib_info_with_pseudo <- pData(cds)

  if (is.null(cds@dim_reduce_type) | is.null(cds@rge_method)){
    stop("Error: dimensionality not reduced or graph is not learned yet. Please call reduce_dimension(), partition_cells() and learn_graph() before calling this function.")
  }


  reduced_dim_coords <- reducedDimK(cds)

  ica_space_df <- Matrix::t(reduced_dim_coords) %>%
    as.data.frame() %>%
    dplyr::select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
    dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))

  dp_mst <- principal_graph(cds)

  if (is.null(dp_mst)){
    stop("You must first call order_cells() before using this function")
  }

  edge_df <- dp_mst %>%
    igraph::as_data_frame() %>%
    dplyr::select_(source = "from", target = "to") %>%
    dplyr::left_join(ica_space_df %>% dplyr::select_(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", source_prin_graph_dim_2="prin_graph_dim_2"), by = "source") %>%
    dplyr::left_join(ica_space_df %>% dplyr::select_(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2"), by = "target")

  S_matrix <-
    reducedDimS(cds)
  data_df <- data.frame(t(S_matrix[c(x,y),]))
  #data_df <- cbind(data_df, sample_state)
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)
  data_df <- merge(data_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")

  return_rotation_mat <- function(theta) {
    theta <- theta / 180 * pi
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
  }
  rot_mat <- return_rotation_mat(theta)

  cn1 <- c("data_dim_1", "data_dim_2")
  cn2 <- c("source_prin_graph_dim_1", "source_prin_graph_dim_2")
  cn3 <- c("target_prin_graph_dim_1", "target_prin_graph_dim_2")
  data_df[, cn1] <- as.matrix(data_df[, cn1]) %*% t(rot_mat)
  edge_df[, cn2] <- as.matrix(edge_df[, cn2]) %*% t(rot_mat)
  edge_df[, cn3] <- as.matrix(edge_df[, cn3]) %*% t(rot_mat)

  markers_exprs <- NULL
  if (is.null(markers) == FALSE) {
    markers_fData <- subset(fData(cds), gene_short_name %in% markers)
    if (nrow(markers_fData) >= 1) {
      markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData),])))
      colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
      markers_exprs <- merge(markers_exprs, markers_fData, by.x = "feature_id", by.y="row.names")
      #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    data_df <- merge(data_df, markers_exprs, by.x="sample_name", by.y="cell_id")
    if(use_color_gradient) {
      if(markers_linear){
        g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + geom_point(aes(color= value), size=I(cell_size), na.rm = TRUE, alpha = alpha) +
          viridis::scale_color_viridis(name = paste0("value"), ...) + facet_wrap(~feature_label)
      } else {
        g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + geom_point(aes(color=log10(value + 0.1)), size=I(cell_size), na.rm = TRUE, alpha = alpha) +
          viridis::scale_color_viridis(name = paste0("log10(value + 0.1)"), ...) + facet_wrap(~feature_label)
      }
    } else {
      if(markers_linear){
        g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2, size= (value * 0.1))) + facet_wrap(~feature_label)
      } else {
        g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2, size=log10(value + 0.1))) + facet_wrap(~feature_label)
      }
    }
  } else {
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2))
  }
  if (show_backbone){
    g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), size=cell_link_size, linetype="solid", na.rm=TRUE, data=edge_df)
  }

  # FIXME: setting size here overrides the marker expression funtionality.
  # Don't do it!
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    if(use_color_gradient) {
      # g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE)
    } else {
      g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE, alpha = alpha)
    }
  }else {
    if(use_color_gradient) {
      # g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE)
    } else {
      g <- g + geom_point(aes_string(color = color_by), size=I(cell_size), na.rm = TRUE, alpha = alpha)
      if (class(pData(cds)[,color_by]) == "numeric"){
        g <- g + viridis::scale_color_viridis(option="C")
      }
    }
  }


  if (show_branch_points && cds@rge_method %in% c('DDRTree', 'SimplePPT', 'L1graph')){
    mst_branch_nodes <- cds@aux_ordering_data[[cds@rge_method]]$branch_points
    branch_point_df <- ica_space_df %>%
      dplyr::slice(match(mst_branch_nodes, sample_name)) %>%
      dplyr::mutate(branch_point_idx = seq_len(n()))

    g <- g +
      geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                 size=5, na.rm=TRUE, branch_point_df, alpha = alpha) +
      geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2", label="branch_point_idx"),
                size=4, color="white", na.rm=TRUE, branch_point_df)
  }
  if (show_cell_names){
    g <- g + geom_text(aes(label=sample_name), size=cell_name_size)
  }
  if (show_state_number){
    g <- g + geom_text(aes(label = sample_state), size = state_number_size)
  }

  g <- g +
    #scale_color_brewer(palette="Set1") +
    monocle_theme_opts() +
    xlab(paste("Component", x)) +
    ylab(paste("Component", y)) +
    theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
    #guides(color = guide_legend(label.position = "top")) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white'))
  g
}





#' Plot a dataset and trajectory in 3 dimensions
#'
#' @param cds cell_data_set for the experiment
#' @param dim the dimensions used to create the 3D plot, by default it is the first three dimensions
#' @param color_by the cell attribute (e.g. the column of pData(cds)) to map to each cell's color
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
#' plot_3d_cell_trajectory(cds, markers=c("Rbfox3, Neurod1", "Sox2"))
#' }
#'
plot_3d_cell_trajectory <- function(cds,
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

  lib_info_with_pseudo <- pData(cds)[cell_sampled, ]
  sample_state <- pData(cds)$State[cell_sampled]

  if (is.null(cds@dim_reduce_type)){
    stop("Error: dimensionality not yet reduced. Please call reduce_dimension() before calling this function.")
  }

  reduced_dim_coords <- reducedDimK(cds)

  if(length(dim) != 3)
    dim <- 1:3

  ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[dim,]))
  colnames(ica_space_df) <- c("prin_graph_dim_1", "prin_graph_dim_2",  "prin_graph_dim_3")

  ica_space_df$sample_name <- row.names(ica_space_df)
  ica_space_df$sample_state <- row.names(ica_space_df)
  #ica_space_with_state_df <- merge(ica_space_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")
  #print(ica_space_with_state_df)
  dp_mst <- principal_graph(cds)

  if (is.null(dp_mst)){
    stop("You must first call orderCells() before using this function")
  }

  edge_list <- as.data.frame(igraph::get.edgelist(dp_mst))
  colnames(edge_list) <- c("source", "target")

  edge_df <- merge(edge_list, ica_space_df, by.x="source", by.y="sample_name", all=F)
  edge_df <- plyr::rename(edge_df, c("prin_graph_dim_1"="source_prin_graph_dim_1", "prin_graph_dim_2"="source_prin_graph_dim_2", "prin_graph_dim_3"="source_prin_graph_dim_3"))
  edge_df <- merge(edge_df, ica_space_df[,c("sample_name", "prin_graph_dim_1", "prin_graph_dim_2", "prin_graph_dim_3")], by.x="target", by.y="sample_name", all=F)
  edge_df <- plyr::rename(edge_df, c("prin_graph_dim_1"="target_prin_graph_dim_1", "prin_graph_dim_2"="target_prin_graph_dim_2", "prin_graph_dim_3"="target_prin_graph_dim_3"))

  S_matrix <- reducedDimS(cds)[, cell_sampled]
  data_df <- data.frame(t(S_matrix[dim,]))
  #data_df <- cbind(data_df, sample_state)
  colnames(data_df) <- c("data_dim_1", "data_dim_2", "data_dim_3")
  data_df$sample_name <- row.names(data_df)
  data_df <- merge(data_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")

  markers_exprs = NULL
  if (is.null(markers) == FALSE){
    markers_fData <- subset(fData(cds), gene_short_name %in% markers)
    if (nrow(markers_fData) >= 1){

      markers_expr_val <- exprs(cds[row.names(markers_fData), cell_sampled])
      markers_expr_val <- Matrix::colSums(markers_expr_val)
      markers_expr_val <- markers_expr_val / pData(cds)$Size_Factor
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
  }else if (is.null(color_by) == FALSE && color_by %in% colnames(pData(cds))){
    point_colors_df = dplyr::data_frame(sample_name=data_df$sample_name,
                                        color_by=as.character(pData(cds)[data_df$sample_name,color_by]))

    if (class(pData(cds)[,color_by]) == "numeric"){
      point_colors_df$point_colors = map2color(pData(cds)[data_df$sample_name,color_by])
      point_colors_df$point_colors[is.na(point_colors_df$point_colors)] = "darkgray"
    }else{
      gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
      }

      if (is.null(palette)){

        if(is.factor(pData(cds)[,color_by])) {
          num_colors = length(levels(pData(cds)[,color_by]))
          point_colors = gg_color_hue(num_colors)
          names(point_colors) = levels(pData(cds)[,color_by])
        } else {
          num_colors = length(unique(pData(cds)[,color_by]))
          point_colors = gg_color_hue(num_colors)
          names(point_colors) = unique(pData(cds)[,color_by])
        }
      }else{
        point_colors = palette
      }

      point_colors_df$point_colors = point_colors[point_colors_df$color_by]
      #point_colors = colors[as.character(pData(cds)[data_df$sample_name,color_by])]
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
    if (show_group_labels && color_by %in% colnames(pData(cds)) && class(pData(cds)[,color_by]) != "numeric"){
      if(!use_plotly) {
        rgl::text3d(x=medoid_df$mean_d1, y=medoid_df$mean_d2, z=medoid_df$mean_d3, texts=as.character(medoid_df$color_by))
      }
    }
  }
  branch_point_df <- NULL
  if (show_branch_points && cds@rge_method %in% c('DDRTree', 'SimplePPT', 'L1graph')){
    mst_branch_nodes <- cds@aux_ordering_data[[cds@rge_method]]$branch_points
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

