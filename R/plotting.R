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
#' @param cds CellDataSet for the experiment
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
    stop("Error: dimensionality not reduced or graph is not learned yet. Please call reduceDimension(), partitionCells() and learnGraph() before calling this function.")
  }


  reduced_dim_coords <- reducedDimK(cds)

  ica_space_df <- Matrix::t(reduced_dim_coords) %>%
    as.data.frame() %>%
    dplyr::select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
    dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))

  dp_mst <- principal_graph(cds)

  if (is.null(dp_mst)){
    stop("You must first call orderCells() before using this function")
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
