
#' Orders cells according to pseudotime.
#'
#' Learns a "trajectory" describing the biological process the cells are
#' going through, and calculates where each cell falls within that trajectory.
#' Monocle learns trajectories in two steps. The first step is reducing the dimensionality
#' of the data with \code{\link{reduce_dimension}()}. The second is this function.
#' function. This function takes as input a cell_data_set and returns it with
#' two new columns: \code{Pseudotime} and \code{State}, which together encode
#' where each cell maps to the trajectory. \code{order_cells()} optionally takes
#' a "root" state, which you can use to specify the start of the trajectory. If
#' you don't provide a root state, one is selected arbitrarily.
#'
#' The \code{reduction_method} argument to \code{\link{reduce_dimension}()}
#' determines which algorithm is used by \code{order_cells()} to learn the trajectory.
#' If \code{reduction_method == "ICA"}, this function uses \emph{polygonal reconstruction}
#' to learn the underlying trajectory. If \code{reduction_method == "DDRTree"},
#' the trajectory is specified by the principal graph learned by the
#' \code{\link[DDRTree]{DDRTree}()} function.
#'
#' Whichever algorithm you use, the trajectory will be composed of segments.
#' The cells from a segment will share the same value of \code{State}. One of
#' these segments will be selected as the root of the trajectory arbitrarily.
#' The most distal cell on that segment will be chosen as the "first" cell in the
#' trajectory, and will have a Pseudotime value of zero. \code{order_cells()} will
#' then "walk" along the trajectory, and as it encounters additional cells, it
#' will assign them increasingly large values of Pseudotime.
#'
#' @param cds the cell_data_set upon which to perform this operation
#' @param root_pr_nodes The starting principal points. We learn a principal graph that passes through the middle of the data points and use it to represent the developmental process.
#' @param root_cells The starting cells. Each cell corresponds to a principal point and multiple cells can correspond to the same principal point.
#' @param reverse Whether to reverse the direction of the trajectory
#' @param orthogonal_proj_tip Whether to perform orthogonal projection for cells corresponding to the tip principal points. Default to be FALSE
#' @param verbose Whether to show running information for order_cells
#'
#' @return an updated cell_data_set object, in which phenoData contains values for State and Pseudotime for each cell
#' @export
order_cells <- function(cds,
                       root_pr_nodes=NULL,
                       reduced_dimension = "UMAP",
                       root_cells=NULL,
                       reverse = FALSE,
                       orthogonal_proj_tip = FALSE,
                       verbose = FALSE){

  assertthat::assert_that(is(cds, "cell_data_set"))
  # if (is.null(cds@dim_reduce_type)){
  #  stop("Error: dimensionality not yet reduced. Please call reduce_dimension() and learnGraph() (for learning principal graph) before calling this function.")
  #}
  #if (is.null(cds@rge_method)){
  #  stop("Error: principal graph has not learned yet. Please call learnGraph() before calling this function.")
  #}
  # reducedDimA, S, and K are not NULL in the cds
  #if (length(cds@reducedDimS) == 0) {
  #  stop("Error: dimension reduction didn't prodvide correct results. Please check your reduce_dimension() step and ensure correct dimension reduction are performed before calling this function.")
  #}
  #if (length(cds@reducedDimK) == 0) {
  #  stop("Error: principal graph learning didn't prodvide correct results. Please check your learnGraph() step and ensure correct principal graph learning are performed before calling this function.")
  #}
  #if(igraph::vcount(principal_graph(cds)) > 10000) {
  #  stop("order_cells doesn't support more than 10k centroids (cells)")
  #}
  #if (is.null(root_pr_nodes) == FALSE & is.null(root_cells) == FALSE){
  #  stop("Error: please specify either root_pr_nodes or root_cells, not both")
  #}
  if (is.null(root_pr_nodes) & is.null(root_cells)){
    if (interactive()){
      root_pr_nodes <- selectTrajectoryRoots(cds)
    }else{
      stop("Error: You must provide one or more root cells (or principal graph nodes) in non-interactive mode")
    }
  }else if(is.null(root_pr_nodes)){
    valid_root_cells <- intersect(root_cells, row.names(colData(cds)))
    if (length(valid_root_cells) == 0){
      stop("Error: no such cell")
    }
    closest_vertex = cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
    root_pr_nodes = closest_vertex[valid_root_cells,]
  }else{
    if (length(intersect(root_pr_nodes, igraph::V(principal_graph(cds)[[reduced_dimension]])$name)) == 0){
      stop("Error: no such principal graph node")
    }
  }

  if (is.null(root_pr_nodes) || length(root_pr_nodes) == 0){
    stop("Error: no valid root principal graph nodes.")
  }


  cds@principal_graph_aux[[reduced_dimension]]$root_pr_nodes <- root_pr_nodes

  cc_ordering <- extract_general_graph_ordering(cds, root_pr_nodes, orthogonal_proj_tip, verbose, reduced_dimension)
  colData(cds)$Pseudotime = cc_ordering[row.names(colData(cds)), ]$pseudo_time
  if(reverse) {
    finite_cells <- is.finite(colData(cds)$Pseudotime)
    colData(cds)$Pseudotime[finite_cells] <- max(colData(cds)$Pseudotime[finite_cells]) - colData(cds)$Pseudotime[finite_cells]
  }

  cds
}

extract_general_graph_ordering <- function(cds, root_cell, orthogonal_proj_tip = FALSE, verbose=T, reduced_dimension)
{
  Z <- t(reducedDims(cds)$UMAP)
  pr_graph <- principal_graph(cds)[[reduced_dimension]]

  res <- list(subtree = pr_graph, root = root_cell)

  parents <- rep(NA, length(igraph::V(pr_graph)))
  states <- rep(NA, length(igraph::V(pr_graph)))

  if(any(is.na(igraph::E(pr_graph)$weight))) {
    igraph::E(pr_graph)$weight <- 1
  }

  # do pseudotime calculation on the cell-wise graph
  # 1. identify nearest cells to the selected principal node
  # 2. build a cell-wise graph for each louvain group
  # 3. run the distance function to assign pseudotime for each cell
  closest_vertex <- findNearestVertex(Z[, root_cell, drop = F], Z)
  closest_vertex_id <- colnames(cds)[closest_vertex]

  cell_wise_graph <- cds@principal_graph_aux[[reduced_dimension]]$pr_graph_cell_proj_tree
  cell_wise_distances <- igraph::distances(cell_wise_graph, v = closest_vertex_id)

  if (length(closest_vertex_id) > 1){
    node_names <- colnames(cell_wise_distances)
    pseudotimes <- apply(cell_wise_distances, 2, min)
  }else{
    node_names <- names(cell_wise_distances)
    pseudotimes <- cell_wise_distances
  }

  names(pseudotimes) <- node_names

  ordering_df <- data.frame(sample_name = igraph::V(cell_wise_graph)$name, # pr_graph
                            # cell_state = states,
                            pseudo_time = as.vector(pseudotimes)
                            # parent = parents
  )
  row.names(ordering_df) <- ordering_df$sample_name
  return(ordering_df)
}





