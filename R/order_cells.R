
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
                       root_cells=NULL,
                       reverse = FALSE,
                       orthogonal_proj_tip = FALSE,
                       verbose = FALSE){

  if(class(cds)[1] != "cell_data_set") {
    stop("Error: cds is not of type 'cell_data_set'")
  }

  if (is.null(cds@dim_reduce_type)){
    stop("Error: dimensionality not yet reduced. Please call reduce_dimension() and learnGraph() (for learning principal graph) before calling this function.")
  }
  if (is.null(cds@rge_method)){
    stop("Error: principal graph has not learned yet. Please call learnGraph() before calling this function.")
  }
  # reducedDimA, S, and K are not NULL in the cds
  if (length(cds@reducedDimS) == 0) {
    stop("Error: dimension reduction didn't prodvide correct results. Please check your reduce_dimension() step and ensure correct dimension reduction are performed before calling this function.")
  }
  if (length(cds@reducedDimK) == 0) {
    stop("Error: principal graph learning didn't prodvide correct results. Please check your learnGraph() step and ensure correct principal graph learning are performed before calling this function.")
  }
  if(igraph::vcount(principal_graph(cds)) > 10000) {
    stop("order_cells doesn't support more than 10k centroids (cells)")
  }
  if (is.null(root_pr_nodes) == FALSE & is.null(root_cells) == FALSE){
    stop("Error: please specify either root_pr_nodes or root_cells, not both")
  }
  if (is.null(root_pr_nodes) & is.null(root_cells)){
    if (interactive()){
      root_pr_nodes <- selectTrajectoryRoots(cds)
    }else{
      stop("Error: You must provide one or more root cells (or principal graph nodes) in non-interactive mode")
    }
  }else if(is.null(root_pr_nodes)){
    valid_root_cells <- intersect(root_cells, row.names(pData(cds)))
    if (length(valid_root_cells) == 0){
      stop("Error: no such cell")
    }
    closest_vertex = cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
    root_pr_nodes = closest_vertex[valid_root_cells,]
  }else{
    if (length(intersect(root_pr_nodes, igraph::V(principal_graph(cds))$name)) == 0){
      stop("Error: no such principal graph node")
    }
  }

  if (is.null(root_pr_nodes) || length(root_pr_nodes) == 0){
    stop("Error: no valid root principal graph nodes.")
  }


  cds@aux_ordering_data[[cds@rge_method]]$root_pr_nodes <- root_pr_nodes

  cc_ordering <- extract_general_graph_ordering(cds, root_pr_nodes, orthogonal_proj_tip, verbose)
  pData(cds)$Pseudotime = cc_ordering[row.names(pData(cds)), ]$pseudo_time
  if(reverse) {
    finite_cells <- is.finite(pData(cds)$Pseudotime)
    pData(cds)$Pseudotime[finite_cells] <- max(pData(cds)$Pseudotime[finite_cells]) - pData(cds)$Pseudotime[finite_cells]
  }

  cds
}

extract_general_graph_ordering <- function(cds, root_cell, orthogonal_proj_tip = FALSE, verbose=T)
{
  Z <- reducedDimS(cds)
  Y <- reducedDimK(cds)
  pr_graph <- principal_graph(cds)

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
  closest_vertex <- findNearestVertex(Y[, root_cell, drop = F], Z)
  closest_vertex_id <- colnames(cds)[closest_vertex]

  cell_wise_graph <- cds@aux_ordering_data[[cds@rge_method]]$pr_graph_cell_proj_tree
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


#' functions to connect the tip points after learning the DDRTree or simplePPT tree
connectTips <- function(pd,
                        R, # kmean cluster
                        stree,
                        reducedDimK_old,
                        reducedDimS_old,
                        k = 25,
                        weight = F,
                        qval_thresh = 0.05,
                        kmean_res,
                        euclidean_distance_ratio = 1,
                        geodestic_distance_ratio = 1/3,
                        medioids,
                        verbose = FALSE,
                        ...) {
  if(is.null(row.names(stree)) & is.null(row.names(stree))) {
    dimnames(stree) <- list(paste0('Y_', 1:ncol(stree)), paste0('Y_', 1:ncol(stree)))
  }

  stree <- as.matrix(stree)
  stree[stree != 0] <- 1
  mst_g_old <- igraph::graph_from_adjacency_matrix(stree, mode = 'undirected')
  if(is.null(kmean_res)) {
    tmp <- matrix(apply(R, 1, which.max))

    row.names(tmp) <- colnames(reducedDimS_old)

    tip_pc_points <- which(igraph::degree(mst_g_old) == 1)

    data <- t(reducedDimS_old[, ])

    louvain_res <- louvain_clustering(data, pd[, ], k = k, weight = weight, verbose = verbose)

    # louvain_res$optim_res$memberships[1, ] <-  tmp[raw_data_tip_pc_points, 1]
    louvain_res$optim_res$membership <- tmp[, 1]
  } else { # use kmean clustering result
    tip_pc_points <- which(igraph::degree(mst_g_old) == 1)
    tip_pc_points_kmean_clusters <- sort(kmean_res$cluster[names(tip_pc_points)])
    # raw_data_tip_pc_points <- which(kmean_res$cluster %in% tip_pc_points_kmean_clusters) # raw_data_tip_pc_points <- which(kmean_res$cluster %in% tip_pc_points)

    data <- t(reducedDimS_old[, ]) # raw_data_tip_pc_points

    louvain_res <- louvain_clustering(data, pd[row.names(data), ], k = k, weight = weight, verbose = verbose)

    # louvain_res$optim_res$memberships[4, ] <-  kmean_res$cluster #[raw_data_tip_pc_points]
    louvain_res$optim_res$membership <- kmean_res$cluster #[raw_data_tip_pc_points]
  }

  # identify edges between only tip cells
  cluster_graph_res <- compute_louvain_connected_components(louvain_res$g, louvain_res$optim_res, qval_thresh=qval_thresh, verbose = verbose)
  dimnames(cluster_graph_res$cluster_mat) <- dimnames(cluster_graph_res$num_links)
  valid_connection <- which(cluster_graph_res$cluster_mat < qval_thresh, arr.ind = T)
  valid_connection <- valid_connection[apply(valid_connection, 1, function(x) all(x %in% tip_pc_points)), ] # only the tip cells

  # prepare the PAGA graph
  G <- cluster_graph_res$cluster_mat
  G[cluster_graph_res$cluster_mat < qval_thresh] <- -1
  G[cluster_graph_res$cluster_mat > 0] <- 0
  G <- - G

  if(all(G == 0, na.rm = T)) { # if no connection based on PAGA (only existed in simulated data), use the kNN graph instead
    return(list(stree = igraph::get.adjacency(mst_g_old), Y = reducedDimK_old, G = G))
    #G <- get_knn(medioids, K = min(5, ncol(medioids)))$G
    #if(!is.null(kmean_res)) {
    #  valid_connection <- which(G > 0, arr.ind = T)
    #  valid_connection <- valid_connection[apply(valid_connection, 1, function(x) all(x %in% tip_pc_points)), ] # only the tip cells
    #}
  }

  if(nrow(valid_connection) == 0) {
    return(list(stree = igraph::get.adjacency(mst_g_old), Y = reducedDimK_old, G = G))
  }

  # calculate length of the MST diameter path
  mst_g <- mst_g_old
  diameter_dis <- igraph::diameter(mst_g_old)
  reducedDimK_df <- reducedDimK_old

  pb4 <- txtProgressBar(max = length(nrow(valid_connection)), file = "", style = 3, min = 0)

  # find the maximum distance between nodes from the MST
  res <- dist(t(reducedDimK_old))
  g <- igraph::graph_from_adjacency_matrix(as.matrix(res), weighted = T, mode = 'undirected')
  mst <- igraph::minimum.spanning.tree(g)
  max_node_dist <- max(igraph::E(mst)$weight)

  # append new edges to close loops in the spanning tree returned from SimplePPT
  for(i in 1:nrow(valid_connection)) {
    # cluster id for the tip point; if kmean_res return valid_connection[i, ] is itself; otherwise the id identified in the tmp file
    edge_vec <- sort(unique(louvain_res$optim_res$membership))[valid_connection[i, ]]
    edge_vec_in_tip_pc_point <- igraph::V(mst_g_old)$name[edge_vec]

    if(length(edge_vec_in_tip_pc_point) == 1) next;
    if(all(edge_vec %in% tip_pc_points) & (igraph::distances(mst_g_old, edge_vec_in_tip_pc_point[1], edge_vec_in_tip_pc_point[2]) >= geodestic_distance_ratio * diameter_dis) &
       (euclidean_distance_ratio * max_node_dist > dist(t(reducedDimK_old[, edge_vec]))) ) {
      if(verbose) message('edge_vec is ', edge_vec[1], '\t', edge_vec[2])
      if(verbose) message('edge_vec_in_tip_pc_point is ', edge_vec_in_tip_pc_point[1], '\t', edge_vec_in_tip_pc_point[2])

      mst_g <- igraph::add_edges(mst_g, edge_vec_in_tip_pc_point)
    }
    setTxtProgressBar(pb = pb4, value = pb4$getVal() + 1)
  }

  close(pb4)

  list(stree = igraph::get.adjacency(mst_g), Y = reducedDimK_df, G = G)
}





