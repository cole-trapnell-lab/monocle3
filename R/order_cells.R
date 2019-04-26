
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
                       reduction_method = "UMAP",
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
      root_pr_nodes <- select_trajectory_roots(cds, reduction_method = reduction_method)
    }else{
      stop("Error: You must provide one or more root cells (or principal graph nodes) in non-interactive mode")
    }
  }else if(is.null(root_pr_nodes)){
    valid_root_cells <- intersect(root_cells, row.names(colData(cds)))
    if (length(valid_root_cells) == 0){
      stop("Error: no such cell")
    }
    closest_vertex = cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex
    root_pr_nodes = closest_vertex[valid_root_cells,]
  }else{
    if (length(intersect(root_pr_nodes, igraph::V(principal_graph(cds)[[reduction_method]])$name)) == 0){
      stop("Error: no such principal graph node")
    }
  }

  if (is.null(root_pr_nodes) || length(root_pr_nodes) == 0){
    stop("Error: no valid root principal graph nodes.")
  }


  cds@principal_graph_aux[[reduction_method]]$root_pr_nodes <- root_pr_nodes

  cc_ordering <- extract_general_graph_ordering(cds, root_pr_nodes, orthogonal_proj_tip, verbose, reduction_method)
  colData(cds)$Pseudotime = cc_ordering[row.names(colData(cds)), ]$pseudo_time
  if(reverse) {
    finite_cells <- is.finite(colData(cds)$Pseudotime)
    colData(cds)$Pseudotime[finite_cells] <- max(colData(cds)$Pseudotime[finite_cells]) - colData(cds)$Pseudotime[finite_cells]
  }

  cds
}

extract_general_graph_ordering <- function(cds, root_cell, orthogonal_proj_tip = FALSE, verbose=T, reduction_method)
{
  Z <- t(reducedDims(cds)[[reduction_method]])
  Y <- cds@principal_graph_aux[[reduction_method]]$dp_mst
  pr_graph <- principal_graph(cds)[[reduction_method]]

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

  cell_wise_graph <- cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_tree
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


#' Select the roots of the principal graph
#' @param cds CellDataSet where roots will be selected from
#' @param x The first dimension to plot
#' @param y The number of dimension to plot
#' @param num_roots Number of roots for the trajectory
#' @param pch Size of the principal graph node
#' @param ... Extra arguments to pass to function
select_trajectory_roots <- function(cds, x=1, y=2,
                                    num_roots = NULL,
                                    pch = 19,
                                    reduction_method,
                                    ...)
{
  #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
  #lib_info_with_pseudo <- colData(cds)

  reduced_dim_coords <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst)

  ica_space_df <- as.data.frame(reduced_dim_coords)
  use_3d <- ncol(ica_space_df) >= 3
  if (use_3d){
    colnames(ica_space_df) = c("prin_graph_dim_1", "prin_graph_dim_2", "prin_graph_dim_3")
  }
  else{
    colnames(ica_space_df) = c("prin_graph_dim_1", "prin_graph_dim_2")
  }

  ica_space_df$sample_name <- row.names(ica_space_df)
  ica_space_df$sample_state <- row.names(ica_space_df)

  dp_mst <- principal_graph(cds)[[reduction_method]]

  if (is.null(dp_mst)){
    stop("You must first call orderCells() before using this function")
  }

  if (use_3d){
    edge_df <- dp_mst %>%
      igraph::as_data_frame() %>%
      select_(source = "from", target = "to") %>%
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
  }else{
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

  if (is.null(num_roots)){
    num_roots = nrow(ica_space_df)
  }
  #xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
  sel <- rep(FALSE, nrow(ica_space_df))

  if (use_3d){
    open3d(windowRect=c(0,0,1024,1024))
    segments3d(matrix(as.matrix(t(edge_df[,c(3,4,5,6,7,8)])), ncol=3, byrow=T),
               lwd=2, col="black", line_antialias=TRUE)
    points3d(Matrix::t(reduced_dim_coords[1:3,]), col="black")
    while(sum(sel) < num_roots) {
      ans <- identify3d(Matrix::t(reduced_dim_coords[1:3,!sel]),
                        labels = which(!sel), n = 1,
                        buttons = c("right", "middle"), ...)
      if(!length(ans)) break
      ans <- which(!sel)[ans]
      #points3d(Matrix::t(reduced_dim_coords[1:3,ans]), col="red")
      sel[ans] <- TRUE
    }
  } else {
    ui <- fluidPage(
      titlePanel("Choose your root nodes"),

      # Sidebar layout with input and output definitions ----
      sidebarLayout(

        # Sidebar panel for inputs ----
        sidebarPanel(
          # done button
          actionButton("choose_toggle", "Choose/unchoose"),
          # clear button
          actionButton("reset", "Clear"),
          # done button
          actionButton("done", "Done")
        ),

        # Main panel for displaying outputs ----
        mainPanel(
          plotOutput("plot1", height = 350,
                     click = "plot1_click",
                     brush = brushOpts(id = "plot1_brush"))
        )
      )
    )

    server <- function(input, output) {

      vals <- reactiveValues(
        keeprows = rep(TRUE, nrow(ica_space_df))
      )

      output$plot1 <- renderPlot({
        # Plot the kept and excluded points as two separate data sets
        keep    <- ica_space_df[ vals$keeprows, , drop = FALSE]
        exclude <- ica_space_df[!vals$keeprows, , drop = FALSE]

        ggplot(keep, aes(prin_graph_dim_1, prin_graph_dim_2)) + geom_point(alpha = .7) +
          geom_point(data = exclude, shape = 21, fill = NA, color = "blue") +
          geom_segment(data = edge_df,  aes(x = source_prin_graph_dim_1,
                                            xend = target_prin_graph_dim_1,
                                            y = source_prin_graph_dim_2,
                                            yend = target_prin_graph_dim_2)) +
          labs(x="Component 1", y="Component 2") + monocle3:::monocle_theme_opts()
      })

      # Toggle points that are clicked
      observeEvent(input$plot1_click, {
        res <- nearPoints(ica_space_df, input$plot1_click, allRows = TRUE)

        vals$keeprows <- xor(vals$keeprows, res$selected_)
      })

      # Toggle points that are brushed, when button is clicked
      observeEvent(input$choose_toggle, {
        res <- brushedPoints(ica_space_df, input$plot1_brush, allRows = TRUE)

        vals$keeprows <- xor(vals$keeprows, res$selected_)
      })

      # Reset all points
      observeEvent(input$reset, {
        vals$keeprows <- rep(TRUE, nrow(ica_space_df))
      })

      observeEvent(input$done, {
        stopApp(vals$keeprows)
      })

    }
    sel <- runApp(shinyApp(ui, server))
  }
  ## return indices of selected points
  as.character(ica_space_df$sample_name[which(sel)])
}





