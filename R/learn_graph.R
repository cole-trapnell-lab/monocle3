#' Learn principal graph from the reduced space using reversed graph embedding
#'
#' @description Monocle3 aims to learn how cells transition through a
#' biological program of gene expression changes in an experiment. Each cell
#' can be viewed as a point in a high-dimensional space, where each dimension
#' describes the expression of a different gene. Identifying the program of
#' gene expression changes is equivalent to learning a \emph{trajectory} that
#' the cells follow through this space. However, the more dimensions there are
#' in the analysis, the harder the trajectory is to learn. Fortunately, many
#' genes typically co-vary with one another, and so the dimensionality of the
#' data can be reduced with a wide variety of different algorithms. Monocle3
#' provides two different algorithms for dimensionality reduction via
#' \code{reduce_dimension}.
#' Both take a cell_data_set object and a number of dimensions allowed for the
#' reduced space. You can also provide a model formula indicating some variables
#' (e.g. batch ID or other technical factors) to "subtract" from the data so it
#' doesn't contribute to the trajectory.
#'
#' @details You can choose two different reduction algorithms: Independent Component
#' Analysis (ICA) and Discriminative Dimensionality Reduction with Trees (DDRTree).
#' The choice impacts numerous downstream analysis steps, including \code{\link{order_cells}}.
#' Choosing ICA will execute the ordering procedure described in Trapnell and Cacchiarelli et al.,
#' which was implemented in Monocle version 1. \code{\link[DDRTree]{DDRTree}} is a more recent manifold
#' learning algorithm developed by Qi Mao and colleages. It is substantially more
#' powerful, accurate, and robust for single-cell trajectory analysis than ICA,
#' and is now the default method.
#'
#' Often, experiments include cells from different batches or treatments. You can
#' reduce the effects of these treatments by transforming the data with a linear
#' model prior to dimensionality reduction. To do so, provide a model formula
#' through \code{residual_model_formula_str}.
#'
#'
#' @param cds the cell_data_set upon which to perform this operation
#' @param rge_method Determines how to transform expression values prior to reducing dimensionality
#' @param scale When this argument is set to TRUE (default), it will scale each gene before running trajectory reconstruction.
#' @param close_loop Whether or not to perform an additional run of loop closing after running DDRTree or SimplePPT to identify potential loop structure in the data space
#' @param euclidean_distance_ratio The maximal ratio between the euclidean distance of two tip nodes in the spanning tree inferred from SimplePPT algorithm and
#' that of the maximum distance between any connecting points on the spanning tree allowed to be connected during the loop closure procedure .
#' @param geodesic_distance_ratio  The minimal ratio between the geodestic distance of two tip nodes in the spanning tree inferred from SimplePPT algorithm and
#' that of the length of the diameter path on the spanning tree allowed to be connected during the loop closure procedure. (Both euclidean_distance_ratio and geodesic_distance_ratio
#' need to be satisfied to introduce the edge for loop closure.)
#' @param prune_graph Whether or not to perform an additional run of graph pruning to remove small insignificant branches
#' @param minimal_branch_len The minimal length of the diameter path for a branch to be preserved during graph pruning procedure
#' @param orthogonal_proj_tip Whether to perform orthogonal projection for cells corresponding to the tip principal points. Default to be FALSE
#' @param verbose Whether to emit verbose output during dimensionality reduction
#' @param ... additional arguments to pass to the dimensionality reduction function
#' @return an updated cell_data_set object
#' @references DDRTree: Qi Mao, Li Wang, Steve Goodison, and Yijun Sun. Dimensionality reduction via graph structure learning. In Proceedings of the 21th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining, pages 765–774. ACM, 2015.
#' @references L1graph (generalized SimplePPT): Qi Mao, Li Wang, Ivor Tsang, and Yijun Sun. Principal graph and structure learning based on reversed graph embedding . IEEE Trans. Pattern Anal. Mach. Intell., 5 December 2016.
#' @references Original SimplePPT: Qi Mao, Le Yang, Li Wang, Steve Goodison, Yijun Sun. SimplePPT: A Simple Principal Tree Algorithm https://epubs.siam.org/doi/10.1137/1.9781611974010.89
#' @export
learn_graph <- function(cds,
                        use_partition = TRUE,
                        close_loop = TRUE,
                        learn_graph_control = NULL,
                        verbose = FALSE) {
  reduction_method <- "UMAP"

  if (!is.null(learn_graph_control)) {
    assertthat::assert_that(is(learn_graph_control, "list"))
    assertthat::assert_that(all(names(learn_graph_control) %in%
                                  c("euclidean_distance_ratio",
                                    "geodesic_distance_ratio",
                                    "minimal_branch_len",
                                    "orthogonal_proj_tip",
                                    "prune_graph",
                                    "scale",
                                    "ncenter",
                                    "maxiter",
                                    "eps",
                                    "L1.gamma",
                                    "L1.sigma")),
                            msg = "Unknown variable in learn_graph_control")
  }

  euclidean_distance_ratio <- ifelse(is.null(learn_graph_control$euclidean_distance_ratio), 1, learn_graph_control$euclidean_distance_ratio)
  geodesic_distance_ratio <- ifelse(is.null(learn_graph_control$geodesic_distance_ratio), 1/3, learn_graph_control$geodesic_distance_ratio)
  minimal_branch_len <- ifelse(is.null(learn_graph_control$minimal_branch_len), 10, learn_graph_control$minimal_branch_len)
  orthogonal_proj_tip <- ifelse(is.null(learn_graph_control$orthogonal_proj_tip), FALSE, learn_graph_control$orthogonal_proj_tip)
  prune_graph <- ifelse(is.null(learn_graph_control$prune_graph), TRUE, learn_graph_control$prune_graph)
  ncenter <- learn_graph_control$ncenter
  scale <- ifelse(is.null(learn_graph_control$scale), FALSE, learn_graph_control$scale)
  maxiter <- ifelse(is.null(learn_graph_control$maxiter), 10, learn_graph_control$maxiter)
  eps <- ifelse(is.null(learn_graph_control$eps), 1e-5, learn_graph_control$eps)
  L1.gamma <- ifelse(is.null(learn_graph_control$L1.gamma), 0.5, learn_graph_control$L1.gamma)
  L1.sigma <- ifelse(is.null(learn_graph_control$L1.sigma), 0.01, learn_graph_control$L1.sigma)

  assertthat::assert_that(is(cds, "cell_data_set"))
  assertthat::assert_that(is.logical(use_partition))
  assertthat::assert_that(is.logical(close_loop))
  assertthat::assert_that(is.logical(verbose))
  assertthat::assert_that(is.logical(orthogonal_proj_tip))
  assertthat::assert_that(is.logical(prune_graph))
  assertthat::assert_that(is.logical(scale))
  assertthat::assert_that(is.numeric(euclidean_distance_ratio))
  assertthat::assert_that(is.numeric(geodesic_distance_ratio))
  assertthat::assert_that(is.numeric(minimal_branch_len))
  if(!is.null(ncenter)) {
    assertthat::assert_that(assertthat::is.count(ncenter))
  }
  assertthat::assert_that(assertthat::is.count(maxiter))
  assertthat::assert_that(is.numeric(eps))
  assertthat::assert_that(is.numeric(L1.sigma))
  assertthat::assert_that(is.numeric(L1.sigma))

  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste("No dimensionality reduction for",
                                      reduction_method, "calculated.",
                                      "Please run reduce_dimensions with",
                                      "reduction_method =", reduction_method,
                                      "and cluster_cells before running",
                                      "learn_graph."))
  assertthat::assert_that(!is.null(cds@clusters[[reduction_method]]),
                          msg = paste("No cell clusters for",
                                      reduction_method, "calculated.",
                                      "Please run cluster_cells with",
                                      "reduction_method =", reduction_method,
                                      "before running learn_graph."))

  if (use_partition) {
    partition_list <- cds@clusters[[reduction_method]]$partitions
  } else {
    partition_list <- rep(1, nrow(colData(cds)))
  }

  multi_tree_DDRTree_res <- multi_component_RGE(cds, scale = scale,
                                                reduction_method = reduction_method,
                                                partition_list = partition_list,
                                                irlba_pca_res = reducedDims(cds)[[reduction_method]],
                                                max_components = max_components,
                                                ncenter = ncenter,
                                                maxiter = maxiter,
                                                eps = eps,
                                                L1.gamma = L1.gamma,
                                                L1.sigma = L1.sigma,
                                                close_loop = close_loop,
                                                euclidean_distance_ratio = euclidean_distance_ratio,
                                                geodesic_distance_ratio = geodesic_distance_ratio,
                                                prune_graph = prune_graph,
                                                minimal_branch_len = minimal_branch_len,
                                                verbose = verbose)

  rge_res_W <- multi_tree_DDRTree_res$ddrtree_res_W
  rge_res_Z <- multi_tree_DDRTree_res$ddrtree_res_Z
  rge_res_Y <- multi_tree_DDRTree_res$ddrtree_res_Y
  cds <- multi_tree_DDRTree_res$cds
  dp_mst <- multi_tree_DDRTree_res$dp_mst


  principal_graph(cds)[[reduction_method]] <- dp_mst
  cds@principal_graph_aux[[reduction_method]]$dp_mst <- rge_res_Y

  cds <- project2MST(cds, project_point_to_line_segment, orthogonal_proj_tip,
                     verbose, reduction_method, rge_res_Y)

  cds
}

multi_component_RGE <- function(cds,
                                scale = FALSE,
                                reduction_method,
                                partition_list,
                                max_components,
                                ncenter,
                                irlba_pca_res,
                                maxiter,
                                eps,
                                L1.gamma,
                                L1.sigma,
                                close_loop = FALSE,
                                euclidean_distance_ratio = 1,
                                geodesic_distance_ratio = 1/3,
                                prune_graph = TRUE,
                                minimal_branch_len = minimal_branch_len,
                                verbose = FALSE) {
  X <- t(irlba_pca_res)

  dp_mst <- NULL
  pr_graph_cell_proj_closest_vertex <- NULL
  cell_name_vec <- NULL
  reducedDimK_coord <- NULL
  merge_rge_res <- NULL
  max_ncenter <- 0

  for(cur_comp in sort(unique(partition_list))) {
    if(verbose) {
      message(paste0('Processing louvain component ', cur_comp))
    }

    X_subset <- X[, partition_list == cur_comp]
    if(verbose) message('Current partition is ', cur_comp)

    #add other parameters...
    if(scale) {
      X_subset <- t(as.matrix(scale(t(X_subset))))
    }

    if(is.null(ncenter)) {
      ncenter <- cal_ncenter(ncol(X_subset))
      if(is.null(ncenter) || ncenter >= ncol(X_subset)) {
        ncenter <- ncol(X_subset) - 1
      }
    } else {
      ncenter <- min(ncol(X_subset) - 1, ncenter)
    }

    kmean_res <- NULL

    centers <- t(X_subset)[seq(1, ncol(X_subset), length.out=ncenter), ,
                           drop = F]
    centers <- centers + matrix(stats::rnorm(length(centers), sd = 1e-10),
                                nrow = nrow(centers)) # add random noise

    kmean_res <- tryCatch({
      stats::kmeans(t(X_subset), centers=centers, iter.max = 100)
    }, error = function(err) {
      stats::kmeans(t(X_subset), centers = ncenter, iter.max = 100)
    })

    if (kmean_res$ifault != 0){
      message(paste("Warning: kmeans returned ifault =", kmean_res$ifault))
    }
    nearest_center <- findNearestVertex(t(kmean_res$centers), X_subset,
                                        process_targets_in_blocks=TRUE)
    medioids <- X_subset[, unique(nearest_center)]
    reduced_dim_res <- t(medioids)
    k <- 25
    mat <- t(X_subset)
    if (is.null(k)) {
      k <- round(sqrt(nrow(mat))/2)
      k <- max(10, k)
    }
    if (verbose)
      message("Finding kNN using RANN with ", k, " neighbors")
    dx <- RANN::nn2(mat, k = min(k, nrow(mat) - 1))
    nn.index <- dx$nn.idx[, -1]
    nn.dist <- dx$nn.dists[, -1]

    if (verbose)
      message("Calculating the local density for each sample based on kNNs ...")

    rho <- exp(-rowMeans(nn.dist))
    mat_df <- as.data.frame(mat)
    tmp <- mat_df %>% tibble::rownames_to_column() %>%
      dplyr::mutate(cluster = kmean_res$cluster, density = rho) %>%
      dplyr::group_by(cluster) %>% dplyr::top_n(n = 1, wt = density) %>%
      dplyr::arrange(-dplyr::desc(cluster))

    # select representative cells by highest density
    medioids <- X_subset[, tmp$rowname]

    reduced_dim_res <- t(medioids)
    graph_args <- list(X = X_subset, C0 = medioids, maxiter = maxiter, eps = eps, L1.gamma = L1.gamma, L1.sigma = L1.sigma,
                         verbose = verbose)

    rge_res <- do.call(calc_principal_graph, graph_args)

    names(rge_res)[c(2, 4, 5)] <- c('Y', 'R','objective_vals')
    stree <- rge_res$W

    stree_ori <- stree
    if(close_loop) {
      connect_tips_res <- connect_tips(colData(cds)[partition_list == cur_comp, ],
                                       R = rge_res$R,
                                       stree = stree,
                                       reducedDimK_old = rge_res$Y,
                                       reducedDimS_old = t(reducedDims(cds)[[reduction_method]])[, partition_list == cur_comp],
                                       kmean_res = kmean_res,
                                       euclidean_distance_ratio = euclidean_distance_ratio,
                                       geodesic_distance_ratio = geodesic_distance_ratio,
                                       medioids = medioids,
                                       verbose = verbose)
      stree <- connect_tips_res$stree
    }
    if(prune_graph) {
      if(verbose) {
        message('Running graph pruning ...')
      }
      stree <- pruneTree_in_learnGraph(stree_ori, as.matrix(stree),
                                       minimal_branch_len = minimal_branch_len)
      # remove the points in Y; mediods, etc.
      rge_res$Y <- rge_res$Y[, match(row.names(stree), row.names(stree_ori))]
      rge_res$R <- rge_res$R[, match(row.names(stree), row.names(stree_ori))]
      medioids <- medioids[, row.names(stree)]
    }

    if(is.null(merge_rge_res)) {
      colnames(rge_res$Y) <- paste0('Y_', 1:ncol(rge_res$Y))
      merge_rge_res <- rge_res
      colnames(merge_rge_res$X) <- colnames(X_subset)
      row.names(merge_rge_res$R) <- colnames(X_subset)
      colnames(merge_rge_res$R) <- paste0('Y_', 1:ncol(merge_rge_res$Y))
      merge_rge_res$R <- list(merge_rge_res$R)
      merge_rge_res$stree <- list(stree)
      merge_rge_res$objective_vals <- list(merge_rge_res$objective_vals)

    } else {
      colnames(rge_res$X) <- colnames(X_subset)
      row.names(rge_res$R) <- colnames(X_subset)
      colnames(rge_res$R) <- paste0('Y_', (ncol(merge_rge_res$Y) + 1):(ncol(merge_rge_res$Y) + ncol(rge_res$Y)), sep = "")
      colnames(rge_res$Y) <- paste("Y_", (ncol(merge_rge_res$Y) + 1):(ncol(merge_rge_res$Y) + ncol(rge_res$Y)), sep = "")
      merge_rge_res$Y <- cbind(merge_rge_res$Y, rge_res$Y)
      merge_rge_res$R <- c(merge_rge_res$R, list(rge_res$R))
      merge_rge_res$stree <- c(merge_rge_res$stree, list(stree))
      merge_rge_res$objective_vals <- c(merge_rge_res$objective_vals, list(rge_res$objective_vals))
    }

    if(is.null(reducedDimK_coord)) {
      curr_cell_names <- paste("Y_", 1:ncol(rge_res$Y), sep = "")
      pr_graph_cell_proj_closest_vertex <- matrix(apply(rge_res$R, 1, which.max))
      cell_name_vec <- colnames(X_subset)
    } else {
      curr_cell_names <- paste("Y_", (ncol(reducedDimK_coord) + 1):(ncol(reducedDimK_coord) + ncol(rge_res$Y)), sep = "")
      pr_graph_cell_proj_closest_vertex <- rbind(pr_graph_cell_proj_closest_vertex, matrix(apply(rge_res$R, 1, which.max) + ncol(reducedDimK_coord)))
      cell_name_vec <- c(cell_name_vec, colnames(X_subset))
    }

    curr_reducedDimK_coord <- rge_res$Y
    dimnames(stree) <- list(curr_cell_names, curr_cell_names)
    cur_dp_mst <- igraph::graph.adjacency(stree, mode = "undirected", weighted = TRUE)

    dp_mst <- igraph::graph.union(dp_mst, cur_dp_mst)
    reducedDimK_coord <- cbind(reducedDimK_coord, curr_reducedDimK_coord)
  }

  row.names(pr_graph_cell_proj_closest_vertex) <- cell_name_vec

  ddrtree_res_W <- as.matrix(rge_res$W)
  ddrtree_res_Z <- reducedDims(cds)[[reduction_method]]
  ddrtree_res_Y <- reducedDimK_coord # ensure the order of column names matches that of the original name ids
  # correctly set up R, stree -- the mapping from each cell to the principal graph points
  R <- Matrix::sparseMatrix(i = 1, j = 1, x = 0, dims = c(ncol(cds), ncol(merge_rge_res$Y))) # use sparse matrix for large datasets
  stree <- Matrix::sparseMatrix(i = 1, j = 1, x = 0, dims = c(ncol(merge_rge_res$Y), ncol(merge_rge_res$Y)))
  curr_row_id <- 1
  curr_col_id <- 1
  R_row_names <- NULL
  for(i in 1:length(merge_rge_res$R)) {
    current_R <- merge_rge_res$R[[i]]

    stree[curr_col_id:(curr_col_id + ncol(current_R) - 1), curr_col_id:(curr_col_id + ncol(current_R) - 1)] <- merge_rge_res$stree[[i]]

    curr_row_id <- curr_row_id + nrow(current_R)
    curr_col_id <- curr_col_id + ncol(current_R)
    R_row_names <- c(R_row_names, row.names(current_R))
  }

  row.names(R) <- R_row_names
  R <- R[colnames(cds), ] # reorder the colnames

  cds@principal_graph_aux[[reduction_method]] <- list(stree = stree, Q = merge_rge_res$Q, R = R, objective_vals = merge_rge_res$objective_vals, history = merge_rge_res$history) # rge_res[c('stree', 'Q', 'R', 'objective_vals', 'history')] #
  cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex <- as.data.frame(pr_graph_cell_proj_closest_vertex)[colnames(cds), , drop = F] # Ensure the row order matches up that of the column order of the cds

  colnames(ddrtree_res_Y) <- paste0("Y_", 1:ncol(ddrtree_res_Y), sep = "")

  return(list(cds = cds,
              ddrtree_res_W = ddrtree_res_W,
              ddrtree_res_Z = ddrtree_res_Z,
              ddrtree_res_Y = ddrtree_res_Y,
              dp_mst = dp_mst))
}

# Function to decide a good number of centers for running DDRTree on big datasets
cal_ncenter <- function(ncells, ncells_limit = 1000){
  if(ncells <= ncells_limit) {
    return(NULL)
  }

  round(5 * ncells_limit * log(ncells)/ (log(ncells) + log(ncells_limit)))
}

#' Finds the nearest principal graph node
#' @param data_matrix the input matrix
#' @param target_points the target points
#' @param block_size the number of input matrix rows to process per bloclk
#' @param process_targets_in_blocks whether to process the targets points in blocks instead
#' @keywords internal
findNearestVertex = function(data_matrix, target_points, block_size=50000, process_targets_in_blocks=FALSE){
  closest_vertex = c()
  if (process_targets_in_blocks == FALSE){
    num_blocks = ceiling(ncol(data_matrix) / block_size)
    for (i in 1:num_blocks){
      if (i < num_blocks){
        block = data_matrix[,((((i-1) * block_size)+1):(i*block_size))]
      }else{
        block = data_matrix[,((((i-1) * block_size)+1):(ncol(data_matrix)))]
      }
      distances_Z_to_Y <- proxy::dist(t(block), t(target_points))
      closest_vertex_for_block <- apply(distances_Z_to_Y, 1, function(z) { which.min(z) } )
      closest_vertex = append(closest_vertex, closest_vertex_for_block)
    }
  }else{
    num_blocks = ceiling(ncol(target_points) / block_size)
    dist_to_closest_vertex = rep(Inf, length(ncol(data_matrix)))
    closest_vertex = rep(NA, length(ncol(data_matrix)))
    for (i in 1:num_blocks){
      if (i < num_blocks){
        block = target_points[,((((i-1) * block_size)+1):(i*block_size))]
      }else{
        block = target_points[,((((i-1) * block_size)+1):(ncol(target_points)))]
      }
      distances_Z_to_Y <- proxy::dist(t(data_matrix), t(block))
      closest_vertex_for_block <- apply(distances_Z_to_Y, 1, function(z) { which.min(z) } )
      new_block_distances = distances_Z_to_Y[cbind(1:nrow(distances_Z_to_Y), closest_vertex_for_block)]
      updated_nearest_idx = which(new_block_distances < dist_to_closest_vertex)
      closest_vertex[updated_nearest_idx] = closest_vertex_for_block[updated_nearest_idx] + (i-1) * block_size
      dist_to_closest_vertex[updated_nearest_idx] = new_block_distances[updated_nearest_idx]
    }
  }
  stopifnot(length(closest_vertex) == ncol(data_matrix))
  return (closest_vertex)
}


#' Function to prune the graph
pruneTree_in_learnGraph <- function(stree_ori, stree_loop_clousre, minimal_branch_len = 10){
  if (ncol(stree_ori) < minimal_branch_len)
    return(stree_loop_clousre);
  dimnames(stree_loop_clousre) <- dimnames(stree_ori)
  stree_ori[stree_ori != 0] <- 1
  stree_ori <- igraph::graph_from_adjacency_matrix(stree_ori, mode = 'undirected', weight = NULL)
  stree_loop_clousre[stree_loop_clousre != 0] <- 1
  stree_loop_clousre <- igraph::graph_from_adjacency_matrix(stree_loop_clousre, mode = 'undirected', weight = NULL)

  # get closed loops:
  added_edges <- igraph::get.edgelist(stree_loop_clousre - stree_ori)
  valid_edges <- matrix(ncol = 2, nrow = 0)
  edges_to_remove_df <- matrix(ncol = 2, nrow = 0)
  vertex_top_keep <- NULL

  if(nrow(added_edges) > 0) {
    edge_dists <- apply(added_edges, 1, function(x) igraph::distances(stree_ori, x[1], x[2]))
    valid_edges <- added_edges[which(edge_dists >= minimal_branch_len), , drop = F]
    edges_to_remove_df <- added_edges[which(edge_dists < minimal_branch_len), , drop = F]
  }
  if(nrow(valid_edges) > 0) {
    vertex_top_keep <- as.character(unlist(apply(valid_edges, 1, function(x) igraph::shortest_paths(stree_ori, x[1], x[2])$vpath[[1]]$name )))
  }

  root_cell <- which(igraph::neighborhood.size(stree_ori) == 2)[1]
  if(is.na(root_cell)) {
    root_cell <- igraph::V(stree_ori)$name[1]
  }

  mst_traversal <- igraph::graph.dfs(stree_ori,
                                     root = root_cell,
                                     neimode = "all",
                                     unreachable=FALSE,
                                     father=TRUE)
  mst_traversal$father <- as.numeric(mst_traversal$father)
  vertex_to_be_deleted <- c()

  # further remove other small branches?
  if (length(mst_traversal$order) > 0){
    for (i in 1:length(mst_traversal$order)){
      curr_node <- tryCatch({mst_traversal$order[i]}, error = function(e) { NA })
      if (is.na(curr_node))
        next;
      curr_node_name <- igraph::V(stree_ori)[curr_node]$name

      if (is.na(mst_traversal$father[curr_node]) == FALSE){
        parent_node <- mst_traversal$father[curr_node]
        parent_node_name <- igraph::V(stree_ori)[parent_node]$name

        if (igraph::degree(stree_ori, v=parent_node_name) > 2){
          parent_neighbors <- igraph::neighbors(stree_ori, v=parent_node_name, mode = 'all')

          parent_neighbors_index <- sort(match(parent_neighbors$name, mst_traversal$order$name)) # follow the order of gene expression
          parent_neighbors <- mst_traversal$order$name[parent_neighbors_index]

          tmp <- igraph::delete.edges(stree_ori, paste0(parent_node_name, "|", parent_neighbors))
          tmp_decomposed <- igraph::decompose.graph(tmp)

          comp_a <- tmp_decomposed[unlist(lapply(tmp_decomposed, function(x) {
            parent_neighbors[2] %in% igraph::V(x)$name
          }))][[1]]

          comp_b <- tmp_decomposed[unlist(lapply(tmp_decomposed, function(x) {
            parent_neighbors[3] %in% igraph::V(x)$name
          }))][[1]]

          diameter_len_a <- igraph::diameter(comp_a) + 1
          diameter_len_b <- igraph::diameter(comp_b) + 1

          if(diameter_len_a < minimal_branch_len) {# if loop closure is not applied to cells on this branch
            vertex_to_be_deleted <- c(vertex_to_be_deleted, igraph::V(comp_a)$name)
          }
          if(diameter_len_b < minimal_branch_len) {# if loop closure is not applied to cells on this branch
            vertex_to_be_deleted <- c(vertex_to_be_deleted, igraph::V(comp_b)$name)
          }
        }
      }
    }
  }

  valid_vertex_to_be_deleted <- setdiff(vertex_to_be_deleted, vertex_top_keep)
  stree_loop_clousre <- igraph::delete_vertices(stree_loop_clousre, valid_vertex_to_be_deleted)

  tmp <- edges_to_remove_df[edges_to_remove_df[, 1] %in% igraph::V(stree_loop_clousre)$name & edges_to_remove_df[, 2] %in% igraph::V(stree_loop_clousre)$name, , drop = FALSE]
  if(nrow(tmp) > 0) {
    edges_to_remove <- paste0(tmp[, 1], '|', tmp[, 2])
    stree_loop_clousre <- igraph::delete.edges(stree_loop_clousre, edges_to_remove)
  }

  return(igraph::get.adjacency(stree_loop_clousre))
}

project2MST <- function(cds, Projection_Method, orthogonal_proj_tip = FALSE, verbose, reduction_method, rge_res_Y){
  dp_mst <- principal_graph(cds)[[reduction_method]]
  Z <- t(reducedDims(cds)[[reduction_method]])
  Y <- rge_res_Y

  cds <- findNearestPointOnMST(cds, reduction_method, rge_res_Y)
  closest_vertex <- cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex

  closest_vertex_names <- colnames(Y)[closest_vertex[, 1]]

  tip_leaves <- names(which(igraph::degree(dp_mst) == 1))

  if(!is.function(Projection_Method)) {
    P <- Y[, closest_vertex]
  } else{ # project cell to the principal graph (each cell will get different coordinates - except certain rare cases)
    P <- matrix(rep(0, length(Z)), nrow = nrow(Z)) #Y
    nearest_edges <- matrix(rep(0, length(Z[1:2, ])), ncol = 2) # nearest principal graph edge for each cell
    row.names(nearest_edges) <-  colnames(cds)
    for(i in 1:length(closest_vertex)) { # This loop is going to be slow
      neighbors <- names(igraph::neighborhood(dp_mst, nodes = closest_vertex_names[i], mode = 'all')[[1]])[-1]
      projection <- NULL
      distance <- NULL
      Z_i <- Z[, i]

      for(neighbor in neighbors) {
        if(closest_vertex_names[i] %in% tip_leaves) {
          if(orthogonal_proj_tip) {
            tmp <- projPointOnLine(Z_i, Y[, c(closest_vertex_names[i], neighbor)]) #projPointOnLine: always perform orthogonal projection to the line
          } else {
            tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], neighbor)]) #projPointOnLine: always perform orthogonal projection to the line
          }
        }
        else {
          tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], neighbor)])
        }
        if(any(is.na(tmp))) { # in case coordinates for nodes closest_vertex_names[i] and neighbor are the same
          tmp <- Y[, neighbor]
        }
        projection <- rbind(projection, tmp)
        distance <- c(distance, dist(rbind(Z_i, tmp)))
      }
      if(class(projection) != 'matrix') {
        projection <- as.matrix(projection)
      }

      which_min <- which.min(distance)

      if(length(which_min) == 0)
        browser()

      P[, i] <- projection[which_min, ] #use only the first index to avoid assignment error
      nearest_edges[i, ] <- c(closest_vertex_names[i], neighbors[which_min])
    }
  }

  colnames(P) <- colnames(Z)

  dp_mst_list <- igraph::decompose.graph(dp_mst)
  dp_mst_df <- NULL
  partitions <- cds@clusters[[reduction_method]]$partitions

  if(length(dp_mst_list) == 1 & length(unique(partitions)) > 1) {
    partitions <- 1
  }

  if(!is.null(partitions)) {
    for(cur_partition in sort(unique(partitions))) {
      data_df <- NULL

      if(verbose) {
        message('\nProjecting cells to principal points for partition: ', cur_partition)
      }

      subset_cds_col_names <- colnames(cds)[cds@clusters[[reduction_method]]$partitions == cur_partition]
      cur_z <- Z[, subset_cds_col_names]
      cur_p <- P[, subset_cds_col_names]

      if (ncol(cur_p) > 0 && nrow(cur_p) > 0){
        cur_centroid_name <- igraph::V(dp_mst_list[[as.numeric(cur_partition)]])$name

        cur_nearest_edges <- nearest_edges[subset_cds_col_names, ] # the nearest edge for each cell
        data_df <- cbind(as.data.frame(t(cur_p)), apply(cur_nearest_edges, 1, sort) %>% t()) # cell by coord + edge (sorted)
        row.names(data_df) <- colnames(cur_p)
        colnames(data_df) <- c(paste0("P_", 1:nrow(cur_p)), 'source', 'target')
        # colnames(data_df)[(ncol(data_df) - (nrow(cur_p) - 1)):ncol(data_df)] <- paste0('S_', 1:nrow(cur_p))

        # sort each cell's distance to the source in each principal edge group
        data_df$distance_2_source <-  sqrt(colSums((cur_p - t(reducedDims(cds)[[reduction_method]])[,data_df[, 'source']])^2))
        data_df <- data_df %>% tibble::rownames_to_column() %>% dplyr::mutate(group = paste(source, target, sep = '_')) %>%
          dplyr::arrange(group, desc(-distance_2_source))

        # add the links from the source to the nearest points belong to the principal edge and also all following connections between those points
        data_df <- data_df %>% dplyr::group_by(group) %>% dplyr::mutate(new_source = dplyr::lag(rowname), new_target = rowname) # NA  1  2  3 -- lag(1:3)
        # use the correct name of the source point
        data_df[is.na(data_df$new_source), "new_source"] <- as.character(as.matrix(data_df[is.na(data_df$new_source), 'source']))

        # add the links from the last point on the principal edge to the target point of the edge
        added_rows <- which(is.na(data_df$new_source) & is.na(data_df$new_target)) # find those rows
        data_df <- as.data.frame(data_df, stringsAsFactors = F) # ????????
        data_df <- as.data.frame(as.matrix(data_df), stringsAsFactors = F)
        data_df[added_rows, c('new_source', 'new_target')] <- data_df[added_rows - 1, c('rowname', 'target')] # assign names for the points

        # calculate distance between each pair
        #cur_p <- cbind(cur_p, cds@reducedDimK[, cur_centroid_name]) # append the coordinates of principal graph points
        aug_P = cbind(cur_p, rge_res_Y)
        data_df$weight <-  sqrt(colSums((aug_P[, data_df$new_source] - aug_P[, data_df$new_target]))^2)
        # add the minimal positive distance between any points to the distance matrix
        data_df$weight <- data_df$weight + min(data_df$weight[data_df$weight > 0])

        # Calculate distance between two connected nodes directly from the original graph
        edge_list <- as.data.frame(igraph::get.edgelist(dp_mst_list[[as.numeric(cur_partition)]]), stringsAsFactors=FALSE)
        dp <- as.matrix(dist(t(rge_res_Y)[cur_centroid_name,]))
        edge_list$weight <- dp[cbind(edge_list[, 1], edge_list[, 2])]
        colnames(edge_list) <- c("new_source", "new_target", 'weight')

        dp_mst_df <- Reduce(rbind, list(dp_mst_df, data_df[, c("new_source", "new_target", 'weight')], edge_list))
        # dp_mst <- graph.union(dp_mst, cur_dp_mst, dp_mst_pc)
      }
    }
  }
  dp_mst <- igraph::graph.data.frame(dp_mst_df, directed = FALSE)
  cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_tree <- dp_mst
  cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_dist <- P

  closest_vertex_df <- as.matrix(closest_vertex)
  row.names(closest_vertex_df) <- row.names(closest_vertex)
  cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df

  cds
}

# Project each point to the nearest on the MST:
findNearestPointOnMST <- function(cds, reduction_method, rge_res_Y){
  dp_mst <- principal_graph(cds)[[reduction_method]]
  dp_mst_list <- igraph::decompose.graph(dp_mst)

  if(length(unique(cds@clusters[[reduction_method]]$partitions)) != length(dp_mst_list)) {
    dp_mst_list <- list(dp_mst)
  }

  closest_vertex_df <- NULL
  cur_start_index <- 0

  for(i in 1:length(dp_mst_list)) {
    cur_dp_mst <- dp_mst_list[[i]]

    if(length(dp_mst_list) == 1) {
      Z <- t(reducedDims(cds)[[reduction_method]])
    } else {
      Z <- t(reducedDims(cds)[[reduction_method]])[, cds@clusters[[reduction_method]]$partitions == i]
    }
    Y <- rge_res_Y[, igraph::V(cur_dp_mst)$name]

    tip_leaves <- names(which(igraph::degree(cur_dp_mst) == 1))

    closest_vertex_ori <- findNearestVertex(Z, Y)
    closest_vertex <- closest_vertex_ori + cur_start_index

    closest_vertex_names <- colnames(Y)[closest_vertex_ori]
    cur_name <- names(closest_vertex)
    closest_vertex <- as.matrix(closest_vertex)
    row.names(closest_vertex) <- cur_name #original cell names for projection
    closest_vertex_df <- rbind(closest_vertex_df, closest_vertex) #index on Z

    cur_start_index <- cur_start_index + igraph::vcount(cur_dp_mst)
  }
  closest_vertex_df <- closest_vertex_df[colnames(cds), , drop = F]
  cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df #as.matrix(closest_vertex)
  cds
}

# Project point to line segment (in >= 2 dimensions)
project_point_to_line_segment <- function(p, df){
  # returns q the closest point to p on the line segment from A to B
  A <- df[, 1]
  B <- df[, 2]
  # vector from A to B
  AB <- (B-A)
  # squared distance from A to B
  AB_squared = sum(AB^2)
  if(AB_squared == 0) {
    # A and B are the same point
    q <- A
  }
  else {
    # vector from A to p
    Ap <- (p-A)
    # from http://stackoverflow.com/questions/849211/
    # Consider the line extending the segment, parameterized as A + t (B - A)
    # We find projection of point p onto the line.
    # It falls where t = [(p-A) . (B-A)] / |B-A|^2
    # t <- max(0, min(1, sum(Ap * AB) / AB_squared))
    t <- sum(Ap * AB) / AB_squared

    if (t < 0.0) {
      # "Before" A on the line, just return A
      q <- A
    }
    else if (t > 1.0) {
      # "After" B on the line, just return B
      q <- B
    }
    else {
      # projection lines "inbetween" A and B on the line
      q <- A + t * AB#
    }
  }
  return(q)
}

#project points to a line (in  >= 2 dimensions)
projPointOnLine <- function(point, line) {
  ap <- point - line[, 1]
  ab <- line[, 2] - line[, 1]

  res <- line[, 1] + c((ap %*% ab) / (ab %*% ab)) * ab
  return(res)
}


#' Function to automatically learn the structure of data by either using L1-graph or the spanning-tree formulization
#' @param X the input data DxN
#' @param C0 the initialization of centroids
#' @param maxiter maximum number of iteraction
#' @param eps relative objective difference
#' @param L1.gamma regularization parameter for k-means (the prefix of 'param' is used to avoid name collision with gamma)
#' @param L1.sigma bandwidth parameter
#' @param verbose emit results from iteraction
#' @return a list of X, C, W, P, objs
#' X is the input data
#' C is the centers for principal graph
#' W is the pricipal graph matrix
#' P is the cluster assignment matrix
#' objs is the objective value for the function
calc_principal_graph <- function(X, C0,
                                 maxiter = 10,
                                 eps = 1e-5,
                                 L1.gamma = 0.5,
                                 L1.sigma = 0.01,
                                 verbose = T) {

  C <- C0;
  K <- ncol(C)
  objs <- c()
  for(iter in 1:maxiter){
    norm_sq <- repmat(t(colSums(C^2)), K, 1) #this part calculates the cost matrix Phi
    Phi <- norm_sq + t(norm_sq) - 2 * t(C) %*% C

    g <- igraph::graph.adjacency(Phi, mode = 'lower', diag = T, weighted = T)
    g_mst <- igraph::mst(g)
    stree <- igraph::get.adjacency(g_mst, attr = 'weight', type = 'lower')
    stree_ori <- stree

    #convert to matrix:
    stree <- as.matrix(stree)
    stree <- stree + t(stree)

    W <- stree != 0
    obj_W <- sum(sum(stree))

    res = soft_assignment(X, C, L1.sigma)
    P <- res$P
    obj_P <- res$obj

    obj <- obj_W + L1.gamma * obj_P
    objs = c(objs, obj)
    if(verbose)
      message('iter = ', iter, ' obj = ', obj)

    if(iter > 1){
      relative_diff = abs( objs[iter-1] - obj) / abs(objs[iter-1]);
      if(relative_diff < eps){
        if(verbose)
          message('eps = ', relative_diff, ', converge.')
        break
      }
      if(iter >= maxiter){
        if(verbose)
          message('eps = ', relative_diff, ' reach maxiter.')
      }
    }

    C <- generate_centers(X, W, P, L1.gamma)

  }

  return(list(X = X, C = C, W = W, P = P, objs = objs))
}


#' function to reproduce the behavior of repmat function in matlab to replicate and tile an matrix
#' @param X matrix for tiling and replicate the data
#' @param m a numeric value for tiling a matrix
#' @param n a numeric value for tiling a matrix
#' @return a matrix
repmat = function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}



#' Function to calculate the third term in the objective function
#' @param X input data
#' @param C center of grap (D * K)
#' @param sigma bandwidth parameter
#' @return a matrix with diagonal element as 1 while other elements as zero (eye matrix)
soft_assignment <- function(X, C, sigma){

  D <- nrow(X); N <- ncol(X)
  K <- ncol(C)
  norm_X_sq <- repmat(t(t(colSums(X^2))), 1, K);
  norm_C_sq <- repmat(t(colSums(C^2)), N, 1);
  dist_XC <- norm_X_sq + norm_C_sq - 2 * t(X) %*% C

  # %% handle numerical problems 0/0 for P
  min_dist <- apply(dist_XC, 1, min) #rowMin(dist_XC)

  dist_XC <- dist_XC - repmat(t(t(min_dist)), 1, K )
  Phi_XC <- exp(- dist_XC / sigma)
  P <- Phi_XC / repmat(t(t(rowSums(Phi_XC))), 1, K)

  obj <- - sigma * sum( log( rowSums( exp(- dist_XC/sigma)) ) #why not \sigma * log (p_{ij})
                        - min_dist/ sigma );

  return(list(P = P, obj = obj))
}

#' Function to reproduce the behavior of eye function in matlab
#' @param X input data
#' @param W the pricipal graph matrix
#' @param P the cluster assignment matrix
#' @param param.gamma regularization parameter for k-means (the prefix of 'param' is used to avoid name collision with gamma)
#' @return A matrix C for the centers for principal graph
#'
generate_centers <- function(X, W, P, param.gamma){
  D <- nrow(X); N <- nrow(X)
  K <- ncol(W)
  # prevent singular
  Q <- 2 *( diag(colSums(W)) - W ) + param.gamma * diag(colSums(P))
  B <-  param.gamma * X %*% P;
  C <- B %*% solve(Q)   #equation 22

  return(C)
}

#' functions to connect the tip points after learning the DDRTree or simplePPT tree
connect_tips <- function(pd,
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
  cluster_graph_res <- compute_partitions(louvain_res$g, louvain_res$optim_res, qval_thresh=qval_thresh, verbose = verbose)
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



