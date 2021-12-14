#' Learn principal graph from the reduced dimension space using reversed graph
#' embedding
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
#' \code{reduce_dimension} (UMAP and tSNE). Both take a cell_data_set object
#' and a number of dimensions allowed for the reduced space. You can also
#' provide a model formula indicating some variables (e.g. batch ID or other
#' technical factors) to "subtract" from the data so it doesn't contribute to
#' the trajectory. The function \code{learn_graph} is the fourth step in the
#' trajectory building process after \code{preprocess_cds},
#' \code{reduce_dimension}, and \code{cluster_cells}. After
#' \code{learn_graph}, \code{order_cells} is typically called.
#'
#' @section Optional \code{learn_graph_control} parameters:
#' \describe{
#'   \item{euclidean_distance_ratio:}{The maximal ratio between the euclidean
#'   distance of two tip nodes in the spanning tree and the maximum distance
#'   between any connecting points on the spanning tree allowed to be connected
#'   during the loop closure procedure. Default is 1.}
#'   \item{geodesic_distance_ratio:}{The minimal ratio between the geodesic
#'   distance of two tip nodes in the spanning tree and the length of the
#'   diameter path on the spanning tree allowed to be connected during the loop
#'   closure procedure. (Both euclidean_distance_ratio and
#'   geodesic_distance_ratio need to be satisfied to introduce the edge for
#'   loop closure). Default is 1/3.}
#'   \item{minimal_branch_len:}{The minimal length of the diameter path for a
#'   branch to be preserved during graph pruning procedure. Default is 10.}
#'   \item{orthogonal_proj_tip:}{ Whether to perform orthogonal projection for
#'   cells corresponding to the tip principal points. Default is FALSE.}
#'   \item{prune_graph:}{Whether or not to perform an additional round of graph
#'   pruning to remove small insignificant branches. Default is TRUE.}
#'   \item{scale:}{}
#'   \item{ncenter:}{}
#'   \item{nn.k:}{Maximum number of nearest neighbors to compute in the
#'   reversed graph embedding. Set k=NULL
#'   to let learn_graph estimate k. Default is 25.}
#'   \item{rann.k:}{nn.k replaces rann.k but rann.k is available for
#'   compatibility with existing code.}
#'   \item{maxiter:}{}
#'   \item{eps:}{}
#'   \item{L1.gamma:}{}
#'   \item{L1.sigma:}{}
#'   \item{nn.method:}{The method to use for finding nearest neighbors.
#'   nn.method can be one of 'nn2', 'annoy', or 'hnsw'.}
#'   \item{nn.metric:}{The distance metric for the annoy or hnsw nearest
#'   neighbor index build. See help(set_nn_control) for more information.}
#'   \item{nn.n_trees:}{The number of trees used to build the annoy nearest
#'   neighbor index. See help(set_nn_control) for more information.}
#'   \item{nn.search_k:}{The number of nodes to search in an annoy index
#'   search. See help(set_nn_control) for more information.}
#'   \item{nn.M:}{Related to internal dimensionality of HNSW index. See
#'   help(set_nn_control) for more information.}
#'   \item{nn.ef_construction:}{Controls the HNSW index build speed/accuracy
#'   tradeoff.}
#'   \item{nn.ef:}{Controls the HNSW index search speed/accuracy tradeoff.
#'   See help(set_nn_control) for more information.}
#'   \item{nn.grain_size:}{Used by annoy and HNSW to set the minimum amount
#'   of work to do per thread. See help(set_nn_control) for more
#'   information.}
#'   \item{nn.cores:}{Used by annoy and HNSW to control the number of
#'   threads used. See help(set_nn_control) for more information.}
#' }
#'
#' @param cds the cell_data_set upon which to perform this operation
#' @param use_partition logical parameter that determines whether to use
#'   partitions calculated during \code{cluster_cells} and therefore to learn
#'   disjoint graph in each partition. When \code{use_partition = FALSE}, a
#'   single graph is learned across all partitions. Default is TRUE.
#' @param close_loop logical parameter that determines whether or not to
#'   perform an additional run of loop closing after estimating the principal
#'   graphs to identify potential loop structure in the data space. Default is
#'   TRUE.
#' @param learn_graph_control NULL or a list of control parameters to be
#'   passed to the reversed graph embedding function. Default is NULL. A list
#'   of potential control parameters is provided in details.
#' @param verbose Whether to emit verbose output during graph learning.
#' @return an updated cell_data_set object
#' @export
learn_graph <- function(cds,
                        use_partition = TRUE,
                        close_loop = TRUE,
                        learn_graph_control = NULL,
                        verbose = FALSE) {
  reduction_method <- "UMAP"
  if (!is.null(learn_graph_control)) {
    assertthat::assert_that(methods::is(learn_graph_control, "list"))
    assertthat::assert_that(all(names(learn_graph_control) %in%
                                  c("euclidean_distance_ratio",
                                    "geodesic_distance_ratio",
                                    "minimal_branch_len",
                                    "orthogonal_proj_tip",
                                    "prune_graph",
                                    "scale",
                                    "ncenter",
                                    "nn.k",
                                    "rann.k",
                                    "maxiter",
                                    "eps",
                                    "L1.gamma",
                                    "L1.sigma",
                                    "nn.method",
                                    "nn.metric",
                                    "nn.n_trees",
                                    "nn.search_k",
                                    "nn.M",
                                    "nn.ef_construction",
                                    "nn.ef",
                                    "nn.grain_size",
                                    "nn.cores")),
                            msg = "Unknown variable in learn_graph_control")
  }

  if(!is.null(learn_graph_control[['rann.k']]) && !is.null(learn_graph_control[['nn.k']])) {
    assertthat::assert_that(learn_graph_control[['rann.k']] == learn_graph_control[['nn.k']],
                            msg=paste0('both learn_graph_control$nn.k and learn_graph_control$rann.k are',
                                       ' defined and are unequal. See help(learn_graph) for more',
                                       ' information.'))
  }

  if(is.null(learn_graph_control[['nn.k']]) && !is.null(learn_graph_control[['rann.k']]))
    learn_graph_control[['nn.k']] <- learn_graph_control[['rann.k']]

  euclidean_distance_ratio <-
    ifelse(is.null(learn_graph_control$euclidean_distance_ratio), 1,
           learn_graph_control$euclidean_distance_ratio)
  geodesic_distance_ratio <-
    ifelse(is.null(learn_graph_control$geodesic_distance_ratio), 1/3,
           learn_graph_control$geodesic_distance_ratio)
  minimal_branch_len <-
    ifelse(is.null(learn_graph_control$minimal_branch_len), 10,
           learn_graph_control$minimal_branch_len)
  orthogonal_proj_tip <-
    ifelse(is.null(learn_graph_control$orthogonal_proj_tip), FALSE,
           learn_graph_control$orthogonal_proj_tip)
  prune_graph <- ifelse(is.null(learn_graph_control$prune_graph), TRUE,
                        learn_graph_control$prune_graph)
  ncenter <- learn_graph_control$ncenter
  scale <- ifelse(is.null(learn_graph_control$scale), FALSE,
                  learn_graph_control$scale)
  nn.k <- ifelse(is.null(learn_graph_control[['nn.k']]), 25,
                   learn_graph_control[['nn.k']])
  maxiter <- ifelse(is.null(learn_graph_control$maxiter), 10,
                    learn_graph_control$maxiter)
  eps <- ifelse(is.null(learn_graph_control$eps), 1e-5,
                learn_graph_control$eps)
  L1.gamma <- ifelse(is.null(learn_graph_control$L1.gamma), 0.5,
                     learn_graph_control$L1.gamma)
  L1.sigma <- ifelse(is.null(learn_graph_control$L1.sigma), 0.01,
                     learn_graph_control$L1.sigma)
  
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(reduction_method %in% c('UMAP'), msg=paste0('unsupported or invalid reduction method \'', reduction_method, '\''))
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
  assertthat::assert_that(assertthat::is.count(nn.k))
  assertthat::assert_that(is.numeric(eps))
  assertthat::assert_that(is.numeric(L1.sigma))
  assertthat::assert_that(is.numeric(L1.sigma))

  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste("No dimensionality reduction for",
                                      reduction_method, "calculated.",
                                      "Please run reduce_dimension with",
                                      "reduction_method =", reduction_method,
                                      "and cluster_cells before running",
                                      "learn_graph."))
  assertthat::assert_that(!is.null(cds@clusters[[reduction_method]]),
                          msg = paste("No cell clusters for",
                                      reduction_method, "calculated.",
                                      "Please run cluster_cells with",
                                      "reduction_method =", reduction_method,
                                      "before running learn_graph."))

  nn_control <- list()
  if(!is.null(learn_graph_control[['nn.method']])) nn_control[['method']] <- learn_graph_control[['nn.method']]
  if(!is.null(learn_graph_control[['nn.metric']])) nn_control[['metric']] <- learn_graph_control[['nn.metric']]
  if(!is.null(learn_graph_control[['nn.n_trees']])) nn_control[['n_trees']] <- learn_graph_control[['nn.n_trees']]
  if(!is.null(learn_graph_control[['nn.search_k']])) nn_control[['search_k']] <- learn_graph_control[['nn.search_k']]
  if(!is.null(learn_graph_control[['nn.M']])) nn_control[['M']] <- learn_graph_control[['nn.M']]
  if(!is.null(learn_graph_control[['nn.ef_construction']])) nn_control[['ef_construction']] <- learn_graph_control[['nn.ef_construction']]
  if(!is.null(learn_graph_control[['nn.ef']])) nn_control[['ef']] <- learn_graph_control[['nn.ef']]
  if(!is.null(learn_graph_control[['nn.grain_size']])) nn_control[['grain_size']] <- learn_graph_control[['nn.grain_size']]
  if(!is.null(learn_graph_control[['nn.cores']])) nn_control[['cores']] <- learn_graph_control[['nn.cores']]

  if(verbose)
    report_nn_control('nn_control: ', nn_control)

  nn_control_default <- get_global_variable('nn_control_annoy_euclidean')
  nn_control <- set_nn_control(mode=3,
                               nn_control=nn_control,
                               nn_control_default=nn_control_default,
                               nn_index=NULL,
                               k=nn.k,
                               verbose=verbose)

  if (use_partition) {
    partition_list <- cds@clusters[[reduction_method]]$partitions
  } else {
    partition_list <- rep(1, nrow(colData(cds)))
  }

  multi_tree_DDRTree_res <-
    multi_component_RGE(cds, scale = scale,
                        reduction_method = reduction_method,
                        partition_list = partition_list,
                        irlba_pca_res = reducedDims(cds)[[reduction_method]],
                        max_components = max_components,
                        ncenter = ncenter,
                        nn.k = nn.k,
                        nn_control = nn_control,
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
                                nn.k=25,
                                nn_control=list(),
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

  for(cur_comp in sort(unique(partition_list))) {  #  for loop 1  start

    if(verbose) {
      message(paste0('Processing partition component ', cur_comp))
    }

    X_subset <- X[, partition_list == cur_comp]

    if(verbose) message('Current partition is ', cur_comp)

    #add other parameters...
    if(scale) {
      X_subset <- t(as.matrix(scale(t(X_subset))))
    }

    if(is.null(ncenter)) {
      num_clusters_in_partition <-
        length(unique(clusters(cds, reduction_method)[colnames(X_subset)]))
      num_cells_in_partition = ncol(X_subset)
      curr_ncenter <- cal_ncenter(num_clusters_in_partition, num_cells_in_partition)
      if(is.null(curr_ncenter) || curr_ncenter >= ncol(X_subset)) {
        curr_ncenter <- ncol(X_subset) - 1
      }
    } else {
      curr_ncenter <- min(ncol(X_subset) - 1, ncenter)
    }
    if (verbose)
      message(paste("Using", curr_ncenter, "nodes for principal graph"))

    kmean_res <- NULL

    centers <- t(X_subset)[seq(1, ncol(X_subset), length.out=curr_ncenter), ,
                           drop = FALSE]
    centers <- centers + matrix(stats::rnorm(length(centers), sd = 1e-10),
                                nrow = nrow(centers)) # add random noise

    kmean_res <- tryCatch({
      stats::kmeans(t(X_subset), centers=centers, iter.max = 100)
    }, error = function(err) {
      stats::kmeans(t(X_subset), centers = curr_ncenter, iter.max = 100)
    })

    if (kmean_res$ifault != 0){
      message(paste("Warning: kmeans returned ifault =", kmean_res$ifault))
    }
    nearest_center <- find_nearest_vertex(t(kmean_res$centers), X_subset,
                                          process_targets_in_blocks=TRUE)
    medioids <- X_subset[, unique(nearest_center)]
    reduced_dim_res <- t(medioids)
    mat <- t(X_subset)
    if (is.null(nn.k)) {
      k <- round(sqrt(nrow(mat))/2)
      k <- max(10, k)
    } else {
      k <- nn.k
    }

    if (verbose)
      message("Finding kNN with ", k, " neighbors")

    nn_method <- nn_control[['method']]
#    dx <- RANN::nn2(mat, k = min(k, nrow(mat) - 1))  # replaced by search_nn_matrix below
    dx <- search_nn_matrix(subject_matrix=mat, query_matrix=mat, k=min(k, nrow(mat) - 1), nn_control=nn_control, verbose=verbose)
    if(nn_method == 'annoy' || nn_method == 'hnsw')
      dx <- swap_nn_row_index_point(nn_res=dx, verbose=verbose)

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
    graph_args <- list(X = X_subset, C0 = medioids, maxiter = maxiter,
                       eps = eps, L1.gamma = L1.gamma, L1.sigma = L1.sigma,
                       verbose = verbose)

    rge_res <- do.call(calc_principal_graph, graph_args)

    names(rge_res)[c(2, 4, 5)] <- c('Y', 'R','objective_vals')
    stree <- rge_res$W

    stree_ori <- stree
    if(close_loop) {
      reduce_dims_old <-
        t(reducedDims(cds)[[reduction_method]])[, partition_list == cur_comp]
      connect_tips_res <-
        connect_tips(cds,
                     pd = colData(cds)[partition_list == cur_comp, ],
                     R = rge_res$R,
                     stree = stree,
                     reducedDimK_old = rge_res$Y,
                     reducedDimS_old = reduce_dims_old,
                     k = 25,
                     nn_control = nn_control,
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
      stree <- prune_tree(stree_ori, as.matrix(stree),
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
      colnames(rge_res$R) <-
        paste0('Y_', (ncol(merge_rge_res$Y) + 1):(ncol(merge_rge_res$Y) +
                                                    ncol(rge_res$Y)), sep = "")
      colnames(rge_res$Y) <-
        paste("Y_", (ncol(merge_rge_res$Y) + 1):(ncol(merge_rge_res$Y) +
                                                   ncol(rge_res$Y)), sep = "")
      merge_rge_res$Y <- cbind(merge_rge_res$Y, rge_res$Y)
      merge_rge_res$R <- c(merge_rge_res$R, list(rge_res$R))
      merge_rge_res$stree <- c(merge_rge_res$stree, list(stree))
      merge_rge_res$objective_vals <- c(merge_rge_res$objective_vals,
                                        list(rge_res$objective_vals))
    }

    if(is.null(reducedDimK_coord)) {
      curr_cell_names <- paste("Y_", 1:ncol(rge_res$Y), sep = "")
      pr_graph_cell_proj_closest_vertex <- matrix(apply(rge_res$R, 1,
                                                        which.max))
      cell_name_vec <- colnames(X_subset)
    } else {
      curr_cell_names <-
        paste("Y_", (ncol(reducedDimK_coord) + 1):(ncol(reducedDimK_coord) +
                                                     ncol(rge_res$Y)),
              sep = "")
      pr_graph_cell_proj_closest_vertex <-
        rbind(pr_graph_cell_proj_closest_vertex,
              matrix(apply(rge_res$R, 1, which.max) + ncol(reducedDimK_coord)))
      cell_name_vec <- c(cell_name_vec, colnames(X_subset))
    }

    curr_reducedDimK_coord <- rge_res$Y
    dimnames(stree) <- list(curr_cell_names, curr_cell_names)
    cur_dp_mst <- igraph::graph.adjacency(stree, mode = "undirected",
                                          weighted = TRUE)

    dp_mst <- igraph::graph.union(dp_mst, cur_dp_mst)
    reducedDimK_coord <- cbind(reducedDimK_coord, curr_reducedDimK_coord)
  }  #  for loop 1  end

  row.names(pr_graph_cell_proj_closest_vertex) <- cell_name_vec

  ddrtree_res_W <- as.matrix(rge_res$W)
  ddrtree_res_Z <- reducedDims(cds)[[reduction_method]]
  ddrtree_res_Y <- reducedDimK_coord

  R <- Matrix::sparseMatrix(i = 1, j = 1, x = 0,
                            dims = c(ncol(cds), ncol(merge_rge_res$Y)))
  stree <- Matrix::sparseMatrix(i = 1, j = 1, x = 0,
                                dims = c(ncol(merge_rge_res$Y),
                                         ncol(merge_rge_res$Y)))
  curr_row_id <- 1
  curr_col_id <- 1
  R_row_names <- NULL
  for(i in 1:length(merge_rge_res$R)) {
    current_R <- merge_rge_res$R[[i]]

    stree[curr_col_id:(curr_col_id + ncol(current_R) - 1),
          curr_col_id:(curr_col_id + ncol(current_R) - 1)] <-
      merge_rge_res$stree[[i]]

    curr_row_id <- curr_row_id + nrow(current_R)
    curr_col_id <- curr_col_id + ncol(current_R)
    R_row_names <- c(R_row_names, row.names(current_R))
  }

  row.names(R) <- R_row_names
  R <- R[colnames(cds), ] # reorder the colnames

  cds@principal_graph_aux[[reduction_method]] <-
    list(stree = stree,
         Q = merge_rge_res$Q,
         R = R,
         objective_vals = merge_rge_res$objective_vals,
         history = merge_rge_res$history)
  cds@principal_graph_aux[[
    reduction_method]]$pr_graph_cell_proj_closest_vertex <-
    as.data.frame(pr_graph_cell_proj_closest_vertex)[colnames(cds), , drop = FALSE]
  colnames(ddrtree_res_Y) <- paste0("Y_", 1:ncol(ddrtree_res_Y), sep = "")

  return(list(cds = cds,
              ddrtree_res_W = ddrtree_res_W,
              ddrtree_res_Z = ddrtree_res_Z,
              ddrtree_res_Y = ddrtree_res_Y,
              dp_mst = dp_mst))
}


cal_ncenter <- function(num_cell_communities, ncells,
                        nodes_per_log10_cells=15) {
  round(num_cell_communities * nodes_per_log10_cells * log10(ncells))
}


#' Finds the nearest principal graph node
#' @param data_matrix the input matrix
#' @param target_points the target points
#' @param block_size the number of input matrix rows to process per block
#' @param process_targets_in_blocks whether to process the targets points in
#'   blocks instead
#' @keywords internal
find_nearest_vertex <- function(data_matrix, target_points, block_size=50000,
                                process_targets_in_blocks=FALSE){
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
      closest_vertex_for_block <- apply(distances_Z_to_Y, 1,
                                        function(z) { which.min(z) } )
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
      closest_vertex_for_block <- apply(distances_Z_to_Y, 1,
                                        function(z) { which.min(z) } )
      new_block_distances <- distances_Z_to_Y[cbind(1:nrow(distances_Z_to_Y),
                                                    closest_vertex_for_block)]
      updated_nearest_idx <- which(new_block_distances < dist_to_closest_vertex)
      closest_vertex[updated_nearest_idx] <-
        closest_vertex_for_block[updated_nearest_idx] + (i-1) * block_size
      dist_to_closest_vertex[updated_nearest_idx] <-
        new_block_distances[updated_nearest_idx]
    }
  }
  stopifnot(length(closest_vertex) == ncol(data_matrix))
  return (closest_vertex)
}


# Function to prune the graph
prune_tree <- function(stree_ori, stree_loop_closure,
                       minimal_branch_len = 10) {
  if (ncol(stree_ori) < minimal_branch_len)
    return(stree_loop_closure);
  dimnames(stree_loop_closure) <- dimnames(stree_ori)
  stree_ori[stree_ori != 0] <- 1
  stree_ori <- igraph::graph_from_adjacency_matrix(stree_ori,
                                                   mode = 'undirected',
                                                   weight = NULL)
  stree_loop_closure[stree_loop_closure != 0] <- 1
  stree_loop_closure <- igraph::graph_from_adjacency_matrix(stree_loop_closure,
                                                            mode = 'undirected',
                                                            weight = NULL)

  # get closed loops:
  added_edges <- igraph::get.edgelist(stree_loop_closure - stree_ori)
  valid_edges <- matrix(ncol = 2, nrow = 0)
  edges_to_remove_df <- matrix(ncol = 2, nrow = 0)
  vertex_top_keep <- NULL

  if(nrow(added_edges) > 0) {
    edge_dists <- apply(added_edges, 1,
                        function(x) igraph::distances(stree_ori, x[1], x[2]))
    valid_edges <- added_edges[which(edge_dists >= minimal_branch_len), ,
                               drop = FALSE]
    edges_to_remove_df <- added_edges[which(edge_dists < minimal_branch_len), ,
                                      drop = FALSE]
  }
  if(nrow(valid_edges) > 0) {
    vertex_top_keep <-
      as.character(unlist(apply(valid_edges, 1,
                                function(x) {
                                  igraph::shortest_paths(stree_ori,
                                                         x[1],
                                                         x[2])$vpath[[1]]$name
                                  })))
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
      curr_node <- tryCatch({mst_traversal$order[i]}, error = function(e) {NA})
      if (is.na(curr_node))
        next;
      curr_node_name <- igraph::V(stree_ori)[curr_node]$name

      if (is.na(mst_traversal$father[curr_node]) == FALSE){
        parent_node <- mst_traversal$father[curr_node]
        parent_node_name <- igraph::V(stree_ori)[parent_node]$name

        if (igraph::degree(stree_ori, v=parent_node_name) > 2){
          parent_neighbors <- igraph::neighbors(stree_ori, v=parent_node_name,
                                                mode = 'all')

          parent_neighbors_index <- sort(match(parent_neighbors$name,
                                               mst_traversal$order$name))
          parent_neighbors <- mst_traversal$order$name[parent_neighbors_index]

          tmp <- igraph::delete.edges(stree_ori, paste0(parent_node_name,
                                                        "|", parent_neighbors))
          tmp_decomposed <- igraph::decompose.graph(tmp)

          comp_a <- tmp_decomposed[unlist(lapply(tmp_decomposed, function(x) {
            parent_neighbors[2] %in% igraph::V(x)$name
          }))][[1]]

          comp_b <- tmp_decomposed[unlist(lapply(tmp_decomposed, function(x) {
            parent_neighbors[3] %in% igraph::V(x)$name
          }))][[1]]

          diameter_len_a <- igraph::diameter(comp_a) + 1
          diameter_len_b <- igraph::diameter(comp_b) + 1

          if(diameter_len_a < minimal_branch_len) {
            vertex_to_be_deleted <- c(vertex_to_be_deleted,
                                      igraph::V(comp_a)$name)
          }
          if(diameter_len_b < minimal_branch_len) {
            vertex_to_be_deleted <- c(vertex_to_be_deleted,
                                      igraph::V(comp_b)$name)
          }
        }
      }
    }
  }

  valid_vertex_to_be_deleted <- setdiff(vertex_to_be_deleted, vertex_top_keep)
  stree_loop_closure <- igraph::delete_vertices(stree_loop_closure,
                                                valid_vertex_to_be_deleted)

  tmp <- edges_to_remove_df[edges_to_remove_df[, 1] %in%
                              igraph::V(stree_loop_closure)$name &
                              edges_to_remove_df[, 2] %in%
                              igraph::V(stree_loop_closure)$name, ,
                            drop = FALSE]
  if(nrow(tmp) > 0) {
    edges_to_remove <- paste0(tmp[, 1], '|', tmp[, 2])
    stree_loop_closure <- igraph::delete.edges(stree_loop_closure,
                                               edges_to_remove)
  }

  return(igraph::get.adjacency(stree_loop_closure))
}

project2MST <- function(cds, Projection_Method, orthogonal_proj_tip = FALSE,
                        verbose, reduction_method, rge_res_Y){
  dp_mst <- principal_graph(cds)[[reduction_method]]
  Z <- t(reducedDims(cds)[[reduction_method]])
  Y <- rge_res_Y

  cds <- findNearestPointOnMST(cds, reduction_method, rge_res_Y)
  closest_vertex <- cds@principal_graph_aux[[
    reduction_method]]$pr_graph_cell_proj_closest_vertex

  closest_vertex_names <- colnames(Y)[closest_vertex[, 1]]

  tip_leaves <- names(which(igraph::degree(dp_mst) == 1))

  if(!is.function(Projection_Method)) {
    P <- Y[, closest_vertex]
  } else {
    P <- matrix(rep(0, length(Z)), nrow = nrow(Z)) #Y
    nearest_edges <- matrix(rep(0, length(Z[1:2, ])), ncol = 2)
    row.names(nearest_edges) <-  colnames(cds)
    for(i in 1:length(closest_vertex)) { # This loop is going to be slow
      neighbors <- names(igraph::neighborhood(dp_mst,
                                              nodes = closest_vertex_names[i],
                                              mode = 'all')[[1]])[-1]
      projection <- NULL
      distance <- NULL
      Z_i <- Z[, i]

      for(neighbor in neighbors) {
        if(closest_vertex_names[i] %in% tip_leaves) {
          if(orthogonal_proj_tip) {
            tmp <- projPointOnLine(Z_i, Y[, c(closest_vertex_names[i],
                                              neighbor)])
          } else {
            tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i],
                                                neighbor)])
          }
        } else {
          tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i],
                                              neighbor)])
        }
        if(any(is.na(tmp))) {
          tmp <- Y[, neighbor]
        }
        projection <- rbind(projection, tmp)
        distance <- c(distance, stats::dist(rbind(Z_i, tmp)))
      }
      if(class(projection)[1] != 'matrix') {
        projection <- as.matrix(projection)
      }

      which_min <- which.min(distance)

      P[, i] <- projection[which_min, ]
      nearest_edges[i, ] <- c(closest_vertex_names[i], neighbors[which_min])
    }
  }

  colnames(P) <- colnames(Z)

  dp_mst_list <- igraph::decompose.graph(dp_mst)
  dp_mst_df <- NULL
  partitions <- cds@clusters[[reduction_method]]$partitions

  # Sanity test.
  # If this fails there may be a problem with getting cell names from
  # cds@clusters[[reduction_method]]$partitions
  # although that seems utterly implausible. I am replacing
  #   subset_cds_col_names <- colnames(cds)[cds@clusters[[reduction_method]]$partitions == cur_partition]
  # in the loop with the partition element names as
  #   subset_cds_col_names <- names(partitions[partitions==cur_partition])
  # so that this functions works when learn_graph() is run with use_partition=FALSE.
  assertthat::assert_that( !is.null( names(cds@clusters[[reduction_method]]$partitions) ), msg='names(cds@clusters[[reduction_method]]$partitions) == NULL' )
  assertthat::assert_that( !( length( colnames(cds) ) != length(names(cds@clusters[[reduction_method]]$partitions))), msg='length( colnames(cds) ) != length(names(cds@clusters[[reduction_method]]$partitions))' )
  assertthat::assert_that( !any( colnames(cds)!=names(cds@clusters[[reduction_method]]$partitions) ), msg='colnames(cds)!=names(cds@clusters[[reduction_method]]$partitions)' )

  if(length(dp_mst_list) == 1 & length(unique(partitions)) > 1) {
    #
    # Adjust for condition learn_graph(use_partition=FALSE).
    #
    partitions[partitions!='1'] <- '1'
  }

  if(!is.null(partitions)) {
    for(cur_partition in sort(unique(partitions))) {
      data_df <- NULL

      if(verbose) {
        message('\nProjecting cells to principal points for partition: ',
                cur_partition)
      }

      subset_cds_col_names <- names(partitions[partitions==cur_partition])

      cur_z <- Z[, subset_cds_col_names]
      cur_p <- P[, subset_cds_col_names]

      if (ncol(cur_p) > 0 && nrow(cur_p) > 0){
        cur_centroid_name <-
          igraph::V(dp_mst_list[[as.numeric(cur_partition)]])$name

        cur_nearest_edges <- nearest_edges[subset_cds_col_names, ]
        data_df <- cbind(as.data.frame(t(cur_p)),
                         apply(cur_nearest_edges, 1, sort) %>% t())
        row.names(data_df) <- colnames(cur_p)
        colnames(data_df) <- c(paste0("P_", 1:nrow(cur_p)), 'source', 'target')

        # sort each cell's distance to the source in each principal edge group
        data_df$distance_2_source <-
          sqrt(colSums((cur_p - rge_res_Y[,as.character(data_df[, 'source'])])^2))
        data_df <- data_df %>% tibble::rownames_to_column() %>%
          dplyr::mutate(group = paste(source, target, sep = '_')) %>%
          dplyr::arrange(group, dplyr::desc(-distance_2_source))

        # add the links from the source to the nearest points belong to the
        # principal edge and also all following connections between those
        # points
        data_df <- data_df %>% dplyr::group_by(group) %>%
          dplyr::mutate(new_source = dplyr::lag(rowname), new_target = rowname)
        # use the correct name of the source point
        data_df[is.na(data_df$new_source), "new_source"] <-
          as.character(as.matrix(data_df[is.na(data_df$new_source), 'source']))

        # add the links from the last point on the principal edge to the target
        # point of the edge
        added_rows <- which(is.na(data_df$new_source) &
                              is.na(data_df$new_target))
        data_df <- as.data.frame(data_df, stringsAsFactors = FALSE)
        data_df <- as.data.frame(as.matrix(data_df), stringsAsFactors = FALSE)
        data_df[added_rows, c('new_source', 'new_target')] <-
          data_df[added_rows - 1, c('rowname', 'target')]

        # calculate distance between each pair
        aug_P = cbind(cur_p, rge_res_Y)
        data_df$weight <-  sqrt(colSums((aug_P[,
                                               data_df$new_source] -
                                           aug_P[, data_df$new_target]))^2)
        # add the minimal positive distance between any points to the distance
        # matrix
        data_df$weight <-
          data_df$weight + min(data_df$weight[data_df$weight > 0])

        # Calculate distance between two connected nodes directly from the
        # original graph
        edge_list <- as.data.frame(igraph::get.edgelist(dp_mst_list[[
          as.numeric(cur_partition)]]), stringsAsFactors=FALSE)
        dp <- as.matrix(stats::dist(t(rge_res_Y)[cur_centroid_name,]))
        edge_list$weight <- dp[cbind(edge_list[, 1], edge_list[, 2])]
        colnames(edge_list) <- c("new_source", "new_target", 'weight')

        dp_mst_df <- Reduce(rbind,
                            list(dp_mst_df, data_df[, c("new_source",
                                                        "new_target",
                                                        'weight')], edge_list))
      }
    }
  }
  dp_mst <- igraph::graph.data.frame(dp_mst_df, directed = FALSE)
  cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_tree <- dp_mst
  cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_dist <- P

  closest_vertex_df <- as.matrix(closest_vertex)
  row.names(closest_vertex_df) <- row.names(closest_vertex)
  cds@principal_graph_aux[[
    reduction_method]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df

  cds
}

# Project each point to the nearest on the MST:
findNearestPointOnMST <- function(cds, reduction_method, rge_res_Y){
  dp_mst <- principal_graph(cds)[[reduction_method]]
  dp_mst_list <- igraph::decompose.graph(dp_mst)

  if(length(unique(cds@clusters[[reduction_method]]$partitions)) !=
     length(dp_mst_list)) {
    dp_mst_list <- list(dp_mst)
  }

  closest_vertex_df <- NULL
  cur_start_index <- 0

  for(i in 1:length(dp_mst_list)) {
    cur_dp_mst <- dp_mst_list[[i]]

    if(length(dp_mst_list) == 1) {
      Z <- t(reducedDims(cds)[[reduction_method]])
    } else {
      Z <- t(reducedDims(cds)[[reduction_method]])[, cds@clusters[[
        reduction_method]]$partitions == i]
    }
    Y <- rge_res_Y[, igraph::V(cur_dp_mst)$name]

    tip_leaves <- names(which(igraph::degree(cur_dp_mst) == 1))

    closest_vertex_ori <- find_nearest_vertex(Z, Y)
    closest_vertex <- closest_vertex_ori + cur_start_index

    closest_vertex_names <- colnames(Y)[closest_vertex_ori]
    cur_name <- names(closest_vertex)
    closest_vertex <- as.matrix(closest_vertex)
    row.names(closest_vertex) <- cur_name #original cell names for projection
    closest_vertex_df <- rbind(closest_vertex_df, closest_vertex) #index on Z

    cur_start_index <- cur_start_index + igraph::vcount(cur_dp_mst)
  }
  closest_vertex_df <- closest_vertex_df[colnames(cds), , drop = FALSE]
  cds@principal_graph_aux[[
    reduction_method]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df
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


#' Function to automatically learn the structure of data by either using
#' L1-graph or the spanning-tree formulization
#' @param X the input data DxN
#' @param C0 the initialization of centroids
#' @param maxiter maximum number of iteration
#' @param eps relative objective difference
#' @param L1.gamma regularization parameter for k-means (the prefix of 'param'
#'   is used to avoid name collision with gamma)
#' @param L1.sigma bandwidth parameter
#' @param verbose emit results from iteration
#' @return a list of X, C, W, P, objs
#' X is the input data
#' C is the centers for principal graph
#' W is the principal graph matrix
#' P is the cluster assignment matrix
#' objs is the objective value for the function
calc_principal_graph <- function(X, C0,
                                 maxiter = 10,
                                 eps = 1e-5,
                                 L1.gamma = 0.5,
                                 L1.sigma = 0.01,
                                 verbose = TRUE) {

  C <- C0;
  K <- ncol(C)
  objs <- c()
  for(iter in 1:maxiter) {
    #this part calculates the cost matrix Phi
    norm_sq <- repmat(t(colSums(C^2)), K, 1)
    Phi <- norm_sq + t(norm_sq) - 2 * t(C) %*% C

    g <- igraph::graph.adjacency(Phi, mode = 'lower', diag = TRUE, weighted = TRUE)
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


#' function to reproduce the behavior of repmat function in matlab to replicate
#' and tile an matrix
#' @param X matrix for tiling and replicate the data
#' @param m a numeric value for tiling a matrix
#' @param n a numeric value for tiling a matrix
#' @return a matrix
repmat = function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=TRUE)
}



#' Function to calculate the third term in the objective function
#' @param X input data
#' @param C center of graph (D * K)
#' @param sigma bandwidth parameter
#' @return a matrix with diagonal element as 1 while other elements as zero
#'   (eye matrix)
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

  obj <- - sigma * sum( log( rowSums( exp(- dist_XC/sigma)) )
                        - min_dist/ sigma );

  return(list(P = P, obj = obj))
}

#' Function to reproduce the behavior of eye function in matlab
#' @param X input data
#' @param W the principal graph matrix
#' @param P the cluster assignment matrix
#' @param param.gamma regularization parameter for k-means (the prefix of
#'   'param' is used to avoid name collision with gamma)
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


connect_tips <- function(cds,
                         pd,
                         R, # kmean cluster
                         stree,
                         reducedDimK_old,
                         reducedDimS_old,
                         k = 25,
                         nn_control = nn_control,
                         weight = FALSE,
                         qval_thresh = 0.05,
                         kmean_res,
                         euclidean_distance_ratio = 1,
                         geodesic_distance_ratio = 1/3,
                         medioids,
                         verbose = FALSE) {
  reduction_method <- 'UMAP'

  if(is.null(row.names(stree)) & is.null(row.names(stree))) {
    dimnames(stree) <- list(paste0('Y_', 1:ncol(stree)),
                            paste0('Y_', 1:ncol(stree)))
  }

  stree <- as.matrix(stree)
  stree[stree != 0] <- 1
  mst_g_old <- igraph::graph_from_adjacency_matrix(stree, mode = 'undirected')
  if(is.null(kmean_res)) {
    tmp <- matrix(apply(R, 1, which.max))

    row.names(tmp) <- colnames(reducedDimS_old)

    tip_pc_points <- which(igraph::degree(mst_g_old) == 1)

    data <- t(reducedDimS_old[, ])

    cluster_result <- louvain_clustering(data=data,
                                         pd=pd[, ],
                                         weight=weight,
                                         nn_index=NULL,
                                         k=k,
                                         nn_control=nn_control,
                                         louvain_iter=1,
                                         random_seed=0L,
                                         verbose=verbose)
    cluster_result$optim_res$membership <- tmp[, 1]
  } else { # use kmean clustering result
    tip_pc_points <- which(igraph::degree(mst_g_old) == 1)
    tip_pc_points_kmean_clusters <-
      sort(kmean_res$cluster[names(tip_pc_points)])

    data <- t(reducedDimS_old[, ]) # raw_data_tip_pc_points

    cluster_result <- louvain_clustering(data=data,
                                         pd=pd[row.names(data), ],
                                         weight=weight,
                                         nn_index=NULL,
                                         k=k,
                                         nn_control=nn_control,
                                         louvain_iter=1,
                                         random_seed=random_seed,
                                         verbose=verbose)
    cluster_result$optim_res$membership <- kmean_res$cluster
  }

  # identify edges between only tip cells
  cluster_graph_res <- compute_partitions(cluster_result$g, cluster_result$optim_res,
                                          qval_thresh=qval_thresh,
                                          verbose = verbose)
  dimnames(cluster_graph_res$cluster_mat) <-
    dimnames(cluster_graph_res$num_links)
  valid_connection <- which(cluster_graph_res$cluster_mat < qval_thresh,
                            arr.ind = TRUE)
  valid_connection <- valid_connection[apply(valid_connection, 1,
                                             function(x) {
                                               all(x %in% tip_pc_points)
                                               }), ] # only the tip cells

  # prepare the PAGA graph
  G <- cluster_graph_res$cluster_mat
  G[cluster_graph_res$cluster_mat < qval_thresh] <- -1
  G[cluster_graph_res$cluster_mat > 0] <- 0
  G <- - G

  # if no connection based on PAGA (only existed in simulated data), use the
  # kNN graph instead
  if(all(G == 0, na.rm = TRUE)) {
    return(list(stree = igraph::get.adjacency(mst_g_old),
                Y = reducedDimK_old, G = G))
  }

  if(nrow(valid_connection) == 0) {
    return(list(stree = igraph::get.adjacency(mst_g_old),
                Y = reducedDimK_old, G = G))
  }

  # calculate length of the MST diameter path
  mst_g <- mst_g_old
  diameter_dis <- igraph::diameter(mst_g_old)
  reducedDimK_df <- reducedDimK_old

  pb4 <- utils::txtProgressBar(max = length(nrow(valid_connection)), file = "",
                        style = 3, min = 0)

  # find the maximum distance between nodes from the MST
  res <- stats::dist(t(reducedDimK_old))
  g <- igraph::graph_from_adjacency_matrix(as.matrix(res), weighted = TRUE,
                                           mode = 'undirected')
  mst <- igraph::minimum.spanning.tree(g)
  max_node_dist <- max(igraph::E(mst)$weight)

  # append new edges to close loops in the spanning tree returned from
  # SimplePPT
  for(i in 1:nrow(valid_connection)) {
    # cluster id for the tip point; if kmean_res return valid_connection[i, ]
    # is itself; otherwise the id identified in the tmp file
    edge_vec <- sort(unique(cluster_result$optim_res$membership))[
      valid_connection[i, ]]
    edge_vec_in_tip_pc_point <- igraph::V(mst_g_old)$name[edge_vec]

    if(length(edge_vec_in_tip_pc_point) == 1) next;
    if(all(edge_vec %in% tip_pc_points) &
       (igraph::distances(mst_g_old, edge_vec_in_tip_pc_point[1],
                          edge_vec_in_tip_pc_point[2]) >=
        geodesic_distance_ratio * diameter_dis) &
       (euclidean_distance_ratio * max_node_dist >
        stats::dist(t(reducedDimK_old[, edge_vec]))) ) {
      if(verbose) message('edge_vec is ', edge_vec[1], '\t', edge_vec[2])
      if(verbose) message('edge_vec_in_tip_pc_point is ',
                          edge_vec_in_tip_pc_point[1], '\t',
                          edge_vec_in_tip_pc_point[2])

      mst_g <- igraph::add_edges(mst_g, edge_vec_in_tip_pc_point)
    }
    utils::setTxtProgressBar(pb = pb4, value = pb4$getVal() + 1)
  }

  close(pb4)

  list(stree = igraph::get.adjacency(mst_g), Y = reducedDimK_df, G = G)
}



