#' This function tries to partition cells into different graphs based on a similar approach proposed by Alex Wolf and colleagues
#'
#' Recently Alex Wolf and colleague first proposed the idea to represent the data with an "abstract partition graph"
#' of clusters identified by Louvain clustering by simply connecting significantly overlapping Louvain clusters (Wolf et al. 2017).
#' Similar methods for "abstract partition graph" are also recently developed and applied in analyzing the zebrafish / frog cell
#' atlas datasets (Wagner et al. 2018; Briggs et al. 2018). This coarse-graining representation of the data address a few limitations
#' of tree-based trajectory inference algorithms, for example, the default principal tree learning algorithm (DDRTree) in Monocle 2.
#' Although the particion graph doesn't learn an explicit simplified principal tree as DDRTree in Monocle 2, it can naturally separate
#' outlier cell groups and potentially also parallel trajectories while DDRTree often requires pre-processing before hand to robustly
#' reconstruct trajectory and cannot handle non-tree like structure. Instead of directly learn a coarse-graining graph of clusters,
#' we instead take advantage of the participation graph and use it as merely a heuristic initial condition for L1-graph algorithm
#' (as explained in the next section) to learn the principal points and principal graph at the same time directly from the reduced
#' UMAP data space. In contrast to the cluster participation method, the principal graph learnt provides an abstraction of the data
#' manifold while also preserves the local information from the original data space as it is directly embedded in the original data space.
#' In Monocle 3, we use the clustering_louvain function from the igraph package to perform community detection and implemented an efficient
#' version of "abstract partition graph" from Alex Wolf. Basically, we first create a design matrix \eqn{X}{} representing the allocation of
#' each cell to a particular louvain cluster. The column of \eqn{X}{} represents a louvain cluster while the row of \eqn{X}{} a particular cell.
#' \eqn{X_{ij} = 1}{} if cell \eqn{i}{} belongs to cluster \eqn{j}{}, otherwise 0. We can further obtain the adjacency matrix \eqn{A}{} of the kNN graph
#' used to perform the louvain clustering where \eqn{A_{ij} = 1}{} if cell \eqn{i}{} connects to \eqn{j}{} in the kNN graph. Then the connection
#' matrix \eqn{M}{} between each cluster is calculated as, \eqn{M \cdot {X'}{A}{X}}{}. Once \eqn{M}{} is constructed, we can then follow
#' Supplemental Note 3.1 from (Wolf et al. 2017) to calculate the significance of the connection between each louvain clustering and
#' consider any clusters with p-value larger than 0.05 by default as not disconnected.
#'
#' @param cds the cell_data_set upon which to perform this operation
#' @param partition_names Which partition groups (column in the pData) should be used to calculate the connectivity between partitions
#' @param use_pca Whether or not to cluster cells based on top PCA component. Default to be FALSE.
#' @param k number of nearest neighbors used for Louvain clustering (pass to louvain_clustering function)
#' @param weight whether or not to calculate the weight for each edge in the kNN graph (pass to louvain_clustering function)
#' @param louvain_iter the number of iteraction for louvain clustering (pass to louvain_clustering function)
#' @param resolution resolution of clustering result, specifiying the granularity of clusters.
#' Default to not use resolution and the standard igraph louvain clustering algorithm will be used.
#' @param louvain_qval The q-val threshold used to determine the partition of cells (pass to compute_louvain_connected_components)
#' @param return_all Whether to return all saved objects from compute_louvain_connected_components function.
#' @param verbose Whether to emit verbose output during louvain clustering
#' @param ... additional arguments to pass to the smoothEmbedding function
#' @return an updated cell_data_set object
#'
#' @export
partition_cells <- function(cds,
                           partition_names = NULL,
                           use_pca = FALSE,
                           k = 20,
                           weight = F,
                           louvain_iter = 1,
                           resolution = NULL,
                           louvain_qval = 0.05,
                           return_all = FALSE,
                           verbose = FALSE, ...){
  extra_arguments <- list(...)
  irlba_pca_res <- cds@normalized_data_projection
  if(nrow(irlba_pca_res) == 0)
    stop("No normalized data projection. Please run preprocess_cds before running partition_cells")

  Y <- reducedDimS(cds)
  reduced_dim_res = Y

  if(verbose)
    message("Running louvain clustering algorithm ...")
  #row.names(umap_res) <- colnames(FM)
  if(nrow(Y) == 0) {
    reduced_dim_res <- t(irlba_pca_res)
  }

  if(use_pca) {
    reduced_dim_res <- t(cds@normalized_data_projection)
  }

  if(is.null(partition_names)) {
    louvain_clustering_args <- c(list(data = t(reduced_dim_res), pd = pData(cds)[row.names(irlba_pca_res), ], k = k,
                                      resolution = resolution, weight = weight, louvain_iter = louvain_iter, verbose = verbose)) # , extra_arguments[names(extra_arguments) %in% c("k", "weight", "louvain_iter")]
    louvain_res <- do.call(louvain_clustering, louvain_clustering_args)

    if(length(unique(louvain_res$optim_res$membership)) == 1) {
      pData(cds)$louvain_component <- 1

      return(cds)
    }

  } else {
    build_asym_kNN_graph_args <- c(list(data = t(reduced_dim_res), k = k, return_graph = T),
                                   extra_arguments[names(extra_arguments) %in% c('dist_type', 'return_graph')])
    louvain_res <- list(g = do.call(build_asym_kNN_graph, build_asym_kNN_graph_args), optim_res = list(membership = NULL))
  }


  if(!is.null(partition_names)) {
    if(!(partition_names %in% colnames(pData(cds)))) {
      stop(paste0('Error: please make sure pData has a column with the name ', partition_names))
    }
    if(partition_names %in% colnames(pData(cds))) {
      louvain_res$optim_res$membership <- pData(cds)[, partition_names]
    }
  }

  cluster_graph_res <- compute_louvain_connected_components(louvain_res$g, louvain_res$optim_res, louvain_qval, verbose)
  louvain_component = igraph::components(cluster_graph_res$cluster_g)$membership[louvain_res$optim_res$membership]
  names(louvain_component) = row.names(irlba_pca_res)
  louvain_component = as.factor(louvain_component)
  pData(cds)$louvain_component <- louvain_component

  cds@aux_clustering_data$partition_cells <- list(cluster_graph_res = cluster_graph_res, louvain_res = louvain_res)

  if(return_all) {
    return(list(cds = cds, cluster_graph_res = cluster_graph_res))
  } else {
    return(cds)
  }
}

#' Function to run louvain clustering algorithm
#'
#' @param data low dimensional space used to perform graph clustering
#' @param pd the dataframe of the phenotype from the cell dataset (pData(cds))
#' @param k number of nearest neighbors used for Louvain clustering
#' @param weight whether or not to calculate the weight for each edge in the kNN graph
#' @param louvain_iter the number of iteraction for louvain clustering
#' @param resolution resolution of clustering result, specifiying the granularity of clusters.
#' Default to not use resolution and the standard igraph louvain clustering algorithm will be used.
#' @param random_seed  the seed used by the random number generator in louvain-igraph package
#' @param verbose Whether to emit verbose output during dimensionality reduction
#' @param ... extra arguments used to run louvain_R
#' @useDynLib monocle3, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @return a list with four elements (g (igraph object for the kNN graph), coord (coordinates of the graph with
#' layout_component, if the number of cells is less than 3000), edge_links (the data frame to plot the edges of
#' the igraph, if the number of cells is less than 3000) and optim_res (the louvain clustering result)).
louvain_clustering <- function(data, pd, k = 20, weight = F, louvain_iter = 1, resolution = NULL, random_seed = 0L, verbose = F, ...) {
  extra_arguments <- list(...)
  cell_names <- row.names(pd)
  if(!identical(cell_names, row.names(pd)))
    stop("phenotype and row name from the data doesn't match")

  if (is.data.frame(data))
    data <- as.matrix(data)
  if (!is.matrix(data))
    stop("Wrong input data, should be a data frame of matrix!")
  if (k < 1) {
    stop("k must be a positive integer!")
  } else if (k > nrow(data) - 2) {
    k <- nrow(data) - 2
    warning("RANN counts the point itself, k must be smaller than\nthe total number of points - 1 (all other points) - 1 (itself)!")
  }
  if (verbose) {
    message("Run kNN based graph clustering starts:", "\n", "  -Input data of ",
            nrow(data), " rows and ", ncol(data), " columns",
            "\n", "  -k is set to ", k)
  }
  if (verbose) {
    cat("  Finding nearest neighbors...")
  }
  t1 <- system.time(tmp <- RANN::nn2(data, data, k +
                                       1, searchtype = "standard"))
  neighborMatrix <- tmp[[1]][, -1]
  distMatrix <- tmp[[2]][, -1]
  if (verbose) {
    cat("DONE ~", t1[3], "s\n", " Compute jaccard coefficient between nearest-neighbor sets ...")
  }
  t2 <- system.time(links <- jaccard_coeff(neighborMatrix,
                                           weight))
  if (verbose) {
    cat("DONE ~", t2[3], "s\n", " Build undirected graph from the weighted links ...")
  }
  links <- links[links[, 1] > 0, ]
  relations <- as.data.frame(links)
  colnames(relations) <- c("from", "to", "weight")

  relations$from <- cell_names[relations$from]
  relations$to <- cell_names[relations$to]
  t3 <- system.time(g <- igraph::graph.data.frame(relations, directed = FALSE))
  if (verbose) {
    cat("DONE ~", t3[3], "s\n", " Run louvain clustering on the graph ...\n")
  }
  t_start <- Sys.time()
  Qp <- -1
  optim_res <- NULL
  best_max_resolution <- 'No resolution'

  if(louvain_iter >= 2) {
    random_seed <- NULL
  }

  for (iter in 1:louvain_iter) {
    if(verbose) {
      cat("Running louvain iteration ", iter, "...\n")
    }
    if(!is.null(resolution)) {
      for(i in 1:length(resolution)) {
        cur_resolution <- resolution[i]
        louvain_args <- c(list(X = igraph::get.adjacency(g), res = as.numeric(cur_resolution), random_seed = random_seed, verbose = verbose),
                          extra_arguments[names(extra_arguments) %in%
                                            c("python_home", "partition_method", "initial_membership", "weights", "node_sizes", 'return_all')])
        Q <- do.call(louvain_R, louvain_args)
        Qt <- max(Q$modularity)
        if(verbose) {
          message('Current iteration is ', iter, '; current resolution is ', cur_resolution, '; Modularity is ', Qt, '; Number of clusters are ', max(Q$membership))
        }
        if (Qt > Qp) {
          optim_res <- Q
          Qp <- Qt
          best_max_resolution <- cur_resolution
        }
      }
    } else {
      Q <- igraph::cluster_louvain(g)
    }
    if (is.null(optim_res)) {
      Qp <- max(Q$modularity)
      optim_res <- Q
    }
    else {
      Qt <- max(Q$modularity)
      if (Qt > Qp) {
        optim_res <- Q
        Qp <- Qt
      }
    }
  }
  if(verbose)
    message('Maximal modularity is ', Qp, '; corresponding resolution is ', best_max_resolution)
  t_end <- Sys.time()
  if (verbose) {
    message("\nRun kNN based graph clustering DONE, totally takes ", t_end -
              t_start, " s.")
    cat("  -Number of clusters:", length(unique(igraph::membership(optim_res))), "\n")
  }

  if(igraph::vcount(g) < 3000) {

    coord <- NULL
    edge_links <- NULL
  } else {
    coord <- NULL
    edge_links <- NULL
  }

  igraph::V(g)$names <- as.character(igraph::V(g))
  return(list(g = g, relations = relations, distMatrix = distMatrix, coord = coord, edge_links = edge_links, optim_res = optim_res))
}

compute_louvain_connected_components <- function(g, optim_res, qval_thresh=0.05, verbose = FALSE){
  cell_membership <- as.factor(igraph::membership(optim_res))
  membership_matrix = sparse.model.matrix( ~ cell_membership + 0)
  num_links = t(membership_matrix) %*% igraph::as_adjacency_matrix(g) %*% membership_matrix
  diag(num_links) = 0
  louvain_modules = levels(cell_membership)

  cluster_mat <- matrix(0, nrow = length(louvain_modules), ncol = length(louvain_modules)) # a matrix storing the overlapping clusters between louvain clusters which is based on the spanning tree
  enrichment_mat <- matrix(0, nrow = length(louvain_modules), ncol = length(louvain_modules)) # a matrix storing the overlapping clusters between louvain clusters which is based on the spanning tree

  overlapping_threshold <- 1e-5

  edges_per_module = rowSums(num_links)
  total_edges = sum(num_links)

  theta <- (as.matrix(edges_per_module) / total_edges) %*% t(edges_per_module / total_edges)
  var_null_num_links <- theta * (1 - theta) / total_edges
  num_links_ij <- num_links / total_edges - theta
  cluster_mat <- pnorm_over_mat(as.matrix(num_links_ij), var_null_num_links) # much faster c++ version

  enrichment_mat <- num_links_ij
  num_links <- num_links_ij / total_edges

  cluster_mat = matrix(p.adjust(cluster_mat), nrow=length(louvain_modules), ncol=length(louvain_modules))

  sig_links <- as.matrix(num_links)
  sig_links[cluster_mat > qval_thresh] = 0
  diag(sig_links) = 0

  cluster_g <- igraph::graph_from_adjacency_matrix(sig_links, weighted = T, mode = 'undirected')
  louvain_modules <- igraph::cluster_louvain(cluster_g)

  # return also the layout coordinates and the edges link for the graph of clusters

  coord <- igraph::layout_components(cluster_g)
  coord <- as.data.frame(coord)
  colnames(coord) <- c('x', 'y')
  row.names(coord) <- 1:nrow(coord)
  coord$Cluster <- 1:nrow(coord)
  coord$louvain_cluster <- as.character(igraph::membership(louvain_modules))

  edge_links <- NULL
  if(length(igraph::E(cluster_g)) > 0) { # run this only when there is edges
    edge <- igraph::get.data.frame(cluster_g)
    edge <- as.data.frame(edge)
    colnames(edge) <- c('start', 'end', 'weight')
    edge_links <- cbind(coord[edge$start, 1:2], coord[edge$end, 1:2])
    edge_links <- as.data.frame(edge_links)
    colnames(edge_links) <- c('x_start', 'x_end', 'y_start', 'y_end')
    edge_links$weight <- edge[, 3]
  }

  list(cluster_g = cluster_g, cluster_optim_res = optim_res, num_links = num_links, cluster_mat = cluster_mat, enrichment_mat = enrichment_mat, cluster_coord = coord, edge_links = edge_links)
}

