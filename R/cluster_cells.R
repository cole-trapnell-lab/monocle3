#' Cluster cells using Louvain/Leiden community detection
#'
#' Unsupervised clustering of cells is a common step in many single-cell
#' expression workflows. In an experiment containing a mixture of cell types,
#' each cluster might correspond to a different cell type. This function takes
#' a cell_data_set as input, clusters the cells using Louvain/Leiden community
#' detection, and returns a cell_data_set with internally stored cluster
#' assignments. In addition to clusters this function calculates partitions,
#' which represent superclusters of the Louvain/Leiden communities that are
#' found using a kNN pruning method. Cluster assignments can be accessed using
#' the \code{\link{clusters}} function and partition assignments can be
#' accessed using the \code{\link{partitions}} function.
#'
#' @param cds The cell_data_set upon which to perform clustering.
#' @param reduction_method The dimensionality reduction method upon which to
#'   base clustering. Options are "UMAP", "tSNE", "PCA" and "LSI".
#' @param k Integer number of nearest neighbors to use when creating the k
#'   nearest neighbor graph for Louvain/Leiden clustering. k is related to the
#'   resolution of the clustering result, a bigger k will result in lower
#'   resolution and vice versa. Default is 20.
#' @param clustering_algorithm String indicating the clustering algorithm to be
#'   used. Options are "louvain", and "leiden". Default is "louvain".
#'   Resoultion parameter ignored if set to "louvain".
#' @param num_iter Integer number of iterations used for Louvain/Leiden
#'   clustering. The clustering result giving the largest modularity score will
#'   be used as the final clustering result. Default is 1. Note that if
#'   num_iter is greater than 1, the random_seed argument will be ignored.
#' @param partition_qval Numeric, the q-value cutoff to determine when to
#'   partition. Default is 0.05.
#' @param weight A logical argument to determine whether or not to use Jaccard
#'   coefficients for two nearest neighbors (based on the overlapping of their
#'   kNN) as the weight used for Louvain clustering. Default is FALSE.
#' @param resolution Parameter that controls the resolution of clustering. If
#'   NULL (Default), the parameter is determined automatically. Only used when
#'   clustering_algorithm is "leiden".
#' @param random_seed The seed used by the random number generator in
#'   louvain-igraph package. This argument will be ignored if num_iter is
#'   larger than 1.
#' @param verbose A logic flag to determine whether or not we should print the
#'   run details.
#' @param ... Extra arguments to pass to clustering algorithm.
#'
#' @return an updated cell_data_set object, with cluster and partition
#'   information stored internally and accessible using
#'   \code{\link{clusters}} and \code{\link{partitions}}
#' @references Rodriguez, A., & Laio, A. (2014). Clustering by fast search and
#'   find of density peaks. Science, 344(6191), 1492-1496.
#'   doi:10.1126/science.1242072
#' @references Vincent D. Blondel, Jean-Loup Guillaume, Renaud Lambiotte,
#'   Etienne Lefebvre: Fast unfolding of communities in large networks.
#'   J. Stat. Mech. (2008) P10008
#' @references Jacob H. Levine and et. al. Data-Driven Phenotypic Dissection of
#'   AML Reveals Progenitor-like Cells that Correlate with Prognosis.
#'   Cell, 2015.
#' @useDynLib monocle3, .registration = TRUE
#' @export

cluster_cells <- function(cds,
                          reduction_method = c("UMAP", "tSNE", "PCA", "LSI"),
                          k = 20,
                          clustering_algorithm = c('louvain','leiden'),
                          num_iter = 1,
                          partition_qval = 0.05,
                          weight = FALSE,
                          resolution = NULL,
                          random_seed = 0L,
                          verbose = F, ...) {

  reduction_method <- match.arg(reduction_method)
  clustering_algorithm <- match.arg(clustering_algorithm)

  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(is.character(reduction_method))
  assertthat::assert_that(assertthat::is.count(k))
  assertthat::assert_that(is.logical(weight))
  assertthat::assert_that(assertthat::is.count(num_iter))
  if(!is.null(resolution) & clustering_algorithm == "louvain") {
    message(paste("Resolution can only be used when clustering_algorithm is",
                  "'leiden', switching to leiden clustering"))
    clustering_algorithm = "leiden"
  }
  if(!is.null(resolution)) {
    assertthat::assert_that(is.numeric(resolution))
  }
  assertthat::assert_that(is.numeric(partition_qval))
  assertthat::assert_that(is.logical(verbose))
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste("No dimensionality reduction for",
                                      reduction_method, "calculated.",
                                      "Please run reduce_dimensions with",
                                      "reduction_method =", reduction_method,
                                      "before running cluster_cells"))

  reduced_dim_res <- reducedDims(cds)[[reduction_method]]

  if(verbose)
    message("Running ", clustering_algorithm, " clustering algorithm ...")

  if(clustering_algorithm=='louvain')
  {
    clust_res <- louvain_clustering(data = reduced_dim_res,
                                    pd = colData(cds),
                                    k = k,
                                    weight = weight,
                                    num_iter = num_iter,
                                    resolution = resolution,
                                    random_seed = random_seed,
                                    verbose = verbose, ...)
  } else if(clustering_algorithm=='leiden'){
    clust_res <- leiden_clustering(data = reduced_dim_res,
                                   pd = colData(cds),
                                   k = k,
                                   weight = weight,
                                   num_iter = num_iter,
                                   resolution = resolution,
                                   random_seed = random_seed,
                                   verbose = verbose, ...)
  }
  if(length(unique(clust_res$optim_res$membership)) > 1) {
    cluster_graph_res <- compute_partitions(clust_res$g,
                                            clust_res$optim_res,
                                            partition_qval, verbose)
    partitions <- igraph::components(
      cluster_graph_res$cluster_g)$membership[clust_res$optim_res$membership]
    names(partitions) <- row.names(reduced_dim_res)
    partitions <- as.factor(partitions)
  } else {
    partitions <- rep(1, nrow(colData(cds)))
  }
  clusters <- factor(igraph::membership(clust_res$optim_res))
  names(clusters) <- row.names(reduced_dim_res)
  cds@clusters[[reduction_method]] <- list(clust_res = clust_res,
                                           partitions = partitions,
                                           clusters = clusters)

  return(cds)
}

louvain_clustering <- function(data,
                               pd,
                               k = 20,
                               weight = F,
                               num_iter = 1,
                               resolution = NULL,
                               random_seed = 0L,
                               verbose = F, ...) {
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
    warning(paste("RANN counts the point itself, k must be smaller than\nthe",
                  "total number of points - 1 (all other points) - 1",
                  "(itself)!"))
  }
  if (verbose) {
    message("Run kNN based graph clustering starts:", "\n",
            "  -Input data of ", nrow(data), " rows and ", ncol(data),
            " columns", "\n", "  -k is set to ", k)
  }
  if (verbose) {
    cat("  Finding nearest neighbors...")
  }
  t1 <- system.time(tmp <- RANN::nn2(data, data, k +
                                       1, searchtype = "standard"))
  neighborMatrix <- tmp[[1]][, -1]
  distMatrix <- tmp[[2]][, -1]
  if (verbose) {
    cat(paste("DONE ~", t1[3], "s\n", " Compute jaccard coefficient between",
              "nearest-neighbor sets ..."))
  }
  t2 <- system.time(links <- jaccard_coeff(neighborMatrix,
                                           weight))
  if (verbose) {
    cat(paste("DONE ~", t2[3], "s\n", " Build undirected graph from the",
              "weighted links ..."))
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

  if(num_iter >= 2) {
    random_seed <- NULL
  }

  for (iter in 1:num_iter) {
    if(verbose) {
      cat("Running louvain iteration ", iter, "...\n")
    }
    if(!is.null(resolution)) {
      for(i in 1:length(resolution)) {
        cur_resolution <- resolution[i]
        louvain_args <- c(list(X = igraph::get.adjacency(g),
                               res = as.numeric(cur_resolution),
                               random_seed = random_seed,
                               verbose = verbose),
                          extra_arguments[names(extra_arguments) %in%
                                            c("python_home",
                                              "partition_method",
                                              "initial_membership", "weights",
                                              "node_sizes", 'return_all')])
        Q <- do.call(louvain_R, louvain_args)
        Qt <- max(Q$modularity)
        if(verbose) {
          message('Current iteration is ', iter, '; current resolution is ',
                  cur_resolution, '; Modularity is ', Qt,
                  '; Number of clusters are ', max(Q$membership))
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
    message('Maximal modularity is ', Qp, '; corresponding resolution is ',
            best_max_resolution)
  t_end <- Sys.time()
  if (verbose) {
    message("\nRun kNN based graph clustering DONE, totally takes ", t_end -
              t_start, " s.")
    cat("  -Number of clusters:",
        length(unique(igraph::membership(optim_res))), "\n")
  }

  if(igraph::vcount(g) < 3000) {

    coord <- NULL
    edge_links <- NULL
  } else {
    coord <- NULL
    edge_links <- NULL
  }

  igraph::V(g)$names <- as.character(igraph::V(g))
  return(list(g = g, relations = relations, distMatrix = distMatrix,
              coord = coord, edge_links = edge_links, optim_res = optim_res))
}


leiden_clustering <- function(data,
                              pd,
                              k = 20,
                              weight = NULL,
                              num_iter = 2,
                              resolution = 0.5,
                              random_seed = NULL,
                              verbose = FALSE, ...) {
  extra_arguments <- list(...)
  if( is.null( resolution ) )
  {
    resolution = 0.5
  }
  if( is.null( num_iter ) )
  {
    num_iter = 2
  }
  if( random_seed == 0L )
  {
    random_seed = NULL
  }

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
    warning(paste("RANN counts the point itself, k must be smaller than\nthe",
                  "total number of points - 1 (all other points) - 1",
                  "(itself)!"))
  }
  if (verbose) {
    message("Run kNN based graph clustering starts:", "\n",
            "  -Input data of ", nrow(data), " rows and ", ncol(data),
            " columns", "\n", "  -k is set to ", k)
  }
  if (verbose) {
    cat("  Finding nearest neighbors...")
  }
  t1 <- system.time(tmp <- RANN::nn2(data, data, k +
                                       1, searchtype = "standard"))
  neighborMatrix <- tmp[[1]][, -1]
  distMatrix <- tmp[[2]][, -1]
  if (verbose) {
    cat(paste("DONE ~", t1[3], "s\n", " Compute jaccard coefficient between",
              "nearest-neighbor sets ..."))
  }
  t2 <- system.time(links <- jaccard_coeff(neighborMatrix,
                                           weight))
  if (verbose) {
    cat(paste("DONE ~", t2[3], "s\n", " Build undirected graph from the",
              "weighted links ..."))
  }
  links <- links[links[, 1] > 0, ]
  relations <- as.data.frame(links)
  colnames(relations) <- c("from", "to", "weight")
  relations$from <- cell_names[relations$from]
  relations$to <- cell_names[relations$to]
  t3 <- system.time(g <- igraph::graph.data.frame(relations, directed = FALSE))
  if (verbose) {
    cat("DONE ~", t3[3], "s\n", " Run leiden clustering on the graph ...\n")
  }
  t_start <- Sys.time()
  Qp <- -1
  optim_res <- NULL
  best_max_resolution <- 'No resolution'

  if( 'partition_method' %in% names( extra_arguments ) )
    partition_method <- extra_arguments[['partition_method']]
  else
    partition_method <- NULL
  if( 'initial_membership' %in% names( extra_arguments ) )
    initial_membership <- extra_arguments[['initial_membership']]
  else
    initial_membership <- NULL
  if( 'weights' %in% names( extra_arguments ) )
    edge_weights <- extra_arguments[['weights']]
  else
    edge_weights <- NULL
  if( 'node_sizes' %in% names( extra_arguments ) )
    node_sizes <- extra_arguments[['node_sizes']]
  else
    node_sizes <- NULL

  for(i in 1:length(resolution)) {
    cur_resolution <- resolution[i]
    leiden_res <-
      leidenbase::leiden_find_partition(g,
                                        partition_type = partition_method,
                                        initial_membership = initial_membership,
                                        edge_weights = edge_weights,
                                        node_sizes = node_sizes,
                                        seed = random_seed,
                                        resolution_parameter = resolution,
                                        num_iter = num_iter,
                                        verbose = verbose )

    Qt <- leiden_res[['quality']]
    if(verbose) {
      message('Current resolution is ',
              cur_resolution, '; Modularity is ', Qt,
              '; Number of clusters are ', max(leiden_res[['membership']]))
    }
    if (Qt > Qp) {
      optim_res <- leiden_res
      Qp <- Qt
    }
  }
  if (is.null(optim_res)) {
    Qp <- leiden_res[['quality']]
    optim_res <- leiden_res
  }
  else {
    Qt <- leiden_res[['quality']]
    if (Qt > Qp) {
      optim_res <- leiden_res
      Qp <- Qt
    }
  }

  if(verbose)
    message('Maximal modularity is ', Qp, '; corresponding resolution is ',
            best_max_resolution)
  t_end <- Sys.time()
  if (verbose) {
    message("\nRun kNN based graph clustering DONE, totally takes ", t_end -
              t_start, " s.")
    cat("  -Number of clusters:",
        max(best_res[['membership']]), "\n")
  }

  if(igraph::vcount(g) < 3000) {
    coord <- NULL
    edge_links <- NULL
  } else {
    coord <- NULL
    edge_links <- NULL
  }
  igraph::V(g)$names <- as.character(igraph::V(g))
  out_res <- list(membership = optim_res[['membership']],
                  modularity = optim_res[['quality']] )
  names(out_res$membership) = cell_names
  return(list(g = g, relations = relations, distMatrix = distMatrix,
              coord = coord, edge_links = edge_links, optim_res = out_res))
}


compute_partitions <- function(g,
                               optim_res,
                               qval_thresh=0.05,
                               verbose = FALSE){
  cell_membership <- as.factor(igraph::membership(optim_res))
  membership_matrix <- Matrix::sparse.model.matrix( ~ cell_membership + 0)
  num_links <- Matrix::t(membership_matrix) %*%
    igraph::as_adjacency_matrix(g) %*% membership_matrix
  diag(num_links) <- 0
  louvain_modules <- levels(cell_membership)

  edges_per_module <- Matrix::rowSums(num_links)
  total_edges <- sum(num_links)

  theta <- (as.matrix(edges_per_module) / total_edges) %*%
    Matrix::t(edges_per_module / total_edges)
  var_null_num_links <- theta * (1 - theta) / total_edges
  num_links_ij <- num_links / total_edges - theta
  cluster_mat <- pnorm_over_mat(as.matrix(num_links_ij), var_null_num_links)

  num_links <- num_links_ij / total_edges

  cluster_mat <- matrix(stats::p.adjust(cluster_mat),
                        nrow=length(louvain_modules),
                        ncol=length(louvain_modules))

  sig_links <- as.matrix(num_links)
  sig_links[cluster_mat > qval_thresh] = 0
  diag(sig_links) <- 0

  cluster_g <- igraph::graph_from_adjacency_matrix(sig_links, weighted = T,
                                                   mode = 'undirected')

  list(cluster_g = cluster_g, num_links = num_links, cluster_mat = cluster_mat)
}
