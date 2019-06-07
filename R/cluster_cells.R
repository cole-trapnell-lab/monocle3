#' Cluster cells using Louvain community detection
#'
#' Unsupervised clustering of cells is a common step in many single-cell
#' expression workflows. In an experiment containing a mixture of cell types,
#' each cluster might correspond to a different cell type. This function takes
#' a cell_data_set as input, clusters the cells using Louvain community
#' detection, and returns a cell_data_set with internally stored cluster
#' assignments. In addition to clusters this function calculates partitions,
#' which represent superclusters of the Louvain communities that are found
#' using a KNN pruning method. Cluster assignments can be accessed using the
#' \code{\link{clusters}} function and partition assignments can be
#' accessed using the \code{\link{partitions}} function.
#'
#' @param cds The cell_data_set upon which to perform clustering.
#' @param reduction_method The dimensionality reduction method upon which to
#'   base clustering. Options are "UMAP", "tSNE" and "PCA".
#' @param k Integer number of nearest neighbors to use when creating the k
#'   nearest neighbor graph for Louvain clustering. k is related to the
#'   resolution of the clustering result, a bigger k will result in lower
#'   resolution and vice versa. Default is 20.
#' @param louvain_iter Integer number of iterations used for Louvain
#'   clustering. The clustering result giving the largest modularity score will
#'   be used as the final clustering result. Default is 1. Note that if
#'   louvain_iter is greater than 1, the random_seed argument will be ignored.
#' @param partition_qval Numeric, the q-value cutoff to determine when to
#'   partition. Default is 0.05.
#' @param weight A logical argument to determine whether or not to use Jaccard
#'   coefficents for two nearest neighbors (based on the overlapping of their
#'   kNN) as the weight used for Louvain clustering. Default is FALSE.
#' @param resolution Parameter that controls the resolution of clustering. If
#'   NULL (Default), the parameter is determined automatically.
#' @param random_seed The seed used by the random number generator in
#'   louvain-igraph package. This argument will be ignored if louvain_iter is
#'   larger than 1.
#' @param verbose A logic flag to determine whether or not we should print the
#'   run details.
#' @param ... Additional arguments passed to louvain-igraph Python package.
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
#' @references Jacob H. Levine and et.al. Data-Driven Phenotypic Dissection of
#'   AML Reveals Progenitor-like Cells that Correlate with Prognosis.
#'   Cell, 2015.
#' @useDynLib monocle3, .registration = TRUE
#' @export

cluster_cells <- function(cds,
                          reduction_method = c("UMAP", "tSNE", "PCA"),
                          k = 20,
                          louvain_iter = 1,
                          partition_qval = 0.05,
                          weight = FALSE,
                          resolution = NULL,
                          random_seed = 0L,
                          verbose = F,
                          ...) {
  method = 'louvain'
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(is(cds, "cell_data_set"))
  assertthat::assert_that(is.character(reduction_method))
  assertthat::assert_that(assertthat::is.count(k))
  assertthat::assert_that(is.logical(weight))
  assertthat::assert_that(assertthat::is.count(louvain_iter))
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
    message("Running louvain clustering algorithm ...")

  louvain_res <- louvain_clustering(data = reduced_dim_res,
                                    pd = colData(cds),
                                    k = k,
                                    weight = weight,
                                    louvain_iter = louvain_iter,
                                    resolution = resolution,
                                    random_seed = random_seed,
                                    verbose = verbose, ...)
  if(length(unique(louvain_res$optim_res$membership)) > 1) {
    cluster_graph_res <- compute_partitions(louvain_res$g,
                                            louvain_res$optim_res,
                                            partition_qval, verbose)
    partitions <- igraph::components(
      cluster_graph_res$cluster_g)$membership[louvain_res$optim_res$membership]
    names(partitions) <- row.names(reduced_dim_res)
    partitions <- as.factor(partitions)
  } else {
    partitions <- rep(1, nrow(colData(cds)))
  }
  clusters <- factor(igraph::membership(louvain_res$optim_res))
  cds@clusters[[reduction_method]] <- list(louvain_res = louvain_res,
                                           partitions = partitions,
                                           clusters = clusters)
  return(cds)
}




#' Cluster cells based on louvain community detection algorithm.
#'
#' @description This function is a wrapper of the louvain function from the
#' python package (louvain-igraph, https://github.com/vtraag/louvain-igraph)
#' The following description is from the original package "This package
#' implements the louvain algorithm in C++ and exposes it to python. It relies
#' on (python-)igraph for it to function. Besides the relative flexibility of
#' the implementation, it also scales well, and can be run on graphs of
#' millions of nodes (as long as they can fit in memory). The core function is
#' find_partition which finds the optimal partition using the louvain algorithm
#' [1] for a number of different methods. The methods currently implemented are
#' (1) modularity [2], (2) Reichardt and Bornholdt's model using the
#' configuration null model and the Erdös-Rényi null model [3], (3) the
#' constant Potts model (CPM) [4], (4) Significance [5], and finally (5)
#' Surprise [6]. In addition, it supports multiplex partition optimisation
#' allowing community detection on for example negative links [7] or multiple
#' time slices [8]. It also provides some support for community detection on
#' bipartite graphs. See the documentation for more information." Please see
#' the github above for the citations. Right now we only support
#' CPMVertexPartition, RBConfigurationVertexPartition, RBERVertexPartition,
#' ModularityVertexPartition SignificanceVertexPartition and
#' SurpriseVertexPartition partition methods.
#'
#' @param X the dataset upon which to perform louvain-igraph
#' @param python_home The python home directory where louvain-igraph is
#'   installed
#' @param partition_method character - either the default "CPMVertexPartition"
#'   or "RBConfigurationVertexPartition" / "RBERVertexPartition".
#' @param initial_membership (list of int) – Initial membership for the
#'   partition. If None then defaults to a singleton partition.
#' @param weights (list of double, or edge attribute) – Weights of edges. Can
#'   be either an iterable or an edge attribute.
#' @param res (double) – Resolution parameter.
#' @param node_sizes  (list of int, or vertex attribute) – Sizes of nodes are
#'   necessary to know the size of communities in aggregate graphs. Usually
#'   this is set to 1 for all nodes, but in specific cases this could be
#'   changed.
#' @param random_seed the seed used by the random number generator in
#'   louvain-igraph package
#' @param verbose bool (optional, default False)
#' @param return_all Whether to return all slots after louvain
#' @return The cluster id if return_all set to be FALSE, otherwise all slots
#'   from the louvain function
#' @encoding UTF-8
#'
louvain_R <- function(X, python_home = system('which python', intern = TRUE),
                      partition_method = 'CPMVertexPartition',
                      initial_membership = NULL,
                      weights = NULL,
                      res = 0.6,
                      node_sizes = NULL,
                      random_seed = 0L,
                      verbose = FALSE,
                      return_all = FALSE) {

  reticulate::use_python(python_home)

  tryCatch({
    reticulate::import("louvain")
  }, warning = function(w) {
  }, error = function(e) {
    print (e)
    stop(paste("Could not find louvain Python package. Please pass the python",
               "home directory where louvain is installed with python_home",
               "argument."))
  }, finally = {
  })

  reticulate::source_python(paste(system.file(package="monocle3"),
                                  "louvain.py", sep="/"))
  # X <- Matrix::t(X)
  if(length(grep('Matrix', class(X))) == 0){
    X <- as(as.matrix(X), 'TsparseMatrix')
  } else {
    X <- as(X, 'TsparseMatrix')
  }

  i <- as.integer(X@i)
  j <- as.integer(X@j)
  val <- X@x

  dim <- as.integer(X@Dim)

  if(is.null(partition_method) == F) {
    partition_method <- as.character(partition_method)
  }
  if(!is.null(random_seed)) {
    random_seed <- as.integer(random_seed)
  }

  louvain_res <- louvain(i, j, val, dim,
                         as.character(partition_method),
                         initial_membership,
                         weights,
                         as.numeric(res),
                         node_sizes,
                         random_seed,
                         as.logical(verbose))
  if(return_all) {
    return(louvain_res)
  } else {
    res = list(membership = louvain_res$membership + 1,
               modularity = louvain_res$modularity)
    names(res$membership) = colnames(X)
    return(res)
  }
}

louvain_clustering <- function(data,
                               pd,
                               k = 20,
                               weight = F,
                               louvain_iter = 1,
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


