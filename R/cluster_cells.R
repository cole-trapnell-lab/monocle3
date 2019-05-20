#' Cluster cells using Louvain community detection
#'
#' Unsupervised clustering of cells is a common step in many single-cell
#' expression workflows. In an experiment containing a mixture of cell types,
#' each cluster might correspond to a different cell type. This method takes a
#' cell_data_set as input and clusters them using Louvain community detection,
#' and returns a cell_data_set with cluster assignments and Louvain partitions.
#' Cluster assignments can then be returned using the \code{get_clusters}
#' function and partition assignments can be returned using the
#' \code{get_partitions} function.
#'
#' cluster assignments stored in the colData table. When number of clusters is set
#' to NULL (num_clusters = NULL), the decision plot as introduced in the reference
#' will be plotted and the users are required to check the decision plot to select
#' the rho and delta to determine the number of clusters to cluster. When the dataset
#' is big, for example > 50 k, we recommend the user to use the Louvain clustering
#' algorithm which is inspired from phenograph paper. Note Louvain doesn't support the
#' num_cluster argument but the k (number of k-nearest neighbors) is relevant to the final
#' clustering number. The deafult implementation of Louvain clustering (when res is set to be NULL)
#' is based on the Rphenograph package but updated based on our requirement (for example,
#' changed the jaccard_coeff function as well as adding louvain_iter argument, etc.). We also
#' support setting resolutioin parameter when performing louvain clustering using the louvain package
#' from python. With different res values, the users can obtain different granularity of the data.
#'
#' @param cds the cell_data_set upon which to perform this operation
#' @param reduction_method
#' @param k number of kNN used in creating the k nearest neighbor graph for Louvain clustering. The number of kNN is related to the resolution of the clustering result, bigger number of kNN gives low resolution and vice versa. Default to be 20
#' @param louvain_iter Integer number of iterations used for Louvain clustering. The clustering result gives the largest modularity score will be used as the final clustering result.  Default to be 1. Note that if louvain_iter is large than 1, the `seed` argument will be ignored.
#' @param louvain_qval
#' @param weight A logic argument to determine whether or not we will use Jaccard coefficent for two nearest neighbors (based on the overlapping of their kNN) as the weight used for Louvain clustering. Default to be FALSE.
#' @param resolution
#' @param random_seed  the seed used by the random number generator in louvain-igraph package. This argument will be ignored if louvain_iter is larger than 1.
#' @param verbose Verbose A logic flag to determine whether or not we should print the running details.
#' @param ... Additional arguments passed to \code{\link{louvain_clustering}()}
#'
#' @return an updated cell_data_set object, in which phenoData contains values for Cluster for each cell
#' @references Rodriguez, A., & Laio, A. (2014). Clustering by fast search and find of density peaks. Science, 344(6191), 1492-1496. doi:10.1126/science.1242072
#' @references Vincent D. Blondel, Jean-Loup Guillaume, Renaud Lambiotte, Etienne Lefebvre: Fast unfolding of communities in large networks. J. Stat. Mech. (2008) P10008
#' @references Jacob H. Levine and et.al. Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like Cells that Correlate with Prognosis. Cell, 2015.
#'
#'
#' @export

cluster_cells <- function(cds,
                          reduction_method = c("UMAP", "tSNE", "PCA"),
                          k = 20,
                          louvain_iter = 1,
                          louvain_qval = 0.05,
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
  ## TO DO what is resolution?
  assertthat::assert_that(is.numeric(louvain_qval))
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
    cluster_graph_res <- compute_louvain_connected_components(louvain_res$g,
                                                              louvain_res$optim_res,
                                                              louvain_qval, verbose)
    louvain_component <- igraph::components(cluster_graph_res$cluster_g)$membership[louvain_res$optim_res$membership]
    names(louvain_component) <- row.names(reduced_dim_res)
    louvain_component <- as.factor(louvain_component)
  } else {
    louvain_component <- rep(1, nrow(colData(cds)))

  }

  cds@clusters[[reduction_method]] <- list(louvain_res = louvain_res,
                                             louvain_component = louvain_component,
                                             clusters = factor(igraph::membership(louvain_res$optim_res)))
  return(cds)
}




#' Cluster cells based on louvain community detection algorithm.
#'
#' @description This function is a wrapper of the louvain function from the python package (louvain-igraph, https://github.com/vtraag/louvain-igraph)
#' The following description is from the original package "This package implements the louvain algorithm in C++ and exposes it to python. It relies on (python-)igraph for it to function.
#' Besides the relative flexibility of the implementation, it also scales well, and can be run on graphs of millions of nodes (as long as they
#' can fit in memory). The core function is find_partition which finds the optimal partition using the louvain algorithm [1] for a number of
#' different methods. The methods currently implemented are (1) modularity [2], (2) Reichardt and Bornholdt's model using the configuration null
#' model and the Erdös-Rényi null model [3], (3) the constant Potts model (CPM) [4], (4) Significance [5], and finally (5) Surprise [6]. In
#' addition, it supports multiplex partition optimisation allowing community detection on for example negative links [7] or multiple time slices
#' [8]. It also provides some support for community detection on bipartite graphs. See the documentation for more information." Please see the github
#' above for the citations. Right now we only support CPMVertexPartition, RBConfigurationVertexPartition, RBERVertexPartition, ModularityVertexPartition
#' SignificanceVertexPartition and SurpriseVertexPartition partition methods.
#'
#' @param X the dataset upon which to perform umap dimension reduction
#' @param python_home The python home directory where umap is installed
#' @param partition_method character - either the default "CPMVertexPartition" or  "RBConfigurationVertexPartition" / "RBERVertexPartition".
#' @param initial_membership (list of int) – Initial membership for the partition. If None then defaults to a singleton partition.
#' @param weights (list of double, or edge attribute) – Weights of edges. Can be either an iterable or an edge attribute.
#' @param res (double) – Resolution parameter.
#' @param node_sizes  (list of int, or vertex attribute) – Sizes of nodes are necessary to know the size of communities in aggregate graphs. Usually this is set to 1 for all nodes, but in specific cases this could be changed.
#' @param random_seed  the seed used by the random number generator in louvain-igraph package
#' @param verbose bool (optional, default False)
#' @param return_all Whether to return all slots after louvain
#' @return The cluster id if return_all set to be FALSE, otherwise all slots from the louvain function
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
    stop('please pass the python home directory where louvain is installed with python_home argument!')
  }, finally = {
  })

  reticulate::source_python(paste(system.file(package="monocle3"), "louvain.py", sep="/"))
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
    res = list(membership = louvain_res$membership + 1, modularity = louvain_res$modularity)
    names(res$membership) = colnames(X)
    return(res)
  }
}


#' Function to run louvain clustering algorithm
#'
#' @param data low dimensional space used to perform graph clustering
#' @param pd the dataframe of the phenotype from the cell dataset (colData(cds))
#' @param k number of nearest neighbors used for Louvain clustering
#' @param weight whether or not to calculate the weight for each edge in the kNN graph
#' @param louvain_iter the number of iteraction for louvain clustering
#' @param resolution resolution of clustering result, specifiying the granularity of clusters.
#' Default to not use resolution and the standard igraph louvain clustering algorithm will be used.
#' @param random_seed  the seed used by the random number generator in louvain-igraph package
#' @param verbose Whether to emit verbose output during dimensionality reduction
#' @param ... extra arguments used to run louvain_R
#' @useDynLib monocle3, .registration = TRUE
#' @return a list with four elements (g (igraph object for the kNN graph), coord (coordinates of the graph with
#' layout_component, if the number of cells is less than 3000), edge_links (the data frame to plot the edges of
#' the igraph, if the number of cells is less than 3000) and optim_res (the louvain clustering result)).
louvain_clustering <- function(data,
                               pd,
                               k = 20,
                               weight = F,
                               louvain_iter = 1,
                               resolution = NULL, random_seed = 0L,
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
  membership_matrix <- Matrix::sparse.model.matrix( ~ cell_membership + 0)
  num_links <- Matrix::t(membership_matrix) %*% igraph::as_adjacency_matrix(g) %*% membership_matrix
  diag(num_links) <- 0
  louvain_modules <- levels(cell_membership)

  edges_per_module <- Matrix::rowSums(num_links)
  total_edges <- sum(num_links)

  theta <- (as.matrix(edges_per_module) / total_edges) %*% Matrix::t(edges_per_module / total_edges)
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


