#' Cluster cells using Louvain/Leiden community detection
#'
#' Unsupervised clustering of cells is a common step in many single-cell
#' expression workflows. In an experiment containing a mixture of cell types,
#' each cluster might correspond to a different cell type. This function takes
#' a cell_data_set as input, clusters the cells using Louvain/Leiden community
#' detection, and returns a cell_data_set with internally stored cluster
#' assignments. In addition to clusters this function calculates partitions,
#' which represent superclusters of the Louvain/Leiden communities that are found
#' using a kNN pruning method. Cluster assignments can be accessed using the
#' \code{\link{clusters}} function and partition assignments can be
#' accessed using the \code{\link{partitions}} function.
#'
#' @param cds The cell_data_set upon which to perform clustering.
#' @param reduction_method The dimensionality reduction method upon which to
#'   base clustering. Options are "UMAP", "tSNE", "PCA" and "LSI".
#' @param k Integer number of nearest neighbors to use when creating the k
#'   nearest neighbor graph for Louvain/Leiden clustering. k is related to the
#'   resolution of the clustering result, a bigger k will result in lower
#'   resolution and vice versa. Default is 20.
#' @param cluster_method String indicating the clustering method to use.
#'   Options are "louvain" or "leiden". Default is "leiden". Resolution
#'   parameter is ignored if set to "louvain".
#' @param num_iter Integer number of iterations used for Louvain/Leiden
#'   clustering. The clustering result giving the largest modularity score will
#'   be used as the final clustering result. Default is 1. Note that if
#'   num_iter is greater than 1, the random_seed argument will be ignored
#'   for the louvain method.
#' @param partition_qval Numeric, the q-value cutoff to determine when to
#'   partition. Default is 0.05.
#' @param weight A logical argument to determine whether or not to use Jaccard
#'   coefficients for two nearest neighbors (based on the overlapping of their
#'   kNN) as the weight used for Louvain clustering. Default is FALSE.
#' @param resolution Parameter that controls the resolution of clustering. If
#'   NULL (Default), the parameter is determined automatically.
#' @param random_seed The seed used by the random number generator in
#'   louvain-igraph package. This argument will be ignored if num_iter is
#'   larger than 1.
#' @param verbose A logic flag to determine whether or not we should print the
#'   run details.
#' @param nn_control A list of parameters used to make the nearest
#'  neighbor index. See the set_nn_control help for detailed information.
#' @param ... Additional arguments passed to the leidenbase package.
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
#' @references V. A. Traag and L. Waltman and N. J. van Eck: From Louvain
#'   to Leiden: guaranteeing well-connected communities. Scientific Reports,
#'   9(1) (2019). doi: 10.1038/s41598-019-41695-z.
#' @references Jacob H. Levine and et. al. Data-Driven Phenotypic Dissection of
#'   AML Reveals Progenitor-like Cells that Correlate with Prognosis.
#'   Cell, 2015.
#' @useDynLib monocle3, .registration = TRUE
#' @export
cluster_cells <- function(cds,
                          reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
                          k = 20,
                          cluster_method = c('leiden', 'louvain'),
                          num_iter = 2,
                          partition_qval = 0.05,
                          weight = FALSE,
                          resolution = NULL,
                          random_seed = NULL,
                          verbose = F,
                          nn_control = list(),
                          ...) {
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'UMAP', 'tSNE', 'PCA', 'LSI', 'Aligned'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(cluster_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "cluster_method must be one of 'leiden', 'louvain'")
  cluster_method <- match.arg(cluster_method)

  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(is.character(reduction_method))
  assertthat::assert_that(is.logical(weight))
  assertthat::assert_that(assertthat::is.count(num_iter))
  assertthat::assert_that(assertthat::is.count(k))

  if (!is.null(resolution) & cluster_method == "louvain") {
    message(paste("Resolution can only be used when cluster_method is",
                  "'leiden'. Switching to leiden clustering."))
    cluster_method <- "leiden"
  }

  if(!is.null(resolution)) {
    assertthat::assert_that(is.numeric(resolution))
  }
  assertthat::assert_that(is.numeric(partition_qval))
  assertthat::assert_that(is.logical(verbose))
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste("No dimensionality reduction for",
                                      reduction_method, "calculated.",
                                      "Please run reduce_dimension with",
                                      "reduction_method =", reduction_method,
                                      "before running cluster_cells"))

  nn_control <- set_nn_control(mode=3,
                               nn_control=nn_control,
                               k=k,
                               nn_control_default=get_global_variable('nn_control_1'),
                               verbose=verbose)
  nn_method <- nn_control[['method']]

  # The nn index is made on the full reducedDims(cds)[[reduction_method]]
  # matrix so use/store the nn index object in the cds. This saves nn index
  # build time if the index is used later in another function. In that case,
  # test for nn index consistency.
# Check for consistency between matrix and index before
# using the following code.
#   if((nn_method == 'annoy' || nn_method == 'hnsw')) {
#      if(!check_cds_nn_index_is_current(cds=cds, reduction_method=reduction_method, nn_control=nn_control, verbose=verbose)) {
#        nn_index <- make_nn_index(subject_matrix=reducedDims(cds)[[reduction_method]],
#                                  nn_control=nn_control,
#                                  verbose=verbose)
#        cds <- set_cds_nn_index(cds=cds,
#                                reduction_method=reduction_method,
#                                nn_index=nn_index,
#                                nn_control=nn_control,
#                                verbose=verbose)
#      }
#      else {
#        nn_index <- get_cds_nn_index(cds=cds,
#                                     reduction_method=reduction_method,
#                                     nn_control=nn_control,
#                                     verbose=verbose)
#      }
#   }
#   else
#   if(nn_method == 'nn2') {
#     nn_index <- NULL
#   }

  # Set nn_index to NULL so that louvain_clustering and
  # leiden_clustering make a new index.
  nn_index <- NULL

  reduced_dim_res <- reducedDims(cds)[[reduction_method]]

  if(is.null(random_seed)) {
    random_seed <- sample.int(.Machine$integer.max, 1)
  }
  if(verbose)
    message("Running ", cluster_method, " clustering algorithm ...")

  if(cluster_method=='louvain') {
    cluster_result <- louvain_clustering(data=reduced_dim_res,
                                         pd=colData(cds),
                                         weight=weight,
                                         nn_index=nn_index,
                                         k=k,
                                         nn_control=nn_control,
                                         louvain_iter=num_iter,
                                         random_seed=random_seed,
                                         verbose=verbose)

    if (length(unique(cluster_result$optim_res$membership)) > 1) {
      cluster_graph_res <- compute_partitions(cluster_result$g,
                                              cluster_result$optim_res,
                                              partition_qval, verbose)
      partitions <- igraph::components(
        cluster_graph_res$cluster_g)$membership[cluster_result$optim_res$membership]
      partitions <- as.factor(partitions)
    } else {
      partitions <- rep(1, nrow(colData(cds)))
    }
    names(partitions) <- row.names(reduced_dim_res)
    clusters <- factor(igraph::membership(cluster_result$optim_res))
    cds@clusters[[reduction_method]] <- list(cluster_result = cluster_result,
                                             partitions = partitions,
                                             clusters = clusters)
  }
  else if(cluster_method=='leiden'){
    cds <- add_citation(cds, "leiden")
    cluster_result <- leiden_clustering(data=reduced_dim_res,
                                        pd=colData(cds),
                                        weight=weight,
                                        nn_index=nn_index,
                                        k=k,
                                        nn_control=nn_control,
                                        num_iter=num_iter,
                                        resolution_parameter=resolution,
                                        random_seed=random_seed,
                                        verbose=verbose, ...)

    if(length(unique(cluster_result$optim_res$membership)) > 1) {
      cluster_graph_res <- compute_partitions(cluster_result$g,
                                              cluster_result$optim_res,
                                              partition_qval, verbose)
      partitions <- igraph::components(
        cluster_graph_res$cluster_g)$membership[cluster_result$optim_res$membership]
      partitions <- as.factor(partitions)
    } else {
      partitions <- rep(1, nrow(colData(cds)))
    }
    names(partitions) <- row.names(reduced_dim_res)
    clusters <- factor(igraph::membership(cluster_result$optim_res))
    cds@clusters[[reduction_method]] <- list(cluster_result = cluster_result,
                                             partitions = partitions,
                                             clusters = clusters)
  }
  cds <- add_citation(cds, "clusters")
  cds <- add_citation(cds, "partitions")
  return(cds)
}


cluster_cells_make_graph <- function(data,
                                     weight,
                                     cell_names,
                                     nn_index=NULL,
                                     k=k,
                                     nn_control=list(),
                                     verbose) {
  if (is.data.frame(data))
    data <- as.matrix(data)
  if (!is.matrix(data))
    stop("Wrong input data, should be a data frame or matrix!")

  if (k < 1) {
    stop("k must be a positive integer!")
  } else
  if (k > nrow(data) - 2) {
    k <- nrow(data) - 2
    warning(paste("The nearest neighbors includes the point itself, k must be smaller than\nthe",
                  "total number of points - 1 (all other points) - 1",
                  "(itself)!",
                  "Total number of points is", nrow(data)))
  }

  if (verbose) {
    message("Run kNN based graph clustering starts:", "\n",
            "  -Input data of ", nrow(data), " rows and ", ncol(data),
            " columns", "\n", "  -k is set to ", k)
    message("  Finding nearest neighbors...")
  }

  nn_method <- nn_control[['method']]
  if(nn_method == 'nn2') {
    t1 <- system.time(tmp <- RANN::nn2(data, data, k+1, searchtype = "standard"))
  }
  else {
    if(is.null(nn_index)) {
      nn_index <- make_nn_index(subject_matrix=data, nn_control=nn_control, verbose=verbose)
    }
    tmp <- search_nn_index(query_matrix=data,
                           nn_index=nn_index,
                           k=k+1,
                           nn_control=nn_control,
                           verbose=verbose)
    if(nn_method == 'annoy' || nn_method == 'hnsw') {
      tmp <- swap_nn_row_index_point(nn_res=tmp, verbose=verbose)
    }
  }

  neighborMatrix <- tmp[['nn.idx']][, -1]
  distMatrix <- tmp[['nn.dists']][, -1]

  if (verbose) {
    if(nn_method == 'nn2') {
      message("DONE. Run time: ", t1[3], "s\n")
    }
    message("Compute jaccard coefficient between nearest-neighbor sets ..." )
  }

  t2 <- system.time(links <- jaccard_coeff(neighborMatrix, weight))

  if (verbose)
    message("DONE. Run time:", t2[3], "s\n", " Build undirected graph from the weighted links ...")

  links <- links[links[, 1] > 0,]
  relations <- as.data.frame(links)
  colnames(relations) <- c("from", "to", "weight")

  relations$from <- cell_names[relations$from]
  relations$to <- cell_names[relations$to]

  t3 <- system.time(g <- igraph::graph.data.frame(relations, directed = FALSE))

  if (verbose)
    message("DONE ~", t3[3], "s\n")

  return(list(g=g, distMatrix=distMatrix, relations=relations))
}


# Notes:
#   o  louvain_clustering does not update the nearest neighbor index
#      stored in the cds because it does not return a cds.
louvain_clustering <- function(data,
                               pd,
                               weight=F,
                               nn_index=NULL,
                               k=20,
                               nn_control=list(),
                               louvain_iter=1,
                               random_seed=0L,
                               verbose=FALSE) {
  assertthat::assert_that(assertthat::is.count(k))

  cell_names <- row.names(pd)

  if(!identical(cell_names, row.names(pd)))
    stop("Phenotype and row name from the data doesn't match")

  graph_result <- cluster_cells_make_graph(data=data,
                                           weight=weight,
                                           cell_names=cell_names,
                                           nn_index,
                                           k=k, 
                                           nn_control=nn_control,
                                           verbose=verbose)

  if(verbose)
    message("  Run louvain clustering ...")

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

    Q <- igraph::cluster_louvain(graph_result[['g']])

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

  if(igraph::vcount(graph_result[['g']]) < 3000) {
    coord <- NULL
    edge_links <- NULL
  } else {
    coord <- NULL
    edge_links <- NULL
  }

  igraph::V(graph_result[['g']])$names <- as.character(igraph::V(graph_result[['g']]))
  return(list(g=graph_result[['g']],
              relations=graph_result[['relations']],
              distMatrix=graph_result[['distMatrix']],
              coord = coord,
              edge_links=edge_links,
              optim_res=optim_res))
}


# Notes:
#   o  leiden_clustering does not update the nearest neighbor index
#      stored in the cds because it does not return a cds.
leiden_clustering <- function(data,
                              pd,
                              weight=NULL,
                              nn_index=NULL,
                              k=20,
                              nn_control=list(),
                              num_iter=2,
                              resolution_parameter=0.0001,
                              random_seed=NULL,
                              verbose=FALSE, ...) {
  extra_arguments <- list(...)
  if( 'partition_type' %in% names( extra_arguments ) )
    partition_type <- extra_arguments[['partition_type']]
  else
    partition_type <- 'CPMVertexPartition'
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

  # Check input parameters.
  assertthat::assert_that(assertthat::is.count(k))

  # The following vertex partitions have no resolution parameter.
  if( partition_type %in% c('ModularityVertexPartition','SignificanceVertexPartition','SurpriseVertexPartition') )
  {
    resolution_parameter = NA
  }
  else if( is.null( resolution_parameter ) )
  {
    resolution_parameter = 0.0001
  }
  if( is.null( num_iter ) )
    num_iter = 2
  if( random_seed == 0L )
    random_seed = NULL
  cell_names <- row.names(pd)
  if(!identical(cell_names, row.names(pd)))
    stop("Phenotype and row name from the data don't match")

  graph_result <- cluster_cells_make_graph(data=data,
                                           weight=weight, 
                                           cell_names=cell_names,
                                           nn_index,
                                           k=k,
                                           nn_control=nn_control,
                                           verbose=verbose)

  if(verbose)
    message("  Run leiden clustering ...")

  t_start <- Sys.time()

  if(verbose)
  {
    table_results <- data.frame(
      resolution_parameter = double(),
      quality              = double(),
      modularity           = double(),
      significance         = double(),
      number_clusters      = integer() )
  }

  best_modularity <- -1
  best_result <- NULL
  best_resolution_parameter <- 'No resolution'
  # These three vertex partition types have a resolution parameter
  # so scan parameter range, if given.
  for(i in 1:length(resolution_parameter)) {
    cur_resolution_parameter <- resolution_parameter[i]
    cluster_result <- leidenbase::leiden_find_partition( graph_result[['g']],
                                                         partition_type = partition_type,
                                                         initial_membership = initial_membership,
                                                         edge_weights = edge_weights,
                                                         node_sizes = node_sizes,
                                                         seed = random_seed,
                                                         resolution_parameter = cur_resolution_parameter,
                                                         num_iter = num_iter,
                                                         verbose = verbose )
    quality      <- cluster_result[['quality']]
    modularity   <- cluster_result[['modularity']]
    significance <- cluster_result[['significance']]

    if(verbose)
      table_results <- rbind( table_results, data.frame(
        resolution_parameter = cur_resolution_parameter,
        quality              = quality,
        modularity           = modularity,
        significance         = significance,
        cluster_count        = max(cluster_result[['membership']]) ) )
    if(verbose)
      message('    Current resolution is ', cur_resolution_parameter,
              '; Modularity is ', modularity,
              '; Quality is ', quality,
              '; Significance is ', significance,
              '; Number of clusters is ', max(cluster_result[['membership']]))
    if(modularity > best_modularity) {
      best_result <- cluster_result
      best_resolution_parameter <- cur_resolution_parameter
      best_modularity <- modularity
    }
    if ( is.null( best_result ) ) {
      best_result <- cluster_result
      best_resolution_parameter <- NULL
      best_modularity <- cluster_result[['modularity']]
    }
  }
  t_end <-
    Sys.time()

  if(verbose)
  {
    message('    Done. Run time: ', t_end - t_start, 's\n')
    message('  Clustering statistics')
    selected <- vector( mode='character',
                        length = length( resolution_parameter ) )
    for( irespar in 1:length( resolution_parameter ) )
    {
      if( identical( table_results[['resolution_parameter']][irespar],
                     best_resolution_parameter ) )
        selected[irespar] <- '*'
      else
        selected[irespar] <- ' '
    }
    print( cbind(' '=' ',table_results, selected ),row.names=FALSE )
    message()
    message('  Cell counts by cluster')
    membership<-best_result[['membership']]
    membership_frequency <- stats::aggregate(data.frame(cell_count = membership),
                                             list(cluster = membership),
                                             length)
    membership_frequency <- cbind(' '=' ',membership_frequency,
                                  cell_fraction=sprintf("%.3f",membership_frequency[['cell_count']]/sum(membership_frequency[['cell_count']])))
    print( membership_frequency,row.names=FALSE)
    message()
    message('  Maximal modularity is ', best_modularity,
            ' for resolution parameter ', best_resolution_parameter)
    message("\n  Run kNN based graph clustering DONE.\n  -Number of clusters: ",
            max(best_result[['membership']]))
  }

  if(igraph::vcount(graph_result[['g']]) < 3000) {
    coord <- NULL
    edge_links <- NULL
  } else {
    coord <- NULL
    edge_links <- NULL
  }

  igraph::V(graph_result[['g']])$names <- as.character(igraph::V(graph_result[['g']]))
  out_result <- list(membership = best_result[['membership']],
                     modularity = best_result[['modularity']] )
  names(out_result$membership) = cell_names

  return(list(g=graph_result[['g']],
              relations=graph_result[['relations']],
              distMatrix=graph_result[['distMatrix']],
              coord=coord,
              edge_links=edge_links,
              optim_res=out_result))
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

