#' Compute a projection of a CellDataSet object into a lower dimensional space with non-linear dimension reduction methods
#'
#' @description Monocle aims to learn how cells transition through a biological program of
#' gene expression changes in an experiment. Each cell can be viewed as a point
#' in a high-dimensional space, where each dimension describes the expression of
#' a different gene in the genome. Identifying the program of gene expression
#' changes is equivalent to learning a \emph{trajectory} that the cells follow
#' through this space. However, the more dimensions there are in the analysis,
#' the harder the trajectory is to learn. Fortunately, many genes typically
#' co-vary with one another, and so the dimensionality of the data can be
#' reduced with a wide variety of different algorithms. Monocle provides two
#' different algorithms for dimensionality reduction via \code{reduce_dimension}.
#' Both take a CellDataSet object and a number of dimensions allowed for the
#' reduced space. You can also provide a model formula indicating some variables
#' (e.g. batch ID or other technical factors) to "subtract" from the data so it
#' doesn't contribute to the trajectory.
#'
#' @details You can choose a few different reduction algorithms: Independent Component
#' Analysis (ICA) and Discriminative Dimensionality Reduction with Trees (DDRTree).
#' The choice impacts numerous downstream analysis steps, including \code{\link{orderCells}}.
#' Choosing ICA will execute the ordering procedure described in Trapnell and Cacchiarelli et al.,
#' which was implemented in Monocle version 1. \code{\link[DDRTree]{DDRTree}} is a more recent manifold
#' learning algorithm developed by Qi Mao and colleages. It is substantially more
#' powerful, accurate, and robust for single-cell trajectory analysis than ICA,
#' and is now the default method.
#'
#' Often, experiments include cells from different batches or treatments. You can
#' reduce the effects of these treatments by transforming the data with a linear
#' model prior to dimensionality reduction. To do so, provide a model formula
#' through \code{residualModelFormulaStr}.
#'
#' Prior to reducing the dimensionality of the data, it usually helps
#' to normalize it so that highly expressed or highly variable genes don't
#' dominate the computation. \code{reduce_dimension()} automatically transforms
#' the data in one of several ways depending on the \code{expression_family} of
#' the CellDataSet object. If the expression_family is \code{negbinomial} or \code{negbinomial.size}, the
#' data are variance-stabilized. If the expression_family is \code{Tobit}, the data
#' are adjusted by adding a pseudocount (of 1 by default) and then log-transformed.
#' If you don't want any transformation at all, set norm_method to "none" and
#' pseudo_expr to 0. This maybe useful for single-cell qPCR data, or data you've
#' already transformed yourself in some way.
#'
#' @param cds the CellDataSet upon which to perform this operation
#' @param max_components the dimensionality of the reduced space
#' @param reduction_method A character string specifying the algorithm to use for dimensionality reduction.
#' @param auto_param_selection when this argument is set to TRUE (default), it will automatically calculate the proper value for the ncenter (number of centroids) parameters which will be passed into DDRTree call.
#' @param scaling When this argument is set to TRUE (default), it will scale each gene before running trajectory reconstruction.
#' @param verbose Whether to emit verbose output during dimensionality reduction
#' @param ... additional arguments to pass to the dimensionality reduction function
#' @return an updated CellDataSet object
#' @references DDRTree: Qi Mao, Li Wang, Steve Goodison, and Yijun Sun. Dimensionality reduction via graph structure learning. In Proceedings of the 21th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining, pages 765–774. ACM, 2015.
#' @references UMAP: McInnes, L, Healy, J, UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction, ArXiv e-prints 1802.03426, 2018
#' @references tSNE: Laurens van der Maaten and Geoffrey Hinton. Visualizing data using t-SNE. J. Mach. Learn. Res., 9(Nov):2579– 2605, 2008.
#' @export
reduce_dimension <- function(cds,
                            max_components=2,
                            reduction_method=c("UMAP", 'tSNE', 'none'),
                            auto_param_selection = TRUE,
                            scaling = TRUE,
                            verbose=FALSE,
                            ...){
  extra_arguments <- list(...)
  set.seed(2016) #ensure results from RNG sensitive algorithms are the same on all calls

  if (verbose)
    message("Retrieving normalized data ...")

  FM <- cds@aux_ordering_data$normalize_expr_data
  irlba_pca_res <- cds@normalized_data_projection

  if(is.null(FM)) {
    message('Warning: The cds has not been pre-processed yet. Running preprocessCDS() with default parameters.')
    cds <- preprocessCDS(cds)
    FM <- cds@aux_ordering_data$normalize_expr_data
    irlba_pca_res <- cds@normalized_data_projection
  }

  #FM <- FM[apply(FM, 1, function(x) all(is.finite(x))), ] #ensure all the expression values are finite values
  if (is.function(reduction_method)) {

    if(scaling){
      FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
      FM <- FM[!is.na(row.names(FM)), ]
    } else FM <- as.matrix(FM)

    reducedDim <- reduction_method(FM, ...)
    colnames(reducedDim) <- colnames(FM)
    reducedDimW(cds) <- as.matrix(reducedDim)
    reducedDimA(cds) <- as.matrix(reducedDim)
    reducedDimS(cds) <- as.matrix(reducedDim)
    reducedDimK(cds) <- as.matrix(reducedDim)
    dp <- as.matrix(dist(reducedDim))
    cellPairwiseDistances(cds) <- dp
    gp <- igraph::graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)
    principal_graph(cds) <- dp_mst
    cds@dim_reduce_type <- "function_passed"
  }
  else{
    reduction_method <- match.arg(reduction_method)
    if (reduction_method == "tSNE") {
      if("num_dim" %in% names(extra_arguments)){ #when you pass pca_dim to the function, the number of dimension used for tSNE dimension reduction is used
        num_dim <- extra_arguments$num_dim #variance_explained
      }
      else{
        num_dim <- 50
      }

      topDim_pca <- irlba_pca_res#[, 1:num_dim]

      #then run tSNE
      if (verbose)
        message("Reduce dimension by tSNE ...")

      tsne_res <- Rtsne::Rtsne(as.matrix(topDim_pca), dims = max_components, pca = F, check_duplicates=FALSE, ...)

      tsne_data <- tsne_res$Y[, 1:max_components]
      row.names(tsne_data) <- colnames(tsne_data)

      reducedDimA(cds) <- t(tsne_data) #this may move to the auxClusteringData environment

      #set the important information from densityClust to certain part of the cds object:
      cds@aux_clustering_data[["tSNE"]]$pca_components_used <- num_dim

      cds@dim_reduce_type <- "tSNE"

      pData(cds)$tsne_1 = reducedDimA(cds)[1,]
      pData(cds)$tsne_2 = reducedDimA(cds)[2,]

    }else if (reduction_method == c("UMAP") ) {
      if (verbose)
        message("Running Uniform Manifold Approximation and Projection")

      umap_args <- c(list(X = irlba_pca_res, log = F, n_component = as.integer(max_components), verbose = verbose, return_all = T),
                     extra_arguments[names(extra_arguments) %in%
                                       c("python_home", "n_neighbors", "metric", "n_epochs", "negative_sample_rate", "learning_rate", "init", "min_dist", "spread",
                                         'set_op_mix_ratio', 'local_connectivity', 'repulsion_strength', 'a', 'b', 'random_state', 'metric_kwds', 'angular_rp_forest', 'verbose')])
      tmp <- do.call(UMAP, umap_args)
      tmp$embedding_ <- (tmp$embedding_ - min(tmp$embedding_)) / max(tmp$embedding_) # normalize UMAP space
      umap_res <- tmp$embedding_;

      adj_mat <- Matrix::sparseMatrix(i = tmp$graph_$indices, p = tmp$graph_$indptr,
                                      x = -as.numeric(tmp$graph_$data), dims = c(ncol(cds), ncol(cds)), index1 = F,
                                      dimnames = list(colnames(cds), colnames(cds)))

      S <- t(umap_res)

      Y <- S
      W <- t(irlba_pca_res)

      principal_graph(cds) <- igraph::graph_from_adjacency_matrix(adj_mat, weighted=TRUE)

      A <- S
      colnames(A) <- colnames(FM)
      reducedDimA(cds) <- A

      colnames(S) <- colnames(FM)
      colnames(Y) <- colnames(FM)
      reducedDimW(cds) <- W
      reducedDimS(cds) <- as.matrix(Y)
      reducedDimK(cds) <- S

      cds@dim_reduce_type <- reduction_method
    } else if(reduction_method == 'none') {
      irlba_pca_res <- t(irlba_pca_res)
      colnames(irlba_pca_res) <- colnames(FM)
      reducedDimS(cds) <- irlba_pca_res
      reducedDimK(cds) <- irlba_pca_res
    }else {
      stop("Error: unrecognized dimensionality reduction method")
    }
  }
  cds
}
