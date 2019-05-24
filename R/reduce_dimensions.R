#' Compute a projection of a cell_data_set object into a lower dimensional
#' space with non-linear dimension reduction methods
#'
#' @description Monocle aims to learn how cells transition through a biological
#' program of gene expression changes in an experiment. Each cell can be viewed
#' as a point in a high-dimensional space, where each dimension describes the
#' expression of a different gene in the genome. Identifying the program of
#' gene expression changes is equivalent to learning a \emph{trajectory} that
#' the cells follow through this space. However, the more dimensions there are
#' in the analysis, the harder the trajectory is to learn. Fortunately, many
#' genes typically co-vary with one another, and so the dimensionality of the
#' data can be reduced with a wide variety of different algorithms. Monocle
#' provides three different algorithms for dimensionality reduction via
#' \code{reduce_dimension}. All take a cell_data_set object and a number of
#' dimensions allowed for the reduced space.
#'
#' @details You can choose a few different reduction algorithms: Independent
#' Component Analysis (ICA) and Discriminative Dimensionality Reduction with
#' Trees (DDRTree). The choice impacts numerous downstream analysis steps,
#' including \code{\link{order_cells}}. Choosing ICA will execute the ordering
#' procedure described in Trapnell and Cacchiarelli et al., which was
#' implemented in Monocle version 1. \code{\link[DDRTree]{DDRTree}} is a more
#' recent manifold learning algorithm developed by Qi Mao and colleages. It is
#' substantially more powerful, accurate, and robust for single-cell trajectory
#' analysis than ICA, and is now the default method.
#'
#' Often, experiments include cells from different batches or treatments. You
#' can reduce the effects of these treatments by transforming the data with a
#' linear model prior to dimensionality reduction. To do so, provide a model
#' formula through \code{residual_model_formula_str}.
#'
#' Prior to reducing the dimensionality of the data, it usually helps to
#' normalize it so that highly expressed or highly variable genes don't
#' dominate the computation. \code{reduce_dimension()} automatically transforms
#' the data in one of several ways depending on the \code{expression_family} of
#' the cell_data_set object. If the expression_family is \code{negbinomial} or
#' \code{negbinomial.size}, the data are variance-stabilized. If the
#' expression_family is \code{Tobit}, the data are adjusted by adding a
#' pseudocount (of 1 by default) and then log-transformed. If you don't want
#' any transformation at all, set norm_method to "none" and pseudo_count to 0.
#' This maybe useful for single-cell qPCR data, or data you've already
#' transformed yourself in some way.
#'
#' @param cds the cell_data_set upon which to perform this operation
#' @param max_components the dimensionality of the reduced space
#' @param reduction_method A character string specifying the algorithm to use
#'   for dimensionality reduction.
#' @param verbose Whether to emit verbose output during dimensionality
#'   reduction
#' @param ... additional arguments to pass to the dimensionality reduction
#'   function
#' @return an updated cell_data_set object
#' @references DDRTree: Qi Mao, Li Wang, Steve Goodison, and Yijun Sun.
#'   Dimensionality reduction via graph structure learning. In Proceedings of
#'   the 21th ACM SIGKDD International Conference on Knowledge Discovery and
#'   Data Mining, pages 765–774. ACM, 2015.
#' @references UMAP: McInnes, L, Healy, J, UMAP: Uniform Manifold Approximation
#'   and Projection for Dimension Reduction, ArXiv e-prints 1802.03426, 2018
#' @references tSNE: Laurens van der Maaten and Geoffrey Hinton. Visualizing
#'   data using t-SNE. J. Mach. Learn. Res., 9(Nov):2579– 2605, 2008.
#' @export
reduce_dimension <- function(cds,
                             max_components=2,
                             reduction_method=c("UMAP", 'tSNE', 'PCA'),
                             preprocess_method=c("PCA", "LSI"),
                             umap.metric = "cosine",
                             umap.min_dist = 0.1,
                             umap.n_neighbors = 15L,
                             umap.fast_sgd = TRUE,
                             umap.nn_method = "annoy",
                             cores=max(1, RcppParallel::defaultNumThreads()/2),
                             verbose=FALSE,
                             ...){
  extra_arguments <- list(...)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'UMAP', 'PCA' or 'tSNE'")

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(preprocess_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "preprocess_method must be one of 'PCA' or 'LSI'")

  reduction_method <- match.arg(reduction_method)
  preprocess_method <- match.arg(preprocess_method)

  assertthat::assert_that(assertthat::is.count(max_components))

  assertthat::assert_that(!is.null(reducedDims(cds)[[preprocess_method]]),
                          msg = paste("Data has not been preprocessed with",
                                      "chosen method:", preprocess_method,
                                      "Please run preprocess_cds with",
                                      "method =", preprocess_method,
                                      "before running reduce_dimension."))
  if(reduction_method == "PCA") {
    assertthat::assert_that(preprocess_method == "PCA",
                            msg = paste("preprocess_method must be 'PCA' when",
                                        "reduction_method = 'PCA'"))
    assertthat::assert_that(!is.null(reducedDims(cds)[["PCA"]]),
                            msg = paste("When reduction_method = 'PCA', the",
                                        "cds must have been preprocessed for",
                                        "PCA. Please run preprocess_cds with",
                                        "method = 'PCA' before running",
                                        "reduce_dimension with",
                                        "reduction_method = 'PCA'."))
  }

  #ensure results from RNG sensitive algorithms are the same on all calls
  set.seed(2016)

  preprocess_mat <- reducedDims(cds)[[preprocess_method]]

  if(reduction_method == "PCA") {
    if (verbose) message("Returning preprocessed PCA matrix")
  } else if (reduction_method == "tSNE") {
    if (verbose) message("Reduce dimension by tSNE ...")

    tsne_res <- Rtsne::Rtsne(as.matrix(preprocess_mat), dims = max_components,
                             pca = F, check_duplicates=FALSE, ...)

    tsne_data <- tsne_res$Y[, 1:max_components]
    row.names(tsne_data) <- colnames(tsne_data)

    reducedDims(cds)$tSNE <- tsne_data

  } else if (reduction_method == c("UMAP")) {
    if (verbose)
      message("Running Uniform Manifold Approximation and Projection")

    umap_res = uwot::umap(as.matrix(preprocess_mat),
                          n_components = max_components,
                          metric = umap.metric,
                          min_dist = umap.min_dist,
                          n_neighbors = umap.n_neighbors,
                          fast_sgd = umap.fast_sgd,
                          n_threads=cores,
                          verbose=verbose,
                          nn_method = umap.nn_method,
                          ...)

    row.names(umap_res) <- colnames(cds)
    reducedDims(cds)$UMAP <- umap_res
  }



  ## Clear out any old graphs:
  cds@principal_graph_aux[[reduction_method]] <- NULL
  cds@principal_graph[[reduction_method]] <- NULL
  cds@clusters[[reduction_method]] <- NULL

  cds
}

