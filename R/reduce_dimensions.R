#' Compute a projection of a cell_data_set object into a lower dimensional
#' space with non-linear dimension reduction methods
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
#' \code{reduce_dimensions} (UMAP and tSNE). The function
#' \code{reduce_dimensions} is the second step in the trajectory building
#' process after \code{preprocess_cds}.
#'
#' UMAP is implemented from the package uwot.
#'
#' @param cds the cell_data_set upon which to perform this operation.
#' @param max_components the dimensionality of the reduced space. Default is 2.
#' @param reduction_method A character string specifying the algorithm to use
#'   for dimensionality reduction. Currently "UMAP", "tSNE", and "PCA" are
#'   supported.
#' @param umap.metric A string indicating the distance metric to be used when
#'   calculating UMAP. Default is "cosine". See uwot package's
#'   \code{\link[umap]{umap}} for details.
#' @param umap.min_dist Numeric indicating the minimum distance to be passed to
#'   UMAP function. Default is 0.1.See uwot package's \code{\link[umap]{umap}}
#'   for details.
#' @param umap.n_neighbors Integer indicating the number of neighbors to use
#'   during KNN graph construction. Default is 15L. See uwot package's
#'   \code{\link[umap]{umap}} for details.
#' @param umap.fast_sgd Logical indicating whether to use fast SGD. Default is
#'   TRUE. See uwot package's \code{\link[umap]{umap}} for details.
#' @param umap.nn_method String indicating the nearest neighbor method to be
#'   used by umap. Default is "annoy". See uwot package's
#'   \code{\link[umap]{umap}} for details.
#' @param verbose Logical, whether to emit verbose output.
#' @param ... additional arguments to pass to the dimensionality reduction
#'   function.
#' @return an updated cell_data_set object
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

