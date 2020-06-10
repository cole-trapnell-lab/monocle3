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
#' \code{reduce_dimension} (UMAP and tSNE). The function
#' \code{reduce_dimension} is the second step in the trajectory building
#' process after \code{preprocess_cds}.
#'
#' UMAP is implemented from the package uwot.
#'
#' @param cds the cell_data_set upon which to perform this operation.
#' @param max_components the dimensionality of the reduced space. Default is 2.
#' @param reduction_method A character string specifying the algorithm to use
#'   for dimensionality reduction. Currently "UMAP", "tSNE", "PCA", "LSI", and "Aligned"
#'   are supported.
#' @param preprocess_method A string indicating the preprocessing method used
#'   on the data. Options are "PCA" and "LSI". Default is "LSI".
#' @param umap.metric A string indicating the distance metric to be used when
#'   calculating UMAP. Default is "cosine". See uwot package's
#'   \code{\link[umap]{umap}} for details.
#' @param umap.min_dist Numeric indicating the minimum distance to be passed to
#'   UMAP function. Default is 0.1.See uwot package's \code{\link[umap]{umap}}
#'   for details.
#' @param umap.n_neighbors Integer indicating the number of neighbors to use
#'   during kNN graph construction. Default is 15L. See uwot package's
#'   \code{\link[umap]{umap}} for details.
#' @param umap.fast_sgd Logical indicating whether to use fast SGD. Default is
#'   TRUE. See uwot package's \code{\link[umap]{umap}} for details.
#' @param umap.nn_method String indicating the nearest neighbor method to be
#'   used by UMAP. Default is "annoy". See uwot package's
#'   \code{\link[umap]{umap}} for details.
#' @param cores Number of compute cores to use.
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
                             reduction_method=c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
                             #preprocess_method=c("PCA", "LSI", "Aligned"),
                             preprocess_method=NULL,
                             umap.metric = "cosine",
                             umap.min_dist = 0.1,
                             umap.n_neighbors = 15L,
                             umap.fast_sgd = FALSE,
                             umap.nn_method = "annoy",
                             cores=1,
                             verbose=FALSE,
                             ...){
  extra_arguments <- list(...)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'UMAP', 'PCA', 'tSNE', 'LSI', 'Aligned'")

  reduction_method <- match.arg(reduction_method)

  if (is.null(preprocess_method)){
    if ("Aligned" %in% names(reducedDims(cds))){
      preprocess_method = "Aligned"
      message(paste("No preprocess_method specified, and aligned coordinates have been computed previously. Using preprocess_method = 'Aligned'"))
    }else{
      preprocess_method = "PCA"
      message(paste("No preprocess_method specified, using preprocess_method = 'PCA'"))
    }
  }else{
    assertthat::assert_that(
      preprocess_method %in% c("PCA", "LSI", "Aligned"),
      msg = "preprocess_method must be one of 'PCA' or 'LSI'")
  }

  #preprocess_method <- match.arg(preprocess_method)


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

  if(reduction_method == "LSI") {
    assertthat::assert_that(preprocess_method == "LSI",
                            msg = paste("preprocess_method must be 'LSI' when",
                                        "reduction_method = 'LSI'"))
    assertthat::assert_that(!is.null(reducedDims(cds)[["LSI"]]),
                            msg = paste("When reduction_method = 'LSI', the",
                                        "cds must have been preprocessed for",
                                        "LSI. Please run preprocess_cds with",
                                        "method = 'LSI' before running",
                                        "reduce_dimension with",
                                        "reduction_method = 'LSI'."))
  }

  if(reduction_method == "Aligned") {
    assertthat::assert_that(preprocess_method == "Aligned",
                            msg = paste("preprocess_method must be 'Aligned' when",
                                        "reduction_method = 'Aligned'"))
    assertthat::assert_that(!is.null(reducedDims(cds)[["Aligned"]]),
                            msg = paste("When reduction_method = 'Aligned', the",
                                        "cds must have been aligned.",
                                        "Please run align_cds before running",
                                        "reduce_dimension with",
                                        "reduction_method = 'Aligned'."))
  }

  #ensure results from RNG sensitive algorithms are the same on all calls
  set.seed(2016)

  if (reduction_method=="UMAP" && (umap.fast_sgd == TRUE || cores > 1)){
    message(paste("Note: reduce_dimension will produce slightly different",
                  "output each time you run it unless you set",
                  "'umap.fast_sgd = FALSE' and 'cores = 1'"))
  }

  preprocess_mat <- reducedDims(cds)[[preprocess_method]]

  if(reduction_method == "PCA") {
    if (verbose) message("Returning preprocessed PCA matrix")
  } else if(reduction_method == "LSI") {
    if (verbose) message("Returning preprocessed LSI matrix")
  } else if (reduction_method == "tSNE") {
    if (verbose) message("Reduce dimension by tSNE ...")

    tsne_res <- Rtsne::Rtsne(as.matrix(preprocess_mat), dims = max_components,
                             pca = F, check_duplicates=FALSE, ...)

    tsne_data <- tsne_res$Y[, 1:max_components]
    row.names(tsne_data) <- colnames(tsne_data)

    reducedDims(cds)$tSNE <- tsne_data

  } else if (reduction_method == c("UMAP")) {
    cds <- add_citation(cds, "UMAP")
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

