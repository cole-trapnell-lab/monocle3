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
#' @param verbose Logical, whether to emit verbose output.
#' @param cores Number of cores to use for computing the UMAP.
#' @param build_nn_index logical When this argument is set to TRUE,
#'   preprocess_cds builds the nearest neighbor index from the
#'   reduced dimension matrix for later use. Default is FALSE.
#' @param nn_control An optional list of parameters used to make the nearest
#'  neighbor index. See the set_nn_control help for detailed information.
#'  The default metric is cosine for reduction_methods PCA, LSI, and Aligned,
#'  and is euclidean for reduction_methods tSNE and UMAP. Note: distances in
#'  tSNE space reflect spatial differences poorly so using nearest neighbors
#'  with it may be meaningless.
#' @param ... additional arguments to pass to the dimensionality reduction
#'   function.
#' @return an updated cell_data_set object
#' @references UMAP: McInnes, L, Healy, J, UMAP: Uniform Manifold Approximation
#'   and Projection for Dimension Reduction, ArXiv e-prints 1802.03426, 2018
#' @references tSNE: Laurens van der Maaten and Geoffrey Hinton. Visualizing
#'   data using t-SNE. J. Mach. Learn. Res., 9(Nov):2579– 2605, 2008.
#'
#' @examples
#'   \donttest{
#'     cds <- load_worm_embryo()
#'     cds <- preprocess_cds(cds)
#'     cds <- reduce_dimension(cds)
#'   }
#'
#' @importFrom methods is
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
                             verbose = FALSE,
                             cores = 1,
                             build_nn_index = FALSE,
                             nn_control = list(),
                             ...){

  extra_arguments <- list(...)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'UMAP', 'PCA', 'tSNE', 'LSI', 'Aligned'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(is.logical(build_nn_index),
                          msg = paste("build_nn_index must be either TRUE or FALSE"))

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

  if(build_nn_index) {
    if(reduction_method == 'tSNE' || reduction_method == 'UMAP')
      nn_control_default <- get_global_variable('nn_control_annoy_euclidean')
    else
      nn_control_default <- get_global_variable('nn_control_annoy_cosine')

    nn_control <- set_nn_control(mode=1,
                                 nn_control=nn_control,
                                 nn_control_default=nn_control_default,
                                 nn_index=NULL,
                                 k=NULL,
                                 verbose=verbose)
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

  #
  # Notes:
  #   o  the functions save_transform_models/load_transform_models
  #      expect that the reduce_dim_aux slot consists of a SimpleList
  #      that stores information about methods with the elements
  #        reduce_dim_aux[[method]][['model']] for the transform elements
  #        reduce_dim_aux[[method]][[nn_method]] for the annoy index
  #      and depends on the elements within model and nn_method.
  #
  if(reduction_method == "PCA") {
    if(build_nn_index && is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']])) {
      nn_index <- make_nn_index(subject_matrix=reducedDims(cds)[[reduction_method]],
                                nn_control=nn_control,
                                verbose=verbose)
      cds <- set_cds_nn_index(cds=cds, reduction_method=reduction_method, nn_index=nn_index, verbose=verbose)
    }
    if (verbose) message("Returning preprocessed PCA matrix")
  } else if(reduction_method == "LSI") {
    if(build_nn_index && is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']])) {
      nn_index <- make_nn_index(subject_matrix=reducedDims(cds)[[reduction_method]],
                                nn_control=nn_control,
                                verbose=verbose)
      cds <- set_cds_nn_index(cds=cds, reduction_method=reduction_method, nn_index=nn_index, verbose=verbose)
    }
    if (verbose) message("Returning preprocessed LSI matrix")
  } else if(reduction_method == "Aligned") {
    if(build_nn_index && is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']])) {
      nn_index <- make_nn_index(subject_matrix=reducedDims(cds)[[reduction_method]],
                                nn_control=nn_control,
                                verbose=verbose)
      cds <- set_cds_nn_index(cds=cds, reduction_method=reduction_method, nn_index=nn_index, verbose=verbose)
    }
    if (verbose) message("Returning preprocessed Aligned matrix")
  } else if (reduction_method == "tSNE") {
    if (verbose) message("Reduce dimension by tSNE ...")

    cds <- initialize_reduce_dim_metadata(cds, 'tSNE')
    cds <- initialize_reduce_dim_model_identity(cds, 'tSNE')

    tsne_res <- Rtsne::Rtsne(as.matrix(preprocess_mat), dims = max_components,
                             pca = FALSE, check_duplicates=FALSE, ...)

    tsne_data <- tsne_res$Y[, 1:max_components]
    row.names(tsne_data) <- colnames(tsne_data)

    reducedDims(cds)$tSNE <- tsne_data


    matrix_id <- get_unique_id(reducedDims(cds)[['tSNE']])
    reduce_dim_matrix_identity <- get_reduce_dim_matrix_identity(cds, preprocess_method)

    set_reduce_dim_matrix_identity(cds, 'tSNE',
                                   'matrix:tSNE',
                                   matrix_id,
                                   reduce_dim_matrix_identity[['matrix_type']],
                                   reduce_dim_matrix_identity[['matrix_id']],
                                   'matrix:tSNE',
                                   matrix_id)
    reduce_dim_model_identity <- get_reduce_dim_model_identity(cds, preprocess_method)
    set_reduce_dim_model_identity(cds, 'tSNE',
                                  'matrix:tSNE',
                                  matrix_id,
                                  reduce_dim_model_identity[['model_type']],
                                  reduce_dim_model_identity[['model_id']])

    # make nearest neighbor index in tSNE space

    if( build_nn_index ) {
      nn_index <- make_nn_index(subject_matrix=reducedDims(cds)[[reduction_method]],
                                nn_control=nn_control,
                                verbose=verbose)
      cds <- set_cds_nn_index(cds=cds, reduction_method=reduction_method, nn_index=nn_index, verbose=verbose)
    }
    else
      cds <- clear_cds_nn_index(cds=cds, reduction_method=reduction_method, nn_method='all')

  }
  else
  if (reduction_method == c("UMAP")) {
    cds <- add_citation(cds, "UMAP")
    if (verbose)
      message("Running Uniform Manifold Approximation and Projection")

    cds <- initialize_reduce_dim_metadata(cds, 'UMAP')
    cds <- initialize_reduce_dim_model_identity(cds, 'UMAP')

    umap_model <- uwot::umap(as.matrix(preprocess_mat),
                             n_components = max_components,
                             metric = umap.metric,
                             min_dist = umap.min_dist,
                             n_neighbors = umap.n_neighbors,
                             fast_sgd = umap.fast_sgd,
                             n_threads=cores,
                             verbose=verbose,
                             nn_method = umap.nn_method,
                             ret_model = TRUE,
                             ...)

    # Notes:
    #   o  uwot::umap_transform() returns a slightly different result in
    #      comparison to uwot::umap() (umap_res$embedding) model, even
    #      when uwot::umap_transform() uses the model from uwot::umap().
    #      However, uwot::umap_transform() gives consistent results using
    #      one model from uwot::umap(). So return the result from
    #      uwot::umap_transform().
    #   o  uwot::umap_transform() depends on the RNG state and we want
    #      consistent results when called here and in *_transform function(s).
    #
    set.seed(2016)
    umap_res <- uwot::umap_transform(X=as.matrix(preprocess_mat), model=umap_model, n_threads=1)
    row.names(umap_res) <- colnames(cds)
    reducedDims(cds)[['UMAP']] <- umap_res

    cds@reduce_dim_aux[['UMAP']][['model']][['umap_preprocess_method']] <- preprocess_method
    cds@reduce_dim_aux[['UMAP']][['model']][['max_components']] <- max_components
    cds@reduce_dim_aux[['UMAP']][['model']][['umap_metric']] <- umap.metric
    cds@reduce_dim_aux[['UMAP']][['model']][['umap_min_dist']] <- umap.min_dist
    cds@reduce_dim_aux[['UMAP']][['model']][['umap_n_neighbors']] <- umap.n_neighbors
    cds@reduce_dim_aux[['UMAP']][['model']][['umap_fast_sgd']] <- umap.fast_sgd
    cds@reduce_dim_aux[['UMAP']][['model']][['umap_model']] <- umap_model

    matrix_id <- get_unique_id(reducedDims(cds)[['UMAP']])
    reduce_dim_matrix_identity <- get_reduce_dim_matrix_identity(cds, preprocess_method)

    cds <- set_reduce_dim_matrix_identity(cds, 'UMAP',
                                          'matrix:UMAP',
                                          matrix_id,
                                          reduce_dim_matrix_identity[['matrix_type']],
                                          reduce_dim_matrix_identity[['matrix_id']],
                                          'matrix:UMAP',
                                          matrix_id)
    reduce_dim_model_identity <- get_reduce_dim_model_identity(cds, preprocess_method)
    cds <- set_reduce_dim_model_identity(cds, 'UMAP',
                                         'matrix:UMAP',
                                         matrix_id,
                                         reduce_dim_model_identity[['model_type']],
                                         reduce_dim_model_identity[['model_id']])

    if( build_nn_index ) {
      nn_index <- make_nn_index(subject_matrix=reducedDims(cds)[[reduction_method]],
                                nn_control=nn_control,
                                verbose=verbose)
      cds <- set_cds_nn_index(cds=cds, reduction_method=reduction_method, nn_index=nn_index, verbose=verbose)
    }
    else
      cds <- clear_cds_nn_index(cds=cds, reduction_method=reduction_method, nn_method='all')
  }

  ## Clear out old graphs:
  cds@principal_graph_aux[[reduction_method]] <- NULL
  cds@principal_graph[[reduction_method]] <- NULL
  cds@clusters[[reduction_method]] <- NULL

  cds
}

