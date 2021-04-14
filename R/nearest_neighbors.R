#' @export
build_annoy_index <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), nn_metric=c('cosine', 'euclidean', 'manhattan', 'hamming')) {
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")

  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(nn_metric) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "nn_metric must be one of 'cosine', 'euclidean', 'manhattan', or 'hamming'")

  nn_metric <- match.arg(nn_metric)

  if(reduction_method == 'PCA') {
    assertthat::assert_that(!is.null(reducedDims(cds)[['PCA']]),
                            msg = paste("When reduction_method = 'PCA', the",
                                        "cds must have been preprocessed for",
                                        "PCA. Please run preprocess_cds with",
                                        "method = 'PCA' before running",
                                        "build_annoy_index with",
                                        "reduction_method = 'PCA'."))
    reduced_matrix <- reducedDims(cds)[['PCA']]
    cds@preprocess_aux[['PCA']][['nn_index']] <- SimpleList()
    annoy_index <- uwot:::annoy_build(X = reduced_matrix, metric = nn_metric)
    cds@preprocess_aux[['PCA']][['nn_index']][['annoy_index']] <- annoy_index
    cds@preprocess_aux[['PCA']][['nn_index']][['annoy_metric']] <- nn_metric
    cds@preprocess_aux[['PCA']][['nn_index']][['annoy_ndim']] <- ncol(reduced_matrix)
  } else
  if(reduction_method == 'LSI') {
    assertthat::assert_that(!is.null(reducedDims(cds)[['LSI']]),
                            msg = paste("When reduction_method = 'LSI', the",
                                        "cds must have been preprocessed for",
                                        "LSI. Please run preprocess_cds with",
                                        "method = 'LSI' before running",
                                        "build_annoy_index with",
                                        "reduction_method = 'LSI'."))

    reduced_matrix <- reducedDims(cds)[['LSI']]
    cds@preprocess_aux[['LSI']][['nn_index']] <- SimpleList()
    annoy_index <- uwot:::annoy_build(X = reduced_matrix, metric = nn_metric)
    cds@preprocess_aux[['LSI']][['nn_index']][['annoy_index']] <- annoy_index
    cds@preprocess_aux[['LSI']][['nn_index']][['annoy_metric']] <- nn_metric
    cds@preprocess_aux[['LSI']][['nn_index']][['annoy_ndim']] <- ncol(reduced_matrix)
  } else
  if(reduction_method == 'Aligned') {
    assertthat::assert_that(!is.null(reducedDims(cds)[['Aligned']]),
                            msg = paste("When reduction_method = 'Aligned', the",
                                        "cds must have been aligned.",
                                        "Please run align_cds before running",
                                        "build_annoy_index with",
                                        "reduction_method = 'Aligned'."))

    reduced_matrix <- reducedDims(cds)[['Aligned']]
    cds@preprocess_aux[['Aligned']][['nn_index']] <- SimpleList()
    annoy_index <- uwot:::annoy_build(X = reduced_matrix, metric = nn_metric)
    cds@preprocess_aux[['Aligned']][['nn_index']][['annoy_index']] <- annoy_index
    cds@preprocess_aux[['Aligned']][['nn_index']][['annoy_metric']] <- nn_metric
    cds@preprocess_aux[['Aligned']][['nn_index']][['annoy_ndim']] <- ncol(reduced_matrix)
  } else
  if(reduction_method == 'tSNE') {
    assertthat::assert_that(!is.null(reducedDims(cds)[['tSNE']]),
                            msg = paste("When reduction_method = 'tSNE', the",
                                        "cds must have been processed with.",
                                        "reduce_dimension.",
                                        "Please run reduce_dimension with",
                                        "method='tSNE' before running",
                                        "build_annoy_index."))
    reduced_matrix <- reducedDims(cds)[['tSNE']]
    cds@reduce_dim_aux[['tSNE']][['nn_index']] <- SimpleList()
    annoy_index <- uwot:::annoy_build(X = reduced_matrix, metric=nn_metric)
    cds@reduce_dim_aux[['tSNE']][['nn_index']][['annoy_index']] <- annoy_index
    cds@reduce_dim_aux[['tSNE']][['nn_index']][['annoy_metric']] <- nn_metric
    cds@reduce_dim_aux[['tSNE']][['nn_index']][['annoy_ndim']] <- ncol(reduced_matrix)
  } else
  if(reduction_method == 'UMAP') {
    assertthat::assert_that(!is.null(reducedDims(cds)[['UMAP']]),
                            msg = paste("When reduction_method = 'UMAP', the",
                                        "cds must have been processed with",
                                        "reduce_dimension.",
                                        "Please run reduce_dimension with",
                                        "method='UMAP' before running",
                                        "build_annoy_index."))

    reduced_matrix <- reducedDims(cds)[['UMAP']]
    cds@reduce_dim_aux[['UMAP']][['nn_index']] <- SimpleList()
    annoy_index <- uwot:::annoy_build(X = reduced_matrix, metric=nn_metric)
    cds@reduce_dim_aux[['UMAP']][['nn_index']][['annoy_index']] <- annoy_index
    cds@reduce_dim_aux[['UMAP']][['nn_index']][['annoy_metric']] <- nn_metric
    cds@reduce_dim_aux[['UMAP']][['nn_index']][['annoy_ndim']] <- ncol(reduced_matrix)
  }
  cds
}


#' Search Annoy index for cells near the cells in cds.
#'
#' Search the Annoy nearest neighbor index for cells near cells in cds.
#'
#' @param cds A cell_data_set.
#' @param stage The stage in which the Annoy index was made.
#' @param method The method for which the Annoy index was made. This
#'   can be method 'PCA' or 'LSI' for stage 'preprocess', method 'Aligned'
#'   for stage 'Align', and method 'tSNE' or 'UMAP' for stage 'reduce_dimension'.
#' @param n An integer for the number of nearest neighbors to return for
#'   each cell.
#' @param search_k An integer used to balance accuracy and speed. Larger
#'   values give more accurate results.
#' @param ... Parameters to pass through to annoy_search.
#'
#' @export
search_nn_index <- function(cds, stage=c('PCA', 'align', 'reduce_dimension'), method=c('PCA', 'LSI', 'Aligned', 'UMAP'), n=5, search_k=100 * n, ...) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(stage) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "stage must be 'preprocess', 'align', or 'reduce_dimension'")

  stage <- match.arg(stage)
  method <- match.arg(method)

  if(stage=='preprocess' ) {
    if(!(method %in% c('PCA', 'LSI'))) {
      stop('Stage \'preprocess\' allows methods \'PCA\' and \'LSI\' but not \'', method)
    }
    stage_aux <- cds@preprocess_aux
  } else
  if(stage == 'align') {
    if(!(method %in% c('Aligned'))) {
      stop('Stage \'align\' allows method \'aligned\' but not \'', method)
    }
    stage_aux <- cds@preprocess_aux
  } else
  if(stage == 'reduce_dimension') {
    if(!(method %in% c('UMAP'))) {
      stop('Stage \'reduce_dimension\' allows method \'UMAP\' but not \'', method)
    }
    stage_aux <- cds@reduce_dim_aux
  } else {
    stop('Unrecognized stage \'', stage, '\'')
  }

  assertthat::assert_that(!is.null(reducedDim(cds, method)),
    msg = paste0("The matrix for ", method,
                " does not exist. Please preprocess the",
                " cds as required."))

  query_matrix <- reducedDim(cds, method)

  # notes:
  #   o  there may be dependency on uwot version
  #   o  ensure that ann_res[[1]] are indices and ann_res[[2]] are distances
  #   o  set list names to nn.idx and nn.dists for compatibility with
  #      RANN::nn2()
  #
  ann_index <- stage_aux[[method]][['nn_index']][['annoy_index']][['ann']]
  tmp <- uwot:::annoy_search(X=query_matrix, ann=ann_index, k=n, search_k=search_k, ...)
  ann_res = list(nn.idx = tmp[['idx']], nn.dists = tmp[['dist']])
}

