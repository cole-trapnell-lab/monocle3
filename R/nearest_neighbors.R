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
