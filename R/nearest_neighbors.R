
# Check whether annoy index exists.
check_nn_index_exists <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), nn_method=c('annoy'), nn_metric=c('euclidean', 'cosine', 'manhattan', 'hamming'), n_trees) {
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(nn_metric) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "nn_metric must be one of 'euclidean', 'cosine', 'manhattan', or 'hamming'")
  nn_metric <- match.arg(nn_metric)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(nn_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "nn_method must be 'annoy'")
  nn_method <- match.arg(nn_method)

  assertthat::assert_that(assertthat::is.count(n_trees))

  reduced_matrix <- reducedDims(cds)[[reduction_method]]

  if(is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']]) ||
     cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy_nn_metric']] != nn_metric ||
     cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy_n_trees']] != n_trees ||
     cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy_matrix_id']] != get_reduce_dim_matrix_identity(cds, reduction_method)[['matrix_id']]) {

     return(FALSE)
  }

  return(TRUE)
}


# Clear annoy index.
clear_nn_index <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), nn_method=c('annoy')) {
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(nn_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "nn_method must be 'annoy'")
  nn_method <- match.arg(nn_method)

  if(!is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']]) &&
     !is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy_index']])) {
     cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy_index']] <- NULL
  }

  return(cds)
}


#' @export
build_nn_index <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), nn_method=c('annoy'), nn_metric=c('euclidean', 'cosine', 'manhattan', 'hamming'), n_trees=50, cores=1) {
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(nn_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "nn_method must be 'annoy'")
  nn_method <- match.arg(nn_method)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(nn_metric) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "nn_metric must be one of 'euclidean', 'cosine', 'manhattan', or 'hamming'")
  nn_metric <- match.arg(nn_metric)

  assertthat::assert_that(assertthat::is.count(n_trees))

  if(reduction_method == 'PCA') {
    assertthat::assert_that(!is.null(reducedDims(cds)[['PCA']]),
                            msg = paste("When reduction_method = 'PCA', the",
                                        "cds must have been preprocessed for",
                                        "PCA. Please run preprocess_cds with",
                                        "method = 'PCA' before running",
                                        "build_nn_index with",
                                        "reduction_method = 'PCA'."))
  } else
  if(reduction_method == 'LSI') {
    assertthat::assert_that(!is.null(reducedDims(cds)[['LSI']]),
                            msg = paste("When reduction_method = 'LSI', the",
                                        "cds must have been preprocessed for",
                                        "LSI. Please run preprocess_cds with",
                                        "method = 'LSI' before running",
                                        "build_nn_index with",
                                        "reduction_method = 'LSI'."))
  } else
  if(reduction_method == 'Aligned') {
    assertthat::assert_that(!is.null(reducedDims(cds)[['Aligned']]),
                            msg = paste("When reduction_method = 'Aligned', the",
                                        "cds must have been aligned.",
                                        "Please run align_cds before running",
                                        "build_nn_index with",
                                        "reduction_method = 'Aligned'."))
  } else
  if(reduction_method == 'tSNE') {
    assertthat::assert_that(!is.null(reducedDims(cds)[['tSNE']]),
                            msg = paste("When reduction_method = 'tSNE', the",
                                        "cds must have been processed with.",
                                        "reduce_dimension.",
                                        "Please run reduce_dimension with",
                                        "method='tSNE' before running",
                                        "build_nn_index."))
  } else
  if(reduction_method == 'UMAP') {
    assertthat::assert_that(!is.null(reducedDims(cds)[['UMAP']]),
                            msg = paste("When reduction_method = 'UMAP', the",
                                        "cds must have been processed with",
                                        "reduce_dimension.",
                                        "Please run reduce_dimension with",
                                        "method='UMAP' before running",
                                        "build_nn_index."))
  }

  reduced_matrix <- reducedDims(cds)[[reduction_method]]

  cds@reduce_dim_aux[[reduction_method]][['nn_index']] <- SimpleList()
  annoy_index <- uwot:::annoy_build(X = reduced_matrix, metric = nn_metric, n_trees = n_trees)
  cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy_index']] <- annoy_index
  cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy_nn_metric']] <- nn_metric
  cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy_n_trees']] <- n_trees
  cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy_matrix_id']] <- get_reduce_dim_matrix_identity(cds, reduction_method)[['matrix_id']]

  return(cds)
}


#' Search index for cells near the cells in cds.
#'
#' Search the nearest neighbor index for cells near cells in cds.
#'
#' @param cds A cell_data_set.
#' @param reduction_method The method for which the index was made. method
#'   can be 'PCA', 'LSI','Aligned', 'tSNE', or 'UMAP'.
#' @param k An integer for the number of nearest neighbors to return for
#'   each cell. Default is 5.
#' @param search_k An integer used to balance accuracy and speed. Larger
#'   values give more accurate results. Default is 100 * k.
#' @param cores Integer The number of threads used for the nearest
#'   neighbor search. Default is 1.
#' @param ... Parameters to pass through to annoy_search.
#'
#' @export
search_nn_index <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), nn_method=c('annoy'), k=5, search_k=100 * k, cores=1, ...) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(!is.null(reducedDim(cds, reduction_method)),
    msg = paste0("The matrix for ", reduction_method,
                " does not exist. Please preprocess the",
                " cds as required."))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(nn_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "nn_method must be 'annoy'")
  nn_method <- match.arg(nn_method)

  assertthat::assert_that(assertthat::is.count(k))
  assertthat::assert_that(assertthat::is.count(search_k))

  query_matrix <- reducedDim(cds, reduction_method)

  # notes:
  #   o  there may be dependency on uwot version
  #   o  ensure that ann_res[[1]] are indices and ann_res[[2]] are distances
  #   o  set list names to nn.idx and nn.dists for compatibility with
  #      RANN::nn2()
  #
  ann_index <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy_index']][['ann']]
  tmp <- uwot:::annoy_search(X=query_matrix, k=k, ann=ann_index, search_k=search_k, n_threads=cores, ...)
  
  ann_res = list(nn.idx = tmp[['idx']], nn.dists = tmp[['dist']])
}


# Given a vector of nearest neighbors from a
# search for which the query and database are
# the same set and the self point is not the
# first element, swap it to the first element.
nn_vector_index_swap <- function(vec) {
  # Most often, the self point is shifted right
  # one column; that is, there is one other point
  # with zero distance.
  if(vec[2] != vec[1]) {
    if(vec[3] == vec[1]) {
      vec[3] <- vec[2]
      vec[2] <- vec[1]
      return(vec)
    } else
    {
      for( i in 2:length(vec)) {
        if(vec[i] == vec[1]) {
          vec[i] <- vec[2]
          vec[2] <- vec[1]
          return(vec)
        }
      }
      # The self point is not in the vector so
      # put it in the second column.
      vec[2] <- vec[1]
      return(vec)
    }
  }
  vec
}


# Sort annoy self-set search results so that self
# points are in the first column.
#
# For use when the query and database sets are the
# same and the self points must be in the first
# column.
#
# This is required because annoy does not sort the
# self points into the first column. That is, in
# cases where non-self points have zero distance,
# the self point may be in a column j > 1.
#
# Warning: do not use this if the search point set
#          differs from the index set. It will
#          make a mess of the result.
#
swap_nn_self_search_index_swap <- function(nn_res) {
  idx <- nn_res[['nn.idx']]
  dists <- nn_res[['nn.dists']]

  iidx <- cbind(1:nrow(idx),idx)
  shift_res <- t(apply(iidx, 1, nn_vector_index_swap))

  list(nn.idx=shift_res[,-1], nn.dists=dists)
}


