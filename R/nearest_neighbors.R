
# Check whether annoy index exists.
check_nn_index_exists <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), nn_method=c('annoy', 'hnsw')) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  if(nn_method == 'annoy') {
    if(is.null(cds@reduce_dim_aux[[reduction_method]]) ||
       is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']]) ||
       is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]])) {
       return(FALSE)
    }
  }

  return(TRUE)
}


# Check wether annoy index exists and is consistent with matrix and parameters.
check_nn_index_current <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), nn_control=list()) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(
    !is.null(reducedDims(cds)[[reduction_method]]),
    msg = paste0('Data has not been processed with',
                 ' chosen reduction_method: ',
                 reduction_method))

  nn_method <- nn_control[['method']]

  if(is.null(cds@reduce_dim_aux[[reduction_method]]) ||
     is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']]) ||
     is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]])) {
    return(FALSE)
  }

  if(nn_method == 'annoy') {
    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_metric']] != nn_control[['metric']]) {
      return(FALSE)
    }

    if(nn_method == 'annoy' && cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['n_trees']] != nn_control[['n_trees']]) {
       return(FALSE)
    }
  }
  else
  if(nn_method == 'hnsw') {
    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_metric']] != nn_control[['metric']]) {
      return(FALSE)
    }

    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_M']] != nn_control[['M']]) { 
       return(FALSE)
    }

    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_ef']] != nn_control[['ef']]) { 
       return(FALSE)
    }
  }

  if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['matrix_id']] != get_reduce_dim_matrix_identity(cds, reduction_method)[['matrix_id']]) {
    return(FALSE)
  }

  return(TRUE)
}


# Clear annoy index.
clear_nn_index <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), nn_method=c('annoy', 'hnsw')) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(nn_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "nn_method must be one of 'annoy' or 'hnsw'")
  nn_method <- match.arg(nn_method)

  if(nn_method == 'annoy') {
    if(!is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']]) &&
       !is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]])) {
      if(!is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']])) {
        cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']] <- NULL
      }
      cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]] <- NULL
    }
  }

  return(cds)
}


#' Verify and set nearest neighbor parameter list.
#'
#' @description Verifies the listed parameter values
#'   that will be passed to the nearest neighbor function
#'   given by nn_control[['method']]. Unspecified
#'   values are set to default values.
#'
#' @section Optional nn_control parameters:
#' \describe{
#'   \item{method}{The method used to find nearest neighbor points.
#'      The available methods are 'nn2', 'annoy', and 'hnsw'. Must
#'      be specified.
#'      Detailed information about each can be found on the WWW sites:
#'      https://cran.r-project.org/web/packages/RANN/,
#'      https://cran.r-project.org/web/packages/RcppAnnoy/index.html,
#'      and https://cran.rstudio.com/web/packages/RcppHNSW/index.html.}
#'   \item{metric}{The distance metric used by the nearest neighbor functions.
#'      Annoy accepts 'euclidean', 'cosine', 'manhattan', and 'hamming'.
#'      HNSW accepts 'euclidean', 'l2', 'cosine', and 'ip'. RANN
#'      uses 'euclidean'.}
#'   \item{n_trees}{The annoy index build parameter that affects the build
#'      time and index size. Larger values give more accurate results,
#'      longer build times, and larger indexes. The default is 50.}
#'   \item{search_k}{The annoy index search parameter that affects the
#'      search accuracy and time. Larger values give more accurate results
#'      and longer search times. Default is 100 * k.}
#'   \item{M}{The HNSW index build parameter that affects the
#'      search accuracy and memory requirements. Larger values give more
#'      accurate search results and increase the index memory use. Default
#'      is 48.}
#'   \item{ef_construction}{The HNSW index build parameter that affects
#'      the search accuracy and index build time. Larger values give more 
#'      accurate search results and longer build times. Default is 200.}
#'   \item{ef}{The HNSW index search parameter that affects the search
#'      accuracy and search time. Larger values give more accurate results
#'      and longer search times. Default is 10.}
#'   \item{grain_size}{The annoy and HNSW parameter the gives the minimum
#'      amount of work to do per thread. Default is 1.}
#'   \item{cores}{The annoy and HNSW parameter that gives the number of
#'      threads to use for the annoy index search and for the HNSW index
#'      build and search. Default is 1.}
#' }
#'
#' @param nn_control list A list of parameters passed
#'   to the nearest neighbor function specified by
#'   nn_control[['method']]. The nn_control list can be empty in which
#'   case the variables are set to their default values.
#'   See below for a description of the valid named values and
#'   their default values.
#' @param k integer The number of desired nearest neighbor
#'   points to return from a search. This value is used only
#'   to set the default annoy search_k parameter. The value
#'   is ignored for index builds. Default is 25.
#'
#' @return An updated nn_control list.
#' @export
set_nn_control <- function(nn_control=list(), k=25, method_default=NULL) {
    assertthat::assert_that(methods::is(nn_control, "list"))
    assertthat::assert_that(all(names(nn_control) %in%
                                  c('method',
                                    'metric',
                                    'n_trees',
                                    'search_k',
                                    'M',
                                    'ef_construction',
                                    'ef',
                                    'grain_size',
                                    'cores')),
                            msg = "Unknown variable in nn_control")

  assertthat::assert_that(!is.null(method_default) && method_default %in% c('nn2', 'annoy', 'hnsw'),
                          msg=paste0("method_default must be one of 'nn2', 'annoy', or 'hnsw'"))

  nn_control[['method']] <- ifelse(is.null(nn_control[['method']]), method_default, nn_control[['method']])

  if(nn_control[['method']] == 'nn2') {
    return(nn_control=nn_control)
  } else
  if(nn_control[['method']] == 'annoy') {
    nn_control[['metric']] <- ifelse(is.null(nn_control[['metric']]), 'euclidean', nn_control[['metric']])
    assertthat::assert_that(nn_control[['metric']] %in% c('euclidean', 'cosine', 'manhattan', 'hamming'),
                            msg=paste0("nearest neighbor method for annoy must be one of 'euclidean', 'cosine', 'manhattan', or 'hamming'"))

    nn_control[['n_trees']] <- ifelse(is.null(nn_control[['n_trees']]), 50, nn_control[['n_trees']])
    nn_control[['search_k']] <- ifelse(is.null(nn_control[['search_k']]), 100 * k, nn_control[['search_k']])
    nn_control[['grain_size']] <- ifelse(is.null(nn_control[['grain_size']]), 1, nn_control[['grain_size']])
    nn_control[['cores']] <- ifelse(is.null(nn_control[['cores']]), 1, nn_control[['cores']])

    assertthat::assert_that(assertthat::is.count(nn_control[['n_trees']]))
    assertthat::assert_that(assertthat::is.count(nn_control[['search_k']]))
    assertthat::assert_that(assertthat::is.count(nn_control[['grain_size']]))
    assertthat::assert_that(assertthat::is.count(nn_control[['cores']]))
  } else
  if(nn_control[['method']] == 'hnsw') {
    nn_control[['metric']] <- ifelse(is.null(nn_control[['metric']]), 'euclidean', nn_control[['metric']])

    assertthat::assert_that(nn_control[['metric']] %in% c('euclidean', 'l2', 'cosine', 'ip'),
                            msg=paste0("nearest neighbor method for HNSW must be one of 'euclidean', 'l2', 'cosine', or 'ip'"))

    nn_control[['M']] <- ifelse(is.null(nn_control[['M']]), 48, nn_control[['M']])
    nn_control[['ef_construction']] <- ifelse(is.null(nn_control[['ef_construction']]), 200, nn_control[['ef_construction']])

    nn_control[['ef']] <- ifelse(is.null(nn_control[['ef']]), 10, nn_control[['ef']])

    nn_control[['grain_size']] <- ifelse(is.null(nn_control[['grain_size']]), 1, nn_control[['grain_size']])
    nn_control[['cores']] <- ifelse(is.null(nn_control[['cores']]), 1, nn_control[['cores']])

    assertthat::assert_that(assertthat::is.count(nn_control[['M']]))
    assertthat::assert_that(assertthat::is.count(nn_control[['ef_construction']]))

    assertthat::assert_that(assertthat::is.count(nn_control[['ef']]))

    assertthat::assert_that(assertthat::is.count(nn_control[['grain_size']]))
    assertthat::assert_that(assertthat::is.count(nn_control[['cores']]))

    assertthat::assert_that(nn_control[['M']] >= 2,
                            msg=paste0('HNSW nearest neighbor M parameter must be >= 2'))

    assertthat::assert_that(nn_control[['ef']] >= k,
                            msg=paste0('HNSW nearest neighbor ef parameter must be >= k'))
  }

  return(nn_control)
}


#' @export
make_nn_index <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), nn_control=list()) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  if(reduction_method == 'PCA') {
    assertthat::assert_that(!is.null(reducedDims(cds)[['PCA']]),
                            msg = paste("When reduction_method = 'PCA', the",
                                        "cds must have been preprocessed for",
                                        "PCA. Please run preprocess_cds with",
                                        "method = 'PCA' before running",
                                        "make_nn_index with",
                                        "reduction_method = 'PCA'."))
  } else
  if(reduction_method == 'LSI') {
    assertthat::assert_that(!is.null(reducedDims(cds)[['LSI']]),
                            msg = paste("When reduction_method = 'LSI', the",
                                        "cds must have been preprocessed for",
                                        "LSI. Please run preprocess_cds with",
                                        "method = 'LSI' before running",
                                        "make_nn_index with",
                                        "reduction_method = 'LSI'."))
  } else
  if(reduction_method == 'Aligned') {
    assertthat::assert_that(!is.null(reducedDims(cds)[['Aligned']]),
                            msg = paste("When reduction_method = 'Aligned', the",
                                        "cds must have been aligned.",
                                        "Please run align_cds before running",
                                        "make_nn_index with",
                                        "reduction_method = 'Aligned'."))
  } else
  if(reduction_method == 'tSNE') {
    assertthat::assert_that(!is.null(reducedDims(cds)[['tSNE']]),
                            msg = paste("When reduction_method = 'tSNE', the",
                                        "cds must have been processed with.",
                                        "reduce_dimension.",
                                        "Please run reduce_dimension with",
                                        "method='tSNE' before running",
                                        "make_nn_index."))
  } else
  if(reduction_method == 'UMAP') {
    assertthat::assert_that(!is.null(reducedDims(cds)[['UMAP']]),
                            msg = paste("When reduction_method = 'UMAP', the",
                                        "cds must have been processed with",
                                        "reduce_dimension.",
                                        "Please run reduce_dimension with",
                                        "method='UMAP' before running",
                                        "make_nn_index."))
  }

  reduced_matrix <- reducedDims(cds)[[reduction_method]]

  if(is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']])) {
    cds@reduce_dim_aux[[reduction_method]][['nn_index']] <- SimpleList()
  }

  nn_method <- nn_control[['method']]

  if(nn_control[['method']] == 'annoy') {
message('make_nn_index: annoy: start')
    annoy_index <- uwot:::annoy_build(X=reduced_matrix,
                                      metric=nn_control[['metric']],
                                      n_trees=nn_control[['n_trees']])

    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]] <- SimpleList()
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']] <- annoy_index
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_metric']] <- nn_control[['metric']]
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['n_trees']] <- nn_control[['n_trees']]
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_ndim']] <- ncol(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['matrix_id']] <- get_reduce_dim_matrix_identity(cds, reduction_method)[['matrix_id']]
message('make_nn_index: annoy: done')
  } else
  if(nn_control[['method']] == 'hnsw') {
    if(nn_control[['M']] > nrow(reduced_matrix))
      nn_control[['M']] <- nrow(reduced_matrix)

message('make_nn_index: hnsw')
    hnsw_index <- RcppHNSW:::hnsw_build(X=reduced_matrix,
                                        distance=nn_control[['metric']],
                                        M=nn_control[['M']],
                                        ef=nn_control[['ef_construction']],
                                        n_threads=nn_control[['cores']],
                                        grain_size=nn_control[['grain_size']])

    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]] <- SimpleList()
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']] <- list()
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['version']] <- packageVersion('RcppHNSW')
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['ann']] <- hnsw_index
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_metric']] <- nn_control[['metric']]
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_M']] <- nn_control[['M']]
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_ef']] <- nn_control[['ef_construction']]
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_ndim']] <- ncol(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['matrix_id']] <- get_reduce_dim_matrix_identity(cds, reduction_method)[['matrix_id']]
  }

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
#' @param ... Parameters to pass through to annoy_search.
#'
#' @export
search_nn_index <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), k, nn_control) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(
    !is.null(reducedDims(cds)[[reduction_method]]),
    msg = paste("Data has not been processed with",
                "chosen reduction_method:",
                reduction_method))

  nn_method <- nn_control[['method']]
  assertthat::assert_that(!is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']]) &&
                          !is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]]),
    msg = paste0('The ', nn_method, ' index for ',
                 reduction_method,
                 ' does not exist. Please process the',
                 ' cds as required.'))

  query_matrix <- reducedDim(cds, reduction_method)
  k <- min(k, nrow(query_matrix))

  if(nn_method == 'annoy') {
    # notes:
    #   o  there may be dependency on uwot version
    #   o  ensure that ann_res[[1]] are indices and ann_res[[2]] are distances
    #   o  set list names to nn.idx and nn.dists for compatibility with
    #      RANN::nn2()
    #
    ann_index <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['ann']]
    tmp <- uwot:::annoy_search(X=query_matrix, k=k, ann=ann_index, nn_control[['search_k']], n_threads=nn_control[['cores']], nn_control[['grain_size']])
    ann_res <- list(nn.idx = tmp[['idx']], nn.dists = tmp[['dist']])
  } else
  if(nn_method == 'hnsw') {
    ann_index <- ann_index <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['ann']]
    tmp <- RcppHNSW:::hnsw_search(X=query_matrix, ann=ann_index, k=k, ef=nn_control[['ef']], n_threads=nn_control[['cores']], grain_size=nn_control[['grain_size']])
    ann_res <- list(nn.idx=tmp[['idx']], nn.dists=tmp[['dist']])
  }

  return(ann_res) 
}


search_nn_matrix <- function(mat, k=10, nn_control=list()) {
  method <- nn_control[['method']]
  k <- min(k, nrow(mat))

  if(method == 'nn2') {
    nn_res <- RANN::nn2(mat, mat, min(k, nrow(mat)), searchtype = "standard")
  } else
  if(method == 'annoy') {
    ann_index <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['ann']]
    tmp <- annoy_nn(X=mat,
                    k=k,
                    metric=nn_control[['metric']],
                    n_trees=nn_control[['n_trees']],
                    n_threads=nn_control[['cores']],
                    grain_size=nn_control[['grain_size']])
message('Warning: fix shuffled/missing self-indices.')
    nn_res <- list(nn.idx = tmp[['idx']], nn.dists = tmp[['dist']])
  } else
  if(method == 'hnsw') {
    tmp <- hnsw_knn(X=mat,
                    k=k,
                    distance=nn_control[['metric']],
                    M=nn_control[['M']],
                    ef_construction=nn_control[['ef_construction']],
                    ef=nn_control[['ef']],
                    n_threads=nn_control[['cores']],
                    grain_size=nn_control[['grain_size']])
message('Warning: fix shuffled/missing self-indices.')
    nn_res <- list(nn.idx = tmp[['idx']], nn.dists = tmp[['dist']])
  }
  return(nn_res)
}


# Check whether the neighbor in the first column == row index.
# The hnsw_search help page states
#   Every item in the dataset is considered to be a neighbor of itself,
#   so the first neighbor of item i should always be i itself. If that
#   isn't the case, then any of M, ef_construction and ef may need
#   increasing.
check_nn_col1 <- function(mat) {
  irow <- seq(ncol(mat))
  if(!all(mat[,1] == irow)) {
    return(FALSE)
  }
  return(TRUE)
}


# Given a vector of nearest neighbors from a
# search for which the query and database are
# the same set and the self point is not the
# first element, swap it to the first element.
nn_vector_index_swap <- function(iidx, idst) {
  # Most often, the self point is shifted right
  # one column; that is, there is one other point
  # with zero distance.
  if(iidx[2] != iidx[1]) {
message('nn_vector_index_swap: swap: ', iidx[1])
    if(iidx[3] == iidx[1]) {
      iidx[3] <- iidx[2]
      iidx[2] <- iidx[1]
      return(iidx)
    } else
    {
      for( i in 2:length(iidx)) {
        if(iidx[i] == iidx[1]) {
          iidx[i] <- iidx[2]
          iidx[2] <- iidx[1]
          return(iidx)
        }
      }
      # The self point is not in the iidx so
      # shift the vector and distances right.
      # Notice that idst is written to the
      # parent frame.
      iidx <- c(rep(NA,times=shift), iidx[1:(length(iidx)-1)])
      idst <<- c(rep(NA,times=shift), idst[1:(length(idst)-1)])
      return(iidx)
    }
  }
  iidx
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
# In case that the idx is missing the self-point
# and the distance in the first column is > 0,
# shift the idx and dists vectors right by 1.
#
# Warning: do not use this if the search point set
#          differs from the index set. It will
#          make a mess of the result.
#
swap_nn_self_search_index <- function(nn_res) {
message('swap_nn_self_search_index: start')
  idx <- nn_res[['nn.idx']]
  dst <- nn_res[['nn.dists']]

  iidx <- cbind(1:nrow(idx),idx)
  idst <- cbind(1:nrow(dst),dst)
  shift_res <- t(apply(iidx, 1, nn_vector_index_swap))

  list(nn.idx=shift_res[,-1], nn.dists=idst[,-1])
}


