# Check whether nearest neighbor index exists.
# Note: if cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]] exists
#       then cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']]
#       exists, which is consistent with the clear_cds_nn_index() behavior.
check_cds_nn_index_exists <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), nn_method=c('annoy', 'hnsw')) {
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
  else
  if(nn_method == 'hnsw') {
    if(is.null(cds@reduce_dim_aux[[reduction_method]]) ||
       is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']]) ||
       is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]])) {
       return(FALSE)
    }
  }
  else
    stop('check_cds_nn_index_exists: unsupported nearest neighbor index type \'', nn_method, '\'')

  return(TRUE)
}


# Check whether annoy index exists and is consistent with matrix and parameters.
check_cds_nn_index_is_current <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), nn_control=list(), verbose=FALSE) {

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

  nn_control <- set_nn_control(nn_control=nn_control, k=1, method_default='annoy', verbose=verbose)

  if(verbose) {
    message('check_cds_nn_index_is_current:')
    report_nn_control('  nn_control: ', nn_control=nn_control)
  }

  nn_method <- nn_control[['method']]

  if(is.null(cds@reduce_dim_aux[[reduction_method]]) ||
     is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']]) ||
     is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]])) {
    if(verbose) {
      message('check_cds_nn_index_is_current: FALSE')
    }
    return(FALSE)
  }

  if(nrow(reducedDims(cds)[[reduction_method]]) != cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nrow']]) {
    if(verbose) {
      message('check_cds_nn_index_is_current: FALSE')
    }
    return(FALSE)
  }

  if(ncol(reducedDims(cds)[[reduction_method]]) != cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['ncol']]) {
    if(verbose) {
      message('check_cds_nn_index_is_current: FALSE')
    }
    return(FALSE)
  }

  if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['matrix_id']] != get_reduce_dim_matrix_identity(cds, reduction_method)[['matrix_id']]) {
    if(verbose) {
      message('check_cds_nn_index_is_current: FALSE')
    }
    return(FALSE)
  }

  if(nn_method == 'annoy') {
    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['metric']] != nn_control[['metric']]) {
      if(verbose) {
        message('check_cds_nn_index_is_current: FALSE')
      }
      return(FALSE)
    }

    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['n_trees']] != nn_control[['n_trees']]) {
        if(verbose) {
          message('check_cds_nn_index_is_current: FALSE')
        }
       return(FALSE)
    }
  }
  else
  if(nn_method == 'hnsw') {
    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['metric']] != nn_control[['metric']]) {
      if(verbose) {
        message('check_cds_nn_index_is_current: FALSE')
      }
      return(FALSE)
    }

    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['M']] != nn_control[['M']]) { 
      if(verbose) {
        message('check_cds_nn_index_is_current: FALSE')
      }
      return(FALSE)
    }

    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['ef_construction']] != nn_control[['ef_construction']]) { 
      if(verbose) {
        message('check_cds_nn_index_is_current: FALSE')
      }
      return(FALSE)
    }
  }
  else
    stop('check_cds_nn_index_is_current: unsupported nearest neighbor index type \'', nn_method, '\'')

  if(verbose) {
    message('check_cds_nn_index_is_current: TRUE')
  }

  return(TRUE)
}


# Clear nearest neighbor index.
# Note: remove cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]]
#       because the objects stored in it are related; that is, if the object
#       cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']]
#       is removed or changed, the other objects in
#       cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]] need to be
#       updated.
clear_cds_nn_index <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), nn_method=c('annoy', 'hnsw', 'all')) {
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
    msg = "nn_method must be one of 'annoy', 'hnsw', or 'all'")

  nn_method <- match.arg(nn_method)

  if(nn_method == 'annoy') {
    if(!is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']]) &&
       !is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]])) {
      cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]] <- NULL
    }
  }
  else
  if(nn_method == 'hnsw') {
    if(!is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']]) &&
       !is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]])) {
      cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]] <- NULL
    }
  }
  else
  if(nn_method == 'all') {
    if(!is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']])) {
      cds@reduce_dim_aux[[reduction_method]][['nn_index']] <- SimpleList()
    }
  }
  else
    stop('clear_cds_nn_index: unsupported nearest neighbor index type \'', nn_method, '\'')

  return(cds)
}


check_cds_nn_search_exists <- function(cds, reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"), search_id=NULL, verbose=FALSE)
{
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'UMAP', 'tSNE', 'PCA', 'LSI', 'Aligned'")
  reduction_method <- match.arg(reduction_method)

  if(is.null(reducedDims(cds)[[reduction_method]])) {
    return(FALSE)
  }

  if(is.null(cds@clusters[[reduction_method]])) {
    return(FALSE)
  }

  if(is.null(cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]])) {
    return(FALSE)
  }

  return(TRUE)
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
#'      The available methods are 'nn2', 'annoy', and 'hnsw'. If not
#'      specified, the method given by the method_default parameter is
#'      used.
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
#' @param method_default The nearest neighbor method to use when not
#'   specified in the nn_control list.
#' @param verbose Whether to emit verbose output.
#'
#' @return An updated nn_control list.
#' @export
set_nn_control <- function(nn_control=list(), k=25, method_default=NULL, verbose=FALSE) {

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
                            msg = "set_nn_control: unknown variable in nn_control")

  assertthat::assert_that(!is.null(method_default) && method_default %in% c('nn2', 'annoy', 'hnsw'),
                          msg=paste0("set_nn_control: method_default must be one of 'nn2', 'annoy', or 'hnsw'"))

  nn_control[['method']] <- ifelse(is.null(nn_control[['method']]), method_default, nn_control[['method']])

  assertthat::assert_that(!is.null(nn_control[['method']]) && nn_control[['method']] %in% c('nn2', 'annoy', 'hnsw'),
                          msg=paste0("set_nn_control: method must be one of 'nn2', 'annoy', or 'hnsw'"))

  assertthat::assert_that(assertthat::is.count(k))

  if(nn_control[['method']] == 'nn2') {
    # Do nothing for nn2. Call report_nn_control when verbose <- TRUE
  } else
  if(nn_control[['method']] == 'annoy') {
    nn_control[['metric']] <- ifelse(is.null(nn_control[['metric']]), 'euclidean', nn_control[['metric']])
    assertthat::assert_that(nn_control[['metric']] %in% c('euclidean', 'cosine', 'manhattan', 'hamming'),
                            msg=paste0("set_nn_control: nearest neighbor method for annoy must be one of 'euclidean', 'cosine', 'manhattan', or 'hamming'"))

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
                            msg=paste0("set_nn_control: nearest neighbor method for HNSW must be one of 'euclidean', 'l2', 'cosine', or 'ip'"))

    nn_control[['M']] <- ifelse(is.null(nn_control[['M']]), 48, nn_control[['M']])
    nn_control[['ef_construction']] <- ifelse(is.null(nn_control[['ef_construction']]), 200, nn_control[['ef_construction']])
    nn_control[['ef']] <- ifelse(is.null(nn_control[['ef']]), 30, nn_control[['ef']])

    nn_control[['grain_size']] <- ifelse(is.null(nn_control[['grain_size']]), 1, nn_control[['grain_size']])
    nn_control[['cores']] <- ifelse(is.null(nn_control[['cores']]), 1, nn_control[['cores']])

    assertthat::assert_that(assertthat::is.count(nn_control[['M']]))
    assertthat::assert_that(assertthat::is.count(nn_control[['ef_construction']]))

    assertthat::assert_that(assertthat::is.count(nn_control[['ef']]))

    assertthat::assert_that(assertthat::is.count(nn_control[['grain_size']]))
    assertthat::assert_that(assertthat::is.count(nn_control[['cores']]))

    assertthat::assert_that(nn_control[['M']] >= 2,
                            msg=paste0('set_nn_control: HNSW nearest neighbor M index build parameter must be >= 2'))

    assertthat::assert_that(nn_control[['ef']] >= k,
                            msg=paste0('set_nn_control: HNSW nearest neighbor ef index search parameter must be >= k (',k,')'))
  }
  else
    stop('set_nn_control: unsupported nearest neighbor method \'', nn_control[['method']], '\'')

  if(verbose) {
    cs <- get_call_stack_as_string()
    message('set_nn_control: call stack: ', cs)
    report_nn_control('  nn_control: ', nn_control)
  }

  return(nn_control)
}


report_nn_control <- function(label=NULL, nn_control=list()) {

  indent <- ''
  if(!is.null(label)) {
    indent <- '  '
  }

  message(ifelse(!is.null(label), label, ''))

  message(indent, '  method: ', ifelse(!is.null(nn_control[['method']]), nn_control[['method']], 'unset'))
  message(indent, '  metric: ', ifelse(!is.null(nn_control[['metric']]), nn_control[['metric']], 'unset'))

  if(is.null(nn_control[['method']])) {
    return()
  }

  if(nn_control[['method']] == 'nn2') {
    message(indent, '  no nn2 parameters')
  }
  else
  if(nn_control[['method']] == 'annoy') {
    message(indent, '  n_trees: ', ifelse(!is.null(nn_control[['n_trees']]), nn_control[['n_trees']], 'unset'))
    message(indent, '  search_k: ', ifelse(!is.null(nn_control[['search_k']]), nn_control[['search_k']], 'unset'))
    message(indent, '  cores: ', ifelse(!is.null(nn_control[['cores']]), nn_control[['cores']], 'unset'))
    message(indent, '  grain_size: ', ifelse(!is.null(nn_control[['grain_size']]), nn_control[['grain_size']], 'unset'))
  }
  else
  if(nn_control[['method']] == 'hnsw') {
    message(indent, '  M: ', ifelse(!is.null(nn_control[['M']]), nn_control[['M']], 'unset'))
    message(indent, '  ef_construction: ', ifelse(!is.null(nn_control[['ef_construction']]), nn_control[['ef_construction']], 'unset'))
    message(indent, '  ef: ', ifelse(!is.null(nn_control[['ef']]), nn_control[['ef']], 'unset'))
    message(indent, '  cores: ', ifelse(!is.null(nn_control[['cores']]), nn_control[['cores']], 'unset'))
    message(indent, '  grain_size: ', ifelse(!is.null(nn_control[['grain_size']]), nn_control[['grain_size']], 'unset'))
  }
  else
    stop('report_nn_control: unsupported nearest neighbor method \'', nn_method, '\'')
}


#' Make and store a nearest neighbor index in the CDS.
#'
#' @description Make the a nearest neighbor index from the specified
#'  reduction_method matrix using either the default nearest neighbor
#'  or the method specified in the nn_control list.
#'
#' @param cds The cell_data_set upon which to perform this operation
#' @param reduction_method A string indicating the reduced matrix that you
#'  want to use to make the nearest neighbor index.
#' @param nn_control A list of parameters used to make the nearest
#'  neighbor index. See the set_nn_control help for detailed information.
#' @param verbose Whether to emit verbose output.
#'
#' @return The cds cell data set updated with the nearest neighbor index.
#' @export
make_nn_index <- function(subject_matrix, nn_control=list(), verbose=FALSE) {

  nn_control <- set_nn_control(nn_control=nn_control, k=1, method_default='annoy', verbose=verbose)

  if(verbose) {
    message('make_nn_index:')
    report_nn_control('  nn_control: ', nn_control)
    tic('make_nn_index: build time')
  }

  nn_method <- nn_control[['method']]

  if(nn_method == 'annoy') {
    nn_index <- uwot:::annoy_build(X=subject_matrix,
                                   metric=nn_control[['metric']],
                                   n_trees=nn_control[['n_trees']],
                                   verbose=verbose)
  }
  else
  if(nn_method == 'hnsw') {
    nn_index <- RcppHNSW:::hnsw_build(X=subject_matrix,
                                      distance=nn_control[['metric']],
                                      M=nn_control[['M']],
                                      ef=nn_control[['ef_construction']],
                                      verbose=verbose,
                                      n_threads=nn_control[['cores']],
                                      grain_size=nn_control[['grain_size']])
  }
  else
    stop('make_nn_index: unsupported nearest neighbor index type \'', nn_method, '\'')

  if(verbose) {
    toc()
  }

  return(nn_index)
}


# Notes:
#   o  a reducedDims(cds)[[reduction_method]] matrix must be used to make the index
#   o  set_cds_nn_index assumes that the matrix currently in reducedDims(cds)[[reduction_method]]
#      was used to make the index so it's safest for make_nn_index to precede directly the call to
#      set_cds_nn_index
set_cds_nn_index <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), nn_index, nn_control=list(), verbose=FALSE) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste0("When reduction_method = '", reduction_method,
                                      "' the cds must have been processed for it.",
                                      " Please run the required processing function",
                                      " before this one."))

  reduced_matrix <- reducedDims(cds)[[reduction_method]]

  nn_control <- set_nn_control(nn_control=nn_control, k=1, method_default='annoy', verbose=verbose)

  nn_method <- nn_control[['method']]

  if(nn_method == 'annoy') {
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]] <- SimpleList()
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']] <- nn_index
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['metric']] <- nn_control[['metric']]
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['n_trees']] <- nn_control[['n_trees']]
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nrow']] <- nrow(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['ncol']] <- ncol(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['index_version']] <- packageVersion('uwot')
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['matrix_id']] <- get_reduce_dim_matrix_identity(cds, reduction_method)[['matrix_id']]
  }
  else
  if(nn_method == 'hnsw') {
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]] <- SimpleList()
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']] <- list()
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['ann']] <- nn_index
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['metric']] <- nn_control[['metric']]
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['M']] <- nn_control[['M']]
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['ef_construction']] <- nn_control[['ef_construction']]
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nrow']] <- nrow(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['ncol']] <- ncol(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['index_version']] <- packageVersion('RcppHNSW')
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['matrix_id']] <- get_reduce_dim_matrix_identity(cds, reduction_method)[['matrix_id']]
  }
  else
    stop('set_cds_nn_index: unsupported nearest neighbor index type \'', nn_method, '\'')

  return(cds)
}


get_cds_nn_index <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), nn_control=list(), verbose=FALSE) {

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

  # set k to a dummy value
  k <- 1
  nn_control <- set_nn_control(nn_control=nn_control, k=k, method_default='annoy', verbose=verbose)

  nn_method <- nn_control[['method']]

  assertthat::assert_that(nn_method %in% c('annoy', 'hnsw'),
    msg = paste0('get_cds_nn_index: unsupported nearest neighbor index type \'',
                nn_method, '\'.'))

  if(nn_method == 'annoy') {
    nn_index <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['ann']]
  }
  else
  if(nn_method == 'hnsw') {
    nn_index <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['ann']]
  }
  else
    stop('get_cds_nn_index: unsupported nearest neighbor index type \'', nn_method, '\'')

return(nn_index)
}


#  Search nearest neighbor index for cells near the cells in the reduced
#  dim matrix.
# 
# 
#  @param query_matrix A matrix used to find the nearest neighbors in the
#   index created from the reduction_method matrix.
#  @param nn_index A nearest_neighbor index.
#  @param k An integer for the number of nearest neighbors to return for
#    each cell. Default is 25.
#  @param nn_control A list of parameters used to make the nearest
#   neighbor index. See the set_nn_control help for detailed information.
#  @param verbose Whether to emit verbose output.
# 
#  @return A list list(nn.idx, nn.dists) where nn.idx is
#   a matrix of nearest neighbor indices and nn.dists is a matrix of
#   the distance between the index given by the row number and the index
#   given in nn.idx. If the same reduced dim matrix is used to make the
#   index and search the index, the index given by the row number should
#   be in the row, usually in the first column.
search_nn_index <- function(query_matrix, nn_index, k=25, nn_control=list(), verbose=FALSE) {

  nn_control <- set_nn_control(nn_control=nn_control, k=k, method_default='annoy', verbose=verbose)

  nn_method <- nn_control[['method']]

  assertthat::assert_that(nn_method %in% c('annoy', 'hnsw'),
    msg = paste0('search_nn_index: unsupported nearest neighbor index type \'',
                nn_method, '\'.'))

  if(verbose) {
    message('search_nn_index:')
    message('  k: ', k)
    report_nn_control('  nn_control: ', nn_control)
    tic('search_nn_index: search time:')
  }

  k <- min(k, nrow(query_matrix))

  if(nn_method == 'annoy') {
    # notes:
    #   o  there may be dependency on uwot version
    #   o  ensure that nn_res[[1]] are indices and nn_res[[2]] are distances
    #   o  set list names to nn.idx and nn.dists for compatibility with
    #      RANN::nn2()
    #
    tmp <- uwot:::annoy_search(X=query_matrix,
                               k=k,
                               ann=nn_index,
                               search_k=nn_control[['search_k']],
                               n_threads=nn_control[['cores']],
                               nn_control[['grain_size']], verbose=verbose)
    nn_res <- list(nn.idx = tmp[['idx']], nn.dists = tmp[['dist']])
  }
  else
  if(nn_method == 'hnsw') {
    assertthat::assert_that(nn_control[['ef']] >= k,
      msg=paste0('search_nn_index: ef must be >= k'))

    tmp <- RcppHNSW:::hnsw_search(X=query_matrix,
                                  ann=nn_index,
                                  k=k,
                                  ef=nn_control[['ef']],
                                  verbose=verbose,
                                  n_threads=nn_control[['cores']],
                                  grain_size=nn_control[['grain_size']])
    nn_res <- list(nn.idx=tmp[['idx']], nn.dists=tmp[['dist']])
  }
  else
    stop('search_nn_index: unsupported nearest neighbor index type \'', nn_method, '\'')

  if(verbose) {
    toc()
  }

  return(nn_res) 

}


# Notice that this function uses the available reducedDims(cds)[[reduction_method]] matrix
# for verifying the query matrix so be certain that the search (and index build) were run on
# this matrix so call this function immediately after running the search.
set_cds_nn_search <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), search_id=NULL, k=NULL, nn_control=NULL, verbose=TRUE) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'  xxx")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste0("When reduction_method = '", reduction_method,
                                      "' the cds must have been processed for it.",
                                      " Please run the required processing function",
                                      " before this one."))
  assertthat::assert_that(!is.null(search_id),
                          msg = paste0('You must give a search_id string.'))

  assertthat::assert_that(!is.null(k),
                          msg = paste0('You must give a k value.'))

  reduced_matrix <- reducedDims(cds)[[reduction_method]]

  nn_control <- set_nn_control(nn_control=nn_control, k=1, method_default='annoy', verbose=verbose)

  nn_method <- nn_control[['method']]

  # Example: cds@reduce_dim_aux[[reduction_method]][['nn_search']][[nn_method]][[search_id]][['search_k']]

  if(is.null(cds@reduce_dim_aux[[reduction_method]][['nn_search']])) {
    cds@reduce_dim_aux[[reduction_method]][['nn_search']] <- SimpleList()
  }

  if(nn_method == 'nn2') {
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]] <- SimpleList()
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['method']] <- nn_method
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['k']] <- nn_control[['k']]
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['nrow']] <- nrow(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['ncol']] <- ncol(reduced_matrix)
  }
  else
  if(nn_method == 'annoy') {
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]] <- SimpleList()
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['method']] <- nn_method
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['search_k']] <- nn_control[['search_k']]
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['k']] <- nn_control[['k']]
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['nrow']] <- nrow(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['ncol']] <- ncol(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['grain_size']] <- nn_control[['grain_size']]
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['cores']] <- nn_control[['cores']]
  }
  else
  if(nn_method == 'hnsw') {
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]] <- SimpleList()
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['method']] <- nn_method
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['ef']] <- nn_control[['ef']]
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['k']] <- nn_control[['k']]
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['nrow']] <- nrow(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['ncol']] <- ncol(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['grain_size']] <- nn_control[['grain_size']]
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['cores']] <- nn_control[['cores']]
  }
  else {
    stop('set_cds_nn_search: unsupported nearest neighbor index type \'', nn_method, '\'')
  }

  return(cds)
}


get_cds_nn_search <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), search_id=NULL, verbose=TRUE) {

  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)
  
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste0("When reduction_method = '", reduction_method,
                                      "' the cds must have been processed for it.",
                                      " Please run the required processing function",
                                      " before this one."))
  assertthat::assert_that(!is.null(search_id),
                          msg = paste0('You must give a search_id string.'))

  assertthat::assert_that(!is.null(cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]]),
                          msg = paste0("Stopping because there is no information in the cds for search_id = '",
                                       search_id, "'."))

  nn_search <- cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]]

  nn_method <- nn_search[['method']]

  if(nn_method == 'nn2') {
    nn_control <- list(method=nn_method,
                       grain_size=nn_search[['grain_size']],
                       cores=nn_search[['cores']])
  }
  else
  if(nn_method == 'annoy') {
    nn_control <- list(method=nn_method,
                       search_k=nn_search[['search_k']],
                       grain_size=nn_search[['grain_size']],
                       cores=nn_search[['cores']])
  }
  else
  if(nn_method == 'hnsw') {
    nn_control <- list(method=nn_method,
                       ef=nn_search[['ef']],
                       grain_size=nn_search[['grain_size']],
                       cores=nn_search[['cores']])
  }
  else {
    stop('get_cds_nn_search: unsupported nearest neighbor search type \'', nn_method, '\'')
  }
                     
  return(nn_control)
}


get_cds_nn_control <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), search_id=NULL, verbose=TRUE) {

  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste0("When reduction_method = '", reduction_method,
                                      "' the cds must have been processed for it.",
                                      " Please run the required processing function",
                                      " before this one."))
  assertthat::assert_that(!is.null(search_id),
                          msg = paste0('You must give a search_id string.'))

  assertthat::assert_that(!is.null(cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]]),
                          msg = paste0("Stopping because there is no information in the cds for search_id = '",
                                       search_id, "'."))

  nn_search <- cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]]

  nn_method <- nn_search[['method']]

  if(nn_method == 'nn2') {
    nn_control <- list(method=nn_method,
                       grain_size=nn_search[['grain_size']],
                       cores=nn_search[['cores']])
  }
  else
  if(nn_method == 'annoy') {
    nn_control <- list(method=nn_method,
                       metric=cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['metric']],
                       n_trees=cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['n_trees']],
                       search_k=nn_search[['search_k']],
                       grain_size=nn_search[['grain_size']],
                       cores=nn_search[['cores']])
  }
  else
  if(nn_method == 'hnsw') {
    nn_control <- list(method=nn_method,
                       metric=cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['metric']],
                       M=cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['M']],
                       ef_construction=cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['ef_construction']],
                       ef=nn_search[['ef']],
                       grain_size=nn_search[['grain_size']],
                       cores=nn_search[['cores']])
  }
  else {
    stop('get_cds_nn_control: unsupported nearest neighbor search type \'', nn_method, '\'')
  }

  if(verbose) {
    cs <- get_call_stack_as_string()
    message('get_cds_nn_control: (',cs,')')
    report_nn_control('  nn_control: ', nn_control=nn_control)
  }

  return(nn_control)
}


search_nn_matrix <- function(subject_matrix, query_matrix, k=25, nn_control=list(), verbose=FALSE) {

  assertthat::assert_that(is.matrix(subject_matrix),
    msg=paste0('search_nn_matrix: the subject_matrix object must be of type matrix'))

  assertthat::assert_that(is.matrix(query_matrix),
    msg=paste0('search_nn_matrix: the query_matrix object must be of type matrix'))

  nn_control <- set_nn_control(nn_control=nn_control, k=k, method_default='annoy', verbose=verbose)

  method <- nn_control[['method']]
  k <- min(k, nrow(subject_matrix))

  if(verbose) {
    message('search_nn_matrix:')
    message('  k: ', k)
    report_nn_control('  nn_control: ', nn_control)
    tic('search_nn_matrix: search_time')
  }

  if(method == 'nn2') {
    nn_res <- RANN::nn2(subject_matrix, query_matrix, min(k, nrow(subject_matrix)), searchtype = "standard")
  } else
  if(method == 'annoy') {
    nn_index <- uwot:::annoy_build(X=subject_matrix,
                                   metric=nn_control[['metric']],
                                   n_trees=nn_control[['n_trees']],
                                   verbose=verbose)

    tmp <- uwot:::annoy_search(X=query_matrix,
                               k=k,
                               ann=nn_index,
                               search_k=nn_control[['search_k']],
                               prep_data=TRUE,
                               tmpdir=tempdir(),
                               n_threads=nn_control[['cores']],
                               grain_size=nn_control[['grain_size']],
                               verbose=verbose)
    nn_res <- list(nn.idx = tmp[['idx']], nn.dists = tmp[['dist']])
    swap_nn_row_index_point(nn_res, verbose=verbose)
  } else
  if(method == 'hnsw') {
    progress <- 'bar'
    nn_index <- RcppHNSW:::hnsw_build(X=subject_matrix,
                                      distance=nn_control[['metric']],
                                      M=nn_control[['M']],
                                      ef=nn_control[['ef_construction']],
                                      verbose=verbose,
                                      progress=progress,
                                      n_threads=nn_control[['cores']],
                                      grain_size=nn_control[['grain_size']])
  
    tmp <- RcppHNSW:::hnsw_search(X=query_matrix,
                                  ann=nn_index,
                                  k=k,
                                  ef=nn_control[['ef']],
                                  verbose=verbose,
                                  progress=progress,
                                  n_threads=nn_control[['cores']],
                                  grain_size=nn_control[['grain_size']])
    nn_res <- list(nn.idx = tmp[['idx']], nn.dists = tmp[['dist']])
    swap_nn_row_index_point(nn_res, verbose=verbose)
  }
  else
    stop('search_nn_matrix: unsupported nearest neighbor method \'', nn_method, '\'')

  if(verbose)
    toc()

  return(nn_res)

}


# Check whether the neighbor in the first column == row index.
# The hnsw_search help page states
#   Every item in the dataset is considered to be a neighbor of itself,
#   so the first neighbor of item i should always be i itself.
check_nn_col1 <- function(mat) {
  irow <- seq(nrow(mat))
  if(any(mat[,1] != irow)) {
    return(FALSE)
  }
  return(TRUE)
}


# When the search query matrix is the same as the
# index matrix, we expect that the row index
# exists in the row. Count the cases for which
# this is not true and the distances are not
# all zero.
count_nn_missing_self_index <- function(nn_res) {
  idx <- nn_res[['nn.idx']]
  dst <- nn_res[['nn.dists']]
  len <- length(idx[[1]])
  num_missing <- 0

  for (irow in 1:nrow(idx)) {
    vidx <- idx[irow,]
    vdst <- dst[irow,]
    if(vidx[[1]] != irow) {
      dself <- FALSE
      dzero <- TRUE
      for(i in 1:len) {
        if(vidx[[i]] == irow) {
          dself = TRUE
          break
        }
        if(vdst[[i]] > .Machine$double.xmin) {
          dzero <- FALSE
        }
      }

      if(dself == FALSE && dzero == FALSE) {
        num_missing <- num_missing + 1
      }
    }
  }

  message('count_nn_missing_self_index:')
  message('  ', num_missing, ' out of ', nrow(nn_res[['nn.idx']]), ' rows are missing the row index')
  message('  \'recall\': ', formatC((as.double(nrow(nn_res[['nn.idx']]) - num_missing) / as.double(nrow(nn_res[['nn.idx']]))) * 100.0), "%")

  return(0)
}


# For searches in which the index and search sets
# are the same, sort nearest neighbor results so
# that row index point is in the first column. If the
# row index point was missed, shift the results to the
# right one place and store the row index point in the
# first element.
#
# This is required because approximate searches may
# not store the row index point in the first column if
# there is more than one point with zero distance,
# or, the function may miss the row index point.
#
# Of course, do not use this if the search point set
# differs from the index set.
#
swap_nn_row_index_point <- function(nn_res, verbose=FALSE) {

  if(verbose) {
    count_nn_missing_self_index(nn_res)
  }

  # Skip if there is no need to swap indices.
  if(check_nn_col1(nn_res$nn.idx))
    return(nn_res)

  idx <- nn_res[['nn.idx']]
  dst <- nn_res[['nn.dists']]

  diagnostics <- FALSE

  for (irow in 1:nrow(idx)) {
    vidx <- idx[irow,]
    vdst <- dst[irow,]
    if(vidx[[1]] != irow) {
      if(diagnostics) {
        message('swap_nn_row_index_point: adjust nn matrix row: ', irow)
        message('swap_nn_row_index_point: idx row pre fix: ', paste(vidx, collapse=' '))
      }

      if(vidx[[2]] == irow) {
        vidx[[2]] <- vidx[[1]]
        vidx[[1]] <- irow
      } else {
        match <- FALSE
        for(i in 2:length(vidx)) {
          if(vidx[[i]] == irow) {
            vidx[[i]] <- vidx[[1]]
            vidx[[1]] <- irow
            match <- TRUE
            break
          }
        }
        if(!match) {
          if(diagnostics) {
            message('swap_nn_row_index_point: dst row pre fix: ', paste(vdst, collapse=' '))
          }
          vidx <- c(irow, vidx[1:(length(vidx)-1)])
          vdst <- c(0, vdst[1:(length(vdst)-1)])
          dst[irow,] <- vdst
          if(diagnostics) {
            message('swap_nn_row_index_point: dst row post fix: ', paste(vdst, collapse=' '))
          }
        }
      }
      idx[irow,] <- vidx
      if(diagnostics) {
        message('swap_nn_row_index_point: idx row post fix: ', paste(vidx, collapse=' '))
      }
    }
  }

  nn_res_out <- list(nn.idx=idx, nn.dists=dst)

  return(nn_res_out)
}

