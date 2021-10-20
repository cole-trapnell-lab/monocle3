# Notes:
#   To add nearest neighbor method, modify
#     o  many of the functions in this file
#     o  learn_graph: learn_graph_control
#     o  io.R: save and load indices
#     o  zzz.R: possibly change defaults
#     o  others?


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

  # We check the index build parameters.
  nn_control <- set_nn_control(mode=1, nn_control=nn_control, k=1, nn_control_default=list(), verbose=verbose)

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


# Check whether nearest neighbor search information exists
# for reduction_method. This is not useful at this time. The
# intention was to store information about the search
# parameters used in cluster_cells() and use them in for
# nearest neighbor searches in later function calls. However,
# this does not seem useful at this time.
check_cds_nn_search_exists <- function(cds, reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"), search_id, verbose=FALSE)
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


#' @title Verify and set nearest neighbor parameter list.
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
#' @param mode the nearest neighbor operation for which the nn_control
#'   list will be used. 1=make index, 2=search index, and 3=both make
#'   and search index. The default mode value is 3.
#' @param nn_control an optional list of parameters passed to
#'   the nearest neighbor function specified by nn_control[['method']].
#'   The nn_control list can be empty in which case defaults are used
#'   as follows: if the nn_control_default list has non-zero length, it
#'   is used; otherwise the annoy nearest neighbor method is used with
#'   metric='euclidean', n_trees=50, and search_k=k * n_trees. If not
#'   all of the values required by a nearest neighbor method are given
#    in nn_control, the missing values are set to the default values
#'   listed below.
#' @param k integer the number of desired nearest neighbor points to
#'   return from a search. This value is used only to set the annoy
#'   search_k parameter when search_k is not given in nn_control. The
#'   value is ignored for index builds. The default is 25.
#' @param nn_control_default an optional nn_control list to use when the
#'   nn_control parameter has length zero.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return an updated nn_control list.
#' @export
set_nn_control <- function(mode=3, nn_control=list(), k=25, nn_control_default=list(), verbose=FALSE) {

  assertthat::assert_that(methods::is(nn_control, "list"))
  assertthat::assert_that(methods::is(nn_control_default, "list"))
  assertthat::assert_that(assertthat::is.count(k))

  if(length(nn_control) > 0 )
    nn_control <- nn_control
  else
  if(length(nn_control_default) > 0)
    nn_control <- nn_control_default
  else {
    nn_control <- list(method = 'annoy',
                       metric = 'euclidean')
  }
 
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

  assertthat::assert_that(assertthat::is.count(mode) && mode >= 1 && mode <= 3,
                          msg = paste0("set_nn_control: invalid mode value. Mode must be an integer with the value 1, 2, or 3."))

  if(nn_control[['method']] == 'nn2') {
    # Do nothing for nn2. Call report_nn_control when verbose <- TRUE
  } else
  if(nn_control[['method']] == 'annoy') {
    nn_control[['metric']] <- ifelse(is.null(nn_control[['metric']]), 'euclidean', nn_control[['metric']])
    assertthat::assert_that(nn_control[['metric']] %in% c('euclidean', 'cosine', 'manhattan', 'hamming'),
                            msg=paste0("set_nn_control: nearest neighbor metric for annoy must be one of 'euclidean', 'cosine', 'manhattan', or 'hamming'"))

    if(bitwAnd(mode, 1)) {
      nn_control[['n_trees']] <- ifelse(is.null(nn_control[['n_trees']]), 50, nn_control[['n_trees']])
      assertthat::assert_that(assertthat::is.count(nn_control[['n_trees']]))
    }

    if(bitwAnd(mode, 2)) {
      nn_control[['n_trees']] <- ifelse(is.null(nn_control[['n_trees']]), 50, nn_control[['n_trees']])
      nn_control[['search_k']] <- ifelse(is.null(nn_control[['search_k']]), nn_control[['n_trees']] * k, nn_control[['search_k']])
      nn_control[['grain_size']] <- ifelse(is.null(nn_control[['grain_size']]), 1, nn_control[['grain_size']])
      nn_control[['cores']] <- ifelse(is.null(nn_control[['cores']]), 1, nn_control[['cores']])
      assertthat::assert_that(assertthat::is.count(nn_control[['search_k']]))
      assertthat::assert_that(assertthat::is.count(nn_control[['grain_size']]))
      assertthat::assert_that(assertthat::is.count(nn_control[['cores']]))
    }

  } else
  if(nn_control[['method']] == 'hnsw') {
    nn_control[['metric']] <- ifelse(is.null(nn_control[['metric']]), 'euclidean', nn_control[['metric']])

    assertthat::assert_that(nn_control[['metric']] %in% c('euclidean', 'l2', 'cosine', 'ip'),
                            msg=paste0("set_nn_control: nearest neighbor metric for HNSW must be one of 'euclidean', 'l2', 'cosine', or 'ip'"))

    nn_control[['grain_size']] <- ifelse(is.null(nn_control[['grain_size']]), 1, nn_control[['grain_size']])
    nn_control[['cores']] <- ifelse(is.null(nn_control[['cores']]), 1, nn_control[['cores']])
  
    assertthat::assert_that(assertthat::is.count(nn_control[['grain_size']]))
    assertthat::assert_that(assertthat::is.count(nn_control[['cores']]))

    if(bitwAnd(mode, 1)) {
      nn_control[['M']] <- ifelse(is.null(nn_control[['M']]), 48, nn_control[['M']])
      nn_control[['ef_construction']] <- ifelse(is.null(nn_control[['ef_construction']]), 200, nn_control[['ef_construction']])
  
      assertthat::assert_that(assertthat::is.count(nn_control[['M']]))
      assertthat::assert_that(assertthat::is.count(nn_control[['ef_construction']]))

      assertthat::assert_that(nn_control[['M']] >= 2,
                              msg=paste0('set_nn_control: HNSW nearest neighbor M index build parameter must be >= 2'))
    }

    if(bitwAnd(mode, 2)) {
      nn_control[['ef']] <- ifelse(is.null(nn_control[['ef']]), 30, nn_control[['ef']])
      assertthat::assert_that(assertthat::is.count(nn_control[['ef']]))
  
      assertthat::assert_that(nn_control[['ef']] >= k,
                              msg=paste0('set_nn_control: HNSW nearest neighbor ef index search parameter must be >= k (',k,')'))
    }
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


# Report nn_control list values. This is used primarily for
# verbose output.
report_nn_control <- function(label=NULL, nn_control) {
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


#' @title Make a nearest neighbor index.
#'
#' @description Make a nearest neighbor index from the specified
#' reduction_method matrix using either the default nearest neighbor
#' method or the method specified in the nn_control list parameter
#' The function returns the index.
#'
#' @param subject_matrix the reduced dimension matrix used to build
#'  the index.
#' @param nn_control a list of parameters used to make the nearest
#'  neighbor index. See the set_nn_control help for details.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return a nearest neighbor index.
#' @export
make_nn_index <- function(subject_matrix, nn_control=list(), verbose=FALSE) {

  # We check the index build parameters.
  nn_control <- set_nn_control(mode=1, nn_control=nn_control, k=1, nn_control_default=list(), verbose=verbose)

  if(verbose) {
    message('make_nn_index:')
    report_nn_control('  nn_control: ', nn_control)
    tic('make_nn_index: build time')
  }

  nn_method <- nn_control[['method']]

  if(nn_method == 'annoy') {
    # set verbose to FALSE because when TRUE and nr is small, uwot::annoy_build can fail in
    #       if (verbose) {
    #         nstars <- 50
    #         progress_for(
    #           nr, nstars,
    #           function(chunk_start, chunk_end) {
    #             for (i in chunk_start:chunk_end) {
    #               ann$addItem(i - 1, X[i, , drop = FALSE])
    #             }
    #           }
    #         )
    #       }
    #       else {
    nn_index <- uwot:::annoy_build(X=subject_matrix,
                                   metric=nn_control[['metric']],
                                   n_trees=nn_control[['n_trees']],
                                   verbose=FALSE)
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


#' @title Set a nearest neighor index in the cell_data_set.
#'
#' @description Store the nearest nn_index neighbor index
#' in the cds. The reduction_method parameter tells
#' set_cds_nn_index where in the cds to store the index.
#'
#' @param cds a cell_data_set in which to store the nearest neighbor
#'   index.
#' @param reduction_method a string giving the reduced dimension matrix used
#'   to make the nn_index nearest neighbor index.
#' @param nn_index a nearest neighbor index to store in cds.
#' @param nn_control a list of parameters used to make the nearest
#'  neighbor index. See the set_nn_control help for details.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return a cell_data_set with the stored index.
#'
#' @export
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

  # We check the index build parameters.
  nn_control <- set_nn_control(mode=1, nn_control=nn_control, k=1, nn_control_default=list(), verbose=verbose)

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


#' @title Make and store a nearest neighbor index in the cell_data_set.
#'
#' @description Make a nearest neighbor index from the specified
#' reduction_method matrix in the cds using either the default
#' nearest neighbor method or the method specified in the nn_control
#' list parameter. This function returns a cell_data_set.
#'
#' @param cds a cell_data_set with the reduced dimension matrix from
#'   which to make the index and in which the nearest neighbor
#'   index is stored.
#' @param reduction_method a string giving the reduced dimension matrix
#'  to use for making the nn_index nearest neighbor index.
#' @param nn_control a list of parameters to use for making the nearest
#'  neighbor index. See the set_nn_control help for details.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return a cell_data_set with the stored index.
#'
#' @export
make_cds_nn_index <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), nn_control=list(), verbose=FALSE) {
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

  nn_control <- set_nn_control(mode=1, nn_control=nn_control, k=1, nn_control_default=list(), verbose=verbose)

  reduced_matrix <- reducedDims(cds)[[reduction_method]]
  nn_index <- make_nn_index(subject_matrix=reduced_matrix, nn_control=nn_control, verbose=verbose)
  cds <- set_cds_nn_index(cds=cds, reduction_method=reduction_method, nn_index=nn_index, nn_control=nn_control, verbose=verbose)

  return(cds)
}


# Retun the nn_index that was made from the reduction_method
# reduced dimension matrix.
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

  # We need the build method.
  nn_control <- set_nn_control(mode=1, nn_control=nn_control, k=1, nn_control_default=list(), verbose=verbose)

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


#' @title Search nearest neighor index.
#'
#' @description Search nn_index nearest neighbor index for cells near
#' those in the query_matrix.
#'
#' @param query_matrix a reduced dimension matrix used to find the
#'   nearest neighbors in the index created from the
#'   reduction_method matrix.
#' @param nn_index a nearest_neighbor index.
#' @param k an integer for the number of nearest neighbors to return for
#'   each cell. Default is 25.
#' @param nn_control a list of parameters used to make the nearest
#'  neighbor index. See the set_nn_control help for details.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return a list list(nn.idx, nn.dists) where nn.idx is
#'  a matrix of nearest neighbor indices and nn.dists is a matrix of
#'  the distance between the index given by the row number and the index
#'  given in nn.idx. If the same reduced dim matrix is used to make the
#'  index and search the index, the index given by the row number should
#'  be in the row, usually in the first column.
#'
#' @export
search_nn_index <- function(query_matrix, nn_index, k=25, nn_control=list(), verbose=FALSE) {

  nn_control <- set_nn_control(mode=2, nn_control=nn_control, k=k, nn_control_default=list(), verbose=verbose)

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
# for verifying the query matrix -- be certain that the search (and index build) were run
#' on this matrix so call this function immediately after running the search.
set_cds_nn_search <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), search_id, k, nn_control=list(), verbose=TRUE) {
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

  nn_control <- set_nn_control(mode=2, nn_control=nn_control, k=1, nn_control_default=list(), verbose=verbose)

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


# Get nearest neighbor search information stored in the cds.
get_cds_nn_search <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), search_id, verbose=TRUE) {

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


# Get index build and search information that's stored in the cds.
get_cds_nn_control <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), search_id, verbose=TRUE) {

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


#' @title Search a subject matrix for nearest neighbors to a query_matrix.
#'
#' @description Make a nearest neighbors index using the subject matrix and
#' search it for nearest neighbors to the query_matrix.
#'
#' @param subject_matrix a reduced dimension matrix used to build a nearest neighbor index.
#' @param query_matrix a reduced dimension matrix used to search the subject_matrix
#'   nearest neighbor index.
#' @param k an integer for the number of nearest neighbors to return for
#'   each cell. Default is 25.
#' @param nn_control a list of parameters used to make the nearest
#'  neighbor index. See the set_nn_control help for details.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return a list list(nn.idx, nn.dists) where nn.idx is
#'  a matrix of nearest neighbor indices and nn.dists is a matrix of
#'  the distance between the index given by the row number and the index
#'  given in nn.idx. If the query_matrix is the same as the subject
#'  matrix, the index given by the row number should be in the row, usually
#' in the first column.
#'
#' @export
search_nn_matrix <- function(subject_matrix, query_matrix, k=25, nn_control=list(), verbose=FALSE) {

  assertthat::assert_that(is.matrix(subject_matrix),
    msg=paste0('search_nn_matrix: the subject_matrix object must be of type matrix'))

  assertthat::assert_that(is.matrix(query_matrix),
    msg=paste0('search_nn_matrix: the query_matrix object must be of type matrix'))

  nn_control <- set_nn_control(mode=3, nn_control=nn_control, k=k, nn_control_default=list(), verbose=verbose)

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
    # set verbose to FALSE because when TRUE and nr is small, uwot::annoy_build can fail in
    #       if (verbose) {
    #         nstars <- 50
    #         progress_for(
    #           nr, nstars,
    #           function(chunk_start, chunk_end) {
    #             for (i in chunk_start:chunk_end) {
    #               ann$addItem(i - 1, X[i, , drop = FALSE])
    #             }
    #           }
    #         )
    #       }
    #       else {
    nn_index <- uwot:::annoy_build(X=subject_matrix,
                                   metric=nn_control[['metric']],
                                   n_trees=nn_control[['n_trees']],
                                   verbose=FALSE)

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
      for(i in seq(1, len, 1)) {
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
        for(i in seq(2, length(vidx), 1)) {
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

