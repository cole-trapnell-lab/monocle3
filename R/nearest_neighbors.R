# Functions that support nearest neighbors use.


# Check whether nn index exists and is consistent with matrix and parameters.
# This function is not in use currently and may fall into disrepair.
check_cds_nn_index_is_current <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), nn_control=list(), verbose=FALSE) {

  assertthat::assert_that(is(cds, 'cell_data_set'),
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
  nn_control <- set_nn_control(mode=1, nn_control=nn_control, nn_control_default=list(), cds=cds, reduction_method=reduction_method, k=NULL, verbose=verbose)

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
  assertthat::assert_that(is(cds, 'cell_data_set'),
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
# This function is not in use currently and may fall into disrepair.
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


select_nn_parameter_value <- function(parameter, nn_control, nn_control_default, default_value) {
  if(!is.null(nn_control[[parameter]])) {
    return(nn_control[[parameter]])
  }
  else
  if(!is.null(nn_control_default[[parameter]])) {
    return(nn_control_default[[parameter]])
  }
  return(default_value)
}


#' set annoy search_k parameter
#' precedence for setting search_k
#'   o  nn_control[['search_k']]
#'   o  mode == 2 and the following values exist/set
#'        o  cds parameter
#'        o  reduction_method parameter
#'        o  k parameter
#'        o  annoy nn_index in cds for reduction_method
#'   o  nn_control[['n_trees']] and k parameter or k default
#'   o  nn_control_default[['search_k']]
#'   o  nn_control_default[['n_trees']] and k parameter or k default
#'   o  n_trees default and k parameter or k default
#' @noRd
select_annoy_search_k <- function(mode, nn_control, nn_control_default, cds, reduction_method, k, default_n_trees, default_k) {
  if(!is.null(k)) {
    use_k <- k
    src_k <- 'parameter'
  } 
  else {
    use_k <- default_k
    src_k <- 'default'
  }

  if(!is.null(nn_control[['search_k']])) {
    return(nn_control[['search_k']])
  }
  else
  if(mode == 2 &&
     !is.null(cds) &&
     !is.null(reduction_method) &&
     !is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy']][['n_trees']]) &&
     !is.null(k)) {
    n_trees <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy']][['n_trees']]
    return(2 * n_trees * k)
  }
  else
  if(!is.null(nn_control[['n_trees']])) {
    return(2 * nn_control[['n_trees']] * use_k)
  }
  else
  if(!is.null(nn_control_default[['search_k']])) {
    return(nn_control_default[['search_k']])
  }
  else
  if(!is.null(nn_control_default[['n_trees']])) {
    return(2 * nn_control_default[['n_trees']] * use_k)
  }

  return(2 * default_n_trees * use_k)
}



#' @title Verify and set nearest neighbor parameter list.
#'
#' @description Verifies the listed parameter values
#'   that will be passed to the nearest neighbor function
#'   given by nn_control\[\['method'\]\]. Unspecified
#'   values are set to default values. To see the default
#'   values, call the function with
#'   nn_control=list(show_values=TRUE).
#'
#' @section Optional nn_control parameters:
#' \describe{
#'   \item{method}{The method used to find nearest neighbor points.
#'      The available methods are 'nn2', 'annoy', and 'hnsw'.
#'      Detailed information about each method can be found on the
#'      WWW sites:
#'      https://cran.r-project.org/web/packages/RANN/,
#'      https://cran.r-project.org/web/packages/RcppAnnoy/index.html,
#'      and https://cran.rstudio.com/web/packages/RcppHNSW/index.html.}
#'   \item{metric}{The distance metric used by the nearest neighbor
#'      functions.  Annoy accepts 'euclidean', 'cosine', 'manhattan',
#'      and 'hamming'. HNSW accepts 'euclidean', 'l2', 'cosine', and
#'      'ip'. RANN uses 'euclidean'.}
#'   \item{n_trees}{The annoy index build parameter that affects the build
#'      time and index size. Larger values give more accurate results,
#'      longer build times, and larger indexes.}
#'   \item{search_k}{The annoy index search parameter that affects the
#'      search accuracy and time. Larger values give more accurate
#'      results and longer search times. Default is 2 * n_trees * k.
#'      In order to set search_k, the following conditions are tested
#'      and the first TRUE condition is used: nn_control\[\['search_k'\]\]
#'      exists; mode=2 and cds, reduction_method, and k are not NULL,
#'      and the expected index (given by cds, reduction_method, and
#'      nn_control\[\['method'\]\]) exists; nn_control\[\['n_trees'\]\] exists;
#'      nn_control_default\[\['search_k'\]\] exists;
#'      nn_control_default\[\['n_trees'\]\] exists. If none of those is
#'      TRUE, the fallback default n_trees value is used. If the
#'      set_nn_control k parameter value is not NULL, it is used;
#'      otherwise, the default is used.}
#'
#'   \item{M}{The HNSW index build parameter that affects the
#'      search accuracy and memory requirements. Larger values give
#'      more accurate search results and increase the index memory
#'      use.}
#'   \item{ef_construction}{The HNSW index build parameter that affects
#'      the search accuracy and index build time. Larger values give more 
#'      accurate search results and longer build times. Default is 200.}
#'   \item{ef}{The HNSW index search parameter that affects the search
#'      accuracy and search time. Larger values give more accurate results
#'      and longer search times. ef must be greater than or equal to k.}
#'   \item{grain_size}{The annoy and HNSW parameter that gives the
#'      minimum amount of work to do per thread.}
#'   \item{cores}{The annoy and HNSW parameter that gives the number of
#'      threads to use for the annoy index search and for the HNSW index
#'      build and search.}
#'   \item{show_values}{A logical value used to show the
#'      nearest neighbor parameters to use, and then exit the function.
#'      When show_values=TRUE is the only nn_control value, the
#'      parameters are the defaults for the function. Each function
#'      that calls set_nn_control may have its own nn_control_default
#'      list.}
#' }
#'
#' @param mode the nearest neighbor operation for which the nn_control
#'   list will be used. 1=make index, 2=search index, and 3=both make
#'   and search index. Required parameter.
#' @param nn_control an optional list of parameters passed to
#'   the nearest neighbor function specified by nn_control\[\['method'\]\].
#'   If a value is not given in nn_control, the value in nn_control_default
#'   is used. If neither is given, a fallback default is assigned.
#' @param nn_control_default an optional nn_control list to use when
#'   a parameter is not given in nn_control.
#' @param cds a cell_data_set where a nearest neighbor index may be
#'   stored. This may be used to look up parameters that were used
#'   to make an index that is stored in the cds. For example, the default
#'   annoy search_k parameter depends on the n_trees values used to
#'   make the index. The default is NULL.
#' @param reduction_method a reduction method where a nearest neighbor
#'   index may be stored in the cell_data_set, cds. This may be used
#'   to look up parameters that were used to make an index that is
#'   stored in the cds. For example, the default search_k parameter
#'   depends on the n_trees values used to make the index. The default
#'   is NULL.
#' @param k integer the number of desired nearest neighbor points to
#'   return from a search. k is used to
#'   * set the annoy search_k parameter when search_k is not given in
#'     nn_control or nn_control_default where search_k is 2 * n_trees * k.
#'   * test the hnsw ef parameter, which must be at least as large
#'     as k.
#'
#'   k is ignored for index builds and does not give the number
#'   of nearest neighbors to return for a search.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return an updated nn_control list.
set_nn_control <- function(mode, nn_control=list(), nn_control_default=list(), cds=NULL, reduction_method=NULL, k=NULL, verbose=FALSE) {

  default_method <- 'annoy'
  default_metric <- 'euclidean'
  default_k <- 25
  default_n_trees <- 50
  default_M <- 48
  default_ef_construction <- 200
  default_ef <- 150
  default_grain_size <- 1
  default_cores <- 1

  assertthat::assert_that(methods::is(nn_control, "list"))
  assertthat::assert_that(methods::is(nn_control_default, "list"))

  allowed_control_parameters <- c('method',
                                  'metric',
                                  'n_trees',
                                  'search_k',
                                  'M',
                                  'ef_construction',
                                  'ef',
                                  'grain_size',
                                  'cores',
                                  'show_values')

  assertthat::assert_that(all(names(nn_control) %in% allowed_control_parameters),
                          msg = "set_nn_control: unknown variable in nn_control")
  assertthat::assert_that(all(names(nn_control_default) %in% allowed_control_parameters),
                          msg = "set_nn_control: unknown variable in nn_control_default")

  assertthat::assert_that(assertthat::is.count(mode) && mode >= 1 && mode <= 3,
                          msg = paste0("set_nn_control: invalid mode value. Mode must be an integer with the value 1, 2, or 3."))

  nn_control_out <- list()

  nn_control_out[['method']] <- select_nn_parameter_value('method', nn_control, nn_control_default, default_method)

  if(nn_control_out[['method']] == 'nn2') {
    # Do nothing for nn2. Call report_nn_control when verbose <- TRUE
  } else
  if(nn_control_out[['method']] == 'annoy') {
    nn_control_out[['metric']] <- select_nn_parameter_value('metric', nn_control, nn_control_default, default_metric)
    assertthat::assert_that(nn_control_out[['metric']] %in% c('euclidean', 'cosine', 'manhattan', 'hamming'),
                            msg=paste0("set_nn_control: nearest neighbor metric for annoy must be one of 'euclidean', 'cosine', 'manhattan', or 'hamming'"))
    if(bitwAnd(mode, 1)) {
      nn_control_out[['n_trees']] <- select_nn_parameter_value('n_trees', nn_control, nn_control_default, default_n_trees)
      assertthat::assert_that(assertthat::is.count(nn_control_out[['n_trees']]))
    }

    if(bitwAnd(mode, 2)) {
      nn_control_out[['search_k']] <- select_annoy_search_k(mode, nn_control, nn_control_default, cds, reduction_method, k, default_n_trees, default_k)
      nn_control_out[['grain_size']] <- select_nn_parameter_value('grain_size', nn_control, nn_control_default, default_grain_size)
      nn_control_out[['cores']] <- select_nn_parameter_value('cores', nn_control, nn_control_default, default_cores)
      assertthat::assert_that(assertthat::is.count(nn_control_out[['search_k']]))
      assertthat::assert_that(assertthat::is.count(nn_control_out[['grain_size']]))
      assertthat::assert_that(assertthat::is.count(nn_control_out[['cores']]))
    }

  } else
  if(nn_control_out[['method']] == 'hnsw') {
    nn_control_out[['metric']] <- select_nn_parameter_value('metric', nn_control, nn_control_default, default_metric)
    assertthat::assert_that(nn_control_out[['metric']] %in% c('euclidean', 'l2', 'cosine', 'ip'),
                            msg=paste0("set_nn_control: nearest neighbor metric for HNSW must be one of 'euclidean', 'l2', 'cosine', or 'ip'"))
    nn_control_out[['grain_size']] <- select_nn_parameter_value('grain_size', nn_control, nn_control_default, default_grain_size)
    nn_control_out[['cores']] <- select_nn_parameter_value('cores', nn_control, nn_control_default, default_cores)
    assertthat::assert_that(assertthat::is.count(nn_control_out[['grain_size']]))
    assertthat::assert_that(assertthat::is.count(nn_control_out[['cores']]))

    if(bitwAnd(mode, 1)) {
      nn_control_out[['M']] <- select_nn_parameter_value('M', nn_control, nn_control_default, default_M)
      nn_control_out[['ef_construction']] <- select_nn_parameter_value('ef_construction', nn_control, nn_control_default, default_ef_construction)
      assertthat::assert_that(assertthat::is.count(nn_control_out[['M']]))
      assertthat::assert_that(assertthat::is.count(nn_control_out[['ef_construction']]))
      assertthat::assert_that(nn_control_out[['M']] >= 2,
                              msg=paste0('set_nn_control: HNSW nearest neighbor M index build parameter must be >= 2'))
    }

    if(bitwAnd(mode, 2)) {
      assertthat::assert_that(assertthat::is.count(k),
                              msg=paste0('set_nn_control: parameter k must be set for method=\'hnsw\' and modes 2 and 3.'))
      nn_control_out[['ef']] <- select_nn_parameter_value('ef', nn_control, nn_control_default, default_ef)
      assertthat::assert_that(assertthat::is.count(nn_control_out[['ef']]))
      assertthat::assert_that(nn_control_out[['ef']] >= k,
                              msg=paste0('set_nn_control: HNSW nearest neighbor ef index search parameter must be >= k (',k,')'))
    }
  }
  else
    stop('set_nn_control: unsupported nearest neighbor method \'', nn_control_out[['method']], '\'')

  if(verbose) {
    cs <- get_call_stack_as_string()
    message('set_nn_control: call stack: ', cs)
    report_nn_control('  nn_control: ', nn_control_out)
  }

  if(!is.null(nn_control[['show_values']]) && nn_control[['show_values']] == TRUE)
  {
    report_nn_control('  nn_control: ', nn_control=nn_control_out)
    stop_no_noise()
  }

  return(nn_control_out)
}


# Report nn_control list values.
report_nn_control <- function(label=NULL, nn_control) {
  indent <- ''
  if(!is.null(label)) {
    indent <- '  '
  }

  message(ifelse(!is.null(label), label, ''))

  message(indent, '  method: ', ifelse(!is.null(nn_control[['method']]), nn_control[['method']], as.character(NA)))
  message(indent, '  metric: ', ifelse(!is.null(nn_control[['metric']]), nn_control[['metric']], as.character(NA)))

  if(is.null(nn_control[['method']])) {
    return()
  }

  if(nn_control[['method']] == 'nn2') {
    message(indent, '  nn2 has no parameters')
  }
  else
  if(nn_control[['method']] == 'annoy') {
    message(indent, '  n_trees: ', ifelse(!is.null(nn_control[['n_trees']]), nn_control[['n_trees']], as.character(NA)))
    message(indent, '  search_k: ', ifelse(!is.null(nn_control[['search_k']]), nn_control[['search_k']], as.character(NA)))
    message(indent, '  cores: ', ifelse(!is.null(nn_control[['cores']]), nn_control[['cores']], as.character(NA)))
    message(indent, '  grain_size: ', ifelse(!is.null(nn_control[['grain_size']]), nn_control[['grain_size']], as.character(NA)))
  }
  else
  if(nn_control[['method']] == 'hnsw') {
    message(indent, '  M: ', ifelse(!is.null(nn_control[['M']]), nn_control[['M']], as.character(NA)))
    message(indent, '  ef_construction: ', ifelse(!is.null(nn_control[['ef_construction']]), nn_control[['ef_construction']], as.character(NA)))
    message(indent, '  ef: ', ifelse(!is.null(nn_control[['ef']]), nn_control[['ef']], as.character(NA)))
    message(indent, '  cores: ', ifelse(!is.null(nn_control[['cores']]), nn_control[['cores']], as.character(NA)))
    message(indent, '  grain_size: ', ifelse(!is.null(nn_control[['grain_size']]), nn_control[['grain_size']], as.character(NA)))
  }
  else
    stop('report_nn_control: unsupported nearest neighbor method \'', nn_method, '\'')
}


new_annoy_index <- function(metric, ndim) {
  nn_class <- switch( metric,
                      cosine = RcppAnnoy::AnnoyAngular,
                      euclidean = RcppAnnoy::AnnoyEuclidean,
                      hamming = RcppAnnoy::AnnoyHamming,
                      manhattan = RcppAnnoy::AnnoyManhattan,
                      stop(paste0('unsupported annoy metric ', metric))
                    )
  nn_index <- new(nn_class, ndim)
  return(nn_index)
}


#' @title Make a nearest neighbor index.
#'
#' @description Make a nearest neighbor index from the
#'   subject_matrix using either the default nearest
#'   neighbor method or the method specified in the
#'   nn_control list parameter. The function returns
#'   the index.
#'
#' @param subject_matrix the matrix used to build the
#'   index.
#' @param nn_control a list of parameters used to make the
#'   nearest neighbor index. See the set_nn_control help
#'   for details.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return a nearest neighbor index.
#' @export
make_nn_index <- function(subject_matrix, nn_control=list(), verbose=FALSE) {
  assertthat::assert_that(is(subject_matrix, 'matrix') ||
                          is_sparse_matrix(subject_matrix),
    msg=paste0('make_nn_matrix: the subject_matrix object must be of type matrix'))

  # We check the index build parameters.
  nn_control <- set_nn_control(mode=1, nn_control=nn_control, nn_control_default=list(), cds=NULL, reduction_method=NULL, k=NULL, verbose=verbose)

  if(verbose) {
    message('make_nn_index:')
    report_nn_control('  nn_control: ', nn_control)
    tictoc::tic('make_nn_index: build time')
  }

  nn_method <- nn_control[['method']]

  if(nn_method == 'nn2') {
    stop('make_nn_index is not valid for method nn2')
  } else
  if(nn_method == 'annoy') {
    num_row <- nrow(subject_matrix)
    num_col <- ncol(subject_matrix)
    nn_index <- new_annoy_index(nn_control[['metric']], num_col)
    n_trees <- nn_control[['n_trees']]
    cores <- nn_control[['cores']]
    if(num_row > 0 ) {
      for(i in 1:num_row)
        nn_index$addItem(i-1, subject_matrix[i,])
      nn_index$build(n_trees)
    }
  }
  else
  if(nn_method == 'hnsw') {
    nn_index <- RcppHNSW::hnsw_build(X=subject_matrix,
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
    tictoc::toc()
  }

  return(nn_index)
}


#' @title Set a nearest neighor index in the cell_data_set.
#'
#' @description Store the given nearest neighbor index
#'   in the cell_data_set. The reduction_method parameter
#'   tells set_cds_nn_index where in the cell_data_set to
#'   store the index.
#'
#' @param cds a cell_data_set in which to store the nearest
#'   neighbor index.
#' @param reduction_method a string giving the reduced
#'   dimension matrix used to make the nn_index nearest
#'   neighbor index, and determines where the index is
#'   stored in the cell_data_set.
#' @param nn_index a nearest neighbor index to store in cds.
#' @param nn_control a list of parameters used to make the nearest
#'  neighbor index. This is used to identify where the index
#'  is stored, based on the nearest neighbor method (annoy or hnsw)
#'  and reduction_method. See the set_nn_control help for details.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return a cell_data_set with the stored index.
#'
#' @export
set_cds_nn_index <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), nn_index, nn_control=list(), verbose=FALSE) {
  assertthat::assert_that(is(cds, 'cell_data_set'),
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
  nn_control <- set_nn_control(mode=1, nn_control=nn_control, nn_control_default=list(), cds=cds, reduction_method=reduction_method, k=NULL, verbose=verbose)

  nn_method <- nn_control[['method']]

  if(nn_method == 'nn2') {
    stop('set_cds_nn_index is not valid for method nn2')
  } else
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
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']] <- nn_index
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
#'   reduction_method matrix in the cell_data_set using either the
#'   default nearest neighbor method or the method specified in the
#'   nn_control list parameter, and store the index in the
#'   cell_data_set. This function returns a cell_data_set.
#'
#' @param cds a cell_data_set with the reduced dimension matrix from
#'   which to make the nearest neighbor index and with which the index
#'   is stored.
#' @param reduction_method a string giving the reduced dimension matrix
#'   to use for making the nn_index nearest neighbor index.
#' @param nn_control a list of parameters to use for making the nearest
#'   neighbor index. See the set_nn_control help for details.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return a cell_data_set with the stored index.
#'
#' @export
make_cds_nn_index <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), nn_control=list(), verbose=FALSE) {
  assertthat::assert_that(is(cds, 'cell_data_set'),
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

  nn_control <- set_nn_control(mode=1, nn_control=nn_control, nn_control_default=list(), cds=cds, reduction_method=reduction_method, k=NULL, verbose=verbose)

  nn_method <- nn_control[['method']]

  if(nn_method == 'nn2') {
    stop('make_cds_nn_index is not valid for method nn2')
  }

  reduced_matrix <- reducedDims(cds)[[reduction_method]]
  nn_index <- make_nn_index(subject_matrix=reduced_matrix, nn_control=nn_control, verbose=verbose)
  cds <- set_cds_nn_index(cds=cds, reduction_method=reduction_method, nn_index=nn_index, nn_control=nn_control, verbose=verbose)

  return(cds)
}


# Return the nn_index that was made from the reduction_method
# reduced dimension matrix.
get_cds_nn_index <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), nn_control=list(), verbose=FALSE) {

  assertthat::assert_that(is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  # We need the build method.
  nn_control <- set_nn_control(mode=1, nn_control=nn_control, nn_control_default=list(), cds=cds, reduction_method=reduction_method, k=NULL, verbose=verbose)

  nn_method <- nn_control[['method']]

  assertthat::assert_that(
    !is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]]),
    msg = paste("There is no nearest neighbor index",
                "for reduction_method",
                reduction_method,
                "and nearest neighbor method",
                nn_method))

  assertthat::assert_that(nn_method %in% c('annoy', 'hnsw'),
    msg = paste0('get_cds_nn_index: unsupported nearest neighbor index type \'',
                nn_method, '\'.'))

  if(nn_method == 'annoy') {
    nn_index <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']]
  }
  else
  if(nn_method == 'hnsw') {
    nn_index <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']]
  }
  else
    stop('get_cds_nn_index: unsupported nearest neighbor index type \'', nn_method, '\'')

  return(nn_index)
}


#' @title Search a nearest neighbor index.
#'
#' @description Search a nearest neighbor index for cells near
#' those in the query_matrix.
#'
#' @param query_matrix a reduced dimension matrix used to find the
#'   nearest neighbors in the index nn_index.
#' @param nn_index a nearest_neighbor index.
#' @param k an integer for the number of nearest neighbors to return for
#'   each cell. Default is 25.
#' @param nn_control a list of parameters used to search the nearest
#'  neighbor index. See the set_nn_control help for details. Note: the
#'  default annoy search_k parameter value is set to the default value
#'  of 2 * n_trees * k. It does not know the value of n_trees that was
#'  used to build the annoy index so if a non-default n_trees value
#'  was used to build the index, you may need to set search_k in
#'  nn_control list when you run search_nn_index.
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
  assertthat::assert_that(is(query_matrix, 'matrix') ||
                          is_sparse_matrix(query_matrix),
    msg=paste0('make_nn_matrix: the query_matrix object must be of type matrix'))

  nn_control <- set_nn_control(mode=2, nn_control=nn_control, nn_control_default=list(), cds=NULL, reduction_method=NULL, k=k, verbose=verbose)

  nn_method <- nn_control[['method']]

  assertthat::assert_that(nn_method %in% c('annoy', 'hnsw'),
    msg = paste0('search_nn_index: unsupported nearest neighbor index type \'',
                nn_method, '\'.'))

  if(verbose) {
    message('search_nn_index:')
    message('  k: ', k)
    report_nn_control('  nn_control: ', nn_control)
    tictoc::tic('search_nn_index: search time:')
  }

  k <- min(k, nrow(query_matrix))

  if(nn_method == 'nn2') {
    stop('search_nn_index is not valid for method nn2')
  } else
  if(nn_method == 'annoy') {
    # notes:
    #   o  set list names to nn.idx and nn.dists for compatibility with
    #      RANN::nn2()
    #
    num_row <- nrow(query_matrix)
    idx <- matrix(nrow=num_row, ncol=k)
    dist <- matrix(nrow=num_row, ncol=k)
    search_k <- nn_control[['search_k']]
    num_bad <- 0
    for(i in 1:num_row) {
      nn_list <- nn_index$getNNsByVectorList(query_matrix[i,], k, search_k, TRUE)
      if(length(nn_list$item) != k)
        num_bad <- num_bad + 1
      idx[i,] <- nn_list$item
      dist[i,] <- nn_list$distance
    }
    if(num_bad)
      stop('annoy was unable to find ', k, ' nearest neighbors for ', num_bad, ' rows. You may need to increase the n_trees and/or search_k parameter values.')
    if(nn_control[['metric']] == 'cosine') {
      dist <- 0.5 * dist * dist
    }
    nn_res <- list(nn.idx=idx+1, nn.dists=dist)
  }
  else
  if(nn_method == 'hnsw') {
    assertthat::assert_that(nn_control[['ef']] >= k,
      msg=paste0('search_nn_index: ef must be >= k'))

    tmp <- RcppHNSW::hnsw_search(X=query_matrix,
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
    tictoc::toc()
  }

  return(nn_res)
}


# Notice that this function uses the available reducedDims(cds)[[reduction_method]] matrix
# for verifying the query matrix -- be certain that the search (and index build) were run
# on this matrix so call this function immediately after running the search.
# Notes:
#   o  this function was written to store nn search information in the cds when
#      we intended to re-use an index for later processing stages but
#      it appears that there is little or no opportunity to re-use
#      indices so we are not storing search information at this time.
# This function is not in use currently and may fall into disrepair.
set_cds_nn_search <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), search_id, k, nn_control=list(), verbose=TRUE) {
  assertthat::assert_that(is(cds, 'cell_data_set'),
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

  nn_control <- set_nn_control(mode=2, nn_control=nn_control, nn_control_default=list(), cds=cds, reduction_method=reduction_method, k=k, verbose=verbose)

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
# This function is not in use currently and may fall into disrepair.
get_cds_nn_search <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), search_id, verbose=TRUE) {

  assertthat::assert_that(is(cds, 'cell_data_set'),
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
# This function is not used at this time. It may fall into disrepair.
get_cds_nn_control <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), search_id, verbose=TRUE) {

  assertthat::assert_that(is(cds, 'cell_data_set'),
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
#' @param subject_matrix a matrix used to build a nearest neighbor index.
#' @param query_matrix a matrix used to search the subject_matrix
#'   nearest neighbor index.
#' @param k an integer for the number of nearest neighbors to return for
#'   each cell. Default is 25.
#' @param nn_control a list of parameters used to make and search the nearest
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
  assertthat::assert_that(is(subject_matrix, 'matrix') ||
                          is_sparse_matrix(subject_matrix),
    msg=paste0('search_nn_matrix: the subject_matrix object must be of type matrix'))

  assertthat::assert_that(is(query_matrix, 'matrix') ||
                          is_sparse_matrix(query_matrix),
    msg=paste0('search_nn_matrix: the query_matrix object must be of type matrix'))

  nn_control <- set_nn_control(mode=3, nn_control=nn_control, nn_control_default=list(), cds=NULL, reduction_method=NULL, k=k, verbose=verbose)

  method <- nn_control[['method']]
  k <- min(k, nrow(subject_matrix))

  if(verbose) {
    message('search_nn_matrix:')
    message('  k: ', k)
    report_nn_control('  nn_control: ', nn_control)
    tictoc::tic('search_nn_matrix: search_time')
  }


  if(method == 'nn2') {
    nn_res <- RANN::nn2(subject_matrix, query_matrix, min(k, nrow(subject_matrix)), searchtype = "standard")
  } else {
    nn_index <- make_nn_index(subject_matrix, nn_control=nn_control, verbose=verbose)
    nn_res <- search_nn_index(query_matrix=query_matrix, nn_index=nn_index, k=k, nn_control=nn_control, verbose=verbose)
  }

  if(verbose)
    tictoc::toc()

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
  len <- length(idx[1,])
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
# Notes:
#   o  we assume that the self index distance is zero
#   o  do not use this if the search point set differs
#      from the index build set.
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
          dzero <- TRUE
          for(i in seq(1, length(vidx), 1)) {
            if(vdst[[i]] > .Machine$double.xmin) {
              dzero <- FALSE
            }
          }
          vidx <- c(irow, vidx[1:(length(vidx)-1)])
          vdst <- c(0, vdst[1:(length(vdst)-1)])
          dst[irow,] <- vdst
          if(diagnostics) {
            message('swap_nn_row_index_point: dst row post fix: ', paste(vdst, collapse=' '))
          }
          if(!dzero)
            message(paste('Warning: at least one row of the nearest neighbor search result is missing\n',
                          'the row number (self). You may need to make the index build and/or search\n',
                          'more sensitive.'))
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

