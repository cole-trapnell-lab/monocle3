# Functions that support nearest neighbors use.


# Check if a cds has an nn_index.
# Indices checked:
#   pca_search_annoy
#   pca_search_hnsw
#   aligned_search_annoy
#   aligned_search_hnsw
#   umap_search_annoy
#   umap_search_hnsw
#   umap_model_hnsw
has_nn_index <- function(cds, nn_index_type) {
  if(nn_index_type == 'pca_search_annoy') {
    # Monocle3 PCA search using annoy.
    nn_index <- cds@reduce_dim_aux[['PCA']][['nn_index']][['annoy']][['nn_index']]
    res <- test_annoy_index(nn_index=nn_index, verbose=FALSE)
  }
  else
  if(nn_index_type == 'pca_search_hnsw') {
    # Monocle3 PCA search using hnsw.
    nn_index <- cds@reduce_dim_aux[['PCA']][['nn_index']][['hnsw']][['nn_index']]
    res <- test_hnsw_index(nn_index=nn_index, verbose=FALSE)
  }
  else
  if(nn_index_type == 'aligned_search_annoy') {
    # Monocle3 Aligned search using annoy.
    nn_index <- cds@reduce_dim_aux[['Aligned']][['nn_index']][['annoy']][['nn_index']]
    res <- test_annoy_index(nn_index=nn_index, verbose=FALSE)
  }
  else
  if(nn_index_type == 'aligned_search_hnsw') {
    # Monocle3 Aligned search using hnsw.
    nn_index <- cds@reduce_dim_aux[['Aligned']][['nn_index']][['hnsw']][['nn_index']]
    res <- test_hnsw_index(nn_index=nn_index, verbose=FALSE)
  }
  else
  if(nn_index_type == 'umap_search_annoy') {
    # Monocle3 UMAP search using annoy.
    nn_index <- cds@reduce_dim_aux[['UMAP']][['nn_index']][['annoy']][['nn_index']]
    res <- test_annoy_index(nn_index=nn_index, verbose=FALSE)
  }
  else
  if(nn_index_type == 'umap_search_hnsw') {
    # Monocle3 UMAP search using hnsw.
    nn_index <- cds@reduce_dim_aux[['UMAP']][['nn_index']][['hnsw']][['nn_index']]
    res <- test_hnsw_index(nn_index=nn_index, verbose=FALSE)
  }
  else
  if(nn_index_type == 'umap_model_annoy') {
    # UWOT UMAP model using annoy.
    nn_index <- cds@reduce_dim_aux[['UMAP']][['model']][['umap_model']][['nn_index']]
    res <- test_annoy_index(nn_index=nn_index, verbose=FALSE)
  }
  else {
    stop('has_nn_index: unrecognized nn_index_type: \'', nn_index_type, '\'')
  }
  return(res)
}


# Check whether nn index exists and is consistent with matrix and parameters.
# This function is not in use currently and may fall into disrepair.
check_cds_nn_index_is_current <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), nn_control=list(), verbose=FALSE) {

  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(
    !is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
    msg = paste0('Data has not been processed with',
                 ' chosen reduction_method: ',
                 reduction_method))

  # We check the index build parameters.
  nn_control <- set_nn_control(mode=1, nn_control=nn_control, nn_control_default=list(), nn_index=NULL, k=NULL, verbose=verbose)

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

  if(nrow(SingleCellExperiment::reducedDims(cds)[[reduction_method]]) != cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['nrow']]) {
    if(verbose) {
      message('check_cds_nn_index_is_current: FALSE')
    }
    return(FALSE)
  }

  if(ncol(SingleCellExperiment::reducedDims(cds)[[reduction_method]]) != cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['ncol']]) {
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
    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['metric']] != nn_control[['metric']]) {
      if(verbose) {
        message('check_cds_nn_index_is_current: FALSE')
      }
      return(FALSE)
    }

    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['n_trees']] != nn_control[['n_trees']]) {
        if(verbose) {
          message('check_cds_nn_index_is_current: FALSE')
        }
       return(FALSE)
    }
  }
  else
  if(nn_method == 'hnsw') {
    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['metric']] != nn_control[['metric']]) {
      if(verbose) {
        message('check_cds_nn_index_is_current: FALSE')
      }
      return(FALSE)
    }

    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['M']] != nn_control[['M']]) { 
      if(verbose) {
        message('check_cds_nn_index_is_current: FALSE')
      }
      return(FALSE)
    }

    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['ef_construction']] != nn_control[['ef_construction']]) { 
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
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
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
      cds@reduce_dim_aux[[reduction_method]][['nn_index']] <- S4Vectors::SimpleList()
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

  if(is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]])) {
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
#'        o  nn_index parameter
#'        o  n_trees value in nn_index parameter object
#'        o  k parameter
#'   o  nn_control[['n_trees']] and k parameter or k default
#'   o  nn_control_default[['search_k']]
#'   o  nn_control_default[['n_trees']] and k parameter or k default
#'   o  n_trees default and k parameter or k default
#' @noRd
select_annoy_search_k <- function(mode, nn_control, nn_control_default, nn_index, k, default_n_trees, default_k) {
  if(!is.null(k)) {
    use_k <- k
    src_k <- 'parameter'
  } 
  else {
    use_k <- default_k
    src_k <- 'default'
  }

  # Sanity test.
  if(mode == 2 &&
     !is.null(nn_index) &&
     is.null(nn_index[['n_trees']])) {
    stop('set_nn_control: unexpected condition: found old version of reduce_dim_aux')
  }

  if(!is.null(nn_control[['search_k']])) {
    return(nn_control[['search_k']])
  }
  else
  if(mode == 2 &&
     !is.null(nn_index) &&
     !is.null(nn_index[['n_trees']]) &&
     !is.null(k)) {
    n_trees <- nn_index[['n_trees']]
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
#'      exists; mode=2 and nn_index and k are not NULL;
#'      nn_control\[\['n_trees'\]\] exists;
#'      nn_control_default\[\['search_k'\]\] exists;
#'      nn_control_default\[\['n_trees'\]\] exists. If none of those is
#'      TRUE, the fallback default n_trees value is used. If the
#'      set_nn_control k parameter value is not NULL, it is used;
#'      otherwise, the default is used.}
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
#'   \item{grain_size}{The HNSW parameter that gives the
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
#' @param nn_index an nn_index. This may be used to look up parameters
#'   that were used to make an index. For example, the default search_k
#'   parameter depends on the n_trees values used to make the index.
#'   The default is NULL.
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
set_nn_control <- function(mode, nn_control=list(), nn_control_default=list(), nn_index=NULL, k=NULL, verbose=FALSE) {

  default_method <- 'annoy'
  default_metric <- 'euclidean'
  default_k <- 25
  default_n_trees <- 50
  default_M <- 48
  default_ef_construction <- 200
  default_ef <- 150
  default_grain_size <- 1
  default_cores <- 1
  default_annoy_random_seed <- 42

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
                                  'show_values',
                                  'annoy_random_seed')

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
      nn_control_out[['annoy_random_seed']] <- select_nn_parameter_value('annoy_random_seed', nn_control, nn_control_default, default_annoy_random_seed)
      assertthat::assert_that(assertthat::is.count(nn_control_out[['n_trees']]))
    }

    if(bitwAnd(mode, 2)) {
      nn_control_out[['search_k']] <- select_annoy_search_k(mode, nn_control, nn_control_default, nn_index, k, default_n_trees, default_k)
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
    stop('report_nn_control: unsupported nearest neighbor method \'', nn_control[['method']], '\'')
}


new_annoy_index <- function(metric, ndim) {
  nn_class <- switch( metric,
                      cosine = RcppAnnoy::AnnoyAngular,
                      euclidean = RcppAnnoy::AnnoyEuclidean,
                      hamming = RcppAnnoy::AnnoyHamming,
                      manhattan = RcppAnnoy::AnnoyManhattan,
                      stop('unsupported annoy metric ', metric)
                    )
  nn_index <- methods::new(nn_class, ndim)
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
#'
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     cds <- preprocess_cds(cds)
#'     nn_index <- make_nn_index(SingleCellExperiment::reducedDims(cds)[['PCA']])
#'   }
#'
#' @importFrom utils packageVersion
#' @export
make_nn_index <- function(subject_matrix, nn_control=list(), verbose=FALSE) {
  assertthat::assert_that(methods::is(subject_matrix, 'matrix') ||
                          is_sparse_matrix(subject_matrix),
    msg=paste0('make_nn_matrix: the subject_matrix object must be of type matrix'))

  # We check the index build parameters.
  nn_control_default <- get_global_variable('nn_control_annoy_euclidean')
  nn_control <- set_nn_control(mode=1, nn_control=nn_control, nn_control_default=nn_control_default, nn_index=NULL, k=NULL, verbose=verbose)

  if(verbose) {
    message('make_nn_index:')
    report_nn_control('  nn_control: ', nn_control)
    tick('make_nn_index: build time')
  }

  nn_method <- nn_control[['method']]
  metric = nn_control[['metric']]

  num_row <- nrow(subject_matrix)
  num_col <- ncol(subject_matrix)
  if(!is.null(rownames(subject_matrix)))
    checksum_rownames <- digest::digest(sort(rownames(subject_matrix)))
  else
    checksum_rownames <- NA_character_

  if(nn_method == 'nn2') {
    stop('make_nn_index is not valid for method nn2')
  } else
  if(nn_method == 'annoy') {
    monocle3_annoy_index_version <- get_global_variable('monocle3_annoy_index_version')
    annoy_index <- new_annoy_index(metric, num_col)
    annoy_random_seed <- nn_control[['annoy_random_seed']]
    annoy_index$setSeed(annoy_random_seed)
    n_trees <- nn_control[['n_trees']]
    if(num_row > 0 ) {
      for(i in 1:num_row)
        annoy_index$addItem(i-1, subject_matrix[i,])
      annoy_index$build(n_trees)
    }
    annoy_index_version <- packageVersion('RcppAnnoy')
    nn_index <- list(method='annoy', annoy_index=annoy_index, version=monocle3_annoy_index_version, annoy_index_version=annoy_index_version, metric=metric, n_trees=n_trees, nrow=num_row, ncol=num_col, checksum_rownames=checksum_rownames, annoy_random_seed=annoy_random_seed)
  }
  else
  if(nn_method == 'hnsw') {
    monocle3_hnsw_index_version <- get_global_variable('monocle3_hnsw_index_version')
    M <- nn_control[['M']]
    ef_construction <- nn_control[['ef_construction']]
    hnsw_index <- RcppHNSW::hnsw_build(X=subject_matrix,
                                       distance=metric,
                                       M=M,
                                       ef=ef_construction,
                                       verbose=verbose,
                                       n_threads=nn_control[['cores']],
                                       grain_size=nn_control[['grain_size']])
    hnsw_index_version <- packageVersion('RcppHNSW')
    nn_index <- list(method='hnsw', hnsw_index=hnsw_index, version=monocle3_hnsw_index_version, hnsw_index_version=hnsw_index_version, metric=metric, M=M, ef_construction=ef_construction, nrow=num_row, ncol=num_col, checksum_rownames=checksum_rownames)
  }
  else
    stop('make_nn_index: unsupported nearest neighbor index type \'', nn_method, '\'')

  if(verbose) {
    tock()
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
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return a cell_data_set with the stored index.
#'
#' @export
set_cds_nn_index <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), nn_index, verbose=FALSE) {
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
                          msg = paste0("When reduction_method = '", reduction_method,
                                      "' the cds must have been processed for it.",
                                      " Please run the required processing function",
                                      " before this one."))

  nn_method <- nn_index[['method']]

  if(nn_method == 'nn2') {
    stop('set_cds_nn_index is not valid for method nn2')
  } else
  if(nn_method == 'annoy') {
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]] <- S4Vectors::SimpleList()
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']] <- nn_index
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['matrix_id']] <- get_reduce_dim_matrix_identity(cds, reduction_method)[['matrix_id']]
  }
  else
  if(nn_method == 'hnsw') {
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]] <- S4Vectors::SimpleList()
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']] <- nn_index
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
#'   Note: distances in tSNE space reflect spatial differences poorly
#'   so using nearest neighbors with it may be meaningless.
#' @param nn_control a list of parameters to use for making the nearest
#'   neighbor index. See the set_nn_control help for details.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return a cell_data_set with the stored index.
#'
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     cds <- preprocess_cds(cds)
#'     cds <- make_cds_nn_index(cds, 'PCA')
#'   }
#'
#' @export
make_cds_nn_index <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), nn_control=list(), verbose=FALSE) {
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
                          msg = paste0("When reduction_method = '", reduction_method,
                                      "' the cds must have been processed for it.",
                                      " Please run the required processing function",
                                      " before this one."))

  nn_control_default <- get_global_variable('nn_control_annoy_euclidean')
  nn_control <- set_nn_control(mode=1, nn_control=nn_control, nn_control_default=nn_control_default, nn_index=NULL, k=NULL, verbose=verbose)

  nn_method <- nn_control[['method']]

  if(nn_method == 'nn2') {
    stop('make_cds_nn_index is not valid for method nn2')
  }

  reduced_matrix <- SingleCellExperiment::reducedDims(cds)[[reduction_method]]
  nn_index <- make_nn_index(subject_matrix=reduced_matrix, nn_control=nn_control, verbose=verbose)
  cds <- set_cds_nn_index(cds=cds, reduction_method=reduction_method, nn_index=nn_index, verbose=verbose)

  return(cds)
}


# Return the nn_index that was made from the reduction_method
# reduced dimension matrix.
get_cds_nn_index <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), nn_method, verbose=FALSE) {

  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

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


# Check that the annoy index exists.
# Returns logical TRUE if the index exists.
test_annoy_index <- function(nn_index, verbose=FALSE) {
  res <- TRUE
  index_obj <- NULL
  if(!is.null(nn_index[['annoy_index']])) {
    index_obj <- nn_index[['annoy_index']]
  }
  else
  if(!is.null(nn_index[['ann']])) {
    index_obj <- nn_index[['ann']]
  }
  else {
    if(verbose) {
      cs <- get_call_stack_as_string()
      message('test_annoy_index: the annoy nearest neighbor does not exist\ncall stack: ', cs)
    }
    return(FALSE)
  }

  tryCatch( {
    dist_res <- index_obj$getDistance(0,1)
  },
  error=function(emsg) {
    if(verbose) {
      cs <- get_call_stack_as_string()
      message('test_annoy_index: the annoy nearest neighbor does not exist\ncall stack: ', cs)
    }
    res <<- FALSE
  } )

  return(res)
}


# Check that the hnsw index exists.
# Returns logical TRUE if the index exists.
test_hnsw_index <- function(nn_index, verbose=FALSE) {
  res <- TRUE
  if(is.null(nn_index[['hnsw_index']])) {
    if(verbose) {
      cs <- get_call_stack_as_string()
      message('test_hnsw_index: the hnsw nearest neighbor does not exist\ncall stack: ', cs)
    }
    return(FALSE)
  }

  tryCatch( {
    size_res <- nn_index[['hnsw_index']]$size()
  },
  error=function(emsg) {
    if(verbose) {
      cs <- get_call_stack_as_string()
      message('test_hnsw_index: the hnsw nearest neighbor does not exist\ncall stack: ', cs)
    }
    res <<- FALSE
  } )

  return(res)
}


search_nn_annoy_index <- function(query_matrix, nn_index, metric, k, search_k, beg_row_index, end_row_index) {
  assertthat::assert_that(beg_row_index <= end_row_index,
                          msg=paste0('search_nn_annoy_index: beg_row_index must be <= end_row_index'))

  if(!is.null(nn_index[['version']]))
    annoy_index <- nn_index[['annoy_index']]
  else
  if(!is.null(nn_index[['type']]))
    annoy_index <- nn_index[['ann']]
  else
    annoy_index <- nn_index
  nrow <- end_row_index - beg_row_index + 1
  idx <- matrix(nrow=nrow, ncol=k)
  dists <- matrix(nrow=nrow, ncol=k)
  num_bad <- 0
  offset <- beg_row_index - 1
  for(i in seq(1, nrow)) {
    nn_list <- annoy_index$getNNsByVectorList(query_matrix[i+offset,], k, search_k, TRUE)
    if(length(nn_list$item) != k)
      num_bad <- num_bad + 1
    idx[i,] <- nn_list$item
    dists[i,] <- nn_list$distance
  }
  if(metric == 'cosine') {
    dists <- 0.5 * dists * dists
  }
  return(list(idx=idx+1, dists=dists, num_bad=num_bad))
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
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     cds <- preprocess_cds(cds)
#'     nn_index <- make_nn_index(SingleCellExperiment::reducedDims(cds)[['PCA']])
#'     nn_res <- search_nn_index(SingleCellExperiment::reducedDims(cds)[['PCA']], nn_index, 10)
#'   }
#'
#' @export
search_nn_index <- function(query_matrix, nn_index, k=25, nn_control=list(), verbose=FALSE) {
  assertthat::assert_that(methods::is(query_matrix, 'matrix') ||
                          is_sparse_matrix(query_matrix),
    msg=paste0('make_nn_matrix: the query_matrix object must be of type matrix'))

  assertthat::assert_that(assertthat::is.count(k))

  nn_control_default <- get_global_variable('nn_control_annoy_euclidean')
  nn_control <- set_nn_control(mode=2, nn_control=nn_control, nn_control_default=nn_control_default, nn_index=nn_index, k=k, verbose=verbose)

  nn_method <- nn_control[['method']]

  assertthat::assert_that(nn_method %in% c('annoy', 'hnsw'),
    msg = paste0('search_nn_index: unsupported nearest neighbor index type \'',
                nn_method, '\'.'))

  if(verbose) {
    message('search_nn_index:')
    message('  k: ', k)
    report_nn_control('  nn_control: ', nn_control)
    tick('search_nn_index: search time:')
  }

  k <- min(k, nrow(query_matrix))

  cores <- nn_control[['cores']]

  if(nn_method == 'nn2') {
    stop('search_nn_index is not valid for method nn2')
  } else
  if(nn_method == 'annoy') {
    if(!test_annoy_index(nn_index=nn_index, verbose=verbose)) {
      stop('search_nn_index: the annoy nearest neighbor does not exist.')
    }

    # notes:
    #   o  set list names to nn.idx and nn.dists for compatibility with
    #      RANN::nn2()
    #
    num_row <- nrow(query_matrix)
    idx <- matrix(nrow=num_row, ncol=k)
    dists <- matrix(nrow=num_row, ncol=k)
    metric <- nn_control[['metric']]
    search_k <- nn_control[['search_k']]

    if(cores <= 1) {
      nn_res <- search_nn_annoy_index(query_matrix=query_matrix, nn_index=nn_index, metric=metric, k=k, search_k=search_k, beg_row_index=1, end_row_index=num_row)
      if(nn_res[['num_bad']])
        stop('annoy was unable to find ', k, ' nearest neighbors for ', nn_res[["num_bad"]], ' rows. You may need to increase the n_trees and/or search_k parameter values.')
      nn_res <- list(nn.idx=nn_res[['idx']], nn.dists=nn_res[['dists']])
    }
    else {
      omp_num_threads <- get_global_variable('omp_num_threads')
      blas_num_threads <- get_global_variable('blas_num_threads')

      RhpcBLASctl::omp_set_num_threads(1L)
      RhpcBLASctl::blas_set_num_threads(1L)

      if(cores > num_row)
        cores <- num_row
      tasks <- tasks_per_block(num_row, cores)
      beg_block <- c(0,cumsum(tasks))[1:cores]+1
      end_block <- cumsum(tasks)
      nn_blocks <- list()
      inplan <- future::plan()
      future::plan(future::multicore, workers=cores)
      on.exit(future::plan(inplan), add=TRUE)
      for(iblock in seq(cores)) {
        nn_blocks[[iblock]] <- future::future( { search_nn_annoy_index(query_matrix=query_matrix, nn_index=nn_index, metric=metric, k=k, search_k=search_k, beg_row_index=beg_block[[iblock]], end_row_index=end_block[[iblock]]) })
      }
      tot_bad <- 0
      for(iblock in seq(cores)) {
        nn_res <- future::value(nn_blocks[[iblock]])
        if(nrow(nn_res[['idx']]) != end_block[[iblock]] - beg_block[[iblock]] + 1) {
          stop('bad row count in nn_res')
        }
        idx[beg_block[[iblock]]:end_block[[iblock]],] <- nn_res[['idx']]
        dists[beg_block[[iblock]]:end_block[[iblock]],] <- nn_res[['dists']]
        tot_bad <- tot_bad + nn_res[['num_bad']]
      }
      if(tot_bad)
        stop('annoy was unable to find ', k, ' nearest neighbors for ', tot_bad, ' rows. You may need to increase the n_trees and/or search_k parameter values.')
      nn_res <- list(nn.idx=idx, nn.dists=dists)

      RhpcBLASctl::omp_set_num_threads(as.integer(omp_num_threads))
      RhpcBLASctl::blas_set_num_threads(as.integer(blas_num_threads))
    }
  }
  else
  if(nn_method == 'hnsw') {
    if(!test_hnsw_index(nn_index=nn_index, verbose=verbose)) {
      stop('search_nn_index: the hnsw nearest neighbor does not exist.')
    }

    assertthat::assert_that(nn_control[['ef']] >= k,
      msg=paste0('search_nn_index: ef must be >= k'))

    tmp <- RcppHNSW::hnsw_search(X=query_matrix,
                                  ann=nn_index[['hnsw_index']],
                                  k=k,
                                  ef=nn_control[['ef']],
                                  verbose=verbose,
                                  n_threads=cores,
                                  grain_size=nn_control[['grain_size']])
    # The RcppHNSW documentation says that the L2 metric is the square
    # of the Euclidean distance but my tests indicate that the L2 metric
    # returns sqrt(sum((SingleCellExperiment::reducedDims(cds)[['UMAP']][i,] - SingleCellExperiment::reducedDims(cds)[['UMAP']][j,])^2)).
    # Additionally, nn2, annoy, and hnsw return the same distance values for the euclidean
    # metric.
    nn_res <- list(nn.idx=tmp[['idx']], nn.dists=tmp[['dist']])
  }
  else
    stop('search_nn_index: unsupported nearest neighbor index type \'', nn_method, '\'')

  if(verbose) {
    tock()
  }

  return(nn_res)
}


#' @title Search a nearest neighbor index that is stored in
#'   the cds.
#'
#' @description Search a nearest neighbor index for cells near
#' those in the query_matrix.
#'
#' @param query_matrix a reduced dimension matrix used to find the
#'   nearest neighbors in the index nn_index.
#' @param cds a cell_data_set in which the nearest neighbor index
#'   is stored.
#' @param reduction_method a string giving the reduced
#'   dimension matrix used to make the nearest neighbor index, and
#'   determines where the index is stored in the cell_data_set.
#'   Note: distances in tSNE space reflect spatial differences poorly
#'   so using nearest neighbors with it may be meaningless.
#' @param k an integer for the number of nearest neighbors to return for
#'   each cell. Default is 25.
#' @param nn_control a list of parameters used to make and search
#'   the nearest neighbors indexes. See the set_nn_control help
#'   for additional details. Note that if nn_control\[\['search_k'\]\]
#'   is not defined, transfer_cell_labels will try to use
#'   search_k <- 2 * n_trees * k where n_trees is the value used
#'   to build the index. The default metric is cosine for
#'   reduction_methods PCA, LSI, and Aligned, and is euclidean for
#'   reduction_methods tSNE and UMAP.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return a list list(nn.idx, nn.dists) where nn.idx is
#'  a matrix of nearest neighbor indices and nn.dists is a matrix of
#'  the distance between the index given by the row number and the index
#'  given in nn.idx. If the same reduced dim matrix is used to make the
#'  index and search the index, the index given by the row number should
#'  be in the row, usually in the first column.
#'
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     cds <- preprocess_cds(cds)
#'     cds <- make_cds_nn_index(cds, 'PCA')
#'     nn_res <- search_cds_nn_index(SingleCellExperiment::reducedDims(cds)[['PCA']], cds, 'PCA', 10)
#'   }
#'
#' @export
search_cds_nn_index <- function(query_matrix, cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), k=25, nn_control=list(), verbose=FALSE) {
  assertthat::assert_that(methods::is(query_matrix, 'matrix') ||
                          is_sparse_matrix(query_matrix),
    msg=paste0('make_nn_matrix: the query_matrix object must be of type matrix'))
  
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
                          msg = paste0("When reduction_method = '", reduction_method,
                                      "' the cds must have been processed for it.",
                                      " Please run the required processing function",
                                      " before this one."))

  assertthat::assert_that(assertthat::is.count(k))

  if(reduction_method == 'tSNE' || reduction_method == 'UMAP')
    nn_control_default <- get_global_variable('nn_control_annoy_euclidean')
  else
    nn_control_default <- get_global_variable('nn_control_annoy_cosine')

  # Use set_nn_control to find nn method, which we need in order to select the correct index,
  # and we may need the index to set the nn_control[['search_k']].
  nn_control_tmp <- set_nn_control(mode=2,
                                   nn_control=nn_control,
                                   nn_control_default=nn_control_default,
                                   nn_index=NULL,
                                   k=k,
                                   verbose=verbose)
  nn_index <- get_cds_nn_index(cds=cds, reduction_method=reduction_method, nn_control_tmp[['method']], verbose=FALSE)

  nn_control <- set_nn_control(mode=2,
                               nn_control=nn_control,
                               nn_control_default=nn_control_default,
                               nn_index=nn_index,
                               k=k,
                               verbose=verbose)
  nn_res <- search_nn_index(query_matrix=query_matrix,
                            nn_index=nn_index,
                            k=k,
                            nn_control=nn_control,
                            verbose=verbose)

  return(nn_res)
}

# Notice that this function uses the available SingleCellExperiment::reducedDims(cds)[[reduction_method]] matrix
# for verifying the query matrix -- be certain that the search (and index build) were run
# on this matrix so call this function immediately after running the search.
# Notes:
#   o  this function was written to store nn search information in the cds when
#      we intended to re-use an index for later processing stages but
#      it appears that there is little or no opportunity to re-use
#      indices so we are not storing search information at this time.
# This function is not in use currently and may fall into disrepair.
set_cds_nn_search <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), search_id, ef, k, nn_control=list(), verbose=TRUE) {
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
                          msg = paste0("When reduction_method = '", reduction_method,
                                      "' the cds must have been processed for it.",
                                      " Please run the required processing function",
                                      " before this one."))
  assertthat::assert_that(!is.null(search_id),
                          msg = paste0('You must give a search_id string.'))

  assertthat::assert_that(!is.null(k),
                          msg = paste0('You must give a k value.'))

  reduced_matrix <- SingleCellExperiment::reducedDims(cds)[[reduction_method]]

  nn_control_default <- get_global_variable('nn_control_annoy_euclidean')
  nn_control <- set_nn_control(mode=2, nn_control=nn_control, nn_control_default=nn_control_default, nn_index=NULL, k=k, verbose=verbose)

  nn_method <- nn_control[['method']]

  if(is.null(cds@reduce_dim_aux[[reduction_method]][['nn_search']])) {
    cds@reduce_dim_aux[[reduction_method]][['nn_search']] <- S4Vectors::SimpleList()
  }

  if(nn_method == 'nn2') {
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]] <- S4Vectors::SimpleList()
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['method']] <- nn_method
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['k']] <- k
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['nrow']] <- nrow(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['ncol']] <- ncol(reduced_matrix)
  }
  else
  if(nn_method == 'annoy') {
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]] <- S4Vectors::SimpleList()
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['method']] <- nn_method
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['search_k']] <- nn_contol[['search_k']]
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['k']] <- k
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['nrow']] <- nrow(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['ncol']] <- ncol(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['grain_size']] <- nn_control[['grain_size']]
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['cores']] <- nn_control[['cores']]
  }
  else
  if(nn_method == 'hnsw') {
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]] <- S4Vectors::SimpleList()
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['method']] <- nn_method
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['ef']] <- ef
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['k']] <- k
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

  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)
  
  assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
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

  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
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
                       metric=cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['metric']],
                       n_trees=cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['n_trees']],
                       search_k=nn_search[['search_k']],
                       grain_size=nn_search[['grain_size']],
                       cores=nn_search[['cores']])
  }
  else
  if(nn_method == 'hnsw') {
    nn_control <- list(method=nn_method,
                       metric=cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['metric']],
                       M=cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['M']],
                       ef_construction=cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['ef_construction']],
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
  assertthat::assert_that(methods::is(subject_matrix, 'matrix') ||
                          is_sparse_matrix(subject_matrix),
    msg=paste0('search_nn_matrix: the subject_matrix object must be of type matrix'))

  assertthat::assert_that(methods::is(query_matrix, 'matrix') ||
                          is_sparse_matrix(query_matrix),
    msg=paste0('search_nn_matrix: the query_matrix object must be of type matrix'))

  nn_control_default <- get_global_variable('nn_control_annoy_euclidean')
  nn_control <- set_nn_control(mode=3, nn_control=nn_control, nn_control_default=nn_control_default, nn_index=NULL, k=k, verbose=verbose)

  method <- nn_control[['method']]
  k <- min(k, nrow(subject_matrix))

  if(verbose) {
    message('search_nn_matrix:')
    message('  k: ', k)
    report_nn_control('  nn_control: ', nn_control)
    tick('search_nn_matrix: search_time')
  }


  if(method == 'nn2') {
    nn_res <- RANN::nn2(subject_matrix, query_matrix, k, searchtype = "standard")
  }
  else {
    nn_index <- make_nn_index(subject_matrix, nn_control=nn_control, verbose=verbose)
    nn_res <- search_nn_index(query_matrix=query_matrix, nn_index=nn_index, k=k, nn_control=nn_control, verbose=verbose)
  }

  if(verbose)
    tock()

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
count_nn_missing_self_index <- function(nn_res, verbose=FALSE) {
  idx <- nn_res[['nn.idx']]
  dst <- nn_res[['nn.dists']]
  len <- length(idx[1,])
  num_missing <- 0

  if(nrow(idx) == 0)
    return(0)

  for (irow in 1:nrow(idx)) {
    vidx <- idx[irow,]
    vdst <- dst[irow,]
    if(vidx[[1]] != irow) {
      dself <- FALSE
      dzero <- TRUE
      if(length(vidx) == 1) {
        num_missing <- num_missing + 1
        next
      }
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

  if(verbose) {
    message('count_nn_missing_self_index:')
    message('  ', num_missing, ' out of ', nrow(nn_res[['nn.idx']]), ' rows are missing the row index')
    message('  \'recall\': ', formatC((as.double(nrow(nn_res[['nn.idx']]) - num_missing) / as.double(nrow(nn_res[['nn.idx']]))) * 100.0), "%")
  }

  return(num_missing)
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
    count_nn_missing_self_index(nn_res, verbose)
  }

  # Skip if there is no need to swap indices.
  if(check_nn_col1(nn_res$nn.idx))
    return(nn_res)

  idx <- nn_res[['nn.idx']]
  dst <- nn_res[['nn.dists']]

  diagnostics <- FALSE

  num_no_recall <- 0

  if(nrow(idx) == 0 )
    return(nn_res) 

  for (irow in 1:nrow(idx)) {
    vidx <- idx[irow,]
    vdst <- dst[irow,]
    if(vidx[[1]] != irow) {
      if(diagnostics) {
        message('swap_nn_row_index_point: adjust nn matrix row: ', irow)
        message('swap_nn_row_index_point: idx row pre fix: ', paste(vidx, collapse=' '))
      }

      if(length(vidx) == 1) {
        num_no_recall <- num_no_recall + 1
        next
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
          for(i in seq(1, length(vidx), 1)) {
            if(vdst[[i]] > .Machine$double.xmin) {
              num_no_recall <- num_no_recall + 1
            }
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

  if(num_no_recall > 0) {
    frac_recall <- (nrow(idx)-num_no_recall) / nrow(idx)
    format_recall <- sprintf('%3.1f', frac_recall * 100.0)
    message('The search result is expected to include the query row value (self)\n',
            'because the NN index includes the query objects; however, this search result\n',
            'is missing ', num_no_recall, ' self values (recall: ', format_recall, '%). Monocle3 has added the self\n',
            'values to the first column of the search result in order to allow further\n',
            'analysis -- but it is missing important nearest neighbors so you need to\n',
            'increase the sensitivity for making and/or searching the index.')
  }
  nn_res_out <- list(nn.idx=idx, nn.dists=dst)

  return(nn_res_out)
}

