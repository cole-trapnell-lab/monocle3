#  set_nn_control
#  make_nn_index
#  set_cds_nn_index
#  make_cds_nn_index
#  search_nn_index
#  search_nn_matrix

context("test-nearest_neighbors")
skip_not_travis <- function ()
{
  if (identical(Sys.getenv("TRAVIS"), "true")) {
    return(invisible(TRUE))
  }
  skip("Not on Travis")
}

cds <- load_a549()
cds <- preprocess_cds(cds)

test_that("nearest_neighbors: set_nn_control", {

  #
  # General (independent) parameters.
  #

  # last resort default value for nn2 method
  lrdn_method <- 'nn2'

  # last resort default values and last resort default annoy method
  # default k
  lrda_method <- 'annoy'
  lrda_metric <- 'euclidean'
  lrda_n_trees <- 50
  lrda_search_k <- 2500
  lrda_grain_size <- 1
  lrda_cores <- 1

  # last resort default hnsw method
  lrdh_method <- 'hnsw'
  lrdh_metric <- 'euclidean'
  lrdh_M <- 48
  lrdh_ef_construction <- 200
  lrdh_ef <- 150
  lrdh_grain_size <- 1
  lrdh_cores <- 1

  # list (nn_control_default) annoy default values
  # search_k is for k=1
  ncda_method <- 'annoy'
  ncda_metric <- 'cosine'
  ncda_n_trees <- 30
  ncda_search_k <- 950
  ncda_grain_size <- 2
  ncda_cores <- 4

  # list (nn_control_default) hnsw default values
  ncdh_method <- 'hnsw'
  ncdh_metric <- 'cosine'
  ncdh_M <- 55
  ncdh_ef_construction <- 250
  ncdh_ef <- 175
  ncdh_grain_size <- 3
  ncdh_cores <- 6

  # list (nn_control) annoy non-default values
  # search_k is for k=1
  ncla_method <- 'annoy'
  ncla_metric <- 'hamming'
  ncla_n_trees <- 40
  ncla_search_k <- 1500
  ncla_grain_size <- 4
  ncla_cores <- 8

  # list (nn_control) hnsw non-default values
  nclh_method <- 'hnsw'
  nclh_metric <- 'l2'
  nclh_M <- 65 
  nclh_ef_construction <- 350
  nclh_ef <- 190
  nclh_grain_size <- 5
  nclh_cores <- 10

  # last resort defaults no method
  nn_control <- set_nn_control(mode=3)
  expect_equal(nn_control[['method']], lrda_method)
  expect_equal(nn_control[['metric']], lrda_metric)
  expect_equal(nn_control[['n_trees']], lrda_n_trees)
  expect_equal(nn_control[['search_k']], lrda_search_k)
  expect_equal(nn_control[['grain_size']], lrda_grain_size)
  expect_equal(nn_control[['cores']], lrda_cores)

  # last resort defaults nn2 method
  nn_control_in <- list(method=lrdn_method)
  nn_control <- set_nn_control(mode=3, nn_control=nn_control_in)
  expect_equal(nn_control[['method']], lrdn_method)

  # last resort defaults annoy method
  nn_control_in <- list(method='annoy')
  nn_control <- set_nn_control(mode=3, nn_control=nn_control_in)
  expect_equal(nn_control[['method']], lrda_method)
  expect_equal(nn_control[['metric']], lrda_metric)
  expect_equal(nn_control[['n_trees']], lrda_n_trees)
  expect_equal(nn_control[['search_k']], lrda_search_k)
  expect_equal(nn_control[['grain_size']], lrda_grain_size)
  expect_equal(nn_control[['cores']], lrda_cores)

  # last resort defaults hnsw method
  nn_control_in <- list(method='hnsw')
  nn_control <- set_nn_control(mode=3, nn_control=nn_control_in, k=2)
  expect_equal(nn_control[['method']], lrdh_method)
  expect_equal(nn_control[['metric']], lrdh_metric)
  expect_equal(nn_control[['M']], lrdh_M)
  expect_equal(nn_control[['ef_construction']], lrdh_ef_construction)
  expect_equal(nn_control[['ef']], lrdh_ef)
  expect_equal(nn_control[['grain_size']], lrdh_grain_size)
  expect_equal(nn_control[['cores']], lrdh_cores)

  # listed defaults (nn_control_default) annoy method
  nn_control_default_in <- list(method=ncda_method, metric=ncda_metric, n_trees=ncda_n_trees, search_k=ncda_search_k, grain_size=ncda_grain_size, cores=ncda_cores)
  nn_control <- set_nn_control(mode=3, nn_control_default=nn_control_default_in)
  expect_equal(nn_control[['method']], ncda_method)
  expect_equal(nn_control[['metric']], ncda_metric)
  expect_equal(nn_control[['n_trees']], ncda_n_trees)
  expect_equal(nn_control[['search_k']], ncda_search_k)
  expect_equal(nn_control[['grain_size']], ncda_grain_size)
  expect_equal(nn_control[['cores']], ncda_cores)

  # listed defaults (nn_control_default) hnsw method
  nn_control_default_in <- list(method=ncdh_method, metric=ncdh_metric, M=ncdh_M, ef_construction=ncdh_ef_construction, ef=ncdh_ef, grain_size=ncdh_grain_size, cores=ncdh_cores)
  nn_control <- set_nn_control(mode=3, nn_control_default=nn_control_default_in, k=2)
  expect_equal(nn_control[['method']], ncdh_method)
  expect_equal(nn_control[['metric']], ncdh_metric)
  expect_equal(nn_control[['M']], ncdh_M)
  expect_equal(nn_control[['ef_construction']], ncdh_ef_construction)
  expect_equal(nn_control[['ef']], ncdh_ef)
  expect_equal(nn_control[['grain_size']], ncdh_grain_size)
  expect_equal(nn_control[['cores']], ncdh_cores)

  # list(nn_control) method annoy and listed defaults (nn_control_default) annoy method
  nn_control_in <- list(method=ncla_method, metric=ncla_metric, n_trees=ncla_n_trees, search_k=ncla_search_k, grain_size=ncla_grain_size, cores=ncla_cores)
  nn_control_default_in <- list(method=ncda_method, metric=ncda_metric, n_trees=ncda_n_trees, search_k=ncda_search_k, grain_size=ncda_grain_size, cores=ncda_cores)
  nn_control <- set_nn_control(mode=3, nn_control=nn_control_in, nn_control_default=nn_control_default_in)
  expect_equal(nn_control[['method']], ncla_method)
  expect_equal(nn_control[['metric']], ncla_metric)
  expect_equal(nn_control[['n_trees']], ncla_n_trees)
  expect_equal(nn_control[['search_k']], ncla_search_k)
  expect_equal(nn_control[['grain_size']], ncla_grain_size)
  expect_equal(nn_control[['cores']], ncla_cores)

  # list(nn_control) method hnsw and listed defaults (nn_control_default) hnsw method
  nn_control_in <- list(method=nclh_method, metric=nclh_metric, M=nclh_M, ef_construction=nclh_ef_construction, ef=nclh_ef, grain_size=nclh_grain_size, cores=nclh_cores)
  nn_control_default_in <- list(method=ncdh_method, metric=ncdh_metric, M=ncdh_M, ef_construction=ncdh_ef_construction, ef=ncdh_ef, grain_size=ncdh_grain_size, cores=ncdh_cores)
  nn_control <- set_nn_control(mode=3, nn_control=nn_control_in, nn_control_default=nn_control_default_in, k=2)
  expect_equal(nn_control[['method']], nclh_method)
  expect_equal(nn_control[['metric']], nclh_metric)
  expect_equal(nn_control[['M']], nclh_M)
  expect_equal(nn_control[['ef_construction']], nclh_ef_construction)
  expect_equal(nn_control[['ef']], nclh_ef)
  expect_equal(nn_control[['grain_size']], nclh_grain_size)
  expect_equal(nn_control[['cores']], nclh_cores)

  # list(nn_control) method annoy and listed defaults (nn_control_default) hnsw method
  nn_control_in <- list(method=ncla_method, metric=ncla_metric, n_trees=ncla_n_trees, search_k=ncla_search_k, grain_size=ncla_grain_size, cores=ncla_cores)
  nn_control_default_in <- list(method=ncdh_method, metric=ncdh_metric, M=ncdh_M, ef_construction=ncdh_ef_construction, ef=ncdh_ef, grain_size=ncdh_grain_size, cores=ncdh_cores)
  nn_control <- set_nn_control(mode=3, nn_control=nn_control_in, nn_control_default=nn_control_default_in)
  expect_equal(nn_control[['method']], ncla_method)
  expect_equal(nn_control[['metric']], ncla_metric)
  expect_equal(nn_control[['n_trees']], ncla_n_trees)
  expect_equal(nn_control[['search_k']], ncla_search_k)
  expect_equal(nn_control[['grain_size']], ncla_grain_size)
  expect_equal(nn_control[['cores']], ncla_cores)

  # list(nn_control) method hnsw and listed defaults (nn_control_default) annoy method
  nn_control_in <- list(method=nclh_method, metric=nclh_metric, M=nclh_M, ef_construction=nclh_ef_construction, ef=nclh_ef, grain_size=nclh_grain_size, cores=nclh_cores)
  nn_control_default_in <- list(method=ncda_method, metric=ncda_metric, n_trees=ncda_n_trees, search_k=ncda_search_k, grain_size=ncda_grain_size, cores=ncda_cores)
  nn_control <- set_nn_control(mode=3, nn_control=nn_control_in, nn_control_default=nn_control_default_in, k=2)
  expect_equal(nn_control[['method']], nclh_method)
  expect_equal(nn_control[['metric']], nclh_metric)
  expect_equal(nn_control[['M']], nclh_M)
  expect_equal(nn_control[['ef_construction']], nclh_ef_construction)
  expect_equal(nn_control[['ef']], nclh_ef)
  expect_equal(nn_control[['grain_size']], nclh_grain_size)
  expect_equal(nn_control[['cores']], nclh_cores)

  #
  # Parameters with dependencies on other parameters (annoy search_k).
  # Note: there are redundant tests.
  #
  cds_sk <- load_a549()
  reduction_method <- 'PCA'
  nn_control <- list(method='annoy', metric='euclidean', n_trees=80)
  cds_sk <- preprocess_cds(cds_sk, method=reduction_method, nn_control=nn_control, build_nn_index=TRUE)


  # mode=3, cds, reduction_method and nn_control[['search_k']] only
  nn_control_in <- list(method='annoy', search_k=56)
  nn_control <- set_nn_control(mode=3, nn_control=nn_control_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 56)
  
  # mode=3, cds, reduction_method and nn_control[['n_trees']] only
  nn_control_in <- list(method='annoy', n_trees=22)
  nn_control <- set_nn_control(mode=3, nn_control=nn_control_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 88)

  # mode=3, cds, reduction_method and nn_control[['search_k']] and nn_control[['n_trees']]
  nn_control_in <- list(method='annoy', search_k=56, n_trees=32)
  nn_control <- set_nn_control(mode=3, nn_control=nn_control_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 56)


  # mode=3, cds, reduction_method and nn_control_default[['search_k']] only
  nn_control_default_in <- list(method='annoy', search_k=58)
  nn_control <- set_nn_control(mode=3, nn_control_default=nn_control_default_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 58)

  # mode=3, cds, reduction_method and nn_control_default[['n_trees']] only
  nn_control_default_in <- list(method='annoy', n_trees=28)
  nn_control <- set_nn_control(mode=3, nn_control_default=nn_control_default_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 112)

  # mode=3, cds, reduction_method and nn_control_default[['search_k']] and nn_control_default[['n_trees']]
  nn_control_default_in <- list(method='annoy', search_k=60, n_trees=40)
  nn_control <- set_nn_control(mode=3, nn_control_default=nn_control_default_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 60)


  # mode=3, cds, reduction_method and nn_control[['search_k']] and nn_control_default[['search_k']]
  nn_control_in <- list(method='annoy', search_k=65)
  nn_control_default_in <- list(method='annoy', search_k=15)
  nn_control <- set_nn_control(mode=3, nn_control=nn_control_in, nn_control_default=nn_control_default_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 65)

  # mode=3, cds, reduction_method and nn_control[['n_trees']] and nn_control_default[['n_trees']]
  nn_control_in <- list(method='annoy', n_trees=15)
  nn_control_default_in <- list(method='annoy', n_trees=25)
  nn_control <- set_nn_control(mode=3, nn_control=nn_control_in, nn_control_default=nn_control_default_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 60)

  # mode=3, cds, reduction_method and nn_control[['search_k']] and nn_control_default[['n_trees']]
  nn_control_in <- list(method='annoy', search_k=60)
  nn_control_default_in <- list(method='annoy', n_trees=40)
  nn_control <- set_nn_control(mode=3, nn_control=nn_control_in, nn_control_default=nn_control_default_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 60)

  # mode=3, cds, reduction_method and nn_control[['n_trees']] and nn_control_default[['search_k']]
  nn_control_in <- list(method='annoy', n_trees=10)
  nn_control_default_in <- list(method='annoy', search_k=48)
  nn_control <- set_nn_control(mode=3, nn_control=nn_control_in, nn_control_default=nn_control_default_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 40)



  # mode=2, cds, reduction_method and nn_control[['search_k']] only
  nn_control_in <- list(method='annoy', search_k=56)
  nn_control <- set_nn_control(mode=2, nn_control=nn_control_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 56)

  # mode=2, cds, reduction_method and nn_control[['n_trees']] only
  nn_control_in <- list(method='annoy', n_trees=22)
  nn_control <- set_nn_control(mode=2, nn_control=nn_control_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 320)

  # mode=2, cds, reduction_method and nn_control[['search_k']] and nn_control[['n_trees']]
  nn_control_in <- list(method='annoy', search_k=56, n_trees=32)
  nn_control <- set_nn_control(mode=2, nn_control=nn_control_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 56)


  # mode=2, cds, reduction_method and nn_control_default[['search_k']] only
  nn_control_default_in <- list(method='annoy', search_k=58)
  nn_control <- set_nn_control(mode=2, nn_control_default=nn_control_default_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 320)

  # mode=2, cds, reduction_method and nn_control_default[['n_trees']] only
  nn_control_default_in <- list(method='annoy', n_trees=28)
  nn_control <- set_nn_control(mode=2, nn_control_default=nn_control_default_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 320)

  # mode=2, cds, reduction_method and nn_control_default[['search_k']] and nn_control_default[['n_trees']]
  nn_control_default_in <- list(method='annoy', search_k=60, n_trees=40)
  nn_control <- set_nn_control(mode=2, nn_control_default=nn_control_default_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 320)


  # mode=2, cds, reduction_method and nn_control[['search_k']] and nn_control_default[['search_k']]
  nn_control_in <- list(method='annoy', search_k=60)
  nn_control_default_in <- list(method='annoy', search_k=15)
  nn_control <- set_nn_control(mode=2, nn_control=nn_control_in, nn_control_default=nn_control_default_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 60)

  # mode=2, cds, reduction_method and nn_control[['n_trees']] and nn_control_default[['n_trees']]
  nn_control_in <- list(method='annoy', n_trees=15)
  nn_control_default_in <- list(method='annoy', n_trees=25)
  nn_control <- set_nn_control(mode=2, nn_control=nn_control_in, nn_control_default=nn_control_default_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 320)

  # mode=2, cds, reduction_method and nn_control[['search_k']] and nn_control_default[['n_trees']]
  nn_control_in <- list(method='annoy', search_k=60)
  nn_control_default_in <- list(method='annoy', n_trees=40)
  nn_control <- set_nn_control(mode=2, nn_control=nn_control_in, nn_control_default=nn_control_default_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 60)

  # mode=2, cds, reduction_method and nn_control[['n_trees']] and nn_control_default[['search_k']]
  nn_control_in <- list(method='annoy', n_trees=10)
  nn_control_default_in <- list(method='annoy', search_k=48)
  nn_control <- set_nn_control(mode=2, nn_control=nn_control_in, nn_control_default=nn_control_default_in, cds=cds_sk, reduction_method=reduction_method, k=2)
  expect_equal(nn_control[['search_k']], 320)

} )


test_that("nearest_neighbors: make and search nn2 index", {
  nn_control=list(method='nn2')

  # The following functions are invalid with nn2 (nn2 does not return an index).
  expect_error(make_nn_index(reducedDims(cds)[['PCA']], nn_control=nn_control))
  expect_error(set_cds_nn_index(cds, reduction_method='PCA', nn_control=nn_control))
  expect_error(make_cds_nn_index(cds, reduction_method='PCA', nn_control=nn_control))
  expect_error(search_nn_index(reducedDims(cds)[['PCA']], nn_control=nn_control))

  nn_res <- search_nn_matrix(reducedDims(cds)[['PCA']], reducedDims(cds)[['PCA']], k=2, nn_control=nn_control)
  expect_equal(nn_res[['nn.idx']][[2,1]], 2)
  expect_equal(nn_res[['nn.idx']][[2,2]], 264)
  expect_equal(nn_res[['nn.dists']][[2,1]], 0, tol=1e-4)
  expect_equal(nn_res[['nn.dists']][[2,2]], 5.402802, tol=1e-3)
} )


test_that("nearest_neighbors: make and search annoy index", {
  nn_control=list(method='annoy')

  # make annoy index
  nn_index <- make_nn_index(reducedDims(cds)[['PCA']], nn_control=nn_control)
  if(!is.null(nn_index[['version']]))
    expect_equal(as.character(class(nn_index[['annoy_index']])), 'Rcpp_AnnoyEuclidean')
  else
  if(!is.null(nn_index[['type']]))
    expect_equal(as.character(class(nn_index[['ann']])), 'Rcpp_AnnoyEuclidean')
  else
    expect_equal(as.character(class(nn_index)), 'Rcpp_AnnoyEuclidean')

  # search annoy index
  nn_res <- search_nn_index(reducedDims(cds)[['PCA']], nn_index=nn_index, k=2, nn_control=nn_control)
  expect_equal(nn_res[['nn.idx']][[2,1]], 2)
  expect_equal(nn_res[['nn.idx']][[2,2]], 264)
  expect_equal(nn_res[['nn.dists']][[2,1]], 0, tol=1e-4)
  expect_equal(nn_res[['nn.dists']][[2,2]], 5.402802, tol=1e-3)

  # store annoy index in cds
  cds2 <- set_cds_nn_index(cds, reduction_method='PCA', nn_index, nn_control=nn_control)
  if(!is.null(cds2@reduce_dim_aux[['PCA']][['nn_index']][['annoy']][['nn_index']][['version']]))
    expect_equal(as.character(class(cds2@reduce_dim_aux[['PCA']][['nn_index']][['annoy']][['nn_index']][['annoy_index']])), 'Rcpp_AnnoyEuclidean')
  else
  if(!is.null(cds2@reduce_dim_aux[['PCA']][['nn_index']][['annoy']][['nn_index']][['type']]))
    expect_equal(as.character(class(cds2@reduce_dim_aux[['PCA']][['nn_index']][['annoy']][['nn_index']][['ann']])), 'Rcpp_AnnoyEuclidean')
  else
    expect_equal(as.character(class(cds2@reduce_dim_aux[['PCA']][['nn_index']][['annoy']][['nn_index']])), 'Rcpp_AnnoyEuclidean')

  nn_index2 <- get_cds_nn_index(cds2, reduction_method='PCA', nn_control=nn_control)
  if(!is.null(nn_index2[['version']]))
    expect_equal(as.character(class(nn_index2[['annoy_index']])), 'Rcpp_AnnoyEuclidean')
  else
  if(!is.null(nn_index2[['type']]))
    expect_equal(as.character(class(nn_index2[['ann']])), 'Rcpp_AnnoyEuclidean')
  else
    expect_equal(as.character(class(nn_index2)), 'Rcpp_AnnoyEuclidean')

  rm(cds2)
  cds2 <- make_cds_nn_index(cds, reduction_method='PCA', nn_control=nn_control)
  if(!is.null(cds2@reduce_dim_aux[['PCA']][['nn_index']][['annoy']][['nn_index']][['version']]))
    expect_equal(as.character(class(cds2@reduce_dim_aux[['PCA']][['nn_index']][['annoy']][['nn_index']][['annoy_index']])), 'Rcpp_AnnoyEuclidean')
  else
  if(!is.null(cds2@reduce_dim_aux[['PCA']][['nn_index']][['annoy']][['nn_index']][['type']]))
    expect_equal(as.character(class(cds2@reduce_dim_aux[['PCA']][['nn_index']][['annoy']][['nn_index']][['ann']])), 'Rcpp_AnnoyEuclidean')
  else
    expect_equal(as.character(class(cds2@reduce_dim_aux[['PCA']][['nn_index']][['annoy']][['nn_index']])), 'Rcpp_AnnoyEuclidean')

  rm(cds2)
  nn_res <- search_nn_matrix(reducedDims(cds)[['PCA']], reducedDims(cds)[['PCA']], k=2, nn_control=nn_control)
  expect_equal(nn_res[['nn.idx']][[2,1]], 2)
  expect_equal(nn_res[['nn.idx']][[2,2]], 264)
  expect_equal(nn_res[['nn.dists']][[2,1]], 0, tol=1e-4)
  expect_equal(nn_res[['nn.dists']][[2,2]], 5.402802, tol=1e-3)
} )


test_that("nearest_neighbors: make and search hnsw index", {
  nn_control=list(method='hnsw')

  # make hnsw index
  nn_index <- make_nn_index(reducedDims(cds)[['PCA']], nn_control=nn_control)
  expect_equal(as.character(class(nn_index)), 'Rcpp_HnswL2')

  # search hnsw index
  nn_res <- search_nn_index(reducedDims(cds)[['PCA']], nn_index=nn_index, k=2, nn_control=nn_control)
  expect_equal(nn_res[['nn.idx']][[2,1]], 2)
  expect_equal(nn_res[['nn.idx']][[2,2]], 264)
  expect_equal(nn_res[['nn.dists']][[2,1]], 0, tol=1e-4)
  expect_equal(nn_res[['nn.dists']][[2,2]], 5.402802, tol=1e-3)

  # store hnsw index in cds
  cds2 <- set_cds_nn_index(cds, reduction_method='PCA', nn_index, nn_control=nn_control)
  expect_equal(as.character(class(cds2@reduce_dim_aux[['PCA']][['nn_index']][['hnsw']][['nn_index']])), 'Rcpp_HnswL2')

  nn_index2 <- get_cds_nn_index(cds2, reduction_method='PCA', nn_control=nn_control)
  expect_equal(as.character(class(nn_index2)), 'Rcpp_HnswL2')

  rm(cds2)
  cds2 <- make_cds_nn_index(cds, reduction_method='PCA', nn_control=nn_control)
  expect_equal(as.character(class(cds2@reduce_dim_aux[['PCA']][['nn_index']][['hnsw']][['nn_index']])), 'Rcpp_HnswL2')

  rm(cds2)
  nn_res <- search_nn_matrix(reducedDims(cds)[['PCA']], reducedDims(cds)[['PCA']], k=2, nn_control=nn_control)
  expect_equal(nn_res[['nn.idx']][[2,1]], 2)
  expect_equal(nn_res[['nn.idx']][[2,2]], 264)
  expect_equal(nn_res[['nn.dists']][[2,1]], 0, tol=1e-4)
  expect_equal(nn_res[['nn.dists']][[2,2]], 5.402802, tol=1e-3)
} )


test_that("nearest_neighbors: functions related to the self indices", {
  # Set up a smallish nn_res.
  vidx <- c(1,2,3,4,5,6,7,8,9,10,1082299,11,7140,8488,938113,944258,1183327,707,22278,1079854,28288,13372,450,2205,929641,900494,510,2379,4976,4977,13560,23868,8942,1070074,1110199,996504,6649,8554,17729,4140,1083654,1193086,6427,37,853726,1044200,520564,8634,4363,55322,1063210,9057,11195,12794,1083333,600384,813713,4349,584,41427,883861,40035,13603,15212,802504,1061316,156137,10925,4531,12911,4883,1214856,511,6872,17652,877827,142679,4896,2216,1023234,18002,2147,4516,2793,1154358,918134,1184998,175,6331,1213981,770782,414630,13327,8662,29335,23502,800610,2371,13629,982962)

  vdists <- c(0,0,0,0,0,0,0,0,0,0,8.18230924837125e-07,0,7.94866505828386e-08,7.70586816980026e-08,6.93755389472248e-08,6.0690216875546e-08,6.97396482463265e-07,7.93544534272153e-08,5.82226468105725e-07,1.75800822288692e-06,9.18039888739163e-07,4.79849942160069e-07,1.5804363350254e-07,4.20261223106574e-07,2.77153922903407e-07,6.09892690893894e-08,7.4050774954391e-07,8.02054843660267e-08,7.13906380813785e-07,2.06462860220845e-06,1.24677772220838e-06,1.36033158839539e-06,1.59849439480597e-07,5.04832893965684e-07,4.85141216647906e-07,6.10340016904626e-08,1.02489688626932e-06,2.45138332555463e-07,7.82324556698744e-07,4.43321337124591e-06,1.78286515626898e-06,1.55374203746194e-06,2.37100294695554e-07,5.35215402581291e-07,5.5445008731204e-07,1.82931271605295e-07,1.5746577552703e-06,3.25164911561968e-07,1.18820082265468e-06,4.80286045137078e-06,3.2795869439129e-06,2.28994249718744e-06,2.40478460585685e-07,7.55950124874463e-07,6.2391381896578e-07,2.42388356784464e-07,1.89901019641466e-06,4.84930242215537e-07,1.34027845492523e-06,5.46432854944834e-06,3.50969080475938e-06,2.8877005600703e-06,3.18955195631812e-07,8.21225531882101e-07,6.93172527322793e-07,2.44938909038538e-07,2.5393893801542e-06,5.54079068339843e-07,1.51186630856411e-06,5.59836316977065e-06,4.36074021362558e-06,4.6215640507369e-06,3.93205726694072e-07,1.06735125148503e-06,6.93332498493279e-07,2.45914671629949e-07,2.87381105201941e-06,6.31157055986449e-07,1.52915833121767e-06,7.61886076409296e-06,4.74314888990581e-06,7.16636542480626e-06,4.02745466835228e-07,1.15852933090012e-06,7.62132054795942e-07,4.28356585894972e-07,4.73208954990671e-06,6.68522468215746e-07,3.13471585617744e-06,8.64934257809827e-06,5.03780002379325e-06,7.2199321173982e-06,5.51644609375164e-07,1.28492204118337e-06,8.31504106999057e-07,5.52185830311335e-07,5.55015313937859e-06,6.71790384510988e-07,3.5507210761803e-06,9.06252936240797e-06)

  nn_res <- list()
  nn_res[['nn.idx']] <- matrix(vidx, ncol=10)
  nn_res[['nn.dists']] <- matrix(vdists, ncol=10)

  # Check that check_nn_col1 returns TRUE when all self indices
  # are in the first column and count_nn_missing_self_index
  # returns zero.
  expect_true(check_nn_col1(nn_res[['nn.idx']]))
  expect_equal(count_nn_missing_self_index(nn_res),0)

  # Don't alter nn_res if the self indices are in the first column.
  nn_res.swap <- swap_nn_row_index_point(nn_res)
  expect_equal(nn_res.swap[['nn.idx']], nn_res[['nn.idx']])
  expect_equal(nn_res.swap[['nn.dists']], nn_res[['nn.dists']])

  # Move two self indices to interior columns and zero dists accordingly.
  nn_res_v1 <- nn_res
  nn_res_v1[['nn.idx']][2,1:3] <- c(11, 13372, 2)
  nn_res_v1[['nn.idx']][6,1:2] <- c(944258, 6)
  nn_res_v1[['nn.dists']][2,1:3] <- 0
  nn_res_v1[['nn.dists']][6,1:2] <- 0

  # Check that check_nn_col1 returns FALSE and
  # count_nn_missing_self_index returns 0.
  expect_false(check_nn_col1(nn_res_v1[['nn.idx']]))
  expect_equal(count_nn_missing_self_index(nn_res_v1),0)

  # Swap indices if they are in the row but not the first column.
  nn_res.swap <- swap_nn_row_index_point(nn_res_v1)
  expect_true(check_nn_col1(nn_res.swap[['nn.idx']]))
  expect_equal(count_nn_missing_self_index(nn_res.swap),0)
  expect_equal(nn_res.swap[['nn.dists']], nn_res_v1[['nn.dists']])

  # Remove self indices from two rows and set
  # all distances in those rows to zero. (We
  # assume that the search found the self
  # indices but they were not in the first k
  # nns.)
  nn_res_v2 <- nn_res
  nn_res_v2[['nn.idx']][2,1] <- 1121
  nn_res_v2[['nn.idx']][6,1] <- 1125
  nn_res_v2[['nn.dists']][2,] <- 0
  nn_res_v2[['nn.dists']][6,] <- 0

  # Check that check_nn_col1 returns FALSE and 
  # count_nn_missing_self_index returns 0.
  expect_false(check_nn_col1(nn_res_v2[['nn.idx']]))
  expect_equal(count_nn_missing_self_index(nn_res_v2),0)

  nn_res.swap <- swap_nn_row_index_point(nn_res_v2)
  expect_true(check_nn_col1(nn_res.swap[['nn.idx']]))
  expect_equal(count_nn_missing_self_index(nn_res.swap),0)
  expect_equal(nn_res.swap[['nn.dists']], nn_res_v2[['nn.dists']])

  # Remove self indices from two rows and set
  # one distance in each those rows to non-zero
  # values. (The search missed the self indices.)
  nn_res_v3 <- nn_res_v2
  nn_res_v3[['nn.dists']][2,10] <- .001
  nn_res_v3[['nn.dists']][6,8] <- 5

  # Check that check_nn_col1 returns FALSE and
  # count_nn_missing_self_index returns 0.
  expect_false(check_nn_col1(nn_res_v3[['nn.idx']]))
  expect_equal(count_nn_missing_self_index(nn_res_v3),2)

  expect_message(swap_nn_row_index_point(nn_res_v3))
  nn_res.swap <- swap_nn_row_index_point(nn_res_v3)
  expect_true(check_nn_col1(nn_res.swap[['nn.idx']]))
  expect_equal(count_nn_missing_self_index(nn_res.swap),0)
  irow <- 2
  vidx <- nn_res_v3[['nn.idx']][irow,]
  expect_equal(nn_res.swap[['nn.idx']][irow,], c(irow, vidx[1:(length(vidx)-1)]))
  vdists <- nn_res_v3[['nn.dists']][irow,]
  expect_equal(nn_res.swap[['nn.dists']][irow,], c(0, vdists[1:(length(vdists)-1)]))
  irow <- 6
  vdists <- nn_res_v3[['nn.dists']][irow,]
  expect_equal(nn_res.swap[['nn.dists']][irow,], c(0, vdists[1:(length(vdists)-1)]))
} )

