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
  expect_equal(as.character(class(nn_index[['ann']])), 'Rcpp_AnnoyEuclidean')

  # search annoy index
  nn_res <- search_nn_index(reducedDims(cds)[['PCA']], nn_index=nn_index, k=2, nn_control=nn_control)
  expect_equal(nn_res[['nn.idx']][[2,1]], 2)
  expect_equal(nn_res[['nn.idx']][[2,2]], 264)
  expect_equal(nn_res[['nn.dists']][[2,1]], 0, tol=1e-4)
  expect_equal(nn_res[['nn.dists']][[2,2]], 5.402802, tol=1e-3)

  # store annoy index in cds
  cds2 <- set_cds_nn_index(cds, reduction_method='PCA', nn_index, nn_control=nn_control)
  expect_equal(as.character(class(cds2@reduce_dim_aux[['PCA']][['nn_index']][['annoy']][['nn_index']][['ann']])), 'Rcpp_AnnoyEuclidean')

  nn_index2 <- get_cds_nn_index(cds2, reduction_method='PCA', nn_control=nn_control)
  expect_equal(as.character(class(nn_index2[['ann']])), 'Rcpp_AnnoyEuclidean')

  rm(cds2)
  cds2 <- make_cds_nn_index(cds, reduction_method='PCA', nn_control=nn_control)
  expect_equal(as.character(class(cds2@reduce_dim_aux[['PCA']][['nn_index']][['annoy']][['nn_index']][['ann']])), 'Rcpp_AnnoyEuclidean')

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

