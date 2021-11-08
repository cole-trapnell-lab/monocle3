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
  nn_control <- set_nn_control()
  expect_equal(nn_control[['method']], 'annoy')
  expect_equal(nn_control[['metric']], 'euclidean')
  expect_equal(nn_control[['n_trees']], 50)
  expect_equal(nn_control[['search_k']], 1250)
  expect_equal(nn_control[['grain_size']], 1)
  expect_equal(nn_control[['cores']], 1)

  nn_control <- set_nn_control(mode=3, k=10)
  expect_equal(nn_control[['method']], 'annoy')
  expect_equal(nn_control[['metric']], 'euclidean')
  expect_equal(nn_control[['n_trees']], 50)
  expect_equal(nn_control[['search_k']], 500)
  expect_equal(nn_control[['grain_size']], 1)
  expect_equal(nn_control[['cores']], 1)

  nn_control <- set_nn_control(mode=3, nn_control=list(grain_size=10, cores=20))
  expect_equal(nn_control[['method']], 'annoy')
  expect_equal(nn_control[['metric']], 'euclidean')
  expect_equal(nn_control[['n_trees']], 50)
  expect_equal(nn_control[['search_k']], 1250)
  expect_equal(nn_control[['grain_size']], 10)
  expect_equal(nn_control[['cores']], 20)

  nn_control <- set_nn_control(mode=3, nn_control=list(metric='cosine', n_trees=10, search_k=30, grain_size=10, cores=20))
  expect_equal(nn_control[['method']], 'annoy')
  expect_equal(nn_control[['metric']], 'cosine')
  expect_equal(nn_control[['n_trees']], 10)
  expect_equal(nn_control[['search_k']], 30)
  expect_equal(nn_control[['grain_size']], 10)
  expect_equal(nn_control[['cores']], 20)

  nn_control <- set_nn_control(mode=3, nn_control=list(method='nn2'))
  expect_equal(nn_control[['method']], 'nn2')

  nn_control <- set_nn_control(mode=3, nn_control=list(cores=1), nn_control_default=list(method='hnsw', metric='cosine'))
  expect_equal(nn_control[['method']], 'annoy')
  expect_equal(nn_control[['metric']], 'euclidean')
  expect_equal(nn_control[['n_trees']], 50)
  expect_equal(nn_control[['search_k']], 1250)
  expect_equal(nn_control[['grain_size']], 1)
  expect_equal(nn_control[['cores']], 1)

  nn_control <- set_nn_control(mode=3, nn_control_default=list(method='hnsw'))
  expect_equal(nn_control[['method']], 'hnsw')
  expect_equal(nn_control[['metric']], 'euclidean')
  expect_equal(nn_control[['M']], 48)
  expect_equal(nn_control[['ef_construction']], 200)
  expect_equal(nn_control[['ef']], 30)
  expect_equal(nn_control[['grain_size']], 1)
  expect_equal(nn_control[['cores']], 1)

  nn_control <- set_nn_control(mode=3, nn_control=list(method='annoy', metric='cosine', n_trees=10, search_k=30, grain_size=10, cores=20))
  expect_equal(nn_control[['method']], 'annoy')
  expect_equal(nn_control[['metric']], 'cosine')
  expect_equal(nn_control[['n_trees']], 10)
  expect_equal(nn_control[['search_k']], 30)
  expect_equal(nn_control[['grain_size']], 10)
  expect_equal(nn_control[['cores']], 20)

  nn_control <- set_nn_control(mode=3, nn_control=list(method='hnsw', metric='cosine', M=35, ef_construction=180, ef=45, grain_size=10, cores=20))
  expect_equal(nn_control[['method']], 'hnsw')
  expect_equal(nn_control[['metric']], 'cosine')
  expect_equal(nn_control[['M']], 35)
  expect_equal(nn_control[['ef_construction']], 180)
  expect_equal(nn_control[['ef']], 45)
  expect_equal(nn_control[['grain_size']], 10)
  expect_equal(nn_control[['cores']], 20)

  nn_control <- set_nn_control(mode=3, nn_control_default=list(method='annoy', metric='cosine', n_trees=10, search_k=30, grain_size=10, cores=20))
  expect_equal(nn_control[['method']], 'annoy')
  expect_equal(nn_control[['metric']], 'cosine')
  expect_equal(nn_control[['n_trees']], 10)
  expect_equal(nn_control[['search_k']], 30)
  expect_equal(nn_control[['grain_size']], 10)
  expect_equal(nn_control[['cores']], 20)

  nn_control <- set_nn_control(mode=3, nn_control_default=list(method='hnsw', metric='cosine', M=35, ef_construction=180, ef=45, grain_size=10, cores=20))
  expect_equal(nn_control[['method']], 'hnsw')
  expect_equal(nn_control[['metric']], 'cosine')
  expect_equal(nn_control[['M']], 35)
  expect_equal(nn_control[['ef_construction']], 180)
  expect_equal(nn_control[['ef']], 45)
  expect_equal(nn_control[['grain_size']], 10)
  expect_equal(nn_control[['cores']], 20)
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
  expect_equal(nn_res[['nn.idx']][[2,2]], 145)
  expect_equal(nn_res[['nn.dists']][[2,1]], 0, tol=1e-4)
  expect_equal(nn_res[['nn.dists']][[2,2]], 5.861312, tol=1e-3)

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
  expect_equal(nn_res[['nn.idx']][[2,2]], 145)
  expect_equal(nn_res[['nn.dists']][[2,1]], 0, tol=1e-4)
  expect_equal(nn_res[['nn.dists']][[2,2]], 5.861312, tol=1e-3)
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

