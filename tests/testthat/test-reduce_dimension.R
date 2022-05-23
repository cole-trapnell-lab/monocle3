context("test-reduce_dimension")
skip_not_travis <- function ()
{
  if (identical(Sys.getenv("TRAVIS"), "true")) {
    return(invisible(TRUE))
  }
  skip("Not on Travis")
}
set.seed(100)
cds <- load_a549()
cds <- estimate_size_factors(cds)


test_that('Nearest neighbors', {
  cds <- preprocess_cds(cds)
  cds <- reduce_dimension(cds, build_nn_index=TRUE)
  expect_equal(cds@reduce_dim_aux[['UMAP']][['nn_index']][['annoy']][['nn_index']][['metric']], 'euclidean')
})


test_that("reduce_dimension runs", {
  skip_on_travis()
  expect_error(cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1),
               "Data has not been preprocessed with chosen method: PCA Please run preprocess_cds with method = PCA before running reduce_dimension.")
  cds <- preprocess_cds(cds, num_dim = 20)

  cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1)
  expect_equal(nrow(reducedDims(cds)$UMAP), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$UMAP), 2)
  expect_equal(as.numeric(reducedDims(cds)$UMAP[1,1]), -2.86,
               tolerance = 1e-2)

  cds <- reduce_dimension(cds, max_components = 3, umap.fast_sgd=FALSE, cores=1, reduction_method = "UMAP")
  expect_equal(nrow(reducedDims(cds)$UMAP), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$UMAP), 3)
  expect_equal(as.numeric(reducedDims(cds)$UMAP[1,1]), 1.69,
               tolerance = 1e-2)

  cds <- reduce_dimension(cds, reduction_method = "tSNE")
  expect_equal(nrow(reducedDims(cds)$tSNE), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$tSNE), 2)
  expect_equal(as.numeric(reducedDims(cds)$tSNE[1,1]), -1.16,
               tolerance = 1e-2)

  cds <- reduce_dimension(cds,  max_components = 3, reduction_method = "tSNE")
  expect_equal(nrow(reducedDims(cds)$tSNE), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$tSNE), 3)
  expect_equal(as.numeric(reducedDims(cds)$tSNE[1,1]), -4.29,
               tolerance = 1e-2)

  cds <- reduce_dimension(cds, reduction_method = "PCA")
  expect_equal(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equal(as.numeric(reducedDims(cds)$PCA["E11_A01_RT_467","PC1"]),
               2.4207391, tolerance = 1e-4)

  cds <- preprocess_cds(cds, num_dim = 20, method = "LSI")
  cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1, preprocess_method = "LSI")
  expect_equal(nrow(reducedDims(cds)$UMAP), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$UMAP), 2)
  expect_equal(as.numeric(reducedDims(cds)$UMAP[1,1]), -0.163,
               tolerance = 1e-3)

  cds <- reduce_dimension(cds, reduction_method = "tSNE", preprocess_method = "LSI")
  expect_equal(nrow(reducedDims(cds)$tSNE), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$tSNE), 2)
  expect_equal(as.numeric(reducedDims(cds)$tSNE[1,1]), -2.15,
               tolerance = 1e-2)

  # check model
  set.seed(100)
  cds <- load_a549()
  cds <- estimate_size_factors(cds)
  cds <- preprocess_cds(cds, num_dim = 20)
  cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1)
  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['umap_preprocess_method']], 'PCA')
  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['max_components']], 2, tol=1e1)
  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['umap_metric']], 'cosine')
  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['umap_min_dist']], 0.1, tol=1e-1)
  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['umap_n_neighbors']], 15, tol=1e1)
  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['umap_fast_sgd']], FALSE)
  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['umap_model']][['embedding']][[1,1]], -1.80, tol=1e-1)
  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['umap_model']][['n_neighbors']][[1]], 15, tol=1e1)


  expect_error(reduce_dimension(cds, reduction_method = "DDRTree"),
               "reduction_method must be one of 'UMAP', 'PCA', 'tSNE', 'LSI', 'Aligned'")
})

test_that("reduce_dimension clears old graphs", {
  skip_on_travis()
  cds <- preprocess_cds(cds)
  cds <- reduce_dimension(cds)
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  cds <- reduce_dimension(cds)
  expect_null(cds@principal_graph_aux[["UMAP"]])
  expect_error(partitions(cds), "No partitions calculated for reduction_method = UMAP. Please first run cluster_cells with reduction_method = UMAP.")
  expect_error(clusters(cds), "No clusters calculated for reduction_method = UMAP. Please first run cluster_cells with reduction_method = UMAP.")
  expect_null(cds@principal_graph[["UMAP"]])

  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  cds <- reduce_dimension(cds, reduction_method = "tSNE")
  expect_is(clusters(cds), "factor")
  testthat::expect_is(cds@principal_graph[["UMAP"]], "igraph")
})


#### TRAVIS ####


set.seed(100)
cds <- load_a549()
cds <- estimate_size_factors(cds)

test_that("reduce_dimension runs", {
  skip_not_travis()
  expect_error(cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1),
               "Data has not been preprocessed with chosen method: PCA Please run preprocess_cds with method = PCA before running reduce_dimension.")
  cds <- preprocess_cds(cds, num_dim = 20)

  cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1)
  expect_equal(nrow(reducedDims(cds)$UMAP), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$UMAP), 2)
  expect_equal(as.numeric(reducedDims(cds)$UMAP[1,1]), -2.14,
               tolerance = 1e-2)

  cds <- reduce_dimension(cds, max_components = 3, umap.fast_sgd=FALSE, cores=1, reduction_method = "UMAP")
  expect_equal(nrow(reducedDims(cds)$UMAP), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$UMAP), 3)
  expect_equal(as.numeric(reducedDims(cds)$UMAP[1,1]), 1.65,
               tolerance = 1e-2)

  cds <- reduce_dimension(cds, reduction_method = "tSNE")
  expect_equal(nrow(reducedDims(cds)$tSNE), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$tSNE), 2)
  expect_equal(as.numeric(reducedDims(cds)$tSNE[1,1]), -2.38,
               tolerance = 1e-2)

  cds <- reduce_dimension(cds,  max_components = 3, reduction_method = "tSNE")
  expect_equal(nrow(reducedDims(cds)$tSNE), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$tSNE), 3)
  expect_equal(as.numeric(reducedDims(cds)$tSNE[1,1]),  0.455,
               tolerance = 1e-1)

  cds <- reduce_dimension(cds, reduction_method = "PCA")
  expect_equal(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equal(as.numeric(reducedDims(cds)$PCA["E11_A01_RT_467","PC1"]),
               2.4207391, tolerance = 1e-4)

  cds <- preprocess_cds(cds, num_dim = 20, method = "LSI")
  cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1, preprocess_method = "LSI")
  expect_equal(nrow(reducedDims(cds)$UMAP), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$UMAP), 2)
  expect_equal(as.numeric(reducedDims(cds)$UMAP[1,1]), -0.0531,
               tolerance = 1e-4)

  cds <- reduce_dimension(cds, reduction_method = "tSNE", preprocess_method = "LSI")
  expect_equal(nrow(reducedDims(cds)$tSNE), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$tSNE), 2)
  expect_equal(as.numeric(reducedDims(cds)$tSNE[1,1]), 0.204,
               tolerance = 1e-2)

  # check model
  set.seed(100)
  cds <- load_a549()
  cds <- estimate_size_factors(cds)
  cds <- preprocess_cds(cds, num_dim = 20)
  cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1)
  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['umap_preprocess_method']], 'PCA')
  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['max_components']], 2, tol=1e1)
  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['umap_metric']], 'cosine')
  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['umap_min_dist']], 0.1, tol=1e-1)
  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['umap_n_neighbors']], 15, tol=1e1)
  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['umap_fast_sgd']], FALSE)
  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['umap_model']][['embedding']][[1,1]], -1.80, tol=1e-1)
  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['umap_model']][['n_neighbors']][[1]], 15, tol=1e1)

  expect_error(reduce_dimension(cds, reduction_method = "DDRTree"),
               "reduction_method must be one of 'UMAP', 'PCA', 'tSNE', 'LSI', 'Aligned'")
})

test_that("reduce_dimension clears old graphs", {
  skip_not_travis()
  cds <- preprocess_cds(cds)
  cds <- reduce_dimension(cds)
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  cds <- reduce_dimension(cds)
  expect_null(cds@principal_graph_aux[["UMAP"]])
  expect_error(partitions(cds), "No partitions calculated for reduction_method = UMAP. Please first run cluster_cells with reduction_method = UMAP.")
  expect_error(clusters(cds), "No clusters calculated for reduction_method = UMAP. Please first run cluster_cells with reduction_method = UMAP.")
  expect_null(cds@principal_graph[["UMAP"]])

  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  cds <- reduce_dimension(cds, reduction_method = "tSNE")
  expect_is(clusters(cds), "factor")
  testthat::expect_is(cds@principal_graph[["UMAP"]], "igraph")
})
