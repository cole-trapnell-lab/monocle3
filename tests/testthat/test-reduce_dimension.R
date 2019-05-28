context("test-reduce_dimension")

cds <- load_a549()
cds <- estimate_size_factors(cds)

test_that("reduce_dimension runs", {
  expect_error(cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1),
               "Data has not been preprocessed with chosen method: PCA Please run preprocess_cds with method = PCA before running reduce_dimension.")
  cds <- preprocess_cds(cds, num_dim = 20)

  cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1)
  expect_equal(nrow(reducedDims(cds)$UMAP), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$UMAP), 2)
  expect_equal(as.numeric(reducedDims(cds)$UMAP[1,1]), -1.815985,
               tolerance = 1e-4)

  cds <- reduce_dimension(cds, max_components = 3, umap.fast_sgd=FALSE, cores=1, reduction_method = "UMAP")
  expect_equal(nrow(reducedDims(cds)$UMAP), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$UMAP), 3)
  expect_equal(as.numeric(reducedDims(cds)$UMAP[1,1]), 0.5965699,
               tolerance = 1e-4)

  cds <- reduce_dimension(cds, reduction_method = "tSNE")
  expect_equal(nrow(reducedDims(cds)$tSNE), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$tSNE), 2)
  expect_equal(as.numeric(reducedDims(cds)$tSNE[1,1]), -4.546885,
               tolerance = 1e-4)

  cds <- reduce_dimension(cds,  max_components = 3, reduction_method = "tSNE")
  expect_equal(nrow(reducedDims(cds)$tSNE), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$tSNE), 3)
  expect_equal(as.numeric(reducedDims(cds)$tSNE[1,1]), -2.388678,
               tolerance = 1e-4)

  cds <- reduce_dimension(cds, reduction_method = "PCA")
  expect_equal(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equal(as.numeric(reducedDims(cds)$PCA["E11_A01_RT_467","PC1"]),
               2.4207391, tolerance = 1e-4)

  cds <- preprocess_cds(cds, num_dim = 20, method = "LSI")
  cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1, preprocess_method = "LSI")
  expect_equal(nrow(reducedDims(cds)$UMAP), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$UMAP), 2)
  expect_equal(as.numeric(reducedDims(cds)$UMAP[1,1]), -0.9881196,
               tolerance = 1e-4)

  cds <- reduce_dimension(cds, reduction_method = "tSNE", preprocess_method = "LSI")
  expect_equal(nrow(reducedDims(cds)$tSNE), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$tSNE), 2)
  expect_equal(as.numeric(reducedDims(cds)$tSNE[1,1]), -0.6852947,
               tolerance = 1e-4)

  expect_error(reduce_dimension(cds, reduction_method = "DDRTree"),
               "reduction_method must be one of 'UMAP', 'PCA' or 'tSNE'")
})

test_that("reduce_dimension clears old graphs", {
  cds <- preprocess_cds(cds)
  cds <- reduce_dimension(cds)
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  cds <- reduce_dimension(cds)
  expect_null(cds@preprocess_aux[["UMAP"]])
  expect_error(partitions(cds), "No partitions calculated for reduction_method = UMAP. Please first run cluster_cells with reduction_method = UMAP.")
  expect_error(clusters(cds), "No clusters calculated for reduction_method = UMAP. Please first run cluster_cells with reduction_method = UMAP.")
  expect_null(cds@principal_graph[["UMAP"]])

  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  cds <- reduce_dimension(cds, reduction_method = "tSNE")
  expect_is(clusters(cds), "factor")
  testthat::expect_is(cds@principal_graph[["UMAP"]], "igraph")
})
