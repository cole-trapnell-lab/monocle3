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

test_that("reduce_dimension runs", {
  skip_on_travis()
  set.seed(100)
  expect_error(cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1, approx_pow=TRUE),
               "Data has not been preprocessed with chosen method: PCA Please run preprocess_cds with method = PCA before running reduce_dimension.")
  cds <- preprocess_cds(cds, num_dim = 20)

  cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1, approx_pow=TRUE)
  expect_equal(nrow(reducedDims(cds)$UMAP), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$UMAP), 2)
  expect_equal(as.numeric(reducedDims(cds)$UMAP[1,1]), -0.9980026,
               tolerance = 1e-4)

  cds <- reduce_dimension(cds, max_components = 3, umap.fast_sgd=FALSE, cores=1, reduction_method = "UMAP", approx_pow=TRUE)
  expect_equal(nrow(reducedDims(cds)$UMAP), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$UMAP), 3)
  expect_equal(as.numeric(reducedDims(cds)$UMAP[1,1]), 2.014246,
               tolerance = 1e-3)

  cds <- reduce_dimension(cds, reduction_method = "tSNE")
  expect_equal(nrow(reducedDims(cds)$tSNE), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$tSNE), 2)
  expect_equal(as.numeric(reducedDims(cds)$tSNE[1,1]), 3.196225,
               tolerance = 1e-3)

  cds <- reduce_dimension(cds,  max_components = 3, reduction_method = "tSNE")
  expect_equal(nrow(reducedDims(cds)$tSNE), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$tSNE), 3)
  expect_equal(as.numeric(reducedDims(cds)$tSNE[1,1]), -1.803658,
               tolerance = 1e-3)

  cds <- reduce_dimension(cds, reduction_method = "PCA")
  expect_equal(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equal(as.numeric(reducedDims(cds)$PCA["E11_A01_RT_467","PC1"]),
               2.4207391, tolerance = 1e-3)

  cds <- preprocess_cds(cds, num_dim = 20, method = "LSI")
  cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1, preprocess_method = "LSI")
  expect_equal(nrow(reducedDims(cds)$UMAP), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$UMAP), 2)
  expect_equal(as.numeric(reducedDims(cds)$UMAP[1,1]), -0.8853563,
               tolerance = 1e-3)

  cds <- reduce_dimension(cds, reduction_method = "tSNE", preprocess_method = "LSI")
  expect_equal(nrow(reducedDims(cds)$tSNE), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$tSNE), 2)
  expect_equal(as.numeric(reducedDims(cds)$tSNE[1,1]), 6.191503,
               tolerance = 1e-3)

  expect_error(reduce_dimension(cds, reduction_method = "DDRTree"),
               "reduction_method must be one of 'UMAP', 'PCA', 'tSNE', 'LSI', 'Aligned'")
})

test_that("reduce_dimension clears old graphs", {
  skip_on_travis()
  set.seed(100)
  cds <- preprocess_cds(cds)
  cds <- reduce_dimension(cds, approx_pow=TRUE)
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  cds <- reduce_dimension(cds, approx_pow=TRUE)
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

set.seed(100)
cds <- load_a549()
cds <- estimate_size_factors(cds)

test_that("reduce_dimension runs", {
#  skip_not_travis()
  set.seed(100)
  expect_error(cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1, approx_pow=TRUE),
               "Data has not been preprocessed with chosen method: PCA Please run preprocess_cds with method = PCA before running reduce_dimension.")
  cds <- preprocess_cds(cds, num_dim = 20)

  cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1, approx_pow=TRUE)
  expect_equal(nrow(reducedDims(cds)$UMAP), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$UMAP), 2)
  expect_equal(as.numeric(reducedDims(cds)$UMAP[1,1]), -0.9980026,
               tolerance = 1e-4)

  cds <- reduce_dimension(cds, max_components = 3, umap.fast_sgd=FALSE, cores=1, reduction_method = "UMAP", approx_pow=TRUE)
  expect_equal(nrow(reducedDims(cds)$UMAP), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$UMAP), 3)
  expect_equal(as.numeric(reducedDims(cds)$UMAP[1,1]), 2.014246,
               tolerance = 1e-3)

  cds <- reduce_dimension(cds, reduction_method = "tSNE")
  expect_equal(nrow(reducedDims(cds)$tSNE), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$tSNE), 2)
  expect_equal(as.numeric(reducedDims(cds)$tSNE[1,1]), 3.196225,
               tolerance = 1e-3)

  cds <- reduce_dimension(cds,  max_components = 3, reduction_method = "tSNE")
  expect_equal(nrow(reducedDims(cds)$tSNE), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$tSNE), 3)
  expect_equal(as.numeric(reducedDims(cds)$tSNE[1,1]),  -1.803658,
               tolerance = 1e-3)

  cds <- reduce_dimension(cds, reduction_method = "PCA")
  expect_equal(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equal(as.numeric(reducedDims(cds)$PCA["E11_A01_RT_467","PC1"]),
               2.4207391, tolerance = 1e-4)

  cds <- preprocess_cds(cds, num_dim = 20, method = "LSI")
  cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1, preprocess_method = "LSI")
  expect_equal(nrow(reducedDims(cds)$UMAP), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$UMAP), 2)
  expect_equal(as.numeric(reducedDims(cds)$UMAP[1,1]), -0.8853563,
               tolerance = 1e-4)

  cds <- reduce_dimension(cds, reduction_method = "tSNE", preprocess_method = "LSI")
  expect_equal(nrow(reducedDims(cds)$tSNE), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$tSNE), 2)
  expect_equal(as.numeric(reducedDims(cds)$tSNE[1,1]), 6.191503,
               tolerance = 1e-3)

  expect_error(reduce_dimension(cds, reduction_method = "DDRTree"),
               "reduction_method must be one of 'UMAP', 'PCA', 'tSNE', 'LSI', 'Aligned'")
})

test_that("reduce_dimension clears old graphs", {
#  skip_not_travis()
  set.seed(100)
  cds <- preprocess_cds(cds)
  cds <- reduce_dimension(cds, approx_pow=TRUE)
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  cds <- reduce_dimension(cds, approx_pow=TRUE)
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
