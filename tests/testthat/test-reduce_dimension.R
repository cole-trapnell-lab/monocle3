context("test-reduce_dimension")

cds <- load_a549()
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = 20)


test_that("reduce_dimension runs", {
  cds <- reduce_dimension(cds)
  expect_equal(nrow(reducedDims(cds)$UMAP), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$UMAP), 2)
  expect_equal(as.numeric(reducedDims(cds)$UMAP[1,1]), 1.018337,
               tolerance = 1e-4)

  cds <- reduce_dimension(cds, max_components = 3, reduction_method = "UMAP")
  expect_equal(nrow(reducedDims(cds)$UMAP), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$UMAP), 3)
  expect_equal(as.numeric(reducedDims(cds)$UMAP[1,1]), 0.140256,
               tolerance = 1e-4)

  cds <- reduce_dimension(cds, reduction_method = "tSNE")
  expect_equal(nrow(reducedDims(cds)$tSNE), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$tSNE), 2)
  expect_equal(as.numeric(reducedDims(cds)$tSNE[1,1]), 2.480872,
               tolerance = 1e-4)

  cds <- reduce_dimension(cds,  max_components = 3, reduction_method = "tSNE")
  expect_equal(nrow(reducedDims(cds)$tSNE), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$tSNE), 3)
  expect_equal(as.numeric(reducedDims(cds)$tSNE[1,1]), 2.195226,
               tolerance = 1e-4)

  cds <- reduce_dimension(cds, reduction_method = "PCA")
  expect_equal(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equal(as.numeric(reducedDims(cds)$PCA["E11_A01_RT_467","PC1"]),
               2.4207391, tolerance = 1e-4)

  cds <- preprocess_cds(cds, num_dim = 20, method = "LSI")
  cds <- reduce_dimension(cds, preprocess_method = "LSI")
  expect_equal(nrow(reducedDims(cds)$UMAP), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$UMAP), 2)
  expect_equal(as.numeric(reducedDims(cds)$UMAP[1,1]), 1.359459,
               tolerance = 1e-4)

  cds <- reduce_dimension(cds, reduction_method = "tSNE", preprocess_method = "LSI")
  expect_equal(nrow(reducedDims(cds)$tSNE), nrow(colData(cds)))
  expect_equal(ncol(reducedDims(cds)$tSNE), 2)
  expect_equal(as.numeric(reducedDims(cds)$tSNE[1,1]), 3.996713,
               tolerance = 1e-4)

  expect_error(reduce_dimension(cds, reduction_method = "DDRTree"),
               "reduction_method must be one of 'UMAP', 'PCA' or 'tSNE'")
})
