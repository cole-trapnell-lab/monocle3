context("test-preprocessing")

cds <- load_a549()
cds <- estimate_size_factors(cds)

test_that("preprocessing stays the same", {
  # norm_method = log
  cds <- preprocess_cds(cds, method = "PCA", num_dim = 20)
  expect_equivalent(ncol(reducedDims(cds)$PCA), 20)
  expect_equivalent(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equivalent(reducedDims(cds)$PCA[1,1], 2.4207391, tol = 1e-5)

  cds <- preprocess_cds(cds, method = "LSI", num_dim = 20)
  expect_equivalent(ncol(reducedDims(cds)$LSI), 20)
  expect_equivalent(nrow(reducedDims(cds)$LSI), nrow(colData(cds)))
  expect_equivalent(reducedDims(cds)$LSI[1,1], 11.112528, tol = 1e-5)

  cds <- preprocess_cds(cds, method = "PCA", norm_method = "size_only", num_dim = 20)
  expect_equivalent(ncol(reducedDims(cds)$PCA), 20)
  expect_equivalent(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equivalent(reducedDims(cds)$PCA[1,1], 2.222207, tol = 1e-5)

  cds <- preprocess_cds(cds, method = "LSI", norm_method = "size_only", num_dim = 20)
  expect_equivalent(ncol(reducedDims(cds)$LSI), 20)
  expect_equivalent(nrow(reducedDims(cds)$LSI), nrow(colData(cds)))
  expect_equivalent(reducedDims(cds)$LSI[1,1], 9.283694, tol = 1e-5)

  # non-standard options
  cds <- preprocess_cds(cds, method = "PCA", scaling=FALSE,
                        verbose = TRUE, norm_method = "size_only",
                        pseudo_count = 1.4, num_dim = 20)
  expect_equal(ncol(reducedDims(cds)$PCA), 20)
  expect_equal(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equal(reducedDims(cds)$PCA[1,1], -24.30544, tol = 1e-5)

  # with reduce dim
  cds <- preprocess_cds(cds, method = "PCA", num_dim = 20,
                        residual_model_formula_str = "~PCR_plate")
  expect_equal(ncol(reducedDims(cds)$PCA), 20)
  expect_equal(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equal(reducedDims(cds)$PCA[2,1], 2.274213, tol = 1e-5)


})



