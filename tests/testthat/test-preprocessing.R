context("test-preprocessing")

cds <- load_a549()
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

test_that("preprocessing stays the same", {
  cds <- preprocess_cds(cds, method = "PCA", num_dim = 20)
  expect_equal(ncol(cds@normalized_data_projection), 20)
  expect_equal(nrow(cds@normalized_data_projection), nrow(pData(cds)))
  expect_equal(cds@normalized_data_projection[1,1], 2.4207391)

  #TODO check this actual behaviour is right
  cds <- preprocess_cds(cds, method = "none", num_dim = 20)
  expect_equal(ncol(cds@normalized_data_projection), 179)
  expect_equal(nrow(cds@normalized_data_projection), nrow(pData(cds)))
  expect_equal(cds@normalized_data_projection[1,1], 2.573678, tol = 1e-5)

  #cds <- preprocess_cds(cds, method = "PCA", norm_method = "none", num_dim = 20)
  #expect_equal(ncol(cds@normalized_data_projection), 20)
  #expect_equal(nrow(cds@normalized_data_projection), nrow(pData(cds)))
  #expect_equal(cds@normalized_data_projection[1,1], 2.4207391)

  cds <- preprocess_cds(cds, method = "none", norm_method = "none", num_dim = 20)
  expect_equal(ncol(cds@normalized_data_projection), 198)
  expect_equal(nrow(cds@normalized_data_projection), nrow(pData(cds)))
  expect_equal(cds@normalized_data_projection[1,1], 1)
})
