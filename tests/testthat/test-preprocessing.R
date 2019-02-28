context("test-preprocessing")

cds <- load_a549()
cds <- estimate_size_factors(cds)

test_that("preprocessing stays the same", {
  # norm_method = log
  cds <- preprocess_cds(cds, method = "PCA", num_dim = 20)
  expect_equal(ncol(reducedDims(cds)$normalized_data_projection), 20)
  expect_equal(nrow(reducedDims(cds)$normalized_data_projection), nrow(colData(cds)))
  expect_equal(reducedDims(cds)$normalized_data_projection[1,1], 2.4207391, tol = 1e-5)

  cds <- preprocess_cds(cds, method = "none", num_dim = 20)
  expect_equal(ncol(reducedDims(cds)$normalized_data_projection), 179)
  expect_equal(nrow(reducedDims(cds)$normalized_data_projection), nrow(colData(cds)))
  expect_equal(reducedDims(cds)$normalized_data_projection[1,1], 2.573678, tol = 1e-5)

  cds <- preprocess_cds(cds, method = "tfidf", num_dim = 20)
  expect_equal(ncol(reducedDims(cds)$normalized_data_projection), 20)
  expect_equal(nrow(reducedDims(cds)$normalized_data_projection), nrow(colData(cds)))
  expect_equal(reducedDims(cds)$normalized_data_projection[1,1], 11.112528, tol = 1e-5)

  # norm_method = none
  cds <- preprocess_cds(cds, method = "PCA", norm_method = "none", num_dim = 20)
  expect_equal(ncol(reducedDims(cds)$normalized_data_projection), 20)
  expect_equal(nrow(reducedDims(cds)$normalized_data_projection), nrow(colData(cds)))
  expect_equal(reducedDims(cds)$normalized_data_projection[1,1], 2.222207, tol = 1e-5)

  cds <- preprocess_cds(cds, method = "none", norm_method = "none", num_dim = 20)
  expect_equal(ncol(reducedDims(cds)$normalized_data_projection), 179)
  expect_equal(nrow(reducedDims(cds)$normalized_data_projection), nrow(colData(cds)))
  expect_equal(reducedDims(cds)$normalized_data_projection["E11_A01_RT_467","ENSG00000228253.1"],
               assays(cds)$exprs["ENSG00000228253.1","E11_A01_RT_467"]/colData(cds)["E11_A01_RT_467",]$Size_Factor)

  cds <- preprocess_cds(cds, method = "tfidf", norm_method = "none", num_dim = 20)
  expect_equal(ncol(reducedDims(cds)$normalized_data_projection), 20)
  expect_equal(nrow(reducedDims(cds)$normalized_data_projection), nrow(colData(cds)))
  expect_equal(reducedDims(cds)$normalized_data_projection[1,1], 9.443847, tol = 1e-5)

  # non-standard options
  cds <- preprocess_cds(cds, method = "PCA", scaling=FALSE,
                        verbose = TRUE, norm_method = "none",
                        pseudo_count = 1.4, num_dim = 20)
  expect_equal(ncol(reducedDims(cds)$normalized_data_projection), 20)
  expect_equal(nrow(reducedDims(cds)$normalized_data_projection), nrow(colData(cds)))
  expect_equal(reducedDims(cds)$normalized_data_projection[1,1], -21.48343, tol = 1e-5)

  # with reduce dim
  cds <- preprocess_cds(cds, method = "PCA", num_dim = 20,
                        residual_model_formula_str = "~PCR_plate")
  expect_equal(ncol(reducedDims(cds)$normalized_data_projection), 20)
  expect_equal(nrow(reducedDims(cds)$normalized_data_projection), nrow(colData(cds)))
  expect_equal(reducedDims(cds)$normalized_data_projection[2,1], 2.659146, tol = 1e-5)


})



