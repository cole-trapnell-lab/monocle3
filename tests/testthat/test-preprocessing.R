context("test-preprocessing")

test_that("preprocessing stays the same", {
  cds <- load_a549()
  cds <- estimate_size_factors(cds)

  # norm_method = log
  cds <- preprocess_cds(cds, method = "PCA", num_dim = 20)
  expect_equivalent(ncol(reducedDims(cds)$PCA), 20)
  expect_equivalent(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equivalent(reducedDims(cds)$PCA[1,1], 2.4207391, tol = 1e-5)

  cds <- preprocess_cds(cds, method = "LSI", num_dim = 20)
  expect_equivalent(ncol(reducedDims(cds)$LSI), 20)
  expect_equivalent(nrow(reducedDims(cds)$LSI), nrow(colData(cds)))
  expect_equivalent(reducedDims(cds)$LSI[1,1], 13.73796, tol = 1e-5)

  cds <- preprocess_cds(cds, method = "PCA", norm_method = "size_only",
                        num_dim = 20)
  expect_equivalent(ncol(reducedDims(cds)$PCA), 20)
  expect_equivalent(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equivalent(reducedDims(cds)$PCA[1,1], 2.222207, tol = 1e-5)

  cds <- preprocess_cds(cds, method = "LSI", norm_method = "size_only",
                        num_dim = 20)
  expect_equivalent(ncol(reducedDims(cds)$LSI), 20)
  expect_equivalent(nrow(reducedDims(cds)$LSI), nrow(colData(cds)))
  expect_equivalent(reducedDims(cds)$LSI[1,1], 13.49733, tol = 1e-5)

  cds <- preprocess_cds(cds, method = "PCA", norm_method = "none",
                        num_dim = 20)
  expect_equivalent(ncol(reducedDims(cds)$PCA), 20)
  expect_equivalent(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equivalent(reducedDims(cds)$PCA[1,1], -2.42836, tol = 1e-5)

  cds <- preprocess_cds(cds, method = "LSI", norm_method = "none",
                        num_dim = 20)
  expect_equivalent(ncol(reducedDims(cds)$LSI), 20)
  expect_equivalent(nrow(reducedDims(cds)$LSI), nrow(colData(cds)))
  expect_equivalent(reducedDims(cds)$LSI[1,1], 13.49733, tol = 1e-5)

  # BPCells counts matrix.
  cds <- load_a549(matrix_control=list(matrix_class='BPCells'))
  cds <- estimate_size_factors(cds)

  expect_true(is(counts(cds), 'IterableMatrix'))

  # norm_method = log
  cds <- preprocess_cds(cds, method = "PCA", num_dim = 20)
  expect_equivalent(ncol(reducedDims(cds)$PCA), 20)
  expect_equivalent(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equivalent(reducedDims(cds)$PCA[1,1], 2.4207391, tol = 1e-5)

  cds <- preprocess_cds(cds, method = "LSI", num_dim = 20)
  expect_equivalent(ncol(reducedDims(cds)$LSI), 20)
  expect_equivalent(nrow(reducedDims(cds)$LSI), nrow(colData(cds)))
  expect_equivalent(reducedDims(cds)$LSI[1,1], 13.73796, tol = 1e-5)

  cds <- preprocess_cds(cds, method = "PCA", norm_method = "size_only",
                        num_dim = 20)
  expect_equivalent(ncol(reducedDims(cds)$PCA), 20)
  expect_equivalent(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equivalent(reducedDims(cds)$PCA[1,1], 2.222207, tol = 1e-5)

  cds <- preprocess_cds(cds, method = "LSI", norm_method = "size_only",
                        num_dim = 20)
  expect_equivalent(ncol(reducedDims(cds)$LSI), 20)
  expect_equivalent(nrow(reducedDims(cds)$LSI), nrow(colData(cds)))
  expect_equivalent(reducedDims(cds)$LSI[1,1], 13.49733, tol = 1e-5)

  cds <- preprocess_cds(cds, method = "PCA", norm_method = "none",
                        num_dim = 20)
  expect_equivalent(ncol(reducedDims(cds)$PCA), 20)
  expect_equivalent(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equivalent(reducedDims(cds)$PCA[1,1], -2.42836, tol = 1e-5)

  cds <- preprocess_cds(cds, method = "LSI", norm_method = "none",
                        num_dim = 20)
  expect_equivalent(ncol(reducedDims(cds)$LSI), 20)
  expect_equivalent(nrow(reducedDims(cds)$LSI), nrow(colData(cds)))
  expect_equivalent(reducedDims(cds)$LSI[1,1], 13.49733, tol = 1e-5)


  # PCA model
  cds <- load_a549()
  cds <- preprocess_cds(cds, method = "PCA", num_dim = 20)
  expect_equal(cds@reduce_dim_aux[['PCA']][['model']][['num_dim']], 20)
  expect_equal(cds@reduce_dim_aux[['PCA']][['model']][['norm_method']], 'log')
  expect_equal(cds@reduce_dim_aux[['PCA']][['model']][['svd_v']][[1,1]], 0.1362, tol=1e-3)
  expect_equal(cds@reduce_dim_aux[['PCA']][['model']][['svd_sdev']][[1]], 2.67, tol=1e-2)
  expect_equal(cds@reduce_dim_aux[['PCA']][['model']][['svd_center']][[1]], 2.22475, tol=1e-4)
  expect_equal(cds@reduce_dim_aux[['PCA']][['model']][['svd_scale']][[1]], 0.9393, tol=1e-3)
  expect_equal(cds@reduce_dim_aux[['PCA']][['model']][['prop_var_expl']][[1]], 0.1493, tol=1e-3)

  # LSI model
  cds <- load_a549()
  cds <- preprocess_cds(cds, method = "LSI", num_dim = 20)
  expect_equal(cds@reduce_dim_aux[['LSI']][['model']][['num_dim']], 20)
  expect_equal(cds@reduce_dim_aux[['LSI']][['model']][['norm_method']], 'log')
  expect_equal(cds@reduce_dim_aux[['LSI']][['model']][['col_sums']][[1]], 18.7, tol=1e-1)
  expect_equal(cds@reduce_dim_aux[['LSI']][['model']][['row_sums']][[1]], 467, tol=1)
  expect_equal(cds@reduce_dim_aux[['LSI']][['model']][['svd_v']][[1,1]], 0.16318, tol=1e-3)
  expect_equal(cds@reduce_dim_aux[['LSI']][['model']][['svd_sdev']][[1]], 32.21, tol=1e-2)

  # non-standard options
  cds <- preprocess_cds(cds, method = "PCA", scaling=FALSE,
                        verbose = TRUE, norm_method = "size_only",
                        pseudo_count = 1.4, num_dim = 20)
  expect_equal(ncol(reducedDims(cds)$PCA), 20)
  expect_equal(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equal(reducedDims(cds)$PCA[1,1], -24.30544, tol = 1e-5)

  cds <- preprocess_cds(cds, method = "PCA", norm_method = "log",
                        pseudo_count = 1, num_dim = 20)
  expect_equal(ncol(reducedDims(cds)$PCA), 20)
  expect_equal(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equal(reducedDims(cds)$PCA[1,1], 2.42, tol = 1e-2)

  # with reduce dim
  cds <- preprocess_cds(cds, method = "PCA", num_dim = 20)
  expect_equal(ncol(reducedDims(cds)$PCA), 20)
  expect_equal(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equal(reducedDims(cds)$PCA[2,1], 1.982232, tol = 1e-5)

  # with use_genes
  cds <- preprocess_cds(cds, method = "PCA", num_dim = 20,
                        use_genes = c(row.names(rowData(cds))[1:100]))
  expect_equal(ncol(reducedDims(cds)$PCA), 20)
  expect_equal(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equal(reducedDims(cds)$PCA[2,1], -0.5347819, tol = 1e-5)

  # nearest neighbor index
  cds <- preprocess_cds(cds, method='PCA', num_dim=20, build_nn_index=TRUE)
  expect_equal(cds@reduce_dim_aux[['PCA']][['nn_index']][['annoy']][['nn_index']][['metric']], 'cosine')
  cds <- preprocess_cds(cds, method='LSI', num_dim=20, build_nn_index=TRUE)
  expect_equal(cds@reduce_dim_aux[['LSI']][['nn_index']][['annoy']][['nn_index']][['metric']], 'cosine')

})


