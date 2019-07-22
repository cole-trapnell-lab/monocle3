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


  # with use_genes
  cds <- preprocess_cds(cds, method = "PCA", num_dim = 20,
                        use_genes = c(row.names(rowData(cds))[1:100]))
  expect_equal(ncol(reducedDims(cds)$PCA), 20)
  expect_equal(nrow(reducedDims(cds)$PCA), nrow(colData(cds)))
  expect_equal(reducedDims(cds)$PCA[2,1], -0.5347819, tol = 1e-5)

})


make_fake_batched_cds <- function(){

  set.seed(42)

  cell_type_A_prob = c(10, 0, 0, 25, 0, 0, 60, 0, 0, 25)
  cell_type_B_prob = c(0, 10, 0, 25, 12, 0, 0, 0, 50, 0)
  cell_type_C_prob = c(0, 0, 10, 25, 0, 0, 0, 10, 0, 0)

  batch_1_cell_type_A = rmultinom(500, 300, cell_type_A_prob)
  batch_1_coldata_type_A = data.frame(row.names = paste("batch_1_A", seq(1, ncol(batch_1_cell_type_A)), sep="_"),
                                      cell_type = rep("A", ncol(batch_1_cell_type_A)),
                                      batch="batch_1")

  batch_1_cell_type_B = rmultinom(1500, 250, cell_type_B_prob)
  batch_1_coldata_type_B = data.frame(row.names = paste("batch_1_B", seq(1, ncol(batch_1_cell_type_B)), sep="_"),
                                      cell_type = rep("B", ncol(batch_1_cell_type_B)),
                                      batch="batch_1")

  batch_1_cell_type_C = rmultinom(250, 500, cell_type_C_prob)
  batch_1_coldata_type_C = data.frame(row.names = paste("batch_1_C", seq(1, ncol(batch_1_cell_type_C)), sep="_"),
                                      cell_type = rep("C", ncol(batch_1_cell_type_C)),
                                      batch="batch_1")


  batch_vec = c(-3, 10, 10, 0, 0, 0, -30, 0, 0, -20)
  batch_2_cell_type_A = rmultinom(500, 300, cell_type_A_prob + batch_vec)
  batch_2_coldata_type_A = data.frame(row.names = paste("batch_2_A", seq(1, ncol(batch_2_cell_type_A)), sep="_"),
                                      cell_type = rep("A", ncol(batch_2_cell_type_A)),
                                      batch="batch_2")

  batch_2_cell_type_B = rmultinom(200, 250, cell_type_B_prob)
  batch_2_coldata_type_B = data.frame(row.names = paste("batch_2_B", seq(1, ncol(batch_2_cell_type_B)), sep="_"),
                                      cell_type = rep("B", ncol(batch_2_cell_type_B)),
                                      batch="batch_2")

  batch_2_cell_type_C = rmultinom(500, 500, cell_type_C_prob)
  batch_2_coldata_type_C = data.frame(row.names = paste("batch_2_C", seq(1, ncol(batch_2_cell_type_C)), sep="_"),
                                      cell_type = rep("C", ncol(batch_2_cell_type_C)),
                                      batch="batch_2")


  expr_data = cbind(batch_1_cell_type_A,
                    batch_1_cell_type_B,
                    batch_1_cell_type_C,
                    batch_2_cell_type_A,
                    batch_2_cell_type_B,
                    batch_2_cell_type_C)
  col_data = rbind(batch_1_coldata_type_A,
                   batch_1_coldata_type_B,
                   batch_1_coldata_type_C,
                   batch_2_coldata_type_A,
                   batch_2_coldata_type_B,
                   batch_2_coldata_type_C)

  batched_cds = suppressWarnings(new_cell_data_set(expr_data, col_data))
  return (batched_cds)
}

test_that("Alignment works on synthetic data", {
  batched_cds = make_fake_batched_cds()
  batched_cds = preprocess_cds(batched_cds, num_dim=2)
  batched_cds = cluster_cells(batched_cds, k=30, reduction_method="PCA", resolution=1e-6)
  #plot_cells(batched_cds, reduction_method="PCA", color_cells_by="cluster")

  expect_equal(length(unique(clusters(batched_cds, reduction_method="PCA"))), 4)

  #batched_cds = preprocess_cds(batched_cds, num_dim=2, residual_model_formula_str="~cell_type")

  #plot_cells(batched_cds, reduction_method="PCA", color_cells_by="cell_type")
  #plot_cells(batched_cds, reduction_method="PCA", color_cells_by="batch")

  batched_cds = preprocess_cds(batched_cds, num_dim=2, residual_model_formula_str="~batch")
  batched_cds = cluster_cells(batched_cds, k=30, reduction_method="PCA", resolution=1e-6)

  expect_equal(length(unique(clusters(batched_cds, reduction_method="PCA"))), 6)

  #plot_cells(batched_cds, reduction_method="PCA", color_cells_by="cluster")

  #plot_cells(batched_cds, reduction_method="PCA", color_cells_by="cell_type")
  #plot_cells(batched_cds, reduction_method="PCA", color_cells_by="batch")
})
