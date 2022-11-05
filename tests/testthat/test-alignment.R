context("test-alignment")


skip_not_travis <- function ()
{
  if (identical(Sys.getenv("TRAVIS"), "true")) {
    return(invisible(TRUE))
  }
  skip("Not on Travis")
}


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


test_that('Nearest neighbor index', {
  batched_cds = make_fake_batched_cds()
  batched_cds = preprocess_cds(batched_cds, num_dim=3)
  batched_cds = align_cds(batched_cds, residual_model_formula_str="~batch", build_nn_index=TRUE)
  expect_equal(batched_cds@reduce_dim_aux[['Aligned']][['nn_index']][['annoy']][['nn_index']][['metric']], 'cosine')
})


#### NOT TRAVIS ####

test_that("Alignment works on synthetic data", {
  skip_on_travis()
  batched_cds = make_fake_batched_cds()
  batched_cds = preprocess_cds(batched_cds, num_dim=3)
  batched_cds = cluster_cells(batched_cds, k=10, reduction_method="PCA", resolution=1e-3)
  plot_cells(batched_cds, reduction_method="PCA", color_cells_by="partition")

  expect_equal(length(unique(partitions(batched_cds, reduction_method="PCA"))), 5)

  #batched_cds = preprocess_cds(batched_cds, num_dim=2, residual_model_formula_str="~cell_type")

  #plot_cells(batched_cds, reduction_method="PCA", color_cells_by="cell_type")
  #plot_cells(batched_cds, reduction_method="PCA", color_cells_by="batch")

  batched_cds = align_cds(batched_cds, residual_model_formula_str="~batch")
  batched_cds = cluster_cells(batched_cds, k=10, reduction_method="Aligned", resolution=1e-3)
  plot_cells(batched_cds, reduction_method="Aligned", color_cells_by="batch")

  expect_equal(batched_cds@reduce_dim_aux[['Aligned']][['model']][['beta']][[1]], 2.071, tol=1e-2)
  expect_equal(batched_cds@reduce_dim_aux[['Aligned']][['model']][['alignment_k']], 20, tol=1e1)
  expect_equal(length(unique(partitions(batched_cds, reduction_method="Aligned"))), 6)

  batched_cds = preprocess_cds(batched_cds, num_dim=3)
  batched_cds = suppressWarnings(align_cds(batched_cds, alignment_group="batch"))
  batched_cds = cluster_cells(batched_cds, k=10, reduction_method="Aligned", resolution=1e-3)
  plot_cells(batched_cds, reduction_method="Aligned", color_cells_by="batch")

  expect_equal(length(unique(partitions(batched_cds, reduction_method="Aligned"))), 10)


    #plot_cells(batched_cds, reduction_method="PCA", color_cells_by="cluster")

  #plot_cells(batched_cds, reduction_method="PCA", color_cells_by="cell_type")
  #plot_cells(batched_cds, reduction_method="PCA", color_cells_by="batch")
})


#### TRAVIS ####

test_that("Alignment works on synthetic data", {
  skip_not_travis()
  batched_cds = make_fake_batched_cds()
  batched_cds = preprocess_cds(batched_cds, num_dim=3)
  batched_cds = cluster_cells(batched_cds, k=10, reduction_method="PCA", resolution=1e-3)
  plot_cells(batched_cds, reduction_method="PCA", color_cells_by="partition")

  expect_equal(length(unique(partitions(batched_cds, reduction_method="PCA"))), 4)

  #batched_cds = preprocess_cds(batched_cds, num_dim=2, residual_model_formula_str="~cell_type")

  #plot_cells(batched_cds, reduction_method="PCA", color_cells_by="cell_type")
  #plot_cells(batched_cds, reduction_method="PCA", color_cells_by="batch")

  batched_cds = align_cds(batched_cds, residual_model_formula_str="~batch")
  batched_cds = cluster_cells(batched_cds, k=10, reduction_method="Aligned", resolution=1e-3)
  plot_cells(batched_cds, reduction_method="Aligned", color_cells_by="batch")

  expect_equal(batched_cds@reduce_dim_aux[['Aligned']][['model']][['beta']][[1]], 2.071, tol=1e-2)
  expect_equal(batched_cds@reduce_dim_aux[['Aligned']][['model']][['alignment_k']], 20, tol=1e1)
  expect_equal(length(unique(partitions(batched_cds, reduction_method="Aligned"))), 6)

  batched_cds = preprocess_cds(batched_cds, num_dim=3)
  batched_cds = suppressWarnings(align_cds(batched_cds, alignment_group="batch"))
  batched_cds = cluster_cells(batched_cds, k=10, reduction_method="Aligned", resolution=1e-3)
  plot_cells(batched_cds, reduction_method="Aligned", color_cells_by="batch")

  expect_equal(length(unique(partitions(batched_cds, reduction_method="Aligned"))), 9)


    #plot_cells(batched_cds, reduction_method="PCA", color_cells_by="cluster")

  #plot_cells(batched_cds, reduction_method="PCA", color_cells_by="cell_type")
  #plot_cells(batched_cds, reduction_method="PCA", color_cells_by="batch")
})


