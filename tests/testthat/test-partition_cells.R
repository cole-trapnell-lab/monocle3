context("test-partition_cells")

cds <- load_a549()

test_that("test partition_cells error messages work", {
  expect_error(cds <- partition_cells(cds),
               "No dimensionality reduction for UMAP calculated. Please run reduce_dimensions with reduction_method = UMAP before running partition_cells.")
  cds <- preprocess_cds(cds)
  expect_error(cds <- partition_cells(cds),
               "No dimensionality reduction for UMAP calculated. Please run reduce_dimensions with reduction_method = UMAP before running partition_cells.")
  cds <- reduce_dimension(cds)
  expect_error(cds <- partition_cells(cds, reduction_method = "tSNE"),
               "No dimensionality reduction for tSNE calculated. Please run reduce_dimensions with reduction_method = tSNE before running partition_cells.")
})

cds <- preprocess_cds(cds)
cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1)
cds <- reduce_dimension(cds, reduction_method = "tSNE")
cds <- cluster_cells(cds)

test_that("partition_cells works", {
  cds <- partition_cells(cds)
  expect_is(cds@partitions[["UMAP"]], "list")
  expect_equal(length(cds@partitions[["UMAP"]]), 2)
  expect_equal(length(cds@partitions[["UMAP"]]$louvain_res$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@partitions[["UMAP"]]$louvain_res$optim_res$membership[1], 9)

  # non-standard opts
  cds <- partition_cells(cds, k=22, weight = T, louvain_iter = 2,
                         louvain_qval = .1)
  expect_is(cds@partitions[["UMAP"]], "list")
  expect_equal(length(cds@partitions[["UMAP"]]), 2)
  expect_equal(length(cds@partitions[["UMAP"]]$louvain_res$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@partitions[["UMAP"]]$louvain_res$optim_res$membership[1], 7)

  cds <- partition_cells(cds, reduction_method = "tSNE")
  expect_is(cds@partitions[["tSNE"]], "list")
  expect_equal(length(cds@partitions[["tSNE"]]), 2)
  expect_equal(length(cds@partitions[["tSNE"]]$louvain_res$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@partitions[["tSNE"]]$louvain_res$optim_res$membership[1], 9)

  # non-standard opts
  cds <- partition_cells(cds, reduction_method = "tSNE", k=22, weight = T,
                         louvain_iter = 2, louvain_qval = .1)
  expect_is(cds@partitions[["tSNE"]], "list")
  expect_equal(length(cds@partitions[["tSNE"]]), 2)
  expect_equal(length(cds@partitions[["tSNE"]]$louvain_res$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@partitions[["tSNE"]]$louvain_res$optim_res$membership[1], 9)

  cds <- partition_cells(cds, reduction_method = "PCA")
  expect_is(cds@partitions[["PCA"]], "list")
  expect_equal(length(cds@partitions[["PCA"]]), 2)
  expect_equal(length(cds@partitions[["PCA"]]$louvain_res$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@partitions[["PCA"]]$louvain_res$optim_res$membership[1], 2)

  # non-standard opts
  cds <- partition_cells(cds, reduction_method = "PCA", k=22, weight = T,
                         louvain_iter = 2, louvain_qval = .1)
  expect_is(cds@partitions[["PCA"]], "list")
  expect_equal(length(cds@partitions[["PCA"]]), 2)
  expect_equal(length(cds@partitions[["PCA"]]$louvain_res$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@partitions[["PCA"]]$louvain_res$optim_res$membership[1], 4)
})

