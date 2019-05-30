context("test-cluster_cells")

cds <- load_a549()

test_that("test cluster_cells error messages work", {
  expect_error(cds <- cluster_cells(cds),
               "No dimensionality reduction for UMAP calculated. Please run reduce_dimensions with reduction_method = UMAP before running cluster_cells")
  cds <- preprocess_cds(cds)
  expect_error(cds <- cluster_cells(cds),
               "No dimensionality reduction for UMAP calculated. Please run reduce_dimensions with reduction_method = UMAP before running cluster_cells")
  cds <- reduce_dimension(cds)
  expect_error(cds <- cluster_cells(cds, reduction_method = "tSNE"),
               "No dimensionality reduction for tSNE calculated. Please run reduce_dimensions with reduction_method = tSNE before running cluster_cells")
})

cds <- preprocess_cds(cds)
cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1)
cds <- reduce_dimension(cds, reduction_method = "tSNE")
cds <- cluster_cells(cds)

test_that("cluster_cells works", {
  cds <- cluster_cells(cds)
  expect_is(cds@clusters[["UMAP"]], "list")
  expect_equal(length(cds@clusters[["UMAP"]]), 3)
  expect_equal(length(cds@clusters[["UMAP"]]$louvain_res$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["UMAP"]]$louvain_res$optim_res$membership[1], 9)

  # non-standard opts
  cds <- cluster_cells(cds, k=22, weight = T, louvain_iter = 2,
                         partition_qval = .1)
  expect_is(cds@clusters[["UMAP"]], "list")
  expect_equal(length(cds@clusters[["UMAP"]]), 3)
  expect_equal(length(cds@clusters[["UMAP"]]$louvain_res$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["UMAP"]]$louvain_res$optim_res$membership[1], 7)

  cds <- cluster_cells(cds, reduction_method = "tSNE")
  expect_is(cds@clusters[["tSNE"]], "list")
  expect_equal(length(cds@clusters[["tSNE"]]), 3)
  expect_equal(length(cds@clusters[["tSNE"]]$louvain_res$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["tSNE"]]$louvain_res$optim_res$membership[1], 9)

  # non-standard opts
  cds <- cluster_cells(cds, reduction_method = "tSNE", k=22, weight = T,
                         louvain_iter = 2, partition_qval = .1)
  expect_is(cds@clusters[["tSNE"]], "list")
  expect_equal(length(cds@clusters[["tSNE"]]), 3)
  expect_equal(length(cds@clusters[["tSNE"]]$louvain_res$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["tSNE"]]$louvain_res$optim_res$membership[1], 9)

  cds <- cluster_cells(cds, reduction_method = "PCA")
  expect_is(cds@clusters[["PCA"]], "list")
  expect_equal(length(cds@clusters[["PCA"]]), 3)
  expect_equal(length(cds@clusters[["PCA"]]$louvain_res$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["PCA"]]$louvain_res$optim_res$membership[1], 2)

  # non-standard opts
  cds <- cluster_cells(cds, reduction_method = "PCA", k=22, weight = T,
                         louvain_iter = 2, partition_qval = .1)
  expect_is(cds@clusters[["PCA"]], "list")
  expect_equal(length(cds@clusters[["PCA"]]), 3)
  expect_equal(length(cds@clusters[["PCA"]]$louvain_res$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["PCA"]]$louvain_res$optim_res$membership[1], 4)

  cds <- cluster_cells(cds, reduction_method = "PCA", k=22, weight = T,
                       louvain_iter = 2, partition_qval = .1, resolution = .1)
  expect_is(cds@clusters[["PCA"]], "list")
  expect_equal(length(cds@clusters[["PCA"]]), 3)
  expect_equal(length(cds@clusters[["PCA"]]$louvain_res$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(as.numeric(cds@clusters[["PCA"]]$louvain_res$optim_res$membership[1]), 1)
  expect_equal(length(unique(clusters(cds, reduction_method = "PCA"))), 5)

  cds <- cluster_cells(cds, reduction_method = "PCA", k=22, weight = T,
                       louvain_iter = 2, partition_qval = .1, resolution = .01)
  expect_equal(length(unique(clusters(cds, reduction_method = "PCA"))), 1)
  cds <- cluster_cells(cds, reduction_method = "PCA", k=22, weight = T,
                       louvain_iter = 2, partition_qval = .1, resolution = 20)
  expect_equal(length(unique(clusters(cds, reduction_method = "PCA"))), 500)
})

cds <- load_a549()
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = 20)
cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1)

test_that("cluster_cells works", {
  cds <- cluster_cells(cds)
  expect_equal(sum(clusters(cds) == 9), 85)
  expect_equal(sum(clusters(cds) == 1), 33)
  expect_equal(as.character(clusters(cds)[1]), "9")
})
