context("test-cluster_cells")

skip_not_travis <- function ()
{
  if (identical(Sys.getenv("TRAVIS"), "true")) {
    return(invisible(TRUE))
  }
  skip("Not on Travis")
}


set.seed(200)
cds <- load_a549()

test_that("test cluster_cells error messages work", {
  skip_on_travis()
  expect_error(cds <- cluster_cells(cds),
               paste("No dimensionality reduction for UMAP calculated. Please",
                     "run reduce_dimension with reduction_method = UMAP",
                     "before running cluster_cells"))
  cds <- preprocess_cds(cds)
  expect_error(cds <- cluster_cells(cds),
               paste("No dimensionality reduction for UMAP calculated. Please",
                     "run reduce_dimension with reduction_method = UMAP",
                     "before running cluster_cells"))
  cds <- reduce_dimension(cds)
  expect_error(cds <- cluster_cells(cds, reduction_method = "tSNE"),
               paste("No dimensionality reduction for tSNE calculated. Please",
                     "run reduce_dimension with reduction_method = tSNE",
                     "before running cluster_cells"))
})

cds <- preprocess_cds(cds)
cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1)
cds <- reduce_dimension(cds, reduction_method = "tSNE")
cds <- cluster_cells(cds)

test_that("cluster_cells works", {
  skip_on_travis()
  ### UMAP
  ## leiden
  cds <- cluster_cells(cds, random_seed = 100)
  expect_is(cds@clusters[["UMAP"]], "list")
  expect_equal(length(cds@clusters[["UMAP"]]), 3)
  expect_equal(length(cds@clusters[["UMAP"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["UMAP"]]$cluster_result$optim_res$membership[[1]],
               1)
  expect_equal(length(unique(clusters(cds, reduction_method = "UMAP"))), 1)

  # non-standard opts
  cds <- cluster_cells(cds, k=22, resolution = c(0.75, 0.3), weight = T,
                       num_iter = 1, partition_qval = .1, verbose = TRUE,
                       random_seed = 100)
  expect_is(cds@clusters[["UMAP"]], "list")
  expect_equal(length(cds@clusters[["UMAP"]]), 3)
  expect_equal(length(cds@clusters[["UMAP"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["UMAP"]]$cluster_result$optim_res$membership[[1]],
               9)
  expect_equal(length(unique(clusters(cds, reduction_method = "UMAP"))), 20)

  ## louvain
  cds <- cluster_cells(cds, cluster_method = "louvain", random_seed = 100)
  expect_is(cds@clusters[["UMAP"]], "list")
  expect_equal(length(cds@clusters[["UMAP"]]), 3)
  expect_equal(length(cds@clusters[["UMAP"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["UMAP"]]$cluster_result$optim_res$membership[[1]],
               3)
  expect_equal(length(unique(clusters(cds, reduction_method = "UMAP"))), 10)

  # non-standard opts
  cds <- cluster_cells(cds, cluster_method = "louvain", k=22, weight = T,
                       num_iter = 1, partition_qval = .1, verbose = TRUE,
                       random_seed = 100)
  expect_is(cds@clusters[["UMAP"]], "list")
  expect_equal(length(cds@clusters[["UMAP"]]), 3)
  expect_equal(length(cds@clusters[["UMAP"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["UMAP"]]$cluster_result$optim_res$membership[[1]],
               9)
  expect_equal(length(unique(clusters(cds, reduction_method = "UMAP"))), 11)

  ### tSNE
  ##leiden
  cds <- cluster_cells(cds, reduction_method = "tSNE", random_seed = 100)
  expect_is(cds@clusters[["tSNE"]], "list")
  expect_equal(length(cds@clusters[["tSNE"]]), 3)
  expect_equal(length(cds@clusters[["tSNE"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["tSNE"]]$cluster_result$optim_res$membership[[1]],
               1)
  expect_equal(length(unique(clusters(cds, reduction_method = "tSNE"))), 1)

  # non-standard opts
  cds <- cluster_cells(cds, reduction_method = "tSNE", k=22,
                       resolution = c(0.75, 0.3), weight = T, num_iter = 1,
                       partition_qval = .1, verbose = TRUE, random_seed = 100)
  expect_is(cds@clusters[["tSNE"]], "list")
  expect_equal(length(cds@clusters[["tSNE"]]), 3)
  expect_equal(length(cds@clusters[["tSNE"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["tSNE"]]$cluster_result$optim_res$membership[[1]],
               10)
  expect_equal(length(unique(clusters(cds, reduction_method = "tSNE"))), 20)

  ### PCA

  cds <- cluster_cells(cds, reduction_method = "PCA", random_seed = 100)
  expect_is(cds@clusters[["PCA"]], "list")
  expect_equal(length(cds@clusters[["PCA"]]), 3)
  expect_equal(length(cds@clusters[["PCA"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["PCA"]]$cluster_result$optim_res$membership[[1]],
               1)
  expect_equal(length(unique(clusters(cds, reduction_method = "PCA"))), 1)

  # non-standard opts
  cds <- cluster_cells(cds, reduction_method = "PCA", k=22,
                       resolution = c(0.75, 0.3), weight = T, num_iter = 2,
                       partition_qval = .1, verbose = TRUE, random_seed = 100)
  expect_is(cds@clusters[["PCA"]], "list")
  expect_equal(length(cds@clusters[["PCA"]]), 3)
  expect_equal(length(cds@clusters[["PCA"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["PCA"]]$cluster_result$optim_res$membership[[1]],
               44)
  expect_equal(length(unique(clusters(cds, reduction_method = "PCA"))), 59)

  cds <- cluster_cells(cds, reduction_method = "PCA", k=22, weight = T,
                       num_iter = 2, partition_qval = .1, resolution = .1,
                       random_seed = 100)
  expect_is(cds@clusters[["PCA"]], "list")
  expect_equal(length(cds@clusters[["PCA"]]), 3)
  expect_equal(length(cds@clusters[["PCA"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(as.numeric(cds@clusters[["PCA"]]$cluster_result$optim_res$membership[1]),
               1)
  expect_equal(length(unique(clusters(cds, reduction_method = "PCA"))), 24)

  cds <- cluster_cells(cds, reduction_method = "PCA", k=22, weight = T,
                       num_iter = 2, partition_qval = .1, resolution = .01,
                       random_seed = 100)
  expect_equal(length(unique(clusters(cds, reduction_method = "PCA"))), 1)
  cds <- cluster_cells(cds, reduction_method = "PCA", k=22, weight = T,
                       num_iter = 2, partition_qval = .1, resolution = 20,
                       random_seed = 100)
  expect_equal(length(unique(clusters(cds, reduction_method = "PCA"))), 500)

  cds <- preprocess_cds(cds, method = "LSI")
  cds <- cluster_cells(cds, reduction_method = "LSI", random_seed = 100)
  expect_is(cds@clusters[["LSI"]], "list")
  expect_equal(length(cds@clusters[["LSI"]]), 3)
  expect_equal(length(cds@clusters[["LSI"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["LSI"]]$cluster_result$optim_res$membership[[1]],
               1)
  expect_equal(length(unique(clusters(cds, reduction_method = "LSI"))), 1)

  # non-standard opts
  cds <- cluster_cells(cds, reduction_method = "LSI", k=22,
                       resolution = c(0.75, 0.3), weight = T, num_iter = 2,
                       partition_qval = .1, verbose = TRUE, random_seed = 100)
  expect_is(cds@clusters[["LSI"]], "list")
  expect_equal(length(cds@clusters[["LSI"]]), 3)
  expect_equal(length(cds@clusters[["LSI"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["LSI"]]$cluster_result$optim_res$membership[[1]],
               10)
  expect_equal(length(unique(clusters(cds, reduction_method = "LSI"))), 61)
})

cds <- load_a549()
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = 20)
cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1)

test_that("cluster_cells works", {
  skip_on_travis()
  cds <- cluster_cells(cds, random_seed = 100)
  expect_equal(sum(clusters(cds) == 1), 500)
  expect_equal(as.character(clusters(cds)[1]), "1")
})


### FOR TRAVIS ###


set.seed(200)
cds <- load_a549()

test_that("test cluster_cells error messages work", {
  skip_not_travis()
  expect_error(cds <- cluster_cells(cds),
               paste("No dimensionality reduction for UMAP calculated. Please",
                     "run reduce_dimension with reduction_method = UMAP",
                     "before running cluster_cells"))
  cds <- preprocess_cds(cds)
  expect_error(cds <- cluster_cells(cds),
               paste("No dimensionality reduction for UMAP calculated. Please",
                     "run reduce_dimension with reduction_method = UMAP",
                     "before running cluster_cells"))
  cds <- reduce_dimension(cds)
  expect_error(cds <- cluster_cells(cds, reduction_method = "tSNE"),
               paste("No dimensionality reduction for tSNE calculated. Please",
                     "run reduce_dimension with reduction_method = tSNE",
                     "before running cluster_cells"))
})

cds <- preprocess_cds(cds)
cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1)
cds <- reduce_dimension(cds, reduction_method = "tSNE")
cds <- cluster_cells(cds)

test_that("cluster_cells works", {
  skip_not_travis()
  ### UMAP
  ## leiden
  cds <- cluster_cells(cds, random_seed = 100)
  expect_is(cds@clusters[["UMAP"]], "list")
  expect_equal(length(cds@clusters[["UMAP"]]), 3)
  expect_equal(length(cds@clusters[["UMAP"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["UMAP"]]$cluster_result$optim_res$membership[[1]],
               1)
  expect_equal(length(unique(clusters(cds, reduction_method = "UMAP"))), 1)

  # non-standard opts
  cds <- cluster_cells(cds, k=22, resolution = c(0.75, 0.3), weight = T,
                       num_iter = 1, partition_qval = .1, verbose = TRUE,
                       random_seed = 100)
  expect_is(cds@clusters[["UMAP"]], "list")
  expect_equal(length(cds@clusters[["UMAP"]]), 3)
  expect_equal(length(cds@clusters[["UMAP"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["UMAP"]]$cluster_result$optim_res$membership[[1]],
               1)
  expect_equal(length(unique(clusters(cds, reduction_method = "UMAP"))), 19)

  ## louvain
  cds <- cluster_cells(cds, cluster_method = "louvain", random_seed = 100)
  expect_is(cds@clusters[["UMAP"]], "list")
  expect_equal(length(cds@clusters[["UMAP"]]), 3)
  expect_equal(length(cds@clusters[["UMAP"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["UMAP"]]$cluster_result$optim_res$membership[[1]],
               9)
  expect_equal(length(unique(clusters(cds, reduction_method = "UMAP"))), 9)

  # non-standard opts
  cds <- cluster_cells(cds, cluster_method = "louvain", k=22, weight = T,
                       num_iter = 1, partition_qval = .1, verbose = TRUE,
                       random_seed = 100)
  expect_is(cds@clusters[["UMAP"]], "list")
  expect_equal(length(cds@clusters[["UMAP"]]), 3)
  expect_equal(length(cds@clusters[["UMAP"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["UMAP"]]$cluster_result$optim_res$membership[[1]],
               10)
  expect_equal(length(unique(clusters(cds, reduction_method = "UMAP"))), 10)

  ### tSNE
  ##leiden
  cds <- cluster_cells(cds, reduction_method = "tSNE", random_seed = 100)
  expect_is(cds@clusters[["tSNE"]], "list")
  expect_equal(length(cds@clusters[["tSNE"]]), 3)
  expect_equal(length(cds@clusters[["tSNE"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["tSNE"]]$cluster_result$optim_res$membership[[1]],
               1)
  expect_equal(length(unique(clusters(cds, reduction_method = "tSNE"))), 1)

  # non-standard opts
  cds <- cluster_cells(cds, reduction_method = "tSNE", k=22,
                       resolution = c(0.75, 0.3), weight = T, num_iter = 1,
                       partition_qval = .1, verbose = TRUE, random_seed = 100)
  expect_is(cds@clusters[["tSNE"]], "list")
  expect_equal(length(cds@clusters[["tSNE"]]), 3)
  expect_equal(length(cds@clusters[["tSNE"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["tSNE"]]$cluster_result$optim_res$membership[[1]],
               14)
  expect_equal(length(unique(clusters(cds, reduction_method = "tSNE"))), 19)

  ### PCA

  cds <- cluster_cells(cds, reduction_method = "PCA", random_seed = 100)
  expect_is(cds@clusters[["PCA"]], "list")
  expect_equal(length(cds@clusters[["PCA"]]), 3)
  expect_equal(length(cds@clusters[["PCA"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["PCA"]]$cluster_result$optim_res$membership[[1]],
               1)
  expect_equal(length(unique(clusters(cds, reduction_method = "PCA"))), 1)

  # non-standard opts
  cds <- cluster_cells(cds, reduction_method = "PCA", k=22,
                       resolution = c(0.75, 0.3), weight = T, num_iter = 2,
                       partition_qval = .1, verbose = TRUE, random_seed = 100)
  expect_is(cds@clusters[["PCA"]], "list")
  expect_equal(length(cds@clusters[["PCA"]]), 3)
  expect_equal(length(cds@clusters[["PCA"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["PCA"]]$cluster_result$optim_res$membership[[1]],
               44)
  expect_equal(length(unique(clusters(cds, reduction_method = "PCA"))), 59)

  cds <- cluster_cells(cds, reduction_method = "PCA", k=22, weight = T,
                       num_iter = 2, partition_qval = .1, resolution = .1,
                       random_seed = 100)
  expect_is(cds@clusters[["PCA"]], "list")
  expect_equal(length(cds@clusters[["PCA"]]), 3)
  expect_equal(length(cds@clusters[["PCA"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(as.numeric(cds@clusters[["PCA"]]$cluster_result$optim_res$membership[1]),
               1)
  expect_equal(length(unique(clusters(cds, reduction_method = "PCA"))), 24)

  cds <- cluster_cells(cds, reduction_method = "PCA", k=22, weight = T,
                       num_iter = 2, partition_qval = .1, resolution = .01,
                       random_seed = 100)
  expect_equal(length(unique(clusters(cds, reduction_method = "PCA"))), 1)
  cds <- cluster_cells(cds, reduction_method = "PCA", k=22, weight = T,
                       num_iter = 2, partition_qval = .1, resolution = 20,
                       random_seed = 100)
  expect_equal(length(unique(clusters(cds, reduction_method = "PCA"))), 500)

  cds <- preprocess_cds(cds, method = "LSI")
  cds <- cluster_cells(cds, reduction_method = "LSI", random_seed = 100)
  expect_is(cds@clusters[["LSI"]], "list")
  expect_equal(length(cds@clusters[["LSI"]]), 3)
  expect_equal(length(cds@clusters[["LSI"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["LSI"]]$cluster_result$optim_res$membership[[1]],
               1)
  expect_equal(length(unique(clusters(cds, reduction_method = "LSI"))), 1)

  # non-standard opts
  cds <- cluster_cells(cds, reduction_method = "LSI", k=22,
                       resolution = c(0.75, 0.3), weight = T, num_iter = 2,
                       partition_qval = .1, verbose = TRUE, random_seed = 100)
  expect_is(cds@clusters[["LSI"]], "list")
  expect_equal(length(cds@clusters[["LSI"]]), 3)
  expect_equal(length(cds@clusters[["LSI"]]$cluster_result$optim_res$membership),
               nrow(colData(cds)))
  expect_equal(cds@clusters[["LSI"]]$cluster_result$optim_res$membership[[1]],
               10)
  expect_equal(length(unique(clusters(cds, reduction_method = "LSI"))), 61)
})

cds <- load_a549()
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = 20)
cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1)

test_that("cluster_cells works", {
  skip_not_travis()
  cds <- cluster_cells(cds, random_seed = 100)
  expect_equal(sum(clusters(cds) == 1), 500)
  expect_equal(as.character(clusters(cds)[1]), "1")
})

