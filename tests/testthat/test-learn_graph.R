context("test-learn_graph")

cds <- load_a549()

test_that("test learn_graph error messages work", {
  expect_error(cds <- learn_graph(cds),
               "No dimensionality reduction for UMAP calculated. Please run reduce_dimensions with reduction_method = UMAP and cluster_cells before running learn_graph.")
  cds <- preprocess_cds(cds)
  expect_error(cds <- learn_graph(cds),
               "No dimensionality reduction for UMAP calculated. Please run reduce_dimensions with reduction_method = UMAP and cluster_cells before running learn_graph.")
  cds <- reduce_dimension(cds)
  expect_error(cds <- learn_graph(cds),
               "No cell clusters for UMAP calculated. Please run cluster_cells with reduction_method = UMAP before running learn_graph.")
  #expect_error(cds <- learn_graph(cds, learn_graph_control = list(FALSE)), "")
  expect_error(cds <- learn_graph(cds, learn_graph_control = list(prune = FALSE)), "Unknown variable in learn_graph_control")
})

set.seed(42)
cds <- preprocess_cds(cds)
cds <- reduce_dimension(cds)
cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1)
cds <- cluster_cells(cds)

test_that("learn_graph stays the same", {
  cds <- learn_graph(cds)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 10)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "96")

  # Force partition
  cds@clusters[["UMAP"]]$partitions <- c(1,2)
  cds <- learn_graph(cds, use_partition = FALSE)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 10)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "96")

  cds <- learn_graph(cds)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 10)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "18")

  cds <- learn_graph(cds, close_loop = TRUE)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 10)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "18")

  cds <- learn_graph(cds, learn_graph_control = list(prune_graph = FALSE))
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 10)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "32")
})

