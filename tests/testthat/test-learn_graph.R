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
temp <- readRDS("../testdata/reduced_dims/a549_umap.RDS")
reducedDims(cds)$UMAP <- temp
cds <- cluster_cells(cds, resolution = .01)

test_that("learn_graph stays the same", {
  cds <- learn_graph(cds)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 10)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "36")

  # Force partition
  temp <- rep(c(1,2), length.out=length(partitions(cds)))
  names(temp) <- names(partitions(cds))
  cds@clusters[["UMAP"]]$partitions <- temp
  cds <- learn_graph(cds, use_partition = FALSE)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 10)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "36")

  cds <- learn_graph(cds)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 10)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "44")

  cds <- learn_graph(cds, close_loop = TRUE)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 10)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "44")

  cds <- learn_graph(cds, learn_graph_control = list(prune_graph = FALSE))
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 10)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "81")
})


