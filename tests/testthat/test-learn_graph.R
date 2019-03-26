context("test-learn_graph")

cds <- load_a549()

test_that("test learn_graph error messages work", {
  expect_error(cds <- learn_graph(cds),
               "No normalized data projection calculated. Please run preprocess_cds, reduce_dimensions, and partition_cells before running learn_graph.")
  cds <- preprocess_cds(cds)
  expect_error(cds <- learn_graph(cds),
               "No dimensionality reduction for UMAP calculated. Please run reduce_dimensions with reduction_method = UMAP and partition_cells before running learn_graph.")
  cds <- reduce_dimension(cds)
  expect_error(cds <- learn_graph(cds),
               "No cell partition for UMAP calculated. Please run partition_cells with reduction_method = UMAP before running learn_graph.")
  #expect_error(cds <- learn_graph(cds, learn_graph_control = list(FALSE)), "")
  expect_error(cds <- learn_graph(cds, learn_graph_control = list(prune = FALSE)), "Unknown variable in learn_graph_control")
})

cds <- preprocess_cds(cds)
cds <- reduce_dimension(cds)
cds <- reduce_dimension(cds)
cds <- partition_cells(cds)

test_that("learn_graph stays the same", {
  cds <- learn_graph(cds)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 10)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "6")

  # Force partition
  colData(cds)$louvain_component <- c(1,2)
  cds <- learn_graph(cds, use_partition = FALSE)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 10)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "6")

  cds <- learn_graph(cds)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 10)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "48")

  cds <- learn_graph(cds, close_loop = TRUE)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 10)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "48")

  cds <- learn_graph(cds, learn_graph_control = list(prune_graph = FALSE))
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 10)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "108")
})

