context("test-learn_graph")
skip_not_travis <- function ()
{
  if (identical(Sys.getenv("TRAVIS"), "true")) {
    return(invisible(TRUE))
  }
  skip("Not on Travis")
}

cds <- load_a549()

test_that("test learn_graph error messages work", {
  skip_on_travis()
  expect_error(cds <- learn_graph(cds),
               "No dimensionality reduction for UMAP calculated. Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.")
  cds <- preprocess_cds(cds)
  expect_error(cds <- learn_graph(cds),
               "No dimensionality reduction for UMAP calculated. Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.")
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
cds <- cluster_cells(cds, resolution = .01)

test_that("learn_graph stays the same", {
  skip_on_travis()
  cds <- learn_graph(cds)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 33)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "24")

  # Force partition
  temp <- rep(c(1,2), length.out=length(partitions(cds)))
  names(temp) <- names(partitions(cds))
  cds@clusters[["UMAP"]]$partitions <- temp
  cds <- learn_graph(cds, use_partition = FALSE)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 33)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "24")

  cds <- learn_graph(cds)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 74)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "36")

  cds <- learn_graph(cds, close_loop = TRUE)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 74)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "36")

  cds <- learn_graph(cds, learn_graph_control = list(prune_graph = FALSE))
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 144)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "8")
})


#### TRAVIS ####


cds <- load_a549()

test_that("test learn_graph error messages work", {
  skip_not_travis()
  expect_error(cds <- learn_graph(cds),
               "No dimensionality reduction for UMAP calculated. Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.")
  cds <- preprocess_cds(cds)
  expect_error(cds <- learn_graph(cds),
               "No dimensionality reduction for UMAP calculated. Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.")
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
cds <- cluster_cells(cds, resolution = .01)

test_that("learn_graph stays the same", {
  skip_not_travis()
  cds <- learn_graph(cds)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 33)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "24")

  # Force partition
  temp <- rep(c(1,2), length.out=length(partitions(cds)))
  names(temp) <- names(partitions(cds))
  cds@clusters[["UMAP"]]$partitions <- temp
  cds <- learn_graph(cds, use_partition = FALSE)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 33)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "24")

  cds <- learn_graph(cds)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 74)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "36")

  cds <- learn_graph(cds, close_loop = TRUE)
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 74)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "36")

  cds <- learn_graph(cds, learn_graph_control = list(prune_graph = FALSE))
  expect_is(principal_graph(cds)[["UMAP"]], "igraph")
  expect_equal(length(principal_graph(cds)[["UMAP"]]), 144)
  expect_equal(as.character(principal_graph(cds)[["UMAP"]][[1]]$Y_1[[1]]), "8")
})


