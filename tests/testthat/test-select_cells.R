context("test-select_cells")

cds <- load_a549()

test_that("test choose_graph_segments error messages work", {
  expect_error(cds <- choose_graph_segments(cds),
               "No dimensionality reduction for UMAP calculated. Please run reduce_dimensions with reduction_method = UMAP, cluster_cells and learn_graph before running choose_graph_segments.")
  cds <- preprocess_cds(cds)
  expect_error(cds <- choose_graph_segments(cds),
               "No dimensionality reduction for UMAP calculated. Please run reduce_dimensions with reduction_method = UMAP, cluster_cells and learn_graph before running choose_graph_segments.")
  cds <- reduce_dimension(cds)
  expect_error(cds <- choose_graph_segments(cds),
               "No cell clusters for UMAP calculated. Please run cluster_cells with reduction_method = UMAP and run learn_graph before running choose_graph_segments.")
  cds <- cluster_cells(cds)
  expect_error(cds <- choose_graph_segments(cds),
               "No principal graph for UMAP calculated. Please run learn_graph with reduction_method = UMAP before running choose_graph_segments.")
  expect_error(cds_sub <- choose_graph_segments(cds, starting_pr_node = "Y_4"),
               "If not using interactive mode, you must provide ending_pr_nodes.")
  expect_error(cds_sub <- choose_graph_segments(cds, ending_pr_nodes = "Y_4"),
               "If not using interactive mode, you must provide a starting_pr_node")
  })

cds <- preprocess_cds(cds)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

test_that("test choose_graph_segments non-interactive works", {
  cds_sub <- choose_graph_segments(cds, starting_pr_node = "Y_1",
                                   ending_pr_nodes = "Y_6")
  assertthat::are_equal(ncol(cds_sub), 281)
  cds_sub <- choose_graph_segments(cds, starting_pr_node = "Y_1",
                                   ending_pr_nodes = c( "Y_3"))
  assertthat::are_equal(ncol(cds_sub), 172)

  cds_sub <- choose_graph_segments(cds, starting_pr_node = "Y_3",
                                   ending_pr_nodes = c( "Y_6", "Y_2"))
  assertthat::are_equal(ncol(cds_sub), 342)
})

