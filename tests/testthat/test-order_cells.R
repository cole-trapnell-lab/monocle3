context("test-order_cells")

temp <- readRDS("../testdata/reduced_dims/a549_umap.RDS")
temp3 <- readRDS("../testdata/reduced_dims/a549_umap_3.RDS")
cds <- load_a549()
set.seed(100)

test_that("order_cells error messages work", {
  expect_error(order_cells(cds),
               "No dimensionality reduction for UMAP calculated. Please run reduce_dimensions with reduction_method = UMAP, cluster_cells, and learn_graph before running order_cells." )
  cds <- estimate_size_factors(cds)
  cds <- preprocess_cds(cds, num_dim = 20)
  reducedDims(cds)$UMAP <- temp
  expect_error(order_cells(cds),
               "No cell clusters for UMAP calculated. Please run cluster_cells with reduction_method = UMAP and run learn_graph before running order_cells.")
  cds <- cluster_cells(cds)
  expect_error(order_cells(cds),
               "No principal graph for UMAP calculated. Please run learn_graph with reduction_method = UMAP before running order_cells.")
  cds <- learn_graph(cds)
  expect_error(order_cells(cds, root_cells = c("G07_B02_RT_587"),
                           root_pr_nodes = c("Y_1")),
               "Please specify either root_pr_nodes or root_cells, not both.")
  expect_error(order_cells(cds, root_cells = c("hannah")),
               "All provided root_cells must be present in the cell data set.")
  expect_error(order_cells(cds, root_pr_nodes = c("hannah")),
               "All provided root_pr_nodes must be present in the principal graph.")
  expect_error(order_cells(cds), paste("When not in interactive mode, either",
                                       "root_pr_nodes or root_cells must be",
                                       "provided."))
  expect_error(order_cells(cds, reduction_method = "tSNE"),
               "Currently only 'UMAP' is accepted as a reduction_method.")
})

cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = 20)
reducedDims(cds)$UMAP <- temp
cds <- cluster_cells(cds, cluster_method = "louvain")
cds <- learn_graph(cds)

test_that("order_cells works", {
  cds <- order_cells(cds, root_pr_nodes = "Y_1")
  expect_equal(max(pseudotime(cds)), 11.41136, tol = 1e-2)
  expect_equal(min(pseudotime(cds)), 0)
  expect_equal(as.numeric(pseudotime(cds)[1]), 2.128202, tol = 1e-2)
  cds <- order_cells(cds, root_pr_nodes = c("Y_1", "Y_10"))
  expect_equal(max(pseudotime(cds)), 8.152965, tol = 1e-2)
  expect_equal(min(pseudotime(cds)), 0)
  expect_equal(as.numeric(pseudotime(cds)[1]), 2.128202, tol = 1e-2)
  cds <- order_cells(cds, root_cells = "G07_B02_RT_587")
  expect_equal(max(pseudotime(cds)), 11.41712, tol = 1e-2)
  expect_equal(min(pseudotime(cds)), 0)
  expect_equal(as.numeric(pseudotime(cds)[1]), 2.133954, tol = 1e-2)
  cds <- order_cells(cds, root_cells = c("G07_B02_RT_587", "F06_A01_RT_598"))
  expect_equal(max(pseudotime(cds)), 11.41712, tol = 1e-2)
  expect_equal(min(pseudotime(cds)), 0)
  expect_equal(as.numeric(pseudotime(cds)[1]), 2.133954, tol = 1e-2)
})

reducedDims(cds)$UMAP <- temp3
cds <- cluster_cells(cds, cluster_method = "louvain")
cds <- learn_graph(cds)

test_that("order_cells works 3d", {
  cds <- order_cells(cds, root_pr_nodes = "Y_1")
  expect_equal(max(pseudotime(cds)), 10.25239, tol = 1e-2)
  expect_equal(min(pseudotime(cds)), 0)
  expect_equal(as.numeric(pseudotime(cds)[1]), 0, tol = 1e-2)
  cds <- order_cells(cds, root_pr_nodes = c("Y_1", "Y_10"))
  expect_equal(max(pseudotime(cds)),  10.25239, tol = 1e-2)
  expect_equal(min(pseudotime(cds)), 0)
  expect_equal(as.numeric(pseudotime(cds)[1]), 0, tol = 1e-2)
  cds <- order_cells(cds, root_cells = "G07_B02_RT_587")
  expect_equal(max(pseudotime(cds)), 12.51549, tol = 1e-2)
  expect_equal(min(pseudotime(cds)), 0)
  expect_equal(as.numeric(pseudotime(cds)[1]), 2.263098, tol = 1e-2)
  cds <- order_cells(cds, root_cells = c("G07_B02_RT_587", "F06_A01_RT_598"))
  expect_equal(max(pseudotime(cds)), 9.960125, tol = 1e-2)
  expect_equal(min(pseudotime(cds)), 0)
  expect_equal(as.numeric(pseudotime(cds)[1]), 2.263098, tol = 1e-2)
})

cds <- cluster_cells(cds, random_seed = 100)
cds <- learn_graph(cds)

test_that("order_cells works leiden", {
  cds <- order_cells(cds, root_pr_nodes = "Y_1")
  expect_equal(max(pseudotime(cds)), 5.255745, tol = 1e-2)
  expect_equal(min(pseudotime(cds)), 0)
  expect_equal(as.numeric(pseudotime(cds)[1]), 1.939774, tol = 1e-2)
  cds <- order_cells(cds, root_pr_nodes = c("Y_1", "Y_2"))
  expect_equal(max(pseudotime(cds)), 5.255745, tol = 1e-2)
  expect_equal(min(pseudotime(cds)), 0)
  expect_equal(as.numeric(pseudotime(cds)[1]), 0.7243531, tol = 1e-2)
  cds <- order_cells(cds, root_cells = "G07_B02_RT_587")
  expect_equal(max(pseudotime(cds)), 7.058048, tol = 1e-2)
  expect_equal(min(pseudotime(cds)), 0)
  expect_equal(as.numeric(pseudotime(cds)[1]), 0.0005813664, tol = 1e-2)
  cds <- order_cells(cds, root_cells = c("G07_B02_RT_587", "F06_A01_RT_598"))
  expect_equal(max(pseudotime(cds)), 3.195603, tol = 1e-2)
  expect_equal(min(pseudotime(cds)), 0)
  expect_equal(as.numeric(pseudotime(cds)[1]), 0.0005813664, tol = 1e-2)
})

reducedDims(cds)$UMAP <- temp3
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

test_that("order_cells works leiden 3d", {
  cds <- order_cells(cds, root_pr_nodes = "Y_1")
  expect_equal(max(pseudotime(cds)), 5.255745, tol = 1e-2)
  expect_equal(min(pseudotime(cds)), 0)
  expect_equal(as.numeric(pseudotime(cds)[1]), 1.939774, tol = 1e-2)
  cds <- order_cells(cds, root_pr_nodes = c("Y_1", "Y_2"))
  expect_equal(max(pseudotime(cds)), 5.255745, tol = 1e-2)
  expect_equal(min(pseudotime(cds)), 0)
  expect_equal(as.numeric(pseudotime(cds)[1]), 0.7243531, tol = 1e-2)
  cds <- order_cells(cds, root_cells = "G07_B02_RT_587")
  expect_equal(max(pseudotime(cds)), 7.058048, tol = 1e-2)
  expect_equal(min(pseudotime(cds)), 0)
  expect_equal(as.numeric(pseudotime(cds)[1]), 0.0005813664, tol = 1e-2)
  cds <- order_cells(cds, root_cells = c("G07_B02_RT_587", "F06_A01_RT_598"))
  expect_equal(max(pseudotime(cds)), 3.195603, tol = 1e-2)
  expect_equal(min(pseudotime(cds)), 0)
  expect_equal(as.numeric(pseudotime(cds)[1]), 0.0005813664, tol = 1e-2)
})


