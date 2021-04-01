context("test-graph_test")
skip_if_offline()

skip_not_travis <- function ()
{
  if (identical(Sys.getenv("TRAVIS"), "true")) {
    return(invisible(TRUE))
  }
  skip("Not on Travis")
}


cds <- monocle3:::load_worm_embryo()
set.seed(42)

cds <- preprocess_cds(cds, num_dim = 100, residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")
cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1)
cds <- cluster_cells(cds)
cds <- learn_graph(cds, learn_graph_control=list(ncenter=1000), close_loop=TRUE)
#plot_cells(cds, color_cells_by="cell.type")

# test_that("test graph_test error messages work", {
#   expect_error(cds <- graph_test(cds),
#                "No dimensionality reduction for UMAP calculated. Please run reduce_dimension with reduction_method = UMAP and partition_cells before running learn_graph.")
#   cds <- preprocess_cds(cds)
#   expect_error(cds <- graph_test(cds),
#                "No dimensionality reduction for UMAP calculated. Please run reduce_dimension with reduction_method = UMAP and partition_cells before running learn_graph.")
#   cds <- reduce_dimension(cds)
#   expect_error(cds <- graph_test(cds),
#                "No cell partition for UMAP calculated. Please run partition_cells with reduction_method = UMAP before running learn_graph.")
#   cds <- reduce_dimension(cds)
#   expect_error(cds <- graph_test(cds),
#                "No cell partition for UMAP calculated. Please run partition_cells with reduction_method = UMAP before running learn_graph.")
#
#   #expect_error(cds <- learn_graph(cds, learn_graph_control = list(FALSE)), "")
#   #expect_error(cds <- graph_test(cds, learn_graph_control = list(prune = FALSE)), "Unknown variable in learn_graph_control")
# })

test_that("test graph_test returns Dex-dependent genes",{
  skip_on_travis()
  test_cds = cds
  pos_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "che-1",]
  pr_test_res = graph_test(pos_ctrl_gene)
  expect_equal(pr_test_res$status[1], "OK")
  expect_equal(pr_test_res$morans_I, 0.65, tolerance=1e-2)
  expect_equal(pr_test_res$morans_test_statistic, 204.72, tolerance=1e-1)
  expect_lt(pr_test_res$p_value[1], 0.05)

  neg_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "R02D3.1",]
  pr_test_res = graph_test(neg_ctrl_gene)
  expect_equal(pr_test_res$status[1], "OK")
  expect_equal(pr_test_res$morans_I, -0.00393, tolerance=1e-4)
  expect_equal(pr_test_res$morans_test_statistic, -1.11, tolerance=1e-1)
  expect_gt(pr_test_res$p_value[1], 0.05)
})

test_that("test graph_test returns few genes under UMAP coordinate randomization",{
  skip_on_travis()
  ciliated_genes = c("che-1",
                     "hlh-17",
                     "nhr-6",
                     "dmd-6",
                     "ceh-36",
                     "ham-1")
  cds_subset = cds[rowData(cds)$gene_short_name %in% ciliated_genes,]

  test_cds = cds_subset
  nr = nrow(reducedDims(test_cds)$UMAP)
  reducedDims(test_cds)$UMAP[,1] = reducedDims(test_cds)$UMAP[sample.int(nr),1]
  reducedDims(test_cds)$UMAP[,2] = reducedDims(test_cds)$UMAP[sample.int(nr),2]
  test_cds <- cluster_cells(test_cds)
  test_cds <- learn_graph(test_cds)
  pr_test_res = graph_test(test_cds, neighbor_graph="principal_graph", k=50)

  num_degs = sum(pr_test_res$q_value < 0.05)
  expect_equal(num_degs, 1)
})



#### TRAVIS ####



test_that("test graph_test returns Dex-dependent genes",{
  skip_not_travis()
  test_cds = cds
  pos_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "che-1",]
  pr_test_res = graph_test(pos_ctrl_gene)
  expect_equal(pr_test_res$status[1], "OK")
  expect_equal(pr_test_res$morans_I, 0.653, tolerance=1e-2)
  expect_equal(pr_test_res$morans_test_statistic, 204.72, tolerance=1e-1)
  expect_lt(pr_test_res$p_value[1], 0.05)

  neg_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "R02D3.1",]
  pr_test_res = graph_test(neg_ctrl_gene)
  expect_equal(pr_test_res$status[1], "OK")
  expect_equal(pr_test_res$morans_I, -0.00153, tolerance=1e-4)
  expect_equal(pr_test_res$morans_test_statistic, -0.402, tolerance=1e-2)
  expect_gt(pr_test_res$p_value[1], 0.05)
})

test_that("test graph_test returns few genes under UMAP coordinate randomization",{
  skip_not_travis()
  ciliated_genes = c("che-1",
                     "hlh-17",
                     "nhr-6",
                     "dmd-6",
                     "ceh-36",
                     "ham-1")
  cds_subset = cds[rowData(cds)$gene_short_name %in% ciliated_genes,]

  test_cds = cds_subset
  nr = nrow(reducedDims(test_cds)$UMAP)
  reducedDims(test_cds)$UMAP[,1] = reducedDims(test_cds)$UMAP[sample.int(nr),1]
  reducedDims(test_cds)$UMAP[,2] = reducedDims(test_cds)$UMAP[sample.int(nr),2]
  test_cds <- cluster_cells(test_cds)
  test_cds <- learn_graph(test_cds)
  pr_test_res = graph_test(test_cds, neighbor_graph="principal_graph", k=50)

  num_degs = sum(pr_test_res$q_value < 0.05)
  expect_equal(num_degs, 0)
})


