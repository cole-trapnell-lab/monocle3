context("test-plotting")

cds <- load_a549()
cds <- estimate_size_factors(cds)

test_that("plot_genes_violin doesn't error", {
  cds_subset <- cds[c("ENSG00000228253.1", "ENSG00000103034.14",
                      "ENSG00000223519.8"),]

  plot_genes_violin(cds_subset, group_cells_by="culture_plate")
  plot_genes_violin(cds_subset, group_cells_by="culture_plate", min_expr = 10)
  plot_genes_violin(cds_subset, group_cells_by="culture_plate", ncol=2)
  # TO DO
#  plot_genes_violin(cds_subset, group_cells_by="culture_plate", color_by = "dose")
  plot_genes_violin(cds_subset, group_cells_by="culture_plate", normalize = FALSE)
  plot_genes_violin(cds_subset, group_cells_by="culture_plate", log_scale = FALSE)
  plot_genes_violin(cds_subset, group_cells_by="culture_plate",
                    label_by_short_name = FALSE)
  expect_equal(0,0)
#  plot_genes_violin(cds_subset, group_cells_by="dose", plot_trend = TRUE,
#                    ncol = 3)
})

test_that("plot_percent_cells_positive doesn't error", {
  cds_subset <- cds[c("ENSG00000228253.1", "ENSG00000103034.14",
                      "ENSG00000223519.8"),]

  plot_percent_cells_positive(cds_subset, group_cells_by="culture_plate")
  plot_percent_cells_positive(cds_subset, group_cells_by="culture_plate",
                              min_expr = 10)
  plot_percent_cells_positive(cds_subset, group_cells_by="culture_plate", ncol=2)
  plot_percent_cells_positive(cds_subset, group_cells_by="culture_plate",
                              normalize = FALSE)
  plot_percent_cells_positive(cds_subset, group_cells_by="culture_plate",
                              plot_as_count = TRUE)
  plot_percent_cells_positive(cds_subset, group_cells_by="culture_plate",
                    label_by_short_name = FALSE)
  plot_percent_cells_positive(cds_subset, group_cells_by="culture_plate",
                              plot_limits = c(0,200))
  expect_equal(1, 1)
})


test_that("plot_pc_variance_explained doesn't error", {
  expect_error(plot_pc_variance_explained(cds),
               "Data has not been preprocessed with PCA. Please run preprocess_cds with method = 'PCA' before running plot_pc_variance_explained.")

  cds <- preprocess_cds(cds)
  plot_pc_variance_explained(cds)

  expect_equal(1, 1)
  })

test_that("plot_genes_in_pseudotime doesn't error", {
  cds_subset <- cds[c("ENSG00000228253.1", "ENSG00000103034.14",
                      "ENSG00000223519.8"),]
  expect_error(plot_genes_in_pseudotime(cds_subset),
               "No pseudotime calculated. Must call order_cells first.")

  cds <- preprocess_cds(cds)
  temp <- readRDS("../testdata/reduced_dims/a549_umap.RDS")
  reducedDims(cds)$UMAP <- temp
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  cds <- order_cells(cds, root_pr_nodes = "Y_1")
  cds_subset <- cds[c("ENSG00000228253.1", "ENSG00000103034.14",
                      "ENSG00000223519.8"),]
  plot_genes_in_pseudotime(cds_subset)
  plot_genes_in_pseudotime(cds_subset, min_expr = 5)
  plot_genes_in_pseudotime(cds_subset, nrow = 5)

  plot_genes_in_pseudotime(cds_subset, ncol = 5,
                           panel_order = c("MT-ATP8", "NDRG4", "KIF28P"))

  plot_genes_in_pseudotime(cds_subset, color_cells_by = "culture_plate")
  plot_genes_in_pseudotime(cds_subset, label_by_short_name = F)
  plot_genes_in_pseudotime(cds_subset, vertical_jitter = 1)
  plot_genes_in_pseudotime(cds_subset, horizontal_jitter = 1)
  expect_equal(1, 1)
})

test_that("plot_cells doesn't error", {
  cds <- preprocess_cds(cds)
  plot_cells(cds, reduction_method = "PCA", x=2, y=5)
  temp <- readRDS("../testdata/reduced_dims/a549_umap.RDS")
  reducedDims(cds)$UMAP <- temp
  plot_cells(cds)
  cds <- cluster_cells(cds)
  plot_cells(cds)
  cds <- learn_graph(cds)
  suppressWarnings(plot_cells(cds))
  cds <- order_cells(cds, root_pr_nodes = "Y_1")
  plot_cells(cds, color_cells_by = "pseudotime")
  plot_cells(cds, color_cells_by = "partition")
  plot_cells(cds, color_cells_by = "cluster")
  plot_cells(cds, group_cells_by="cluster", color_cells_by="partition")
  #plot_cells(cds, group_cells_by="culture_plate", color_cells_by="partition")
  plot_cells(cds, genes = c("NDRG4"), min_expr = 1, rasterize = TRUE)
  plot_cells(cds, genes = c("NDRG4", "MT-ATP8", "ENSG00000240216.7"))
  plot_cells(cds, show_trajectory_graph = FALSE)
  plot_cells(cds, trajectory_graph_color = "blue",
             trajectory_graph_segment_size = 3)
  plot_cells(cds, label_cell_groups = FALSE)
  plot_cells(cds, label_groups_by_cluster = FALSE, group_cells_by="partition")
  plot_cells(cds, group_label_size = 10)
  plot_cells(cds, labels_per_group = 4)
  plot_cells(cds, label_branch_points = FALSE, label_roots = FALSE,
             label_leaves = FALSE, graph_label_size = 5)
  plot_cells(cds, cell_size = 1, alpha = .1)

})



