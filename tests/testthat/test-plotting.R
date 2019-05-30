context("test-plotting")

cds <- load_a549()
cds <- estimate_size_factors(cds)

test_that("plot_genes_violin works", {
  cds_subset <- cds[c("ENSG00000228253.1", "ENSG00000103034.14",
                      "ENSG00000223519.8"),]
  plot_genes_violin(cds_subset, grouping="culture_plate")
  plot_genes_violin(cds_subset, grouping="culture_plate", min_expr = 10)
  plot_genes_violin(cds_subset, grouping="culture_plate", ncol=2)
  # TO DO
#  plot_genes_violin(cds_subset, grouping="culture_plate", color_by = "dose")
  plot_genes_violin(cds_subset, grouping="culture_plate", normalize = FALSE)
  plot_genes_violin(cds_subset, grouping="culture_plate", log_scale = FALSE)
  plot_genes_violin(cds_subset, grouping="culture_plate",
                    label_by_short_name = FALSE)
#  plot_genes_violin(cds_subset, grouping="dose", plot_trend = TRUE,
#                    ncol = 3)
})

test_that("plot_percent_cells_positive works", {
  cds_subset <- cds[c("ENSG00000228253.1", "ENSG00000103034.14",
                      "ENSG00000223519.8"),]
  plot_percent_cells_positive(cds_subset, grouping="culture_plate")
  plot_percent_cells_positive(cds_subset, grouping="culture_plate",
                              min_expr = 10)
  plot_percent_cells_positive(cds_subset, grouping="culture_plate", ncol=2)
  plot_percent_cells_positive(cds_subset, grouping="culture_plate",
                              normalize = FALSE)
  plot_percent_cells_positive(cds_subset, grouping="culture_plate",
                              plot_as_count = TRUE)
  plot_percent_cells_positive(cds_subset, grouping="culture_plate",
                    label_by_short_name = FALSE)
  plot_percent_cells_positive(cds_subset, grouping="culture_plate",
                              plot_limits = c(0,200))
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
               "Pseudotime must be a column in colData. Please run order_cells before running plot_genes_in_pseudotime.")

  cds <- preprocess_cds(cds)
  cds <- reduce_dimension(cds)
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

test_that("plot_cells_3d doesn't error", {
  cds <- preprocess_cds(cds)
  cds <- reduce_dimension(cds, max_components = 3)
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  cds <- order_cells(cds, root_pr_nodes = "Y_1")
  expect_equal(1, 1)
})


