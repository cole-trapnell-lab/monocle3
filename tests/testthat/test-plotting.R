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

cds <- load_a549()
cds <- estimate_size_factors(cds)


#test_that("plot_cell_trajectory works", {
#  expect_error(plot_cell_trajectory(cds,
#                                    color_by = "cell_type2"), "No normalized data projection calculated. Please run preprocess_cds before running plot_cell_trajectory")
#  cds <- preprocess_cds(cds, num_dim = 20)
#  expect_error(plot_cell_trajectory(cds,
#                                    color_by = "cell_type2"), "No dimensionality reduction for UMAP calculated. Please run reduce_dimensions with reduction_method = UMAP before running partition_cells.")
#})
