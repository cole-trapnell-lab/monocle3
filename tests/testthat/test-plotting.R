context("test-plotting")

cds <- load_a549()
cds <- estimate_size_factors(cds)

test_that("plot_genes_violin doesn't error", {
  cds_subset <- cds[c("ENSG00000228253.1", "ENSG00000103034.14",
                      "ENSG00000223519.8"),]
  plot_genes_violin(cds_subset)
  plot_genes_violin(cds_subset, group_cells_by="culture_plate")
  plot_genes_violin(cds_subset[1,], group_cells_by="culture_plate")
  plot_genes_violin(cds_subset, group_cells_by="culture_plate", min_expr = 10)
  plot_genes_violin(cds_subset, group_cells_by="culture_plate", nrow=5)
  plot_genes_violin(cds_subset, group_cells_by="culture_plate", nrow=2,
                    ncol = 2)
  plot_genes_violin(cds_subset, group_cells_by="culture_plate", ncol=2)
  plot_genes_violin(cds_subset, group_cells_by="culture_plate",
                    panel_order = c("MT-ATP8", "NDRG4", "KIF28P"))
  plot_genes_violin(cds_subset, group_cells_by="culture_plate",
                    normalize = FALSE)
  plot_genes_violin(cds_subset, group_cells_by="culture_plate",
                    log_scale = FALSE)
  plot_genes_violin(cds_subset, group_cells_by="culture_plate",
                    label_by_short_name = FALSE)
  expect_equal(1, 1)
})

test_that("plot_percent_cells_positive doesn't error", {
  cds_subset <- cds[c("ENSG00000228253.1", "ENSG00000103034.14",
                      "ENSG00000223519.8"),]
  plot_percent_cells_positive(cds_subset)
  plot_percent_cells_positive(cds_subset, group_cells_by="culture_plate")
  plot_percent_cells_positive(cds_subset, group_cells_by="culture_plate",
                              min_expr = 10)
  plot_percent_cells_positive(cds_subset, group_cells_by="culture_plate",
                              ncol=2)
  plot_percent_cells_positive(cds_subset, group_cells_by="culture_plate",
                              normalize = FALSE,
                              panel_order = c("MT-ATP8", "NDRG4", "KIF28P"))
  plot_percent_cells_positive(cds_subset, group_cells_by="culture_plate",
                              plot_as_count = TRUE)
  plot_percent_cells_positive(cds_subset, group_cells_by="culture_plate",
                    label_by_short_name = FALSE)
  plot_percent_cells_positive(cds_subset, group_cells_by="culture_plate",
                              plot_limits = c(0,200))
  expect_equal(1, 1)
})

cds <- load_a549()
cds <- estimate_size_factors(cds)

test_that("plot_percent_cells_positive doesn't error", {




  expect_equal(1, 1)
  })


