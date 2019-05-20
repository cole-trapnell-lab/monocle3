
  small_a549_colData_df <- readRDS(system.file("extdata", "small_a549_dex_pdata.rda", package = "monocle3"))
  small_a549_rowData_df <- readRDS(system.file("extdata", "small_a549_dex_fdata.rda", package = "monocle3"))
  small_a549_exprs <- readRDS(system.file("extdata", "small_a549_dex_exprs.rda", package = "monocle3"))
  small_a549_exprs <- small_a549_exprs[,row.names(small_a549_colData_df)]


test_that("new_cell_data_set works" ,{
  expect_error(cds <- new_cell_data_set(
    expression_data = as.data.frame(as.matrix(small_a549_exprs)),
    cell_metadata = small_a549_colData_df,
    gene_metadata = small_a549_rowData_df),
    "Argument expression_data must be a matrix - either sparse from the Matrix package or dense")

  expect_warning(cds <- new_cell_data_set(expression_data = small_a549_exprs,
                                          cell_metadata = small_a549_colData_df))

  expect_error(cds <- new_cell_data_set(expression_data = small_a549_exprs,
                                        cell_metadata = small_a549_colData_df[1:100,],
                                        gene_metadata = small_a549_rowData_df),
               "cell_metadata must be NULL or have the same number of rows as columns in expression_data")

  expect_error(cds <- new_cell_data_set(expression_data = small_a549_exprs,
                                        cell_metadata = small_a549_colData_df,
                                        gene_metadata = small_a549_rowData_df[1:100,]),
               "gene_metadata must be NULL or have the same number of rows as rows in expression_data")
  temp <- small_a549_colData_df
  row.names(temp)[1] <- "HP"
  expect_error(cds <- new_cell_data_set(expression_data = small_a549_exprs,
                                        cell_metadata = temp,
                                        gene_metadata = small_a549_rowData_df),
               "row.names of cell_metadata must be equal to colnames of expression_data")

  temp <- small_a549_rowData_df
  row.names(temp)[1] <- "HP"
  expect_error(cds <- new_cell_data_set(expression_data = small_a549_exprs,
                                        cell_metadata = small_a549_colData_df,
                                        gene_metadata = temp),
               "row.names of gene_metadata must be equal to row.names of expression_data")
  cds <- new_cell_data_set(expression_data = small_a549_exprs,
                           cell_metadata = small_a549_colData_df,
                           gene_metadata = small_a549_rowData_df)
  expect_is(cds, "cell_data_set")

})

