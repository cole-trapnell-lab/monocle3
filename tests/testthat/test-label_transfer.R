context("test-label_transfer")
skip_not_travis <- function ()
{
  if (identical(Sys.getenv("TRAVIS"), "true")) {
    return(invisible(TRUE))
  }
  skip("Not on Travis")
}

# use tst_table1 and tst_table2 with as.vector() and as.vector(names()) to
# make the ref_vectors* used below.
# ref_vector1 <- as.integer(c(162,80,485,134,311,205,137,32,40,165,40,343,211,154,313,218,133,73,214,184,291,6,38,111,104,115,142,1747))
# names(ref_vector1) <- c("ADF","ADF_AWB","ADL","ADL_parent","AFD","ASE","ASE_parent","ASEL","ASER","ASG","ASG_AWA","ASH","ASI","ASI_parent","ASJ","ASK","ASK_parent","AUA","AWA","AWB","AWC","AWC_ON","Neuroblast_ADF_AWB","Neuroblast_AFD_RMD","Neuroblast_ASE_ASJ_AUA","Neuroblast_ASG_AWA","Neuroblast_ASJ_AUA",NA)
# ref_table1 <- as.table(ref_vector1)
#
# ref_vector2 <- as.integer(c(176,147,506,199,313,205,143,32,40,166,52,368,211,187,318,219,145,73,216,184,293,6,76,863,332,544,155,19))
# names(ref_vector2) <- c("ADF","ADF_AWB","ADL","ADL_parent","AFD","ASE","ASE_parent","ASEL","ASER","ASG","ASG_AWA","ASH","ASI","ASI_parent","ASJ","ASK","ASK_parent","AUA","AWA","AWB","AWC","AWC_ON","Neuroblast_ADF_AWB","Neuroblast_AFD_RMD","Neuroblast_ASE_ASJ_AUA","Neuroblast_ASG_AWA","Neuroblast_ASJ_AUA",NA)
# ref_table2 <- as.table(ref_vector2)

#
# Test label transfer functions.
#
test_that("label transfer", {
  cds <- load_worm_embryo()
  set.seed(2016)
  cds <- preprocess_cds(cds)
  cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")
  set.seed(2016)
  cds <- reduce_dimension(cds, build_nn_index=TRUE)
  cds <- transfer_cell_labels(cds, reduction_method='UMAP', colData(cds), 'cell.type', 'cell.type.xfr')
  cds <- fix_missing_cell_labels(cds, reduction_method='UMAP', from_column_name='cell.type.xfr', to_column_name='cell.type.fix')

#  tst_table1 <- table(colData(cds)[['cell.type.xfr']], useNA='ifany')
#  tst_1_1 <- identical(as.vector(tst_table1), as.vector(ref_table1))
#  tst_1_2 <- identical(as.vector(names(tst_table1)), as.vector(names(ref_table1)))


#  tst_table2 <- table(colData(cds)[['cell.type.fix']], useNA='ifany')
#  tst_2_1 <- identical(as.vector(tst_table2), as.vector(ref_table2))
#  tst_2_2 <- identical(as.vector(names(tst_table2)), as.vector(names(ref_table2)))

#  These tests don't work because the returned nearest
#  neighbors differ on different machines.
#  expect_equivalent(as.vector(tst_table1), as.vector(ref_table1))
#  expect_equivalent(as.vector(names(tst_table1)), as.vector(names(ref_table1)))

#  expect_equivalent(as.vector(tst_table2), as.vector(ref_table2))
#  expect_equivalent(as.vector(names(tst_table2)), as.vector(names(ref_table2)))

  nmatch1 <- nrow(colData(cds)[!is.na(colData(cds)[['cell.type']]) & !is.na(colData(cds)[['cell.type.xfr']]) & colData(cds)[['cell.type']] == colData(cds)[['cell.type.xfr']],])
  nmatch2 <- nrow(colData(cds)[!is.na(colData(cds)[['cell.type.xfr']]) & !is.na(colData(cds)[['cell.type.fix']]) & colData(cds)[['cell.type.xfr']] == colData(cds)[['cell.type.fix']],])
  num_row <- nrow(colData(cds))
#
  # expect_gt apparently allows only integer comparisons.
  expect_gt(100.0 * as.double(nmatch1) / as.double(num_row), 60.0)
  expect_gt(nmatch2, nmatch1)
} )

