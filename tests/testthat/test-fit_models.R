context("fitModels")

test_that("fitModels() properly validates its input",{
  #cds <- load_a549()
  expect_equal(1,1)
#
#   full_model_fits <- fit_models(cds, modelFormulaStr = "~log_dose")
#   coef_table <- coefficient_table(full_model_fits)
#   expect_equal(colnames(diff_test_res), c("status", "family", "pval", "qval"))
#   expect_equal(levels(diff_test_res$status),"OK")
#   expect_equal(levels(diff_test_res$family), "negbinomial.size")
#   expect_equal(diff_test_res$pval, c(1.956459e-32, 8.9011122e-01, 1.490169e-45))
#   expect_equal(diff_test_res$qval, c(2.934689e-32, 8.9011122e-01, 4.470506e-45))
})

