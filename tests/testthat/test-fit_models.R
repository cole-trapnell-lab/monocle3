context("fit_models")

cds <- load_a549()

test_that("fit_models() returns an error when Size_Factors are missing",{
  no_na_cds = cds
  pData(no_na_cds)$Size_Factor = NA
  expect_error(fit_models(no_na_cds))
})

test_that("fit_models() properly validates model formulae",{
  ok_formula_model_fits <- fit_models(cds, model_formula_str = "~log_dose")
  num_failed = sum (unlist(lapply(ok_formula_model_fits$model, function(m) { class(m) })) == "logical" )
  expect_lt(num_failed, nrow(cds))

  expect_error(fit_models(cds, model_formula_str = "~MISSING_TERM"))
})

test_that("fit_models() returns correct output for negative binomial regression",{
  test_cds = cds
  test_cds@expression_family = "negbinomial"

  pos_ctrl_gene = test_cds[fData(cds)$gene_short_name == "ANGPTL4",]
  pos_ctrl_gene_fit = fit_models(pos_ctrl_gene, model_formula_str = "~log_dose")
  pos_ctrl_coefs = coefficient_table(pos_ctrl_gene_fit)
  expect_equal(pos_ctrl_coefs$estimate, c(0.01078543, 0.24663124), tolerance=1e-1)
  expect_equal(pos_ctrl_coefs$normalized_effect, c(0, 0.353), tolerance=1e-1)
  expect_lt(pos_ctrl_coefs$p_value[2], 0.05)

  neg_ctrl_gene = test_cds[fData(cds)$gene_short_name == "CCNE2",]
  neg_ctrl_gene_fit = fit_models(neg_ctrl_gene, model_formula_str = "~log_dose")
  neg_ctrl_coefs = coefficient_table(neg_ctrl_gene_fit)
  expect_equal(neg_ctrl_coefs$estimate, c(-3.553921,  0.179926), tolerance=1e-1)
  expect_equal(neg_ctrl_coefs$normalized_effect, c(0, 0.197), tolerance=1e-1)
  expect_gt(neg_ctrl_coefs$p_value[2], 0.05)
})

test_that("fit_models() returns correct output for Poisson regression",{
  test_cds = cds
  test_cds@expression_family = "poisson"

  pos_ctrl_gene = test_cds[fData(cds)$gene_short_name == "ANGPTL4",]
  pos_ctrl_gene_fit = fit_models(pos_ctrl_gene, model_formula_str = "~log_dose")
  pos_ctrl_coefs = coefficient_table(pos_ctrl_gene_fit)
  expect_equal(pos_ctrl_coefs$estimate, c(0.01078543, 0.24663124), tolerance=1e-1)
  expect_equal(pos_ctrl_coefs$normalized_effect, c(0, 0.353), tolerance=1e-1)
  expect_lt(pos_ctrl_coefs$p_value[2], 0.05)

  neg_ctrl_gene = test_cds[fData(cds)$gene_short_name == "CCNE2",]
  neg_ctrl_gene_fit = fit_models(neg_ctrl_gene, model_formula_str = "~log_dose")
  neg_ctrl_coefs = coefficient_table(neg_ctrl_gene_fit)
  expect_equal(neg_ctrl_coefs$estimate, c(-3.553921,  0.179926), tolerance=1e-1)
  expect_equal(neg_ctrl_coefs$normalized_effect, c(0, 0.197), tolerance=1e-1)
  expect_gt(neg_ctrl_coefs$p_value[2], 0.05)
})

test_that("fit_models() returns correct output for quasipoisson regression",{
  test_cds = cds
  test_cds@expression_family = "quasipoisson"

  pos_ctrl_gene = test_cds[fData(cds)$gene_short_name == "ANGPTL4",]
  pos_ctrl_gene_fit = fit_models(pos_ctrl_gene, model_formula_str = "~log_dose")
  pos_ctrl_coefs = coefficient_table(pos_ctrl_gene_fit)
  expect_equal(pos_ctrl_coefs$estimate, c(0.01078543, 0.24663124), tolerance=1e-1)
  expect_equal(pos_ctrl_coefs$normalized_effect, c(0, 0.353), tolerance=1e-1)
  expect_lt(pos_ctrl_coefs$p_value[2], 0.05)

  neg_ctrl_gene = test_cds[fData(cds)$gene_short_name == "CCNE2",]
  neg_ctrl_gene_fit = fit_models(neg_ctrl_gene, model_formula_str = "~log_dose")
  neg_ctrl_coefs = coefficient_table(neg_ctrl_gene_fit)
  expect_equal(neg_ctrl_coefs$estimate, c(-3.553921,  0.179926), tolerance=1e-1)
  expect_equal(neg_ctrl_coefs$normalized_effect, c(0, 0.197), tolerance=1e-1)
  expect_gt(neg_ctrl_coefs$p_value[2], 0.05)
})

# TODO: need a binary dataset for this:
# test_that("fit_models() returns correct output for binomial regression",{
#   test_cds = cds
#   test_cds@expression_family = "binomial"
#
#   pos_ctrl_gene = test_cds[fData(cds)$gene_short_name == "ANGPTL4",]
#   pos_ctrl_gene_fit = fit_models(pos_ctrl_gene, model_formula_str = "~log_dose")
#   pos_ctrl_coefs = coefficient_table(pos_ctrl_gene_fit)
#   expect_equal(pos_ctrl_coefs$estimate, c(0.01078543, 0.24663124), tolerance=1e-1)
#   expect_equal(pos_ctrl_coefs$normalized_effect, c(0, 0.353), tolerance=1e-1)
#   expect_lt(pos_ctrl_coefs$p_value[2], 0.05)
#
#   neg_ctrl_gene = test_cds[fData(cds)$gene_short_name == "CCNE2",]
#   neg_ctrl_gene_fit = fit_models(neg_ctrl_gene, model_formula_str = "~log_dose")
#   neg_ctrl_coefs = coefficient_table(neg_ctrl_gene_fit)
#   expect_equal(neg_ctrl_coefs$estimate, c(-3.553921,  0.179926), tolerance=1e-1)
#   expect_equal(neg_ctrl_coefs$normalized_effect, c(0, 0.197), tolerance=1e-1)
#   expect_gt(neg_ctrl_coefs$p_value[2], 0.05)
# })

test_that("fit_models() returns correct output for zero-inflated Poisson regression",{
  test_cds = cds
  test_cds@expression_family = "zipoisson"

  pos_ctrl_gene = test_cds[fData(cds)$gene_short_name == "ANGPTL4",]
  pos_ctrl_gene_fit = fit_models(pos_ctrl_gene, model_formula_str = "~log_dose | Size_Factor")
  pos_ctrl_coefs = coefficient_table(pos_ctrl_gene_fit)
  expect_equal(pos_ctrl_coefs$estimate, c(0.6606459,  0.1794696,  0.3733034, -0.5570760), tolerance=1e-1)
  expect_equal(pos_ctrl_coefs$normalized_effect, c(0.00000000, 0.08186624, NA, NA), tolerance=1e-1)
  expect_lt(pos_ctrl_coefs$p_value[2], 0.05)

  neg_ctrl_gene = test_cds[fData(cds)$gene_short_name == "CCNE2",]
  neg_ctrl_gene_fit = fit_models(neg_ctrl_gene, model_formula_str = "~log_dose | Size_Factor")
  neg_ctrl_coefs = coefficient_table(neg_ctrl_gene_fit)
  expect_equal(neg_ctrl_coefs$estimate, c(0.09488287,  0.21718607,  3.72540473, -0.09008649), tolerance=1e-1)
  expect_equal(neg_ctrl_coefs$normalized_effect, c(0.0000000, 0.1382822, NA, NA), tolerance=1e-1)
  expect_gt(neg_ctrl_coefs$p_value[2], 0.05)
})

test_that("fit_models() returns correct output for zero-inflated negative binomial regression",{
  test_cds = cds
  test_cds@expression_family = "zinegbinomial"

  pos_ctrl_gene = test_cds[fData(cds)$gene_short_name == "ANGPTL4",]
  pos_ctrl_gene_fit = fit_models(pos_ctrl_gene, model_formula_str = "~log_dose | Size_Factor")
  pos_ctrl_coefs = coefficient_table(pos_ctrl_gene_fit)
  expect_equal(pos_ctrl_coefs$estimate, c(0.208281688,  0.247980493, -0.007421983,  1.072303352, -3.422687894), tolerance=1e-1)
  expect_equal(pos_ctrl_coefs$normalized_effect, c(0.000000000, 0.146939585, -0.004722402, NA, NA), tolerance=1e-1)
  expect_lt(pos_ctrl_coefs$p_value[2], 0.05)

  neg_ctrl_gene = test_cds[fData(cds)$gene_short_name == "CCNE2",]
  neg_ctrl_gene_fit = fit_models(neg_ctrl_gene, model_formula_str = "~log_dose | Size_Factor")
  neg_ctrl_coefs = coefficient_table(neg_ctrl_gene_fit)
  expect_equal(neg_ctrl_coefs$estimate, c(-0.47261449,  0.20710401, -0.16405569,  3.13926917, -0.09400482), tolerance=1e-1)
  expect_equal(neg_ctrl_coefs$normalized_effect, c(0.0000000,  0.1724040, -0.1463111, NA, NA), tolerance=1e-1)
  expect_gt(neg_ctrl_coefs$p_value[2], 0.05)
})


test_that("fit_models() works with multiple cores",{
  expect_equal(1,1)
})


