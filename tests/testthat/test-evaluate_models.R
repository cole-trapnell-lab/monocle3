context("evaluate_models")

cds <- load_a549()
cds <- estimate_size_factors(cds)
test_that("evaluate_models() returns correct output for poisson models",{
  test_cds = cds[rowData(cds)$gene_short_name == "ANGPTL4",]
  fit_m = fit_models(test_cds, model_formula_str = "~log_dose", expression_family = "poisson")
  evaluated_fit = suppressWarnings(evaluate_fits(fit_m))
  expect_equal(evaluated_fit$null_deviance, 1135, tolerance=1e-3)
  expect_equal(evaluated_fit$df_null, 499)
  expect_equal(evaluated_fit$logLik, -803, tolerance=1e-3)
  expect_equal(evaluated_fit$AIC, 1610, tolerance=1e-3)
  expect_equal(evaluated_fit$BIC, NA_real_)
  expect_equal(evaluated_fit$deviance, 999, tolerance=1e-3)
  expect_equal(evaluated_fit$df_residual, 498)
})

test_that("evaluate_models() returns correct output for quasipoisson models",{
  test_cds = cds[rowData(cds)$gene_short_name == "ANGPTL4",]
  fit_m = fit_models(test_cds, model_formula_str = "~log_dose", expression_family = "quasipoisson")
  evaluated_fit = suppressWarnings(evaluate_fits(fit_m))
  expect_equal(evaluated_fit$null_deviance, 1135, tolerance=1e-3)
  expect_equal(evaluated_fit$df_null, 499)
  expect_equal(evaluated_fit$logLik, NA)
  expect_equal(evaluated_fit$AIC, NA_real_)
  expect_equal(evaluated_fit$BIC, NA_real_)
  expect_equal(evaluated_fit$deviance, 999, tolerance=1e-3)
  expect_equal(evaluated_fit$df_residual, 498)
})

### FIXME: Need cases for binomial (binary) data

test_that("evaluate_models() returns correct output for negative binomial models",{
  test_cds = cds[rowData(cds)$gene_short_name == "ANGPTL4",]
  fit = fit_models(test_cds, model_formula_str = "~log_dose", expression_family = "negbinomial")
  evaluated_fit = suppressWarnings(evaluate_fits(fit))
  expect_equal(evaluated_fit$null_deviance, 497, tolerance=1e-3)
  expect_equal(evaluated_fit$df_null, 499)
  expect_equal(evaluated_fit$logLik, -692, tolerance=1e-3)
  expect_equal(evaluated_fit$AIC, 1389, tolerance=1e-3)
  expect_equal(evaluated_fit$BIC, 1389, tolerance=1e-3)
  expect_equal(evaluated_fit$deviance, 467, tolerance=1e-3)
  expect_equal(evaluated_fit$df_residual, 498)
})

test_that("evaluate_models() returns correct output for zero-inflated poisson models",{
  test_cds = cds[rowData(cds)$gene_short_name == "ANGPTL4",]
  fit = fit_models(test_cds, model_formula_str = "~log_dose", expression_family = "zipoisson")
  evaluated_fit = suppressWarnings(evaluate_fits(fit))
  expect_equal(evaluated_fit$null_deviance, NA_real_)
  expect_equal(evaluated_fit$df_null, 498)
  expect_equal(evaluated_fit$logLik, -748, tolerance=1e-3)
  expect_equal(evaluated_fit$AIC, 1504, tolerance=1e-3)
  expect_equal(evaluated_fit$BIC, 1504, tolerance=1e-3)
  expect_equal(evaluated_fit$deviance, NA_real_)
  expect_equal(evaluated_fit$df_residual, 496)
})


test_that("evaluate_models() returns correct output for zero-inflated negative binomial models",{
  test_cds = cds[rowData(cds)$gene_short_name == "ANGPTL4",]
  fit = fit_models(test_cds, model_formula_str = "~log_dose", expression_family = "zinegbinomial")
  evaluated_fit = suppressWarnings(evaluate_fits(fit))
  expect_equal(evaluated_fit$null_deviance, NA_real_)
  expect_equal(evaluated_fit$df_null, 498)
  expect_equal(evaluated_fit$logLik, -689, tolerance=1e-3)
  expect_equal(evaluated_fit$AIC, 1387, tolerance=1e-3)
  expect_equal(evaluated_fit$BIC, 1387, tolerance=1e-3)
  expect_equal(evaluated_fit$deviance, NA_real_)
  expect_equal(evaluated_fit$df_residual, 495)
})
