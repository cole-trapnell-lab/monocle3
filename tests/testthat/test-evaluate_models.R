context("evaluate_models")

cds <- load_a549()

test_that("evaluate_models() returns correct output for negative binomial models",{
  test_cds = cds[fData(cds)$gene_short_name == "ANGPTL4",]
  test_cds@expression_family = "negbinomial"
  nb_fit = fit_models(test_cds, model_formula_str = "~log_dose")
  nb_reduced_fit = fit_models(test_cds, model_formula_str = "~1")
  evaluated_fit = evaluate_fits(nb_fit)
  expect_equal(evaluated_fit$null.deviance, 486, tolerance=1e-3)
  expect_equal(evaluated_fit$df.null, 499, tolerance=1e-3)
  expect_equal(evaluated_fit$logLik, -735., tolerance=1e-3)
  expect_equal(evaluated_fit$AIC, 1477, tolerance=1e-3)
  expect_equal(evaluated_fit$BIC, 1489, tolerance=1e-3)
  expect_equal(evaluated_fit$deviance, 468, tolerance=1e-3)
  expect_equal(evaluated_fit$df.residual, 498, tolerance=1e-3)
})

