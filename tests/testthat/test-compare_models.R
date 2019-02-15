context("compare_models")
cds <- load_a549()

test_that("compare_models() correctly deems NB better than Poisson",{
  test_cds = cds[rowData(cds)$gene_short_name == "ANGPTL4",]
  zinb_cds = test_cds
  zinb_cds@expression_family = "zinegbinomial"
  zinb_fit = fit_models(zinb_cds, model_formula_str = "~log_dose")

  zipoisson_cds = test_cds
  zipoisson_cds@expression_family = "zipoisson"
  zipoisson_fit = fit_models(zipoisson_cds, model_formula_str = "~log_dose")

  nb_cds = test_cds
  nb_cds@expression_family = "negbinomial"
  nb_fit = fit_models(nb_cds, model_formula_str = "~log_dose")
  nb_reduced_fit = fit_models(nb_cds, model_formula_str = "~1")
  nb_comparison = compare_models(nb_fit, nb_reduced_fit)
  expect_equal(nb_comparison$p_value[1], 0.0000228, tolerance=1e-3)

  nb_fit = fit_models(nb_cds, model_formula_str = "~log_dose", clean_model = FALSE)
  nb_reduced_fit = fit_models(nb_cds, model_formula_str = "~1", clean_model = FALSE)
  require("lmtest") # TODO: skip the test below if lmtest isn't available
  lmtest_lrt_pval = lmtest::lrtest(nb_fit$model[[1]], nb_reduced_fit$model[[1]])[2,5]
  expect_equal(nb_comparison$p_value[1], lmtest_lrt_pval)

  skip("currently failing")
  nb_vs_zipoisson_comparison = compare_models(nb_fit, zipoisson_fit)
  # The function below should return NA, as you can't compare a zipoisson to a negbinomial via lrt:
  expect_equal(nb_vs_zipoisson_comparison$p_value[1], NA_real_)

  zinb_vs_zipoisson_comparison = compare_models(zinb_fit, zipoisson_fit)
  # The function below should return NA, as you can't compare a zipoisson to a zinegbinomial via lrt:
  expect_equal(zinb_vs_zipoisson_comparison$p_value[1], NA_real_)
})

