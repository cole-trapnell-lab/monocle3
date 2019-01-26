context("compare_models")

cds <- load_a549()

test_that("compare_models() correctly deems NB better than Poisson",{
  test_cds = cds[fData(cds)$gene_short_name == "ANGPTL4",]
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
  expect_equal(nb_comparison$lr_test_p_value[1], 0.0000228, tolerance=1e-3)

  nb_vs_zipoisson_comparison = compare_models(nb_fit, zipoisson_fit)
  # The test below should fail, as you can't compare a zipoisson to a negbinomial via lrt:
  expect_equal(nb_vs_zipoisson_comparison$lr_test_p_value[1], NA_real_)

  zinb_vs_zipoisson_comparison = compare_models(zinb_fit, zipoisson_fit)
  # The test below should fail, as you can't compare a zipoisson to a zinegbinomial via lrt:
  expect_equal(zinb_vs_zipoisson_comparison$lr_test_p_value[1], NA_real_)
})
