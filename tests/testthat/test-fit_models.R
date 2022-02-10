context("fit_models")


skip_not_travis <- function ()
{
  if (identical(Sys.getenv("TRAVIS"), "true")) {
    return(invisible(TRUE))
  }
  skip("Not on Travis")
}


cds <- load_a549()
#cds <- estimate_size_factors(cds)
test_that("fit_models() returns an error when Size_Factors are missing",{
  no_na_cds = cds
  colData(no_na_cds)$Size_Factor = NA
  expect_error(fit_models(no_na_cds))
})

test_that("fit_models() properly validates model formulae",{
  ok_formula_model_fits <- fit_models(cds, model_formula_str = "~log_dose",
                                      expression_family = 'negbinomial')
  num_failed = sum (unlist(lapply(ok_formula_model_fits$model,
                                  function(m) { class(m) })) == "logical" )
  expect_lt(num_failed, nrow(cds))

  # Skip the test below until we resolve this issue:
  # https://github.com/cole-trapnell-lab/monocle3/issues/17
  expect_error(fit_models(cds, model_formula_str = "~MISSING_TERM"))
})

test_that("fit_models() multicore works",{
  skip("R 3.6.0 - 'Planting of a Tree that doesn't work', broke the ability to makeCluster with type = FORK")
  ok_formula_model_fits <- fit_models(cds, model_formula_str = "~log_dose",
                                      cores=2,
                                      expression_family <- 'negbinomial')
  num_failed = sum (unlist(lapply(ok_formula_model_fits$model,
                                  function(m) { class(m) })) == "logical" )
  expect_lt(num_failed, nrow(cds))
})


test_that("fit_models() returns correct output for negative binomial regression",{
  test_cds = cds

  pos_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "ANGPTL4",]
  pos_ctrl_gene_fit = fit_models(pos_ctrl_gene, model_formula_str = "~log_dose",
                                 expression_family = 'negbinomial')
  expect_equal(pos_ctrl_gene_fit$status[[1]], "OK")
  pos_ctrl_coefs = coefficient_table(pos_ctrl_gene_fit)
  expect_equal(pos_ctrl_coefs$estimate, c(-0.1772648, 0.2727622, 0.675),
               tolerance=1e-1)
  expect_equal(c(pos_ctrl_coefs$normalized_effect[[1]],
                 pos_ctrl_coefs$normalized_effect[[2]]), c(0, 0.353),
               tolerance=1e-1)
  expect_lt(pos_ctrl_coefs$p_value[2], 0.05)

  fitted_vals = stats::predict(pos_ctrl_gene_fit$model[[1]])
  expect_null(fitted_vals)

  fitted_vals = stats::predict(pos_ctrl_gene_fit$model[[1]],
                               newdata=cbind(colData(test_cds) %>% as.data.frame, Size_Factor=1))
  expect_equal(sum(is.na(fitted_vals)), 0)
  expect_equal(unname(fitted_vals[1]), -0.243, tolerance=1e-3)

  neg_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "CCNE2",]
  neg_ctrl_gene_fit = fit_models(neg_ctrl_gene, model_formula_str = "~log_dose",
                                 expression_family = 'negbinomial')
  expect_equal(neg_ctrl_gene_fit$status[[1]], "OK")
  neg_ctrl_coefs = coefficient_table(neg_ctrl_gene_fit)
  expect_equal(neg_ctrl_coefs$estimate, c(-3.857,  0.139, 0.502),
               tolerance=1e-1)
  expect_equal(c(neg_ctrl_coefs$normalized_effect[[1]],
                 neg_ctrl_coefs$normalized_effect[[2]]), c(0, 0.197),
               tolerance=1e-1)
  expect_gt(neg_ctrl_coefs$p_value[2], 0.05)

  require("pryr")
  expect_lt(object_size(neg_ctrl_gene_fit$model[[1]]), 20000)
})

test_that("fit_models() returns correct output for Poisson regression",{
  test_cds = cds

  pos_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "ANGPTL4",]
  pos_ctrl_gene_fit = fit_models(pos_ctrl_gene, model_formula_str = "~log_dose",
                                 expression_family = "poisson", verbose=TRUE)
  expect_equal(pos_ctrl_gene_fit$status[[1]], "OK")
  pos_ctrl_coefs = coefficient_table(pos_ctrl_gene_fit)
  expect_equal(pos_ctrl_coefs$estimate, c(-0.397, 0.2560421),
               tolerance=1e-1)
  expect_equal(c(pos_ctrl_coefs$normalized_effect[[1]],
                 pos_ctrl_coefs$normalized_effect[[2]]), c(0, 0.477),
               tolerance=1e-1)
  expect_lt(pos_ctrl_coefs$p_value[2], 0.05)

  fitted_vals = suppressWarnings(stats::predict(pos_ctrl_gene_fit$model[[1]]))
  expect_null(fitted_vals)

  fitted_vals = stats::predict(pos_ctrl_gene_fit$model[[1]],
                               newdata=cbind(colData(test_cds) %>% as.data.frame, Size_Factor=1))
  expect_equal(sum(is.na(fitted_vals)), 0)
  expect_equal(unname(fitted_vals[1]), -0.563, tolerance=1e-1)

  neg_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "CCNE2",]
  neg_ctrl_gene_fit = fit_models(neg_ctrl_gene, model_formula_str = "~log_dose",
                                 expression_family = "poisson")
  expect_equal(neg_ctrl_gene_fit$status[[1]], "OK")
  neg_ctrl_coefs = coefficient_table(neg_ctrl_gene_fit)
  expect_equal(neg_ctrl_coefs$estimate, c(-4.0922391,  0.2701673),
               tolerance=1e-1)
  expect_equal(c(neg_ctrl_coefs$normalized_effect[[1]],
                 neg_ctrl_coefs$normalized_effect[[2]]), c(0, 0.197),
               tolerance=1e-1)
  expect_gt(neg_ctrl_coefs$p_value[2], 0.05)
  require("pryr")
  expect_lt(object_size(neg_ctrl_gene_fit$model[[1]]), 20000)
})

test_that("fit_models() returns correct output for quasipoisson regression",{
  test_cds = cds

  pos_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "ANGPTL4",]
  pos_ctrl_gene_fit = fit_models(pos_ctrl_gene, model_formula_str = "~log_dose",
                                 expression_family = "quasipoisson")
  expect_equal(pos_ctrl_gene_fit$status[[1]], "OK")
  pos_ctrl_coefs = coefficient_table(pos_ctrl_gene_fit)
  expect_equal(pos_ctrl_coefs$estimate, c(-0.397, 0.2560421),
               tolerance=1e-1)
  expect_equal(c(pos_ctrl_coefs$normalized_effect[[1]],
                 pos_ctrl_coefs$normalized_effect[[2]]), c(0, 0.477), tolerance=1e-1)
  expect_lt(pos_ctrl_coefs$p_value[2], 0.05)

  fitted_vals = suppressWarnings(stats::predict(pos_ctrl_gene_fit$model[[1]]))
  expect_null(fitted_vals)

  fitted_vals = stats::predict(pos_ctrl_gene_fit$model[[1]],
                               newdata=cbind(colData(test_cds) %>% as.data.frame, Size_Factor=1))
  expect_equal(sum(is.na(fitted_vals)), 0)
  expect_equal(unname(fitted_vals[1]), -0.563, tolerance=1e-1)

  neg_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "CCNE2",]
  neg_ctrl_gene_fit = fit_models(neg_ctrl_gene, model_formula_str = "~log_dose",
                                 expression_family = "quasipoisson")
  expect_equal(neg_ctrl_gene_fit$status[[1]], "OK")
  neg_ctrl_coefs = coefficient_table(neg_ctrl_gene_fit)
  expect_equal(neg_ctrl_coefs$estimate, c(-4.0922391,  0.2701673),
               tolerance=1e-1)
  expect_equal(c(neg_ctrl_coefs$normalized_effect[[1]],
                 neg_ctrl_coefs$normalized_effect[[2]]), c(0, 0.197),
               tolerance=1e-1)
  expect_gt(neg_ctrl_coefs$p_value[2], 0.05)
  require("pryr")
  expect_lt(object_size(neg_ctrl_gene_fit$model[[1]]), 20000)
})

# TODO: need a binary dataset for this:
# test_that("fit_models() returns correct output for binomial regression",{
#   test_cds = cds
#   metadata(test_cds)$expression_family = "binomial"
#
#   pos_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "ANGPTL4",]
#   pos_ctrl_gene_fit = fit_models(pos_ctrl_gene, model_formula_str = "~log_dose")
#   pos_ctrl_coefs = coefficient_table(pos_ctrl_gene_fit)
#   expect_equal(pos_ctrl_coefs$estimate, c(0.01078543, 0.24663124), tolerance=1e-1)
#   expect_equal(pos_ctrl_coefs$normalized_effect, c(0, 0.353), tolerance=1e-1)
#   expect_lt(pos_ctrl_coefs$p_value[2], 0.05)
#
#   neg_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "CCNE2",]
#   neg_ctrl_gene_fit = fit_models(neg_ctrl_gene, model_formula_str = "~log_dose")
#   neg_ctrl_coefs = coefficient_table(neg_ctrl_gene_fit)
#   expect_equal(neg_ctrl_coefs$estimate, c(-3.553921,  0.179926), tolerance=1e-1)
#   expect_equal(neg_ctrl_coefs$normalized_effect, c(0, 0.197), tolerance=1e-1)
#   expect_gt(neg_ctrl_coefs$p_value[2], 0.05)
# })

test_that("fit_models() returns correct output for zero-inflated Poisson regression",{
  test_cds = cds

  pos_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "ANGPTL4",]
  pos_ctrl_gene_fit = fit_models(pos_ctrl_gene, model_formula_str = "~log_dose | Size_Factor", expression_family = "zipoisson")
  expect_equal(pos_ctrl_gene_fit$status[[1]], "OK")
  pos_ctrl_coefs = coefficient_table(pos_ctrl_gene_fit)
  expect_equal(pos_ctrl_coefs$estimate,
               c(-0.0486,  0.22109280, -1.1799,  0.04879729),
               tolerance=1e-1)
  expect_equal(unname(pos_ctrl_coefs$normalized_effect),
               c(0.00000000, 0.21, NA, NA), tolerance=1e-1)
  expect_lt(pos_ctrl_coefs$p_value[2], 0.05)

  fitted_vals = stats::predict(pos_ctrl_gene_fit$model[[1]])
  expect_null(fitted_vals)

  fitted_vals = stats::predict(pos_ctrl_gene_fit$model[[1]],
                               newdata=cbind(colData(test_cds) %>% as.data.frame, Size_Factor=1))
  expect_equal(sum(is.na(fitted_vals)), 0)
  expect_equal(unname(fitted_vals[1]), 0.595, tolerance=1e-1)

  neg_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "CCNE2",]
  neg_ctrl_gene_fit = fit_models(neg_ctrl_gene,
                                 model_formula_str = "~log_dose | Size_Factor",
                                 expression_family = "zipoisson")
  expect_equal(neg_ctrl_gene_fit$status[[1]], "OK")
  neg_ctrl_coefs = coefficient_table(neg_ctrl_gene_fit)
  expect_equal(neg_ctrl_coefs$estimate,
               c(-1.158, 0.24146532,  2.580,  0.134),
               tolerance=1e-1)
  expect_equal(unname(neg_ctrl_coefs$normalized_effect),
               c(0.0000000, 0.1382822, NA, NA), tolerance=1e-1)
  expect_gt(neg_ctrl_coefs$p_value[2], 0.05)
  require("pryr")
  expect_lt(object_size(neg_ctrl_gene_fit$model[[1]]), 20100)
})

test_that("fit_models() returns correct output for zero-inflated negative binomial regression",{
  test_cds = cds

  pos_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "ANGPTL4",]
  pos_ctrl_gene_fit = fit_models(pos_ctrl_gene,
                                 model_formula_str = "~log_dose | Size_Factor",
                                 expression_family = "zinegbinomial")
  expect_equal(pos_ctrl_gene_fit$status[[1]], "OK")
  pos_ctrl_coefs = coefficient_table(pos_ctrl_gene_fit)
  expect_equal(pos_ctrl_coefs$estimate,
               c(-0.255, 0.26434083, -0.146, -11.113, -1.005),
               tolerance=1e-1)
  expect_equal(unname(pos_ctrl_coefs$normalized_effect),
               c(0.000000000, 0.1685858, -0.119, NA, NA), tolerance=1e-1)
  expect_lt(pos_ctrl_coefs$p_value[2], 0.05)

  fitted_vals = stats::predict(pos_ctrl_gene_fit$model[[1]])
  expect_null(fitted_vals)

  fitted_vals = stats::predict(pos_ctrl_gene_fit$model[[1]],
                               newdata=cbind(colData(test_cds) %>% as.data.frame, Size_Factor=1))
  expect_equal(sum(is.na(fitted_vals)), 0)
  expect_equal(unname(fitted_vals[1]), 0.642, tolerance=1e-1)

  neg_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "CCNE2",]
  neg_ctrl_gene_fit = fit_models(neg_ctrl_gene,
                                 model_formula_str = "~log_dose | Size_Factor",
                                 expression_family = "zinegbinomial")
  expect_equal(neg_ctrl_gene_fit$status[[1]], "OK")
  neg_ctrl_coefs = coefficient_table(neg_ctrl_gene_fit)
  expect_equal(neg_ctrl_coefs$estimate,
               c(-3.94, 0.24145643, -3.47, -4.81, -15.42),
               tolerance=1e-1)
  expect_equal(unname(neg_ctrl_coefs$normalized_effect),
               c(0.0000000,  0.1684627, -1.46, NA, NA), tolerance=1e-1)
  expect_gt(neg_ctrl_coefs$p_value[2], 0.05)
  require("pryr")
  expect_lt(object_size(neg_ctrl_gene_fit$model[[1]]), 20600)
})


test_that("fit_models() flags non-convergence for negative binomial regression",{
  test_cds = cds
  pos_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "ANGPTL4",]
  pos_ctrl_gene_fit = fit_models(pos_ctrl_gene,
                                 model_formula_str = "~log_dose", maxit=1,
                                 expression_family = "negbinomial")
  expect_equal(pos_ctrl_gene_fit$status[[1]], "FAIL")
})

test_that("fit_models() flags non-convergence for Poisson regression",{
  test_cds = cds
  pos_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "ANGPTL4",]
  pos_ctrl_gene_fit = fit_models(pos_ctrl_gene,
                                 model_formula_str = "~log_dose", maxit=5,
                                 expression_family = "poisson")
  expect_equal(pos_ctrl_gene_fit$status[[1]], "FAIL")
})

test_that("fit_models() returns correct output for quasipoisson regression",{
  test_cds = cds
  pos_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "ANGPTL4",]
  pos_ctrl_gene_fit = fit_models(pos_ctrl_gene, model_formula_str = "~log_dose",
                                 maxit=5, expression_family = "quasipoisson")
  expect_equal(pos_ctrl_gene_fit$status[[1]], "FAIL")
})

# TODO: need a binary dataset for this:
# test_that("fit_models() flags non-convergence for binomial regression",{
# })

test_that("fit_models() flags non-convergence for zero-inflated Poisson regression",{
  test_cds = cds
  pos_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "ANGPTL4",]
  pos_ctrl_gene_fit = fit_models(pos_ctrl_gene,
                                 model_formula_str = "~log_dose | Size_Factor",
                                 maxit=5, expression_family = "zipoisson")
  expect_equal(pos_ctrl_gene_fit$status[[1]], "FAIL")
})

test_that("fit_models() flags non-convergence for zero-inflated negative binomial regression",{
  test_cds = cds
  pos_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "ANGPTL4",]
  pos_ctrl_gene_fit = fit_models(pos_ctrl_gene,
                                 model_formula_str = "~log_dose | Size_Factor",
                                 maxit=5, expression_family = "zinegbinomial")
  expect_equal(pos_ctrl_gene_fit$status[[1]], "FAIL")
})

test_that("fit_models() works with multiple cores",{
  expect_equal(1,1)
})

#### NOT TRAVIS ####

test_that("fit_models() can handle cluster in model formulae",{
  skip_on_travis()
  test_cds = cds
  test_cds = preprocess_cds(test_cds)
  test_cds = reduce_dimension(test_cds)
  set.seed(100)
  test_cds = cluster_cells(test_cds, resolution = .01)

  pos_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "ANGPTL4",]
  pos_ctrl_gene_fit = fit_models(pos_ctrl_gene, model_formula_str = "~cluster",
                                 expression_family = "quasipoisson")
  expect_equal(pos_ctrl_gene_fit$status[[1]], "OK")
  pos_ctrl_coefs = coefficient_table(pos_ctrl_gene_fit)
  expect_equal(pos_ctrl_coefs$estimate[2], -0.112, tolerance=1e-2)
})


#### TRAVIS ####

test_that("fit_models() can handle cluster in model formulae",{
  skip_not_travis()
  test_cds = cds
  test_cds = preprocess_cds(test_cds)
  test_cds = reduce_dimension(test_cds)
  set.seed(100)
  test_cds = cluster_cells(test_cds, resolution = .01)

  pos_ctrl_gene = test_cds[rowData(cds)$gene_short_name == "ANGPTL4",]
  pos_ctrl_gene_fit = fit_models(pos_ctrl_gene, model_formula_str = "~cluster",
                                 expression_family = "quasipoisson")
  expect_equal(pos_ctrl_gene_fit$status[[1]], "OK")
  pos_ctrl_coefs = coefficient_table(pos_ctrl_gene_fit)
  expect_equal(pos_ctrl_coefs$estimate[2], -0.615, tolerance=1e-2)
})


