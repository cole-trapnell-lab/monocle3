# Removes a bunch of auxiliary data from VGAM objects that we usually don't need.
clean_vgam_model_object = function(model) {
  attributes(model@terms$terms)$.Environment = NULL
  model@misc$formula = NULL
  model@x = matrix()
  model@y = matrix()
  model@fitted.values = matrix()
  model@residuals = matrix()
  model@qr = list()
  return(model)
}

#' @title Helper function for model fitting
#' @description test
#' @param x test
#' @param modelFormulaStr a formula string specifying the model to fit for the genes.
#' @param expression_family specifies the family function used for expression responses
#' @param relative_expr Whether to transform expression into relative values
#' @param disp_func test
#' @param verbose Whether to show VGAM errors and warnings. Only valid for cores = 1.
#' @param ... test
#' @name fit_model_helper
#' @importFrom stats formula
#' @importFrom speedglm speedglm
#' @importFrom pscl zeroinfl
#' @importFrom MASS glm.nb
#' @keywords internal
fit_model_helper <- function(x,
                             modelFormulaStr,
                             expression_family,
                             relative_expr,
                             disp_func = NULL,
                             clean_model = TRUE,
                             verbose = FALSE,
                             ...) {
  modelFormulaStr <- paste("f_expression", modelFormulaStr,
                           sep = "")
  orig_x <- x
  # FIXME: should we be using this here?
  # x <- x + pseudocount
  if (expression_family %in% c("negbinomial", "poisson", "zinegbinomial", "zipoisson", "quasipoisson")) {
    if (relative_expr) {
      x <- x / Size_Factor
    }
    f_expression <- round(x)
    # if (is.null(disp_func) == FALSE) {
    #   disp_guess <- calculate_NB_dispersion_hint(disp_func,
    #                                              round(orig_x))
    #   if (is.null(disp_guess) == FALSE && disp_guess >
    #       0 && is.na(disp_guess) == FALSE) {
    #     size_guess <- 1 / disp_guess
    #     if (expression_family== "negbinomial")
    #       expression_family <- negbinomial(isize = 1 / disp_guess, ...)
    #     else
    #       expression_family <-
    #         negbinomial.size(size = 1 / disp_guess, ...)
    #   }
    # }
  }
  else if (expression_family %in% c("gaussian", "binomial")) {
    f_expression <- x
  }
  else {
    # FIXME: maybe emit a warning or error here instead.
    f_expression <- log10(x)
  }
  f_expression = as.numeric(f_expression)
  model_formula = as.formula(modelFormulaStr)

  tryCatch({
    FM_fit = switch(expression_family,
                    "negbinomial" = MASS::glm.nb(model_formula, epsilon=1e-3),
                    "poisson" = speedglm::speedglm(model_formula, family = poisson(), acc=1e-3),
                    "quasipoisson" = speedglm::speedglm(model_formula, family = quasipoisson(), acc=1e-3),
                    "zipoisson" = pscl::zeroinfl(model_formula, dist="poisson"),
                    "zinegbinomial" = pscl::zeroinfl(model_formula, dist="negbin")
                    )
    # if (clean_model)
    #   FM_fit = clean_vgam_model_object(FM_fit)
    FM_fit
  }, error = function(e) {
    print (e)
    NA
  })
}

#' Fits a model for each gene in a CellDataSet object.
#'
#' This function fits a vector generalized additive model (VGAM) from the VGAM package for each gene in a CellDataSet.
#' By default, expression levels are modeled as smooth functions of the Pseudotime value of each
#' cell. That is, expression is a function of progress through the biological process.  More complicated formulae can be provided to account for
#' additional covariates (e.g. day collected, genotype of cells, media conditions, etc).
#'
#' This function fits a vector generalized additive model (VGAM) from the VGAM package for each gene in a CellDataSet.
#' By default, expression levels are modeled as smooth functions of the Pseudotime value of each
#' cell. That is, expression is a function of progress through the biological process.  More complicated formulae can be provided to account for
#' additional covariates (e.g. day collected, genotype of cells, media conditions, etc).
#'
#' @param cds the CellDataSet upon which to perform this operation
#' @param modelFormulaStr a formula string specifying the model to fit for the genes.
#' @param relative_expr Whether to fit a model to relative or absolute expression. Only meaningful for count-based expression data. If TRUE, counts are normalized by Size_Factor prior to fitting.
#' @param cores the number of processor cores to be used during fitting.
#' @return a tibble containing VGAM model objects
#' @importFrom qlcMatrix rowMax
#' @export
fit_models <- function(cds,
                     modelFormulaStr = "~sm.ns(Pseudotime, df=3)",
                     relative_expr = TRUE,
                     cores = 1,
                     clean_model = TRUE,
                     verbose = FALSE) {
  if (cores > 1) {
    f <-
      mc_es_apply(
        cds,
        1,
        fit_model_helper,
        required_packages = c("BiocGenerics", "Biobase", "pscl", "speedglm", "plyr", "Matrix"),
        cores = cores,
        modelFormulaStr = modelFormulaStr,
        expression_family = cds@expression_family,
        relative_expr = relative_expr,
        disp_func = cds@dispFitInfo[["blind"]]$disp_func,
        clean_model = clean_model,
        verbose = verbose
      )
    f
  } else{
    f <- smart_es_apply(
      cds,
      1,
      fit_model_helper,
      convert_to_dense = TRUE,
      modelFormulaStr = modelFormulaStr,
      expression_family = cds@expression_family,
      relative_expr = relative_expr,
      disp_func = cds@dispFitInfo[["blind"]]$disp_func,
      clean_model = clean_model,
      verbose = verbose
    )
    f
  }
  term_labels =  unlist(head(lapply(lapply(f[which(is.na(f) == FALSE)], coef), names), n =
                               1))

  M_f = as_tibble(fData(cds))
  M_f$model = f
  M_f
}

extract_coefficient_helper = function(model, term_labels, pseudo_expr =
                                        0.01) {
  if (length(intersect(class(model), c("glm", "speedglm"))) >= 1) {
    SM = summary(model)
    coef_mat = SM$coefficients # first row is intercept
    coef_mat = apply(coef_mat, 2, function(x) {as.numeric(as.character(x)) }) # We need this because some summary methods "format" the coefficients into a factor...
    row.names(coef_mat) = row.names(SM$coefficients)
     log_eff_over_int = log2((model$family$linkinv(coef_mat[, 1] + coef_mat[1, 1]) + pseudo_expr) /
                            rep(model$family$linkinv(coef_mat[1, 1]) + pseudo_expr, times = nrow(coef_mat)))
    log_eff_over_int[1] = 0
    coef_mat = as_tibble(coef_mat, rownames = "term")
    coef_mat$normalized_effect = log_eff_over_int
    coef_mat$model_component = "count"
    return (coef_mat)
  } else if (class(model) == "zeroinfl"){
    SM = summary(model)
    count_coef_mat = SM$coefficients$count # first row is intercept
    log_eff_over_int = log2((model$linkinv(count_coef_mat[, 1] + count_coef_mat[1, 1]) + pseudo_expr) /
                              rep(model$linkinv(count_coef_mat[1, 1]) + pseudo_expr, times = nrow(count_coef_mat)))
    log_eff_over_int[1] = 0
    count_coef_mat = as_tibble(count_coef_mat, rownames = "term")
    count_coef_mat$normalized_effect = log_eff_over_int
    count_coef_mat$model_component = "count"

    zero_coef_mat = SM$coefficients$zero # first row is intercept
    zero_coef_mat = as_tibble(zero_coef_mat, rownames = "term")
    zero_coef_mat$normalized_effect = NA
    zero_coef_mat$model_component = "zero"
    coef_mat = dplyr::bind_rows(count_coef_mat, zero_coef_mat)
    return (count_coef_mat)
  }else {
    coef_mat = matrix(NA, nrow = 1, ncol = 5)
    colnames(coef_mat) = c('Estimate',
                           'Std. Error',
                           'z value',
                           'Pr(>|z|)',
                           'normalized_effect')
    coef_mat = as_tibble(coef_mat)
    coef_mat$term = NA
    coef_mat$model_component = NA
    return(coef_mat)
  }
}

#' Extracts a table of coefficients from a tibble containing model objects

#' @importFrom purrr map
#' @importFrom tidyr unnest
#' @importFrom tibble as_tibble
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @export
coefficient_table <- function(model_tbl) {
  M_f = model_tbl %>%
    dplyr::mutate(terms = purrr::map(.f = purrr::possibly(extract_coefficient_helper, NA_real_), .x = model, term_labels)) %>% tidyr::unnest(terms)
  return(M_f)
}
