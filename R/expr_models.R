# Removes a bunch of auxiliary data from VGAM objects that we usually don't
# need.
clean_vgam_model_object <- function(model) {
  attributes(model@terms$terms)$.Environment <- NULL
  model@misc$formula = NULL
  model@x <- matrix()
  model@y <- matrix()
  model@fitted.values <- matrix()
  model@residuals <- matrix()
  model@qr <- list()
  return(model)
}

clean_mass_model_object <- function(cm) {
  cm$y <- c()
  #cm$model = c()
  cm$residuals <- c()
  cm$fitted.values <- c()
  cm$effects <- c()
  cm$qr$qr <- c()
  cm$linear.predictors <- c()
  cm$weights <- c()
  cm$prior.weights <- c()
  cm$data <- c()
  
  cm$family$variance <- c()
  cm$family$dev.resids <- c()
  cm$family$aic <- c()
  cm$family$validmu <- c()
  cm$family$simulate <- c()
  cm$family$control <- c()
  attr(cm$terms,".Environment") <- c()
  attr(cm$formula,".Environment") <- c()
  return(cm)
}

clean_speedglm_model_object = function(cm) {
  cm$y = c()
  #cm$model = c()
  
  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()
  
  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$aic = c()
  cm$family$validmu = c()
  cm$family$simulate = c()
  cm$family$control = c()
  attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()
  return(cm)
}


clean_zeroinfl_model_object = function(cm) {
  cm$y = c()
  cm$model = c()
  
  cm$residuals = c()
  cm$fitted.values = c()
  #cm$effects = c()
  #cm$qr$qr = c()
  #cm$linear.predictors = c()
  cm$weights = c()
  #cm$prior.weights = c()
  cm$data = c()
  
  
  # cm$family$variance = c()
  # cm$family$dev.resids = c()
  # cm$family$aic = c()
  # cm$family$validmu = c()
  # cm$family$simulate = c()
  # cm$family$control = c()
  attr(cm$terms$count,".Environment") = c()
  attr(cm$terms$zero,".Environment") = c()
  attr(cm$terms$full,".Environment") = c()
  attr(cm$formula,".Environment") = c()
  return(cm)
}

clean_glmerMod_model_object <- function(model) {
  rcl="glmResp"
  #trivial.y=FALSE
  #fac <- as.numeric(rcl != "nlsResp")
  model@devcomp$cmp <- c(model@devcomp$cmp, 
                         logLik=logLik(model),
                         AIC=AIC(logLik(model)),
                         BIC=BIC(logLik(model)),
                         df_residual = df.residual(model))
  cm <- new(switch(rcl, lmerResp="lmerMod", glmResp="glmerMod", nlsResp="nlmerMod"),
            call=model@call, 
            #frame=model@frame, # takes up quite a bit of space, but messes up evaluate_fits if we delete
            flist=model@flist, 
            cnms=model@cnms,
            Gp=model@Gp, 
            theta=model@theta, 
            beta=model@beta,
            #u=if (trivial.y) rep(NA_real_,nrow(model@pp$Zt)) else model@pp$u(fac),
            #lower=model@lower, 
            devcomp=model@devcomp,
            pp=model@pp, 
            resp=model@resp
            
  )
  return(cm)
}


clean_glmerMod_summary_object <- function(model_summary) {
  model_summary$residuals <- c()
  return(model_summary)
}


clean_model_object = function(model) {
  if (class(model)[1] == "negbin") {
    model = clean_mass_model_object(model)
  } else if (length(intersect(class(model), c("speedglm"))) >= 1) {
    model = clean_speedglm_model_object(model)
  } else if (class(model)[1] == "zeroinfl"){
    model = clean_zeroinfl_model_object(model)
  } else if (class(model)[1] == "glmerMod") {
    model = clean_glmerMod_model_object(model) 
  } else {
    stop("Unrecognized model class")
  }
}

#' @title Helper function for model fitting
#' @description test
#' @param x test
#' @param model_formula_str a formula string specifying the model to fit for
#'   the genes.
#' @param expression_family specifies the family function used for expression
#'   responses
#' @param disp_func test
#' @param verbose Whether to show VGAM errors and warnings. Only valid for
#'   cores = 1.
#' @param ... test
#' @name fit_model_helper
#' @keywords internal
fit_model_helper <- function(x,
                             model_formula_str,
                             expression_family,
                             disp_func = NULL,
                             clean_model = TRUE,
                             verbose = FALSE,
                             ...) {
  model_formula_str <- paste("f_expression", model_formula_str,
                             sep = "")
  orig_x <- x
  # FIXME: should we be using this here?
  # x <- x + pseudocount
  if (expression_family %in% c("negbinomial", "poisson", "zinegbinomial",
                               "zipoisson", "quasipoisson", "mixed-negbinomial")) {
    x <- x / Size_Factor
    f_expression <- round(x)
  }
  else if (expression_family %in% c("binomial", "gaussian")) {
    f_expression <- x
  }
  else {
    # FIXME: maybe emit a warning or error here instead.
    f_expression <- log10(x)
  }
  f_expression = as.numeric(f_expression)
  model_formula = stats::as.formula(model_formula_str)
  
  tryCatch({
    if (verbose) messageWrapper = function(expr) { expr }
    else messageWrapper = suppressWarnings
    
    FM_fit = messageWrapper(switch(expression_family,
                                   "negbinomial" = MASS::glm.nb(model_formula, epsilon=1e-3,
                                                                model=FALSE, y=FALSE, ...),
                                   "poisson" = speedglm::speedglm(model_formula,
                                                                  family = stats::poisson(),
                                                                  acc=1e-3, model=FALSE,
                                                                  y=FALSE, ...),
                                   "quasipoisson" = speedglm::speedglm(model_formula,
                                                                       family = stats::quasipoisson(),
                                                                       acc=1e-3, model=FALSE,
                                                                       y=FALSE, ...),
                                   "binomial" = speedglm::speedglm(model_formula,
                                                                   family = stats::binomial(),
                                                                   acc=1e-3, model=FALSE,
                                                                   y=FALSE, ...),
                                   "gaussian" = speedglm::speedglm(model_formula,
                                                                   family = stats::gaussian(),
                                                                   acc=1e-3, model=FALSE,
                                                                   y=FALSE, ...),
                                   "zipoisson" = pscl::zeroinfl(model_formula,
                                                                dist="poisson", ...),
                                   "zinegbinomial" = pscl::zeroinfl(model_formula,
                                                                    dist="negbin", ...),
                                   "mixed-negbinomial" = glmer.nb(model_formula, 
                                                                    nAGQ=0,
                                                                    control=glmerControl(optimizer = "nloptwrap"),
                                                                    ...)
    ))
    FM_summary = summary(FM_fit)
    if (clean_model){
      FM_fit = clean_model_object(FM_fit) 
      if (class(FM_fit)[1] == "glmerMod"){
        FM_summary = clean_glmerMod_summary_object(FM_summary)
      }
    }
    df = list(model=FM_fit, model_summary=FM_summary)
    df
  }, error = function(e) {
    if (verbose) { print (e) }
    list(model=NA, model_summary=NA)
  })
}

#' Fits a model for each gene in a cell_data_set object.
#'
#' This function fits a generalized linear model for each gene in a
#' cell_data_set. Formulae can be provided to account for additional covariates
#' (e.g. day collected, genotype of cells, media conditions, etc).
#'
#' @param cds The cell_data_set upon which to perform this operation.
#' @param model_formula_str A formula string specifying the model to fit for
#'   the genes.
#' @param expression_family Specifies the family function used for expression
#'   responses. Can be one of "quasipoisson", "negbinomial", "poisson",
#'   "binomial", "gaussian", "zipoisson", "zinegbinomial", or "mixed-negbinomial".
#'   Default is "quasipoisson".
#' @param reduction_method Which method to use with clusters() and
#'   partitions(). Default is "UMAP".
#' @param cores The number of processor cores to use during fitting.
#' @param clean_model Logical indicating whether to clean the model. Default is
#'   TRUE.
#' @param verbose Logical indicating whether to emit progress messages.
#' @param ... Additional arguments passed to model fitting functions.
#'
#' @return a tibble where the rows are genes and columns are
#'   * id character vector from `rowData(cds)$id`
#'   * gene_short_names character vector from `rowData(cds)$gene_short_names`
#'   * num_cells_expressed int vector from `rowData(cds)$num_cells_expressed`
#'   * gene_id character vector from row.names(rowData(cds))`
#'   * model GLM model list returned by speedglm
#'   * model_summary model summary list returned by `summary(model)`
#'   * status character vector of model fitting status: OK when model converged, otherwise FAIL
#'
#' @export
fit_models <- function(cds,
                       model_formula_str,
                       expression_family = "quasipoisson",
                       reduction_method="UMAP",
                       cores = 1,
                       clean_model = TRUE,
                       verbose = FALSE,
                       ...) {
  
  model_form <- stats::as.formula(model_formula_str)
  if (!"num_cells_expressed" %in% names(rowData(cds))) {
    cds <- detect_genes(cds)
  }
  coldata_df = colData(cds)
  tryCatch({
    coldata_df$cluster = clusters(cds, reduction_method)[colnames(cds)]
    coldata_df$partition = partitions(cds, reduction_method)[colnames(cds)]
    coldata_df$pseudotime = pseudotime(cds, reduction_method)
  }, error = function(e) {} )
  
  # Test model formula validity.
  # Notes:
  #  o  allow for formulas in model formula: e.g.,  ~ splines::ns(pseudotime, df=3)
  #  o  watch for NA, NaN, Inf in terms
  #  o  is.na counts NA and NaN
  #  o  splines::ns(pseudotime, df=3) fails if there are Inf values, at least,
  #     which causes model.frame( model_formula, ...) to fail
  #  o  model.frame catches mis-spelled functions
  err_msg <- NULL
  mf_terms <- all.vars(lme4::subbars(model_form))
  for( mf_term in mf_terms )
  {
    if(!( mf_term %in% names(coldata_df)))
    {
      err_msg <- paste0(err_msg,'  \'', mf_term, '\': not in cds\n')
      next
    }
  }
  if(length(err_msg) > 0)
    stop( '\n-- bad fit_models terms --\n', err_msg )
  tryCatch({
    stats::model.frame(lme4::subbars(model_form), data=coldata_df)
  }, error = function( cnd ) {
    info_msg <- ''
    for( mf_term in mf_terms )
    {
      mf_length  <- length(coldata_df[[mf_term]])
      mf_num_inf <- sum(is.infinite(coldata_df[[mf_term]]))
      mf_num_nan <- sum(is.nan(coldata_df[[mf_term]]))
      mf_num_na  <- sum(is.na(coldata_df[[mf_term]]))
      if( mf_num_inf > 0 )
        info_msg <- paste0(info_msg, '  \'', mf_term, '\': ' , mf_num_inf, ' of ', mf_length, ' values are Inf\n')
      if( mf_num_nan > 0 )
        info_msg <- paste0(info_msg, '  \'', mf_term, '\': ' , mf_num_nan, ' of ', mf_length, ' values are NaN\n')
      if( mf_num_na - mf_num_nan > 0 )
        info_msg <- paste0(info_msg, '  \'', mf_term, '\': ' , mf_num_na - mf_num_nan, ' of ', mf_length, ' values are NA\n')
    }
    rm( mf_term, mf_length, mf_num_inf, mf_num_nan, mf_num_na )
    stop (paste0( 'Error in model formula: ', conditionMessage( cnd ), '\n', info_msg ) )
  })
  
  disp_func <- NULL
  
  if (cores > 1) {
    fits <-
      mc_es_apply(
        cds,
        1,
        fit_model_helper,
        required_packages = c("BiocGenerics", "Biobase", "MASS", "purrr",
                              "pscl", "speedglm", "dplyr", "Matrix", "lme4"),
        cores = cores,
        reduction_method = reduction_method,
        model_formula_str = model_formula_str,
        expression_family = expression_family ,
        disp_func = disp_func,
        clean_model = clean_model,
        verbose = verbose,
        ...
      )
    fits
  } else{
    fits <- smart_es_apply(
      cds,
      1,
      fit_model_helper,
      convert_to_dense = TRUE,
      model_formula_str = model_formula_str,
      expression_family = expression_family,
      reduction_method = reduction_method,
      disp_func = disp_func,
      clean_model = clean_model,
      verbose = verbose,
      ...
    )
    fits
  }
  
  rowData(cds)$gene_id <- row.names(rowData(cds))
  fits <- tibble::as_tibble(purrr::transpose(fits))
  M_f <- tibble::as_tibble(rowData(cds))
  M_f <- dplyr::bind_cols(M_f, fits)
  M_f <- M_f %>%
    dplyr::mutate(status = purrr::map(.f = purrr::possibly(
      extract_model_status_helper, NA_real_), .x = model)) %>%
    tidyr::unnest(status)
  return(M_f)
}

extract_model_status_helper <- function(model){
  if (class(model)[1] == "speedglm") {
    status_str <- ifelse(model$convergence, "OK", "FAIL")
    return (status_str)
    
  } else if (class(model)[1] == "negbin"){
    status_str <- ifelse(model$converged, "OK", "FAIL")
    return (status_str)
  } else if (class(model)[1] == "zeroinfl"){
    status_str <- ifelse(model$converged, "OK", "FAIL")
    return (status_str)
  }else if (class(model)[1] == "glmerMod"){
    status_str <- ifelse(any( grepl("failed to converge", model@optinfo$conv$lme4$messages) ), "FAIL", "OK") # check to see if any warning messages
    return (status_str)
  } else {
    return("FAIL")
  }
}

extract_coefficient_helper = function(model, model_summary,
                                      pseudo_count = 0.01) {
  if (class(model)[1] == "speedglm") {
    coef_mat <- model_summary$coefficients # first row is intercept
    # We need this because some summary methods "format" the coefficients into
    # a factor...
    coef_mat <- apply(coef_mat, 2, function(x) {as.numeric(as.character(x)) })
    row.names(coef_mat) = row.names(model_summary$coefficients)
    colnames(coef_mat) = c('estimate',
                           'std_err',
                           'test_val',
                           'p_value')
    log_eff_over_int = log2((model$family$linkinv(coef_mat[, 1] +
                                                    coef_mat[1, 1]) +
                               pseudo_count) /
                              rep(model$family$linkinv(coef_mat[1, 1]) +
                                    pseudo_count, times = nrow(coef_mat)))
    log_eff_over_int[1] = 0
    coef_mat = tibble::as_tibble(coef_mat, rownames = "term")
    coef_mat$normalized_effect = log_eff_over_int
    coef_mat$model_component = "count"
    return (coef_mat)
    
  } else if (class(model)[1] == "negbin"){
    coef_mat = model_summary$coefficients # first row is intercept
    # We need this because some summary methods "format" the coefficients into
    # a factor...
    coef_mat = apply(coef_mat, 2, function(x) {as.numeric(as.character(x)) })
    row.names(coef_mat) = row.names(model_summary$coefficients)
    colnames(coef_mat) = c('estimate',
                           'std_err',
                           'test_val',
                           'p_value')
    log_eff_over_int = log2((model$family$linkinv(coef_mat[, 1] +
                                                    coef_mat[1, 1]) +
                               pseudo_count) /
                              rep(model$family$linkinv(coef_mat[1, 1]) +
                                    pseudo_count, times = nrow(coef_mat)))
    log_eff_over_int[1] = 0
    coef_mat = tibble::as_tibble(coef_mat, rownames = "term")
    coef_mat$normalized_effect = log_eff_over_int
    
    coef_mat$model_component = "count"
    return (coef_mat)
  } else if (class(model)[1] == "zeroinfl"){
    count_coef_mat = model_summary$coefficients$count # first row is intercept
    colnames(count_coef_mat) = c('estimate',
                                 'std_err',
                                 'test_val',
                                 'p_value')
    log_eff_over_int = log2((model$linkinv(count_coef_mat[, 1] +
                                             count_coef_mat[1, 1]) +
                               pseudo_count) /
                              rep(model$linkinv(count_coef_mat[1, 1]) +
                                    pseudo_count,
                                  times = nrow(count_coef_mat)))
    log_eff_over_int[1] = 0
    count_coef_mat = tibble::as_tibble(count_coef_mat, rownames = "term")
    count_coef_mat$normalized_effect = log_eff_over_int
    count_coef_mat$model_component = "count"
    
    zero_coef_mat = model_summary$coefficients$zero # first row is intercept
    colnames(zero_coef_mat) = c('estimate',
                                'std_err',
                                'test_val',
                                'p_value')
    zero_coef_mat = tibble::as_tibble(zero_coef_mat, rownames = "term")
    zero_coef_mat$normalized_effect = NA
    zero_coef_mat$model_component = "zero"
    coef_mat = dplyr::bind_rows(count_coef_mat, zero_coef_mat)
    return (coef_mat)
  } else if (class(model)[1] == "glmerMod"){
    
    coef_mat = model_summary$coefficients # first row is intercept
    # We need this because some summary methods "format" the coefficients into
    # a factor...
    coef_mat = apply(coef_mat, 2, function(x) {as.numeric(as.character(x)) })
    row.names(coef_mat) = row.names(model_summary$coefficients)
    colnames(coef_mat) = c('estimate',
                           'std_err',
                           'test_val',
                           'p_value')
    log_eff_over_int = log2((exp(coef_mat[, 1] + # inverse-link function for neg-binomial is exp, done manually bc no inv-link function in merMod
                                   coef_mat[1, 1]) +
                               pseudo_count) /
                              rep(exp(coef_mat[1, 1]) +
                                    pseudo_count, times = nrow(coef_mat)))
    log_eff_over_int[1] = 0
    coef_mat = tibble::as_tibble(coef_mat, rownames = "term")
    coef_mat$normalized_effect = log_eff_over_int
    
    coef_mat$model_component = "count"
    return (coef_mat)
    
  } else {
    coef_mat = matrix(NA_real_, nrow = 1, ncol = 5)
    colnames(coef_mat) = c('estimate',
                           'std_err',
                           'test_val',
                           'p_value',
                           'normalized_effect')
    coef_mat = tibble::as_tibble(coef_mat)
    coef_mat$term = NA_character_
    coef_mat$model_component = NA_character_
    return(coef_mat)
  }
}

#' Extracts a table of coefficients from a tibble containing model objects
#'
#' @param model_tbl A tibble of model objects, generally the output of
#'   \code{\link{fit_models}}.
#' @importFrom dplyr %>%
#' @export
coefficient_table <- function(model_tbl) {
  M_f = model_tbl %>%
    dplyr::mutate(terms = purrr::map2(.f = purrr::possibly(
      extract_coefficient_helper, NA_real_), .x = model,
      .y = model_summary)) %>%
    tidyr::unnest(terms)
  M_f = M_f %>% dplyr::group_by(model_component, term) %>%
    dplyr::mutate(q_value = stats::p.adjust(p_value)) %>% dplyr::ungroup()
  return(M_f)
}

#' Compares goodness of fit for two ways of fitting a set of genes' expression
#'
#' @param model_tbl_full A tibble of model objects, generally output of
#'   \code{\link{fit_models}}, to be compared with \code{model_tbl_reduced}
#' @param model_tbl_reduced A tibble of model objects, generally output of
#'   \code{\link{fit_models}}, to be compared with \code{model_tbl_full}.
#'
#' @export
compare_models <- function(model_tbl_full, model_tbl_reduced){
  model_x_eval <- evaluate_fits(model_tbl_full)
  model_y_eval <- evaluate_fits(model_tbl_reduced)
  if ("gene_short_name" %in% names(model_x_eval) &
      "gene_short_name" %in% names(model_y_eval)) {
    joined_fits <- dplyr::full_join(model_x_eval, model_y_eval,
                                    by=c("gene_id", "gene_short_name",
                                         "num_cells_expressed"))
    
  } else {
    joined_fits <- dplyr::full_join(model_x_eval, model_y_eval,
                                    by=c("gene_id", "num_cells_expressed"))
  }
  
  joined_fits <- joined_fits %>% dplyr::mutate(
    dfs = round(abs(df_residual.x  - df_residual.y)),
    LLR = 2 * abs(logLik.x  - logLik.y),
    p_value = stats::pchisq(LLR, dfs, lower.tail = FALSE)
  )
  
  joined_fits <- joined_fits %>%
    dplyr::select_if(names(.) %in% c("gene_id", "gene_short_name",
                                     "num_cells_expressed", "p_value"))
  joined_fits$q_value <- stats::p.adjust(joined_fits$p_value)
  # joined_fits = joined_fits %>%
  #               dplyr::mutate(lr_test_p_value = purrr::map2_dbl(
  #                             model_summary.x, model_summary.x,
  #                   function(x, y) {
  #                     if (identical(class(x), class(y))){
  #                       likelihood_ratio_test_pval(x,y)
  #                     } else {
  #                       NA
  #                     }
  #
  #                     } ) )
  return (joined_fits)
}


#' Evaluate the fits of model objects.
#'
#' @importFrom dplyr %>%
#' @param model_tbl A tibble of model objects, generally output of
#'   \code{\link{fit_models}}.
#' @export
evaluate_fits <- function(model_tbl){
  private_glance <- function(m){
    
    
    mass_glance <- function(m) {
      tibble::tibble(null_deviance = m$null.deviance,
                     df_null = m$df.null,
                     logLik = as.numeric(stats::logLik(m)),
                     AIC = stats::AIC(m),
                     BIC = stats::AIC(m),
                     deviance = m$deviance,
                     df_residual = m$df.residual)
    }
    
    zeroinfl_glance <- function(m) {
      tibble::tibble(null_deviance = NA_real_,
                     df_null = m$df.null,
                     logLik = as.numeric(stats::logLik(m)),
                     AIC = stats::AIC(m),
                     BIC = stats::AIC(m),
                     deviance = NA_real_,
                     df_residual = m$df.residual)
    }
    speedglm_glance <- function(m) {
      tibble::tibble(null_deviance = m$nulldev,
                     df_null = m$nulldf,
                     logLik = m$logLik,
                     AIC = m$aic,
                     BIC = NA_real_,
                     deviance = m$deviance,
                     df_residual = m$df)
    }
    default_glance <- function(m) {
      tibble::tibble(null_deviance = NA_real_,
                     df_null = NA_integer_,
                     logLik = NA_real_,
                     AIC = NA_real_,
                     BIC = NA_real_,
                     deviance = NA_real_,
                     df_residual = NA_integer_)
    }
    glmerMod_glance <- function(m) {
      tibble::tibble(null_deviance = NA_real_,
                     df_null = NA_integer_,
                     logLik = m@devcomp$cmp["logLik"],
                     AIC = m@devcomp$cmp["AIC"],
                     BIC = m@devcomp$cmp["BIC"],
                     deviance = m@devcomp$cmp["dev"],
                     df_residual = m@devcomp$cmp["df_residual"])
    }
    
    tryCatch({switch(class(m)[1],
                     "negbin" = mass_glance(m),
                     "zeroinfl" = zeroinfl_glance(m),
                     "speedglm" = speedglm_glance(m),
                     "glmerMod" = glmerMod_glance(m),
                     default_glance(m)
    )
    })
    
  }
  model_tbl %>% dplyr::mutate(glanced = purrr::map(model, private_glance)) %>%
    tidyr::unnest(glanced, .drop=TRUE)
}

likelihood_ratio_test_pval <- function(model_summary_x, model_summary_y) {
  dfs = round(abs(model_summary_x$df.residual  - model_summary_y$df.residual))
  LLR = 2 * abs(model_summary_x$logLik  - model_summary_y$logLik)
  p_val = stats::pchisq(LLR, dfs, lower.tail = FALSE)
}

#' Predict output of fitted models and return as a matrix
#'
#' @param model_tbl A tibble of model objects, generally output of
#'   \code{\link{fit_models}}.
#' @param new_data A data frame of new data to be passed to predict for
#'   prediction.
#' @param type String of type to pass to predict. Default is "response".
#'
#' @export
model_predictions <- function(model_tbl, new_data, type="response") {
  predict_helper <- function(model, new_data){
    tryCatch({
      stats::predict(model, newdata=new_data, type=type)
    }, error = function(e){
      retval = rep_len(NA, nrow(new_data))
      names(retval) = row.names(new_data)
      return(retval)
    })
  }
  model_tbl <- model_tbl %>%
    dplyr::mutate(predictions = purrr::map(model, predict_helper, new_data))
  pred_matrix = t(do.call(cbind, model_tbl$predictions))
  return(pred_matrix)
}

