

#' @param modelFormulaStr A model formula, passed as a string, specifying how to group the cells prior to estimated dispersion.
#' The default groups all cells together.
#' @param relative_expr Whether to transform expression into relative values
#' @param min_cells_detected Only include genes detected above lowerDetectionLimit in at least this many cells in the dispersion calculation
#' @param remove_outliers Whether to remove outliers (using Cook's distance) when estimating dispersions
#' @param cores The number of cores to use for computing dispersions
#' @export
estimate_dispersions <- function(cds, modelFormulaStr="~ 1",
                                 relative_expr=TRUE, min_cells_detected=1,
                                 remove_outliers=TRUE, cores=1,...) {
  dispModelName="blind"
  stopifnot( is( cds, "cell_data_set" ) )

  if(!(identical("negbinomial.size", cds@expression_family) || identical("negbinomial", cds@expression_family))){
    stop("Error: estimate_dispersions only works, and is only needed, when you're using a cell_data_set with a negbinomial or negbinomial.size expression family")
  }

  if( any( is.na( size_factors(cds) ) ) )
    stop( "NAs found in size factors. Have you called 'estimate_size_factors'?" )

  if( length(list(...)) != 0 )
    warning( "in estimate_dispersions: Ignoring extra argument(s)." )

  # Remove results from previous fits
  cds@disp_fit_info = new.env( hash=TRUE )
  if(!(('negbinomial' == cds@expression_family) ||
       ('negbinomial.size' == cds@expression_family))) {
    stop("Error: estimate_dispersions only works, and is only needed, when you're using a cell_data_set with a negbinomial or negbinomial.size expression family")
  }

  mu <- NA
  model_terms <- unlist(lapply(stringr::str_split(modelFormulaStr, "~|\\+|\\*"),
                               stringr::str_trim))
  model_terms <- model_terms[model_terms != ""]
  progress_opts <- options()$dplyr.show_progress
  options(dplyr.show_progress = T)

  # FIXME: this needs refactoring, badly.
  if (cds@expression_family %in% c("negbinomial", "negbinomial.size")){
    if (length(model_terms) > 1 || (length(model_terms) == 1 && model_terms[1] != "1")){
      cds_pdata <- dplyr::group_by_(dplyr::select_(tibble::rownames_to_column(pData(cds)), "rowname", .dots=model_terms), .dots=model_terms)
      disp_table <- as.data.frame(cds_pdata %>% do(disp_calc_helper_NB(cds[,.$rowname], min_cells_detected)))
    }else{
      cds_pdata <- dplyr::group_by_(dplyr::select_(tibble::rownames_to_column(pData(cds)), "rowname"))
      disp_table <- as.data.frame(disp_calc_helper_NB(cds, min_cells_detected))
      #disp_table <- data.frame(rowname = row.names(type_res), CellType = type_res)
    }

    #message("fitting disersion curves")
    #print (disp_table)
    if(!is.list(disp_table))
      stop("Parametric dispersion fitting failed, please set a different lowerDetectionLimit")
    #disp_table <- do.call(rbind.data.frame, disp_table)
    disp_table <- subset(disp_table, is.na(mu) == FALSE)
    res <- parametricDispersionFit(disp_table, verbose = FALSE)
    fit <- res[[1]]
    coefs <- res[[2]]
    if (remove_outliers){
      CD <- cooks.distance(fit)
      cooksCutoff <- 4/nrow(disp_table)
      message (paste("Removing", length(CD[CD > cooksCutoff]), "outliers"))
      outliers <- union(names(CD[CD > cooksCutoff]), setdiff(row.names(disp_table), names(CD)))
      res <- parametricDispersionFit(disp_table[row.names(disp_table) %in% outliers == FALSE,], verbose = FALSE)
      fit <- res[[1]]
      coefs <- res[[2]]
    }
    names( coefs ) <- c( "asymptDisp", "extraPois" )
    ans <- function( q )
      coefs[1] + coefs[2] / q
    attr( ans, "coefficients" ) <- coefs

  }

  cds@disp_fit_info[[dispModelName]] <- list(disp_table = disp_table, disp_func = ans)

  validObject( cds )
  cds
}


disp_calc_helper_NB <- function(cds, min_cells_detected) {
  rounded <- round(exprs(cds))
  nzGenes <- Matrix::rowSums(rounded > cds@lower_detection_limit)
  nzGenes <- names(nzGenes[nzGenes > min_cells_detected])

  # Note: we do these operations as DelayedArray ops because standard operations will trigger a conversion
  # to an in-memory dense matrix. DelayedArray uses block processing. Control block size with:
  # options(DelayedArray.block.size=100e6)
  # We should make this clear in the documentation, and possibly
  # emit a message to users on calling this (and possibly other) functions.

  #Progress Bar here
  x <- DelayedArray::DelayedArray(Matrix::t(Matrix::t(rounded[nzGenes,]) / pData(cds[nzGenes,])$Size_Factor))

  xim <- mean(1/ pData(cds[nzGenes,])$Size_Factor)

  if (is_sparse_matrix(exprs(cds))){
    f_expression_mean <- methods::as(DelayedMatrixStats::rowMeans2(x), "sparseVector")
  }else{
    f_expression_mean <- DelayedMatrixStats::rowMeans2(x)
  }

  # For NB: Var(Y)=mu*(1+mu/k)
  f_expression_var <- DelayedMatrixStats::rowVars(x)

  disp_guess_meth_moments <- f_expression_var - xim * f_expression_mean

  disp_guess_meth_moments <- disp_guess_meth_moments / (f_expression_mean^2) #fix the calculation of k


  res <- data.frame(mu=as.vector(f_expression_mean), disp=as.vector(disp_guess_meth_moments))
  res[res$mu == 0]$mu = NA
  res[res$mu == 0]$disp = NA
  res$disp[res$disp < 0] <- 0

  res <- cbind(gene_id=row.names(fData(cds[nzGenes,])), res)
  res
}

## This function was swiped from DESeq (Anders and Huber) and modified for our purposes
parametricDispersionFit <- function( disp_table, verbose = FALSE, initial_coefs=c(1e-6, 1) ) {
  coefs <- initial_coefs
  iter <- 0
  while(TRUE) {
    residuals <- disp_table$disp / ( coefs[1] + coefs[2] / disp_table$mu )
    good <- disp_table[which(disp_table$disp > 0 & (residuals < 10000) ),]

    if(verbose)
      fit <- stats::glm( disp ~ I(1/mu), data=good,
                  family=stats::Gamma(link="identity"), start=coefs )
    else
      suppressWarnings(fit <- stats::glm( disp ~ I(1/mu), data=good,
                                   family=stats::Gamma(link="identity"), start=coefs ))

    oldcoefs <- coefs
    coefs <- stats::coefficients(fit)
    if (coefs[1] < initial_coefs[1]){
      coefs[1] <- initial_coefs[1]
    }
    if (coefs[2] < 0){
      stop( "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimate_dispersions')" )
    }
    if( sum( log( coefs / oldcoefs )^2 ) < coefs[1] )
      break
    iter <- iter + 1

    if( iter > 10 ) {
      warning( "Dispersion fit did not converge." )
      break
    }
  }

  if( !all( coefs > 0 ) ){
    stop( "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimate_dispersions')" )
  }

  rm(disp_table)

  list(fit, coefs)
}


