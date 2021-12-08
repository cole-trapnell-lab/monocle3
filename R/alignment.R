#' Align cells from different groups within a cds
#'
#' @description Data sets that contain cells from different groups often
#' benefit from alignment to subtract differences between them. Alignment
#' can be used to remove batch effects, subtract the effects of treatments,
#' or even potentially compare across species.
#' \code{align_cds} executes alignment and stores these adjusted coordinates.
#'
#' This function can be used to subtract both continuous and discrete batch
#' effects. For continuous effects, \code{align_cds} fits a linear model to the
#' cells' PCA or LSI coordinates and subtracts them using Limma. For discrete
#' effects, you must provide a grouping of the cells, and then these groups are
#' aligned using Batchelor, a "mutual nearest neighbor" algorithm described in:
#'
#' Haghverdi L, Lun ATL, Morgan MD, Marioni JC (2018). "Batch effects in
#' single-cell RNA-sequencing data are corrected by matching mutual nearest
#' neighbors." Nat. Biotechnol., 36(5), 421-427. doi: 10.1038/nbt.4091
#'
#' @param cds the cell_data_set upon which to perform this operation
#' @param preprocess_method a string specifying the low-dimensional space
#'   in which to perform alignment, currently either PCA or LSI. Default is
#'   "PCA".
#' @param residual_model_formula_str NULL or a string model formula specifying
#'   any effects to subtract from the data before dimensionality reduction.
#'   Uses a linear model to subtract effects. For non-linear effects, use
#'   alignment_group. Default is NULL.
#' @param alignment_group String specifying a column of colData to use for
#'  aligning groups of cells. The column specified must be a factor.
#'  Alignment can be used to subtract batch effects in a non-linear way.
#'  For correcting continuous effects, use residual_model_formula_str.
#'  Default is NULL.
#' @param alignment_k The value of k used in mutual nearest neighbor alignment
#' @param verbose Whether to emit verbose output during dimensionality
#'   reduction
#' @param build_nn_index logical When this argument is set to TRUE,
#'   align_cds builds the nearest neighbor index from the
#'   aligned reduced matrix for later use. Default is FALSE.
#' @param nn_control A list of parameters used to make the nearest
#'  neighbor index. See the set_nn_control help for detailed information.
#' @param ... additional arguments to pass to limma::lmFit if
#'   residual_model_formula is not NULL
#' @return an updated cell_data_set object
#' @export
align_cds <- function(cds,
                      preprocess_method = c("PCA", "LSI"),
                      alignment_group=NULL,
                      alignment_k=20,
                      residual_model_formula_str=NULL,
                      verbose=FALSE,
                      build_nn_index=FALSE,
                      nn_control=list(),
                      ...){
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(preprocess_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "preprocess_method must be one of 'PCA' or 'LSI'")
  preprocess_method <- match.arg(preprocess_method)

  preproc_res <- reducedDims(cds)[[preprocess_method]]
  assertthat::assert_that(!is.null(preproc_res),
                          msg = paste0("Preprocessing for '",
                                      preprocess_method, "' does not exist. ",
                                      "Please make sure you have run ",
                                      "preprocess_cds with",
                                      "preprocess_method = '",
                                      preprocess_method,
                                      "' before calling align_cds."))

  if(build_nn_index) {
    nn_control <- set_nn_control(mode=1,
                                 nn_control=nn_control,
                                 nn_control_default=get_global_variable('nn_control_annoy_cosine'),
                                 cds=NULL,
                                 reduction_method=NULL,
                                 k=NULL,
                                 verbose=verbose)
  }

  set.seed(2016)

  cds <- initialize_reduce_dim_metadata(cds, 'Aligned')
  cds <- initialize_reduce_dim_model_identity(cds, 'Aligned')

  if (!is.null(residual_model_formula_str)) {
    if (verbose) message("Removing residual effects")
    X.model_mat <- Matrix::sparse.model.matrix(
      stats::as.formula(residual_model_formula_str),
      data = colData(cds),
      drop.unused.levels = TRUE)

    fit <- limma::lmFit(Matrix::t(preproc_res), X.model_mat, ...)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    preproc_res <- Matrix::t(as.matrix(Matrix::t(preproc_res)) -
                               beta %*% Matrix::t(X.model_mat[, -1]))
    cds@reduce_dim_aux[['Aligned']][['model']][['beta']] <- beta
  }

  if(!is.null(alignment_group)) {
    message(paste("Aligning cells from different batches using Batchelor.",
                  "\nPlease remember to cite:\n\t Haghverdi L, Lun ATL,",
                  "Morgan MD, Marioni JC (2018). 'Batch effects in",
                  "single-cell RNA-sequencing data are corrected by matching",
                  "mutual nearest neighbors.' Nat. Biotechnol., 36(5),",
                  "421-427. doi: 10.1038/nbt.4091"))
    corrected_PCA = batchelor::reducedMNN(as.matrix(preproc_res),
                                          batch=colData(cds)[,alignment_group],
                                          k=alignment_k)
    preproc_res = corrected_PCA$corrected
    cds <- add_citation(cds, "MNN_correct")
  }
  reducedDims(cds)[["Aligned"]] <- as.matrix(preproc_res)

  #
  # Notes:
  #   o  the functions save_transform_models/load_transform_models
  #      expect that the reduce_dim_aux slot consists of a SimpleList
  #      that stores information about methods with the elements
  #        reduce_dim_aux[[method]][['model']] for the transform elements
  #        reduce_dim_aux[[method]][[nn_method]] for the annoy index
  #      and depends on the elements within model and nn_method.
  #
  cds@reduce_dim_aux[['Aligned']][['model']][['preprocess_method']] <- preprocess_method
  cds@reduce_dim_aux[['Aligned']][['model']][['alignment_group']] <- alignment_group
  cds@reduce_dim_aux[['Aligned']][['model']][['alignment_k']] <- alignment_k
  cds@reduce_dim_aux[['Aligned']][['model']][['residual_model_formula_str']] <- residual_model_formula_str

  matrix_id <- get_unique_id(reducedDims(cds)[["Aligned"]])
  reduce_dim_matrix_identity <- get_reduce_dim_matrix_identity(cds, preprocess_method) 

  cds <- set_reduce_dim_matrix_identity(cds, 'Aligned',
                                        'matrix:Aligned',
                                        matrix_id,
                                        reduce_dim_matrix_identity[['matrix_type']],
                                        reduce_dim_matrix_identity[['matrix_id']],
                                        'matrix:Aligned',
                                        matrix_id)
  reduce_dim_model_identity <- get_reduce_dim_model_identity(cds, preprocess_method)
  cds <- set_reduce_dim_model_identity(cds, 'Aligned',
                                       'matrix:Aligned',
                                       matrix_id,
                                       reduce_dim_model_identity[['model_type']],
                                       reduce_dim_model_identity[['model_id']])

  if( build_nn_index ) {
    nn_index <- make_nn_index(subject_matrix=reducedDims(cds)[['Aligned']],
                              nn_control=nn_control,
                              verbose=verbose)
    cds <- set_cds_nn_index(cds=cds,
                            reduction_method='Aligned',
                            nn_index,
                            nn_control=nn_control,
                            verbose=verbose)
  }
  else
    cds <- clear_cds_nn_index(cds=cds, reduction_method='Aligned', 'all')

  cds
}

