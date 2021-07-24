#' Apply a preprocess transform model to a cell_data_set.
#'
#' Applies a previously calculated preprocess transform
#' model to a new count matrix. For more information read
#' the help information for save_transform_models.
#'
#' @param cds A cell_data_set to be transformed.
#' @param method A previously loaded transform model that
#'   is used to reduce the dimensions of the count matrix
#'   in cds.
#' @param block_size A numeric value for the DelayedArray
#'   block size used only in this function. Default is
#'   NULL, which does not affect the current block size.
#' @param cores The number of cores to use for the matrix multiplication.
#' @return A cell_data_set with a preprocess reduced count
#'   matrix.
#'
#' @export
#'
preprocess_transform <- function(cds, method=c('PCA'), block_size=NULL, cores=1) {
  #
  # Need to add processing for LSI. TF-IDF transform etc.
  #
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "method must be 'PCA'")

  method <- match.arg(method)

  assertthat::assert_that(!is.null(size_factors(cds)),
             msg = paste("You must call estimate_size_factors before calling",
                         "preprocess_cds."))
  assertthat::assert_that(sum(is.na(size_factors(cds))) == 0,
                          msg=paste("One or more cells has a size factor of",
                                    "NA."))
  assertthat::assert_that(!is.null(cds@reduce_dim_aux[[method]]),
                          msg=paste0("Method '", method, "' is not in the model",
                                    " object."))

  if(!is.null(block_size)) {
    block_size0 <- DelayedArray::getAutoBlockSize()
    DelayedArray::setAutoBlockSize(block_size)
  }

  set.seed(2016)
  norm_method <- cds@reduce_dim_aux[[method]][['model']][['norm_method']]
  pseudo_count <- cds@reduce_dim_aux[[method]][['model']][['pseudo_count']]
  rotation_matrix <- cds@reduce_dim_aux[[method]][['model']]$svd_v
  vcenter <- cds@reduce_dim_aux[[method]][['model']]$svd_center
  vscale <- cds@reduce_dim_aux[[method]][['model']]$svd_scale

  FM <- normalize_expr_data(cds, norm_method=norm_method, pseudo_count=pseudo_count)
  if (nrow(FM) == 0) {
    stop("Error: all rows have standard deviation zero")
  }

  # Don't select matrix rows by use_genes because intersect() does
  # it implicitly through the rotation matrix.

  fm_rowsums = Matrix::rowSums(FM)
  FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]
  
  # Thank you Maddy.
  intersect_genes <- intersect(rownames(rotation_matrix), rownames(FM))
  intersect_indices <- match(intersect_genes, rownames(rotation_matrix))

  # [intersect_genes,] orders FM rows by intersect_genes
  FM <- FM[intersect_genes,]
  
  xt <- Matrix::t(FM)
  xtda <- DelayedArray::DelayedArray(xt)

  vcenter <- vcenter[intersect_indices]
  vscale <- vscale[intersect_indices]

  xtdasc <- t(xtda) - vcenter
  xtdasc <- t(xtdasc / vscale)

  irlba_res <- list()

#  irlba_res$x <- xtdasc %*% rotation_matrix[intersect_indices,]
  irlba_res$x <- matrix_multiply_multicore(mat_a=xtdasc,
                                         mat_b=rotation_matrix[intersect_indices,],
                                         cores)

  irlba_res$x <- as.matrix(irlba_res$x)
  class(irlba_res) <- c('irlba_prcomp', 'prcomp')

  # 'reference' gene names are in the cds@preproc
  reducedDims(cds)[[method]] <- irlba_res$x

  if(!is.null(block_size)) {
    DelayedArray::setAutoBlockSize(block_size0)
  }

  matrix_id <- get_unique_id()
  counts_identity <- get_counts_identity(cds)
  reduce_dim_model_identity <- get_reduce_dim_model_identity(cds, method)
  cds <- initialize_reduce_dim_metadata(cds, method)
  cds <- set_reduce_dim_matrix_identity(cds, method,
                                        paste0('matrix:',method),
                                        matrix_id,
                                        counts_identity[['matrix_type']],
                                        counts_identity[['matrix_id']],
                                        reduce_dim_model_identity[['model_type']],
                                        reduce_dim_model_identity[['model_id']])
  # Keep the existing, loaded, reduce_dim_model_identity information.

  return(cds)
}


#' Apply an align_cds transform model to a cell_data_set.
#'
#' Applies a previously calculated align_cds transform model
#' to a new preprocess transformed matrix. For more
#' information read the help information for
#' save_transform_models.
#'
#' Note: this function is a place holder. It does not
#' map the transformed count matrix to aligned space
#' at this time because I don't know how to do so.
#'
#' @param cds A cell_data_set to be transformed.
#' @param method A previously loaded transform model that
#'   is used to reduce the dimensions of the preprocessed
#'   count matrix in cds.
#'
#' @return A cell_data_set with an align_cds transformed
#'   reduced count matrix.
#'
#' @export
#'
align_transform <- function(cds, method=c('Aligned')) {
  #
  # Need to add transformation code.
  #
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "method must be 'Aligned'")

  method <- match.arg(method)

  preprocess_method <- cds@reduce_dim_aux[['Aligned']][['model']][['preprocess_method']]
  preproc_res <- reducedDims(cds)[[preprocess_method]]
  assertthat::assert_that(!is.null(preproc_res),
                          msg=paste("Preprocessing for '",
                                    preprocess_method,
                                    "' does not exist.",
                                    " Please preprocess the matrix before",
                                    " calling align_transform using preprocess_transform."))

#  stop('This function is a place holder. It does not map the transformed count matrix to aligned space at this time because I don\'t know how to make it do so.')

  set.seed(2016)
  alignment_group <- cds@reduce_dim_aux[['Aligned']][['model']][['alignment_group']]
  alignment_k <- cds@reduce_dim_aux[['Aligned']][['model']][['alignment_k']]
  residual_model_formula_str <- cds@reduce_dim_aux[['Aligned']][['model']][['residual_model_formula_str']]
  nn_metric <- cds@reduce_dim_aux[['Aligned']][['nn_index']][['annoy_metric']]

  X.model_mat <- Matrix::sparse.model.matrix( stats::as.formula(residual_model_formula_str), data = colData(cds), drop.unused.levels = TRUE)
  fit <- limma::lmFit(Matrix::t(preproc_res), X.model_mat)
  beta <- fit$coefficients[, -1, drop = FALSE]
  cds@reduce_dim_aux[['Aligned']][['model']][['beta']] <- beta
  preproc_res <- Matrix::t(as.matrix(Matrix::t(preproc_res)) - beta %*% Matrix::t(X.model_mat[, -1]))
  corrected_PCA = batchelor::reducedMNN(as.matrix(preproc_res), batch=colData(cds)[,alignment_group], k=alignment_k)
  preproc_res = corrected_PCA$corrected
  reducedDims(cds)[['Aligned']] <- as.matrix(preproc_res)

  matrix_id <- get_unique_id()
  reduce_dim_matrix_identity <- get_reduce_dim_matrix_identity(cds, preprocess_method)
  reduce_dim_model_identity <- get_reduce_dim_model_identity(cds, method)
  cds <- initialize_reduce_dim_metadata(cds, method)
  cds <- set_reduce_dim_matrix_identity(cds, method,
                                        paste0('matrix:', method),
                                        matrix_id,
                                        reduce_dim_matrix_identity[['matrix_type']],
                                        reduce_dim_matrix_identity[['matrix_id']],
                                        reduce_dim_model_identity[['model_type']],
                                        reduce_dim_model_identity[['model_id']])
  # Keep the existing, loaded, reduce_dim_model_identity information.

  return(cds)
}


#' Apply a reduce_transform transform model to a cell_data_set.
#'
#' Applies a previously calculated reduce_dimension transform
#' model to a new preprocess transformed matrix. For more
#' information read the help information for
#' save_transform_models.
#'
#' @param cds A cell_data_set to be transformed.
#' @param preprocess_method A previously loaded preprocess method.
#'   The default is NULL, which uses the preprocess_method that
#'   was used when the reduce_dimension model was built.
#' @param method A previously loaded reduce_dimension transform
#'   model that is used to reduce the dimensions of the
#'   preprocessed count matrix in cds.
#'
#' @return A cell_data_set with a reduce_dimension transformed
#'   reduced count matrix.
#'
#' @export
#'
reduce_dimension_transform <- function(cds, preprocess_method=NULL, method=c('UMAP')) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "method must be 'UMAP'")

  if(is.null(preprocess_method)) {
    preprocess_method <- cds@reduce_dim_aux[[method]][['model']][['umap_preprocess_method']]
  } else
  if(!is.null(preprocess_method) && !(preprocess_method %in% c('PCA', 'LSI', 'Aligned'))) {
    stop('Preprocess_method must be one of \'PCA\', \'LSI\', or \'Aligned\'.')
  }

  if(is.null(cds@reduce_dim_aux[[preprocess_method]])) {
    stop('There is no transform model for preprocess_method \'', preprocess_method, '\'.')
  }

  method <- match.arg(method)

  preproc_res <- reducedDims(cds)[[preprocess_method]]
  assertthat::assert_that(!is.null(preproc_res),
                          msg=paste("Preprocessing for '",
                                    preprocess_method,
                                    "' does not exist.",
                                    " Please preprocess the matrix before",
                                    " calling reduce_dimension_transform using preprocess_transform",
                                    " or align_transform."))

  # Notes:
  #   o  uwot::umap_transform() depends on the RNG state and we want
  #      consistent results when called here and in reduce_dimensions()
  #      function(s).
  #
  set.seed(2016)
  umap_model <- cds@reduce_dim_aux[[method]][['model']][['umap_model']]
  reducedDims(cds)[[method]] <- uwot:::umap_transform(X=preproc_res, model=umap_model, init='weighted', n_sgd_threads=1)

  matrix_id <- get_unique_id()
  reduce_dim_matrix_identity <- get_reduce_dim_matrix_identity(cds, preprocess_method)
  reduce_dim_model_identity <- get_reduce_dim_model_identity(cds, method)
  cds <- initialize_reduce_dim_metadata(cds, method)
  cds <- set_reduce_dim_matrix_identity(cds, method,
                                        paste0('matrix:', method),
                                        matrix_id,
                                        reduce_dim_matrix_identity[['matrix_type']],
                                        reduce_dim_matrix_identity[['matrix_id']],
                                        reduce_dim_model_identity[['model_type']],
                                        reduce_dim_model_identity[['model_id']])
  # Keep the existing, loaded, reduce_dim_model_identity information.

  return(cds)
}


#
# Note: this transformation is not generally effective.
#
align_beta_transform <- function(cds, preprocess_method = 'PCA') {
  method <- 'Aligned'
  preproc_res <- reducedDims(cds)[[preprocess_method]]
  beta <- cds@reduce_dim_aux[['Aligned']][['model']][['beta']]
  residual_model_formula_str <- cds@reduce_dim_aux[['Aligned']][['model']][['residual_model_formula_str']]
  X.model_mat <- Matrix::sparse.model.matrix(
    stats::as.formula(residual_model_formula_str),
    data = colData(cds),
    drop.unused.levels = TRUE)
  reducedDims(cds)[['Aligned']] <- Matrix::t(as.matrix(Matrix::t(preproc_res)) -
                                   beta %*% Matrix::t(X.model_mat[, -1]))

  matrix_id <- get_unique_id()
  cds <- initialize_reduce_dim_metadata(cds, method)
  reduce_dim_matrix_identity <- get_reduce_dim_matrix_identity(cds, preprocess_method)
  reduce_dim_model_identity <- get_reduce_dim_model_identity(cds, method)
  cds <- set_reduce_dim_matrix_identity(cds, method,
                                        paste0('matrix:', method),
                                        matrix_id,
                                        reduce_dim_matrix_identity[['matrix_type']],
                                        reduce_dim_matrix_identity[['matrix_id']],
                                        reduce_dim_model_identity[['model_type']],
                                        reduce_dim_model_identity[['model_id']])
  # Keep the existing, loaded, reduce_dim_model_identity information.

  return(cds)
}

