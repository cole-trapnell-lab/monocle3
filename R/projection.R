#' @title Apply a preprocess transform model to a cell_data_set.
#'
#' @description Applies a previously calculated preprocess
#' transform model to a new count matrix. For more
#' information read the help information for save_transform_models.
#'
#' @param cds a cell_data_set to be transformed.
#' @param reduction_method a previously loaded transform
#'   model that is used to reduce the dimensions of the
#'   count matrix in the cell_data_set. The "PCA" and "LSI"
#'   transforms are supported. The default is "PCA".
#' @param block_size a numeric value for the DelayedArray
#'   block size used only in this function. Default is
#'   NULL, which does not affect the current block size.
#' @param cores the number of cores to use for the matrix
#'   multiplication. The default is 1.
#'
#' @return a cell_data_set with a preprocess reduced count
#'   matrix.
#'
#' @export
preprocess_transform <- function(cds, reduction_method=c('PCA', 'LSI'), block_size=NULL, cores=1) {
  #
  # Need to add processing for LSI. TF-IDF transform etc.
  #
  assertthat::assert_that(is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be 'PCA' or 'LSI'")

  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(!is.null(size_factors(cds)),
             msg = paste("You must call estimate_size_factors before calling",
                         "preprocess_cds."))
  assertthat::assert_that(sum(is.na(size_factors(cds))) == 0,
                          msg=paste("One or more cells has a size factor of",
                                    "NA."))
  assertthat::assert_that(!is.null(cds@reduce_dim_aux[[reduction_method]]),
                          msg=paste0("Reduction method '", reduction_method, "' is not in the model",
                                    " object."))

  if(!is.null(block_size)) {
    block_size0 <- DelayedArray::getAutoBlockSize()
    DelayedArray::setAutoBlockSize(block_size)
  }

  set.seed(2016)

  if(reduction_method == 'PCA') {
    norm_method <- cds@reduce_dim_aux[[reduction_method]][['model']][['norm_method']]
    pseudo_count <- cds@reduce_dim_aux[[reduction_method]][['model']][['pseudo_count']]
    rotation_matrix <- cds@reduce_dim_aux[[reduction_method]][['model']]$svd_v
    vcenter <- cds@reduce_dim_aux[[reduction_method]][['model']]$svd_center
    vscale <- cds@reduce_dim_aux[[reduction_method]][['model']]$svd_scale
  
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

    if(any(is.na(intersect_indices))) {
      stop('gene sets differ: genes in the subject matrix are not in the reference matrix')
    }

    if((length(intersect_indices)/length(rownames(FM))) < 0.5) {
      warning('fewer than half the genes in the subject matrix are also in the reference matrix: are the matrices prepared using the same gene set?')
    } 

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
    reducedDims(cds)[[reduction_method]] <- irlba_res$x

  }
  else if(reduction_method == 'LSI') {
    norm_method <- cds@reduce_dim_aux[[reduction_method]][['model']][['norm_method']]
    pseudo_count <- cds@reduce_dim_aux[[reduction_method]][['model']][['pseudo_count']]
    rotation_matrix <- cds@reduce_dim_aux[[reduction_method]][['model']][['svd_v']]
    log_scale_tf <- cds@reduce_dim_aux[[reduction_method]][['model']][['log_scale_tf']]
    frequencies <- cds@reduce_dim_aux[[reduction_method]][['model']][['frequencies']]
    scale_factor <- cds@reduce_dim_aux[[reduction_method]][['model']][['scale_factor']]
    col_sums <- cds@reduce_dim_aux[[reduction_method]][['model']][['col_sums']]
    row_sums <- cds@reduce_dim_aux[[reduction_method]][['model']][['row_sums']]
    num_cols <- cds@reduce_dim_aux[[reduction_method]][['model']][['num_cols']]

    FM <- normalize_expr_data(cds, norm_method=norm_method, pseudo_count=pseudo_count)
    if (nrow(FM) == 0) {
      stop("Error: all rows have standard deviation zero")
    }
 
    fm_rowsums = Matrix::rowSums(FM)
    FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]
   
    intersect_features <- intersect(rownames(rotation_matrix), rownames(FM))
    intersect_indices <- match(intersect_features, rownames(rotation_matrix))

    if(any(is.na(intersect_indices))) {
      stop('feature sets differ: features in the subject matrix are not in the reference matrix')
    }

    if((length(intersect_indices)/length(rownames(FM))) < 0.5) {
      warning('fewer than half the features in the subject matrix are also in the reference matrix: are the matrices prepared using the same feature set?')
    }

    # [intersect_features,] orders FM rows by intersect_features
    FM <- FM[intersect_features,]

    # Andrew's tfidf with modifications
    if (frequencies) {
      # "term frequency" method
      tf <- Matrix::t(Matrix::t(FM) / Matrix::colSums(FM))
    } else {
      # "raw count" method
      tf <- FM
    }

    # Either TF method can optionally be log scaled
    if (log_scale_tf) {
      if (frequencies) {
        tf@x = log1p(tf@x * scale_factor)
      } else {
        tf@x = log1p(tf@x * 1)
      }
    }

    idf <- log(1 + num_cols / row_sums)
    idf <- idf[intersect_indices]

    tf_idf_counts = tryCatch({
      tf_idf_counts = tf * idf
    }, error = function(e) {
      print(paste("TF*IDF multiplication too large for in-memory, falling back",
                  "on DelayedArray."))
      options(DelayedArray.block.size=block_size)
      DelayedArray:::set_verbose_block_processing(TRUE)
  
      tf = DelayedArray::DelayedArray(tf)
      idf = as.matrix(idf)
  
      tf_idf_counts = tf * idf
    })

    rownames(tf_idf_counts) = rownames(FM)
    colnames(tf_idf_counts) = colnames(FM)
    tf_idf_counts = methods::as(tf_idf_counts, "sparseMatrix")

    xt <- Matrix::t(tf_idf_counts)
    xtda <- DelayedArray::DelayedArray(xt)

    irlba_res <- list()
    irlba_res$x <- matrix_multiply_multicore(mat_a=xtda,
                                           mat_b=rotation_matrix[intersect_indices,],
                                           cores)
    irlba_res$x <- as.matrix(irlba_res$x)
    class(irlba_res) <- c('irlba_prcomp', 'prcomp')

    # 'reference' gene names are in the cds@preproc
    reducedDims(cds)[[reduction_method]] <- irlba_res$x

  }

  if(!is.null(block_size)) {
    DelayedArray::setAutoBlockSize(block_size0)
  }

  matrix_id <- get_unique_id(reducedDims(cds)[[reduction_method]])
  counts_identity <- get_counts_identity(cds)
  reduce_dim_model_identity <- get_reduce_dim_model_identity(cds, reduction_method)
  cds <- initialize_reduce_dim_metadata(cds, reduction_method)
  cds <- set_reduce_dim_matrix_identity(cds, reduction_method,
                                        paste0('matrix:',reduction_method),
                                        matrix_id,
                                        counts_identity[['matrix_type']],
                                        counts_identity[['matrix_id']],
                                        reduce_dim_model_identity[['model_type']],
                                        reduce_dim_model_identity[['model_id']])
  # Keep the existing, loaded, reduce_dim_model_identity information.

  return(cds)
}


#' @title Apply an alignment transform model to a cell_data_set.
#'
#' @description align_transform is not supported. Co-embed your
#'  data sets if you need batch correction.
#'
#' @param cds a cell_data_set to be transformed.
#' @param reduction_method a previously loaded transform
#'   model that is used to reduce the dimensions of the
#'   count matrix in the cell_data_set. The "Aligned"
#'   transform is not supported.
#'
#' @return The cds is returned without processing.
#'
#' @export
align_transform <- function(cds, reduction_method=c('Aligned')) {
  stop(paste('align_transform is not supported. If you need batch correction',
              'you will need to co-embed the data sets.'))
  return(cds)
}


#  #' Apply an align_cds transform model to a cell_data_set.
#  #'
#  #' Applies a previously calculated align_cds transform model
#  #' to a new preprocess transformed matrix. For more
#  #' information read the help information for
#  #' save_transform_models.
#  #'
#  #' Note: this function is a place holder. It does not
#  #' map the transformed count matrix to aligned space
#  #' at this time because I don't know how to do so.
#  #'
#  #' @param cds A cell_data_set to be transformed.
#  #' @param reduction_method A previously loaded transform model that
#  #'   is used to reduce the dimensions of the preprocessed
#  #'   count matrix in cds.
#  #'
#  #' @return A cell_data_set with an align_cds transformed
#  #'   reduced count matrix.
#  #'
#  #' @export
# align_transform <- function(cds, reduction_method=c('Aligned')) {
#   #
#   # Need to add transformation code.
#   #
#   assertthat::assert_that(is(cds, 'cell_data_set'),
#                           msg=paste('cds parameter is not a cell_data_set'))
#   assertthat::assert_that(
#     tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
#              error = function(e) FALSE),
#     msg = "reduction_method must be 'Aligned'")
# 
#   reduction_method <- match.arg(reduction_method)
# 
#   preprocess_method <- cds@reduce_dim_aux[['Aligned']][['model']][['preprocess_method']]
#   preproc_res <- reducedDims(cds)[[preprocess_method]]
#   assertthat::assert_that(!is.null(preproc_res),
#                           msg=paste("Preprocessing for '",
#                                     preprocess_method,
#                                     "' does not exist.",
#                                     " Please preprocess the matrix before",
#                                     " calling align_transform using preprocess_transform."))
# 
#   stop('This function is a place holder. It does not map the transformed count matrix to aligned space at this time because we don\'t know how to it yet.')
# 
#   set.seed(2016)
#   alignment_group <- cds@reduce_dim_aux[['Aligned']][['model']][['alignment_group']]
#   alignment_k <- cds@reduce_dim_aux[['Aligned']][['model']][['alignment_k']]
#   residual_model_formula_str <- cds@reduce_dim_aux[['Aligned']][['model']][['residual_model_formula_str']]
#   X.model_mat <- Matrix::sparse.model.matrix( stats::as.formula(residual_model_formula_str), data = colData(cds), drop.unused.levels = TRUE)
#   fit <- limma::lmFit(Matrix::t(preproc_res), X.model_mat)
#   beta <- fit$coefficients[, -1, drop = FALSE]
#   cds@reduce_dim_aux[['Aligned']][['model']][['beta']] <- beta
#   preproc_res <- Matrix::t(as.matrix(Matrix::t(preproc_res)) - beta %*% Matrix::t(X.model_mat[, -1]))
#   corrected_PCA = batchelor::reducedMNN(as.matrix(preproc_res), batch=colData(cds)[,alignment_group], k=alignment_k)
#   preproc_res = corrected_PCA$corrected
#   reducedDims(cds)[['Aligned']] <- as.matrix(preproc_res)
# 
#   matrix_id <- get_unique_id()
#   reduce_dim_matrix_identity <- get_reduce_dim_matrix_identity(cds, preprocess_method)
#   reduce_dim_model_identity <- get_reduce_dim_model_identity(cds, reduction_method)
#   cds <- initialize_reduce_dim_metadata(cds, reduction_method)
#   cds <- set_reduce_dim_matrix_identity(cds, reduction_method,
#                                         paste0('matrix:', reduction_method),
#                                         matrix_id,
#                                         reduce_dim_matrix_identity[['matrix_type']],
#                                         reduce_dim_matrix_identity[['matrix_id']],
#                                         reduce_dim_model_identity[['model_type']],
#                                         reduce_dim_model_identity[['model_id']])
#   # Keep the existing, loaded, reduce_dim_model_identity information.
# 
#   return(cds)
# }


#' @title Apply a reduce_dimension transform model to a cell_data_set.
#'
#' @description Applies a previously calculated reduce_dimension
#'   transform model to a new preprocess transformed matrix. For
#'   more information read the help information for
#'   save_transform_models.
#'
#' @param cds a cell_data_set to be transformed.
#' @param preprocess_method the reduced dimension matrix to be
#'   transformed using the reduction_method transform model.
#'   The default is NULL, which uses the preprocess_method that
#'   was used when the reduce_dimension model was built.
#' @param reduction_method a previously loaded reduce_dimension transform
#'   model that is used to reduce the dimensions of the preprocessed
#'   matrix in the cell_data_set. Only "UMAP" is supported.
#'
#' @return a cell_data_set with a transformed
#'   reduced count matrix.
#'
#' @importFrom methods is
#' @export
#'
reduce_dimension_transform <- function(cds, preprocess_method=NULL, reduction_method=c('UMAP')) {
  assertthat::assert_that(is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be 'UMAP'")

  if(is.null(preprocess_method)) {
    preprocess_method <- cds@reduce_dim_aux[[reduction_method]][['model']][['umap_preprocess_method']]
  } else
  if(!is.null(preprocess_method) && !(preprocess_method %in% c('PCA', 'LSI', 'Aligned'))) {
    stop('Preprocess_method must be one of \'PCA\', \'LSI\', or \'Aligned\'.')
  }

  if(is.null(cds@reduce_dim_aux[[preprocess_method]])) {
    stop('There is no transform model for preprocess_method \'', preprocess_method, '\'.')
  }

  reduction_method <- match.arg(reduction_method)

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
  umap_model <- cds@reduce_dim_aux[[reduction_method]][['model']][['umap_model']]
  reducedDims(cds)[[reduction_method]] <- uwot::umap_transform(X=preproc_res, model=umap_model, init='weighted', n_sgd_threads=1)

  matrix_id <- get_unique_id(reducedDims(cds)[[reduction_method]])
  reduce_dim_matrix_identity <- get_reduce_dim_matrix_identity(cds, preprocess_method)
  reduce_dim_model_identity <- get_reduce_dim_model_identity(cds, reduction_method)
  cds <- initialize_reduce_dim_metadata(cds, reduction_method)
  cds <- set_reduce_dim_matrix_identity(cds, reduction_method,
                                        paste0('matrix:', reduction_method),
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
  reduction_method <- 'Aligned'
  preproc_res <- reducedDims(cds)[[preprocess_method]]
  beta <- cds@reduce_dim_aux[['Aligned']][['model']][['beta']]
  residual_model_formula_str <- cds@reduce_dim_aux[['Aligned']][['model']][['residual_model_formula_str']]
  X.model_mat <- Matrix::sparse.model.matrix(
    stats::as.formula(residual_model_formula_str),
    data = colData(cds),
    drop.unused.levels = TRUE)
  reducedDims(cds)[['Aligned']] <- Matrix::t(as.matrix(Matrix::t(preproc_res)) -
                                   beta %*% Matrix::t(X.model_mat[, -1]))

  matrix_id <- get_unique_id(reducedDims(cds)[['Aligned']])
  cds <- initialize_reduce_dim_metadata(cds, reduction_method)
  reduce_dim_matrix_identity <- get_reduce_dim_matrix_identity(cds, preprocess_method)
  reduce_dim_model_identity <- get_reduce_dim_model_identity(cds, reduction_method)
  cds <- set_reduce_dim_matrix_identity(cds, reduction_method,
                                        paste0('matrix:', reduction_method),
                                        matrix_id,
                                        reduce_dim_matrix_identity[['matrix_type']],
                                        reduce_dim_matrix_identity[['matrix_id']],
                                        reduce_dim_model_identity[['model_type']],
                                        reduce_dim_model_identity[['model_id']])
  # Keep the existing, loaded, reduce_dim_model_identity information.

  return(cds)
}

