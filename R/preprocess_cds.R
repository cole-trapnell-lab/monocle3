#' Preprocess a cds to prepare for trajectory inference
#'
#' @description Most analyses (including trajectory inference, and clustering)
#' in Monocle3, require various normalization and preprocessing steps.
#' \code{preprocess_cds} executes and stores these preprocessing steps.
#'
#' Specifically, depending on the options selected, \code{preprocess_cds} first
#' normalizes the data by log and size factor to address depth differences, or
#' by size factor only. Next, \code{preprocess_cds} calculates a lower
#' dimensional space that will be used as the input for further dimensionality
#' reduction like tSNE and UMAP.
#'
#' @param cds the cell_data_set upon which to perform this operation
#' @param method a string specifying the initial dimension method to use,
#'   currently either "PCA" or "LSI". For "LSI" (latent semantic indexing), it
#'   converts the (sparse) expression matrix into a tf-idf matrix and then
#'   performs SVD to decompose the gene expression / cells into certain
#'   modules / topics. Default is "PCA".
#' @param num_dim the dimensionality of the reduced space.
#' @param norm_method Determines how to transform expression values prior to
#'   reducing dimensionality. Options are "log", "size_only", and "none".
#'   Default is "log". Users should only use "none" if they are confident that
#'   their data is already normalized.
#' @param use_genes NULL or a list of gene IDs. If a list of gene IDs, only
#'   this subset of genes is used for dimensionality reduction. Default is
#'   NULL.
#' @param pseudo_count NULL or the amount to increase expression values before
#'   normalization and dimensionality reduction. If NULL (default), a
#'   pseudo_count of 1 is added for log normalization and 0 is added for size
#'   factor only normalization.
#' @param scaling When this argument is set to TRUE (default), it will scale
#'   each gene before running trajectory reconstruction. Relevant for
#'   method = PCA only.
#' @param verbose Whether to emit verbose output during dimensionality
#'   reduction
#' @param build_nn_index logical When this argument is set to TRUE,
#'   preprocess_cds builds and stores the nearest neighbor index from the
#'   reduced dimension matrix for later use. Default is FALSE.
#' @param nn_control An optional list of parameters used to make the nearest
#'  neighbor index. See the set_nn_control help for detailed information.
#' @param matrix_control An optional list of parameters that control the
#'  storage of intermediate matrices that are required for cds
#'  preprocessing. By default, matrices are stored in memory as dgCMatrix
#'  class (compressed sparse matrix) objects, which can be set explicitly
#'  using the matrix_class="dgCMatrix" list value. A very large matrix can
#'  be stored in a file and accessed by Monocle3 as if it were in memory.
#'  For this, Monocle3 uses the BPCells R package. Here the matrix_control
#'  list values are set to matrix_class="BPCells" and matrix_mode="dir".
#'  Then the count matrix is stored in a directory, on-disk, that's created
#'  by Monocle3 in the directory where you run Monocle3. This directory has
#'  a name with the form "monocle.bpcells.*.tmp" where the asterisk is a
#'  random string that makes the name unique. Do not remove this directory
#'  while Monocle3 is running! Monocle3 tries to remove the BPCells matrix
#'  directory when the preprocess_cds function finishes running; however,
#'  sometimes a matrix directory may persist after preprocess_cds finishes.
#'  In this case, the user must remove the directory after the session
#'  ends. For additional information about the matrix_control list, see the
#'  set_matrix_control help.
#' @return an updated cell_data_set object
#'
#' @examples
#'   \donttest{
#'     cell_metadata <- readRDS(system.file('extdata',
#'                                          'worm_embryo/worm_embryo_coldata.rds',
#'                                          package='monocle3'))
#'     gene_metadata <- readRDS(system.file('extdata',
#'                                          'worm_embryo/worm_embryo_rowdata.rds',
#'                                          package='monocle3'))
#'     expression_matrix <- readRDS(system.file('extdata',
#'                                              'worm_embryo/worm_embryo_expression_matrix.rds',
#'                                              package='monocle3'))
#'     cds <- new_cell_data_set(expression_data=expression_matrix,
#'                              cell_metadata=cell_metadata,
#'                              gene_metadata=gene_metadata)
#'     cds <- preprocess_cds(cds)
#'
#'     # For typical count matrices with a small to large number of cells,
#'     # we suggest that you use the default matrix_control list by not
#'     # not setting the matrix_control parameter. In this case, the
#'     # intermediate matrices are stored in memory as sparse matrices
#'     # in the dgCMatrix format, as they have in the past. It is also
#'     # possible to set the matrix_control list explicitly to use this
#'     # in-memory dgCMatrix format by setting the matrix_control parameter
#'     # list to
#'     #
#'       preprocess_cds(..., matrix_control=list(matrix_class='dgCMatrix'))
#'     #
#'     # For larger count matrices, we suggest that you try storing the
#'     # intermediate matrices as on-disk BPCells class objects by setting
#'     # the matrix_control parameter list as follows
#'     #
#'       preprocess_cds(..., matrix_control=list(matrix_class='BPCells', matrix_mode='dir'))
#'     #
#'   }
#'
#' @export
preprocess_cds <- function(cds,
                           method = c('PCA', "LSI"),
                           num_dim = 50,
                           norm_method = c("log", "size_only", "none"),
                           use_genes = NULL,
                           pseudo_count = NULL,
                           scaling = TRUE,
                           verbose = FALSE,
                           build_nn_index = FALSE,
                           nn_control = list(),
                           matrix_control = list()) {

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "method must be one of 'PCA' or 'LSI'")
  method <- match.arg(method)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(norm_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "norm_method must be one of 'log', 'size_only' or 'none'")
  norm_method <- match.arg(norm_method)

  assertthat::assert_that(assertthat::is.count(num_dim))

  if(!is.null(use_genes)) {
    assertthat::assert_that(is.character(use_genes))
    assertthat::assert_that(all(use_genes %in% row.names(rowData(cds))),
                            msg = paste("use_genes must be NULL, or all must",
                            "be present in the row.names of rowData(cds)"))
  }
  assertthat::assert_that(!is.null(size_factors(cds)),
             msg = paste("You must call estimate_size_factors before calling",
                         "preprocess_cds."))
  assertthat::assert_that(sum(is.na(size_factors(cds))) == 0,
                          msg = paste("One or more cells has a size factor of",
                                      "NA."))

  if(build_nn_index) {
    nn_control <- set_nn_control(mode=1,
                                 nn_control=nn_control,
                                 nn_control_default=get_global_variable('nn_control_annoy_cosine'),
                                 nn_index=NULL,
                                 k=NULL,
                                 verbose=verbose)
  }

  #ensure results from RNG sensitive algorithms are the same on all calls
  set.seed(2016)

  #
  # Notes:
  #   o  set_matrix_class commits BPCells queued
  #      operations to make FM but does not commit
  #      operations in counts(cds). Commit additional
  #      FM queued operations before submitting to
  #      the SVD function.
#message('\n==== preprocess_cds: set_matrix_class ====')
#  matrix_control_res <- set_pca_matrix_control(mat=SingleCellExperiment::counts(cds), matrix_control=matrix_control)
#  FM <- set_matrix_class(mat=SingleCellExperiment::counts(cds), matrix_control=matrix_control_res)

message('\n==== preprocess_cds: make FM matrix ====')
  FM <- SingleCellExperiment::counts(cds)

  FM <- normalize_expr_data(FM=FM, size_factors=size_factors(cds), norm_method=norm_method, pseudo_count=pseudo_count)

  if (nrow(FM) == 0) {
    stop("all rows have standard deviation zero")
  }

  if (!is.null(use_genes)) {
    FM <- FM[use_genes, ]
  }

  #
  # Notes:
  #   o  the functions save_transform_models/load_transform_models
  #      expect that the reduce_dim_aux slot consists of a S4Vectors::SimpleList
  #      that stores information about methods with the elements
  #        reduce_dim_aux[[method]][['model']] for the transform elements
  #        reduce_dim_aux[[method]][[nn_method]] for the nn index
  #      and depends on the elements within model and nn_method.
  #
  if(method == 'PCA') {
    cds <- initialize_reduce_dim_metadata(cds, 'PCA')
    cds <- initialize_reduce_dim_model_identity(cds, 'PCA')

    if (verbose) message("Remove noise by PCA ...")

    if(is(FM, 'dgCMatrix') || is(FM, 'dgeMatrix')) {

      if(verbose) {
        message('preprocess_cds: FM matrix class: ', class(FM))
      }

      fm_rowsums = Matrix::rowSums(FM)
      FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]
      irlba_res <- sparse_prcomp_irlba(Matrix::t(FM),
                                       n = min(num_dim,min(dim(FM)) - 1),
                                       center = scaling, scale. = scaling,
                                       verbose = verbose)
    }
    else
    if(is(FM, 'IterableMatrix')) {

      if(verbose) {
        message('preprocess_cds: FM matrix info:')
        message(show_matrix_info(get_matrix_info(FM), '  '), appendLF=FALSE)
      }

      fm_rowsums = BPCells::rowSums(FM)
      FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]
      irlba_res <- bpcells_prcomp_irlba(BPCells::t(FM),
                                        n = min(num_dim,min(dim(FM)) - 1),
                                        center = scaling, scale. = scaling,
                                        matrix_control=matrix_control,
                                        verbose = verbose)
    }
    else {
      stop('Unrecognized expression matrix class')
    }

    preproc_res <- irlba_res$x
    row.names(preproc_res) <- colnames(cds)
    SingleCellExperiment::reducedDims(cds)[[method]] <- as.matrix(preproc_res)

    irlba_rotation <- irlba_res$rotation
    row.names(irlba_rotation) <- rownames(FM)
    # we need svd_v downstream so
    # calculate gene_loadings in cluster_cells.R
    cds@reduce_dim_aux[['PCA']][['model']][['num_dim']] <- num_dim
    cds@reduce_dim_aux[['PCA']][['model']][['norm_method']] <- norm_method
    cds@reduce_dim_aux[['PCA']][['model']][['use_genes']] <- use_genes
    cds@reduce_dim_aux[['PCA']][['model']][['pseudo_count']] <- pseudo_count
    cds@reduce_dim_aux[['PCA']][['model']][['svd_v']] <- irlba_rotation
    cds@reduce_dim_aux[['PCA']][['model']][['svd_sdev']] <- irlba_res$sdev
    cds@reduce_dim_aux[['PCA']][['model']][['svd_center']] <- irlba_res$center
    cds@reduce_dim_aux[['PCA']][['model']][['svd_scale']] <- irlba_res$svd_scale
    # Note that prop_var_expl is the fraction of variance explained by the retained
    # PCs, not the fraction of total variance.
    cds@reduce_dim_aux[['PCA']][['model']][['prop_var_expl']] <- irlba_res$sdev^2 / sum(irlba_res$sdev^2)

    if(is(FM, 'IterableMatrix')) {
      rm_bpcells_dir(FM)
    }

    matrix_id <- get_unique_id(SingleCellExperiment::reducedDims(cds)[['PCA']])
    counts_identity <- get_counts_identity(cds)

    cds <- set_reduce_dim_matrix_identity(cds, 'PCA',
                                          'matrix:PCA',
                                          matrix_id,
                                          counts_identity[['matrix_type']],
                                          counts_identity[['matrix_id']],
                                          'matrix:PCA',
                                          matrix_id)
    cds <- set_reduce_dim_model_identity(cds, 'PCA',
                                         'matrix:PCA',
                                         matrix_id,
                                         'none',
                                         'none')

    if( build_nn_index ) {
      nn_index <- make_nn_index(subject_matrix=SingleCellExperiment::reducedDims(cds)[[method]], nn_control=nn_control, verbose=verbose)
      cds <- set_cds_nn_index(cds=cds, reduction_method=method, nn_index=nn_index, verbose=verbose)
    }
    else
      cds <- clear_cds_nn_index(cds=cds, reduction_method=method, nn_method='all')

  }
  else
  if(method == "LSI") {
    cds <- initialize_reduce_dim_metadata(cds, 'LSI')
    cds <- initialize_reduce_dim_model_identity(cds, 'LSI')

    fm_rowsums = Matrix::rowSums(FM)
    FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]

#    preproc_res <- tfidf(FM)
    tfidf_res <- tfidf(FM)
    preproc_res <- tfidf_res[['tf_idf_counts']]

    num_col <- ncol(preproc_res)
    irlba_res <- irlba::irlba(Matrix::t(preproc_res),
                              nv = min(num_dim,min(dim(FM)) - 1))

    preproc_res <- irlba_res$u %*% diag(irlba_res$d)
    row.names(preproc_res) <- colnames(cds)
    SingleCellExperiment::reducedDims(cds)[[method]] <- as.matrix(preproc_res)

    irlba_rotation = irlba_res$v
    row.names(irlba_rotation) = rownames(FM)
    cds@reduce_dim_aux[['LSI']][['model']][['num_dim']] <- num_dim
    cds@reduce_dim_aux[['LSI']][['model']][['norm_method']] <- norm_method
    cds@reduce_dim_aux[['LSI']][['model']][['use_genes']] <- use_genes
    cds@reduce_dim_aux[['LSI']][['model']][['pseudo_count']] <- pseudo_count
    cds@reduce_dim_aux[['LSI']][['model']][['log_scale_tf']] <- tfidf_res[['log_scale_tf']]
    cds@reduce_dim_aux[['LSI']][['model']][['frequencies']] <- tfidf_res[['frequencies']]
    cds@reduce_dim_aux[['LSI']][['model']][['scale_factor']] <- tfidf_res[['scale_factor']]
    cds@reduce_dim_aux[['LSI']][['model']][['col_sums']] <- tfidf_res[['col_sums']]
    cds@reduce_dim_aux[['LSI']][['model']][['row_sums']] <- tfidf_res[['row_sums']]
    cds@reduce_dim_aux[['LSI']][['model']][['num_cols']] <- tfidf_res[['num_cols']]
    cds@reduce_dim_aux[['LSI']][['model']][['svd_v']] <- irlba_rotation
    cds@reduce_dim_aux[['LSI']][['model']][['svd_sdev']] <- irlba_res$d/sqrt(max(1, num_col - 1))

    # we need svd_v downstream so
    # calculate gene_loadings in cluster_cells.R

    matrix_id <- get_unique_id(SingleCellExperiment::reducedDims(cds)[['LSI']])
    counts_identity <- get_counts_identity(cds)

    cds <- set_reduce_dim_matrix_identity(cds, 'LSI',
                                          'matrix:LSI',
                                          matrix_id,
                                          counts_identity[['matrix_type']],
                                          counts_identity[['matrix_id']],
                                          'matrix:LSI',
                                          matrix_id)
    cds <- set_reduce_dim_model_identity(cds, 'LSI',
                                         'matrix:LSI',
                                         matrix_id,
                                         'none',
                                         'none')

    if( build_nn_index ) {
      nn_index <- make_nn_index(subject_matrix=SingleCellExperiment::reducedDims(cds)[[method]], nn_control=nn_control, verbose=verbose)
      cds <- set_cds_nn_index(cds=cds, reduction_method=method, nn_index=nn_index, verbose=verbose)
    }
    else
      cds <- clear_cds_nn_index(cds=cds, reduction_method=method, nn_method='all')
  }

  if(!is.null(cds@reduce_dim_aux[['Aligned']]) && !is.null(cds@reduce_dim_aux[['Aligned']][['model']][['beta']])) {
    cds@reduce_dim_aux[['Aligned']][['model']][['beta']] <- NULL
  }

  cds
}


# Helper function to normalize the expression data prior to dimensionality
# reduction
normalize_expr_data <- function(FM, size_factors=NULL,
                                norm_method = c("log", "size_only", "none"),
                                pseudo_count = NULL) {

  assertthat::assert_that(!is.null(size_factors))
  assertthat::assert_that(length(size_factors) == ncol(FM))

  norm_method <- match.arg(norm_method)


  # If we're going to be using log, and the user hasn't given us a
  # pseudocount set it to 1 by default.
  if (is.null(pseudo_count)){
    if(norm_method == "log")
      pseudo_count <- 1
    else
      pseudo_count <- 0
  }

  if (is(FM, 'IterableMatrix')) {
    if(norm_method != 'log' || pseudo_count != 1) {
      stop('BPCells count matrix requires norm_method = \'log\' and pseudo_count = 1.')
    }
    FM <- BPCells::t(BPCells::t(FM)/size_factors)
    FM <- log1p(FM) / log(2)
  }
  else
  if (norm_method == "log") {
    # If we are using log, normalize by size factor before log-transforming

    FM <- Matrix::t(Matrix::t(FM)/size_factors)

    if (pseudo_count != 1 || is_sparse_matrix(FM) == FALSE){
      FM <- FM + pseudo_count
      FM <- log2(FM)
    }
    else {
      FM@x = log2(FM@x + 1)
    }

  }
  else if (norm_method == "size_only") {
    FM <- Matrix::t(Matrix::t(FM)/size_factors)
    FM <- FM + pseudo_count
  }
  return (FM)
}


# Andrew's tfidf
tfidf <- function(count_matrix, frequencies=TRUE, log_scale_tf=TRUE,
                  scale_factor=100000, block_size=2000e6) {
  # Use either raw counts or divide by total counts in each cell
  if (frequencies) {
    # "term frequency" method
    col_sums <- Matrix::colSums(count_matrix)
    tf <- Matrix::t(Matrix::t(count_matrix) / col_sums)
  } else {
    # "raw count" method
    col_sums <- NA
    tf <- count_matrix
  }

  # Either TF method can optionally be log scaled
  if (log_scale_tf) {
    if (frequencies) {
      tf@x = log1p(tf@x * scale_factor)
    } else {
      tf@x = log1p(tf@x * 1)
    }
  }

  # IDF w/ "inverse document frequency smooth" method
  num_cols <- ncol(count_matrix)
  row_sums <- Matrix::rowSums(count_matrix > 0)
  idf = log(1 + num_cols / row_sums)

  # Try to just to the multiplication and fall back on delayed array
  # TODO hopefully this actually falls back and not get jobs killed in SGE
  tf_idf_counts = tryCatch({
    tf_idf_counts = tf * idf
    tf_idf_counts
  }, error = function(e) {
    print(paste("TF*IDF multiplication too large for in-memory, falling back",
                "on DelayedArray."))
    options(DelayedArray.block.size=block_size)
    DelayedArray:::set_verbose_block_processing(TRUE)

    tf = DelayedArray::DelayedArray(tf)
    idf = as.matrix(idf)

    tf_idf_counts = tf * idf
    tf_idf_counts
  })

  rownames(tf_idf_counts) = rownames(count_matrix)
  colnames(tf_idf_counts) = colnames(count_matrix)
  tf_idf_counts = methods::as(tf_idf_counts, "sparseMatrix")
  return(list(tf_idf_counts=tf_idf_counts, frequencies=frequencies, log_scale_tf=log_scale_tf, scale_factor=scale_factor, col_sums=col_sums, row_sums=row_sums, num_cols=num_cols))
}


