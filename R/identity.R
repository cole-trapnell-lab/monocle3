
initialize_counts_metadata <- function(cds) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  if(is.null(int_metadata(cds))) {
    int_metadata(cds) <- list()
  }
  int_metadata(cds)[['counts_metadata']] <- list()

  return(cds)
}


set_counts_identity <- function(cds, matrix_type, matrix_id) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(!is.null(int_metadata(cds)[['counts_metadata']]),
    msg = paste0('call initialize_counts_metadata() before set_counts_identity()'))                  

  int_metadata(cds)[['counts_metadata']][['identity']][['matrix_type']] <- matrix_type
  int_metadata(cds)[['counts_metadata']][['identity']][['matrix_id']] <- matrix_id

  return(cds)
}


get_counts_identity <- function(cds) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  if(is.null(int_metadata(cds)[['counts_metadata']])) {
    initialize_counts_metadata(cds)
  }

  if(!is.null(int_metadata(cds)[['counts_metadata']][['identity']][['matrix_id']])) {
    matrix_id <- int_metadata(cds)[['counts_metadata']][['identity']][['matrix_id']]
  } else {
    matrix_id <- 'none'
  }
  if(!is.null(int_metadata(cds)[['counts_metadata']][['identity']][['matrix_type']])) {
    matrix_type <- int_metadata(cds)[['counts_metadata']][['identity']][['matrix_type']]
  } else {
    matrix_type <- 'matrix:counts'
  }

  return(list(matrix_id=matrix_id, matrix_type=matrix_type))
}


# Note: int_metadata(cds) requires a list
initialize_reduce_dim_metadata <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  if(is.null(int_metadata(cds))) {
    int_metadata(cds) <- list()
  }
  if(is.null(int_metadata(cds)[['reduce_dim_metadata']])) {
    int_metadata(cds)[['reduce_dim_metadata']] <- list()
  }
  int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]] <- list()
  int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']] <- list()

  return(cds)
}


set_reduce_dim_matrix_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'),
                                           matrix_type, matrix_id,
                                           prev_matrix_type, prev_matrix_id,
                                           model_type, model_id) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  if(is.null(int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]])) {
    cds <- initialize_reduce_dim_metadata(cds=cds, reduction_method=reduction_method)
  }

  int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['matrix_type']] <- matrix_type
  int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['matrix_id']] <- matrix_id
  int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['prev_matrix_type']] <- prev_matrix_type
  int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['prev_matrix_id']] <- prev_matrix_id
  int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['model_type']] <- model_type
  int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['model_id']] <- model_id

  return(cds)
}


get_reduce_dim_matrix_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  if(is.null(int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]])) {
    cds <- initialize_reduce_dim_metadata(cds=cds, reduction_method=reduction_method)
    return(list(identity_exists=FALSE))
  }

  return(list(identity_exists=TRUE,
              matrix_type=int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['matrix_type']],
              matrix_id=int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['matrix_id']],
              prev_matrix_type=int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['prev_matrix_type']],
              prev_matrix_id=int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['prev_matrix_id']],
              model_type=int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['model_type']],
              model_id=int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['model_id']]))
}


initialize_reduce_dim_model_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  if(is.null(cds@reduce_dim_aux[[reduction_method]])) {
    cds@reduce_dim_aux[[reduction_method]] <- SimpleList()
  }
  cds@reduce_dim_aux[[reduction_method]][['model']] <- SimpleList()
  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']] <- SimpleList()

  return(cds)
}


set_reduce_dim_model_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'),
                                          model_type, model_id,
                                          prev_model_type, prev_model_id,
                                          model_path='none') {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  if(is.null(cds@reduce_dim_aux[[reduction_method]][['model']][['identity']])) {
    cds <- initialize_reduce_dim_model_identity(cds=cds, reduction_method=reduction_method)
  }

  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_type']] <- model_type
  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_id']] <- model_id
  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['prev_model_type']] <- prev_model_type
  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['prev_model_id']] <- prev_model_id
  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_path']] <- model_path

  global_variable_name <- list(PCA='reduce_dim_pca_model_version',
                               LSI='reduce_dim_lsi_model_version',
                               Aligned='reduce_dim_aligned_model_version',
                               tSNE='reduce_dim_tsne_model_version',
                               UMAP='reduce_dim_umap_model_version')

  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_version']] <- get_global_variable(global_variable_name[[reduction_method]])

  return(cds)
}


get_reduce_dim_model_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  if(is.null(cds@reduce_dim_aux[[reduction_method]][['model']][['identity']])) {
    cds <- initialize_reduce_dim_model_identity(cds=cds, reduction_method=reduction_method)
    return(list(identity_exists=FALSE))
  }

  return(list(identity_exists=TRUE,
              model_type=cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_type']],
              model_id=cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_id']],
              prev_model_type=cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['prev_model_type']],
              prev_model_id=cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['prev_model_id']],
              model_path=cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_path']],
              version=cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_version']]))
}


set_model_identity_path <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), model_path='none') {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  if(is.null(cds@reduce_dim_aux[[reduction_method]][['model']][['identity']])) {
    cds <- initialize_reduce_dim_model_identity(cds=cds, reduction_method=reduction_method)
  }

  if(!is.null(cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_path']])) {
    cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_path']] <- model_path
  } else {
    cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_path']] <- 'none'
  }

  return(cds)
}


#' @title Report the matrix and model identity information.
#'
#' @description Write the cell_data_set matrix and model identity
#'   information to stdout.
#'
#' @details A matrix identity is a time stamp that is made when
#'   a reduced dimension matrix is created. The time stamp is an
#'   MD5 checksum string made from the time and a counter, which
#'   make it unique. Some functions that read count
#'   matrices also make a matrix identity, for example,
#'   load_mm_data(). In addition to the matrix identity for the
#'   reduced dimension matrix, Monocle3 stores the identity of the
#'   matrix from which the matrix was made. For example, in the
#'   case of a UMAP reduced dimension matrix made from a PCA reduced
#'   dimension matrix, the cds has the identity of both the UMAP
#'   and the PCA matrices. The identity_table() function stores
#'   the UMAP identity as 'matrix_id' and the PCA identity as
#'   'prev_matrix_id'. This gives a way to trace a matrix to its
#'   origin. The ability to trace the matrix may be helpful when
#'   a cds is saved and restored, and then partially reprocessed.
#'
#'   The matrix identity string is stored in the internal metadata
#'   slot of the cell_data_set.
#'
#'   Note that the same matrix created at different times will
#'   have different matrix identities. We chose to use a time
#'   stamp because the alternative of making a checksum from a
#'   very large matrix can take a long time.
#'
#'
#' @param cds the cell_data_set to use.
#'
#' @return Write identity information to stdout.
#'
#' @export
identity_table <- function(cds) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))

  reduction_methods <- c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')
  counts_identity <- get_counts_identity(cds)
  write(sprintf('Count matrix identity'), stdout())
  write(sprintf('  matrix_id: %s', counts_identity[['matrix_id']]), stdout())
  write(sprintf('  matrix_type: %s', counts_identity[['matrix_type']]), stdout())
  write('', stdout())

  write(sprintf('Reduced dimension matrix identity'), stdout())
  for(reduction_method in reduction_methods) {
    matrix_identity <- get_reduce_dim_matrix_identity(cds, reduction_method)
    if(!is.null(reducedDims(cds)[[reduction_method]]) && matrix_identity[['identity_exists']]) {
      write(sprintf('  %s', reduction_method), stdout())
      write(sprintf('    %s\t%s', 'matrix_type', matrix_identity[['matrix_type']]), stdout())
      write(sprintf('    %s\t%s', 'matrix_id', matrix_identity[['matrix_id']]), stdout())
      write(sprintf('    %s\t%s', 'prev_matrix_type', matrix_identity[['prev_matrix_type']]), stdout())
      write(sprintf('    %s\t%s', 'prev_matrix_id', matrix_identity[['prev_matrix_id']]), stdout())
      write(sprintf('    %s\t%s', 'model_type', matrix_identity[['model_type']]), stdout())
      write(sprintf('    %s\t%s', 'model_id', matrix_identity[['model_id']]), stdout())
    } else {
       write(sprintf('  %s\n    no identity information', reduction_method), stdout())
    }
    write('', stdout())
  }

  write(sprintf('Reduced dimension model identity'), stdout())
  for(reduction_method in reduction_methods) {
    model_identity <- get_reduce_dim_model_identity(cds, reduction_method)
    if(model_identity[['identity_exists']]) {
      write(sprintf('  %s', reduction_method), stdout())
      write(sprintf('    %s\t%s', 'model_type', model_identity[['model_type']]), stdout())
      write(sprintf('    %s\t%s', 'model_id', model_identity[['model_id']]), stdout())
      write(sprintf('    %s\t%s', 'prev_model_type', model_identity[['prev_model_type']]), stdout())
      write(sprintf('    %s\t%s', 'prev_model_id', model_identity[['prev_model_id']]), stdout())
      write(sprintf('    %s\t%s', 'model_path', model_identity[['model_path']]), stdout())
      write(sprintf('    %s\t%s', 'model_version', model_identity[['model_version']]), stdout())
    } else {
       write(sprintf('  %s\n    no identity information', reduction_method), stdout())
    }
    write('', stdout())
  }
}


