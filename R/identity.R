# Functions for matrix and model identities.


#' @importFrom methods is
#' @importFrom SingleCellExperiment int_metadata
#' @importFrom SingleCellExperiment "int_metadata<-"
initialize_counts_metadata <- function(cds) {
  assertthat::assert_that(is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  if(is.null(int_metadata(cds))) {
    int_metadata(cds) <- list()
  }
  int_metadata(cds)[['counts_metadata']] <- list()

  return(cds)
}


#' @importFrom methods is
#' @importFrom SingleCellExperiment int_metadata
#' @importFrom SingleCellExperiment "int_metadata<-"
set_counts_identity <- function(cds, matrix_type, matrix_id) {
  assertthat::assert_that(is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(!is.null(int_metadata(cds)[['counts_metadata']]),
    msg = paste0('call initialize_counts_metadata() before set_counts_identity()'))                  

  int_metadata(cds)[['counts_metadata']][['identity']][['matrix_type']] <- matrix_type
  int_metadata(cds)[['counts_metadata']][['identity']][['matrix_id']] <- matrix_id

  return(cds)
}


#' @importFrom methods is
#' @importFrom SingleCellExperiment int_metadata
get_counts_identity <- function(cds) {
  assertthat::assert_that(is(cds, 'cell_data_set'),
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
#' @importFrom methods is
#' @importFrom SingleCellExperiment int_metadata
#' @importFrom SingleCellExperiment "int_metadata<-"
initialize_reduce_dim_metadata <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')) {
  assertthat::assert_that(is(cds, 'cell_data_set'),
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


#' @importFrom methods is
#' @importFrom SingleCellExperiment int_metadata
#' @importFrom SingleCellExperiment "int_metadata<-"
set_reduce_dim_matrix_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'),
                                           matrix_type, matrix_id,
                                           prev_matrix_type, prev_matrix_id,
                                           model_type, model_id) {
  assertthat::assert_that(is(cds, 'cell_data_set'),
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


#' @importFrom methods is
#' @importFrom SingleCellExperiment int_metadata
get_reduce_dim_matrix_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')) {
  assertthat::assert_that(is(cds, 'cell_data_set'),
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


#' @importFrom methods is
initialize_reduce_dim_model_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')) {
  assertthat::assert_that(is(cds, 'cell_data_set'),
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


#' @importFrom methods is
set_reduce_dim_model_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'),
                                          model_type, model_id,
                                          prev_model_type, prev_model_id,
                                          model_path='none') {
  assertthat::assert_that(is(cds, 'cell_data_set'),
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


#' @importFrom methods is
get_reduce_dim_model_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')) {
  assertthat::assert_that(is(cds, 'cell_data_set'),
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


#' @importFrom methods is
set_model_identity_path <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), model_path='none') {
  assertthat::assert_that(is(cds, 'cell_data_set'),
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


identity_as_string <- function( object_id ) {
  if(is.character(object_id)) {
    object_id_string <- object_id
  }
  else
  if(!is.null(object_id[['checksum']])) {
    object_id_string <- sprintf('%s  dim: %s', object_id[['checksum']], paste(object_id[['dim']], collapse=':'))
  }
  else {
    object_id_string <- object_id
  }

  return( object_id_string)
}


#' @title Report matrix and model identity information.
#'
#' @description Write the cell_data_set matrix and model 
#'   identity information to stdout.
#'
#' @details A matrix identity is a checksum that is stored 
#'   in the cell_data_set when a reduced dimension matrix is 
#'   created and when certain functions read count matrices
#'   into the cell_data_set, such as load_mm_data(). At the 
#'   same time, the same checksum is stored as the model 
#'   identity in order to link the model to its matrix.
#'   
#'   Additionally, Monocle3 stores the identity of the matrix
#'   from which the matrix was made. For example, in the case
#'   of a UMAP reduced dimension matrix made from a PCA
#'   reduced dimension matrix, the cell_data_set has the
#'   identities of both the UMAP and the PCA matrices. The
#'   UMAP identity is stored as 'matrix_id' and the PCA as
#'   'prev_matrix_id'. Similarly, the model and the previous 
#'   model identities are stored as 'model_id' and 
#'   'prev_model_id'. This allows one to trace a matrix
#'   to its origin, which may be helpful when a cell_data_set
#'   is partially reprocessed; for example, if preprocess_cds()
#'   is re-run but reduce_dimension() is not. Also, it may be
#'   helpful when transform models are loaded with the
#'   load_transform_models() function, in which case
#'   the matrix and model identities will differ.
#'   
#'   The identity of the model used to transform a matrix
#'   is stored with the matrix identity information as 
#'   'model_id'. Ordinarily, the matrix 'matrix_id' and
#'   'model_id' and the corresponding model 'model_id' will
#'   have the same string value. However, they differ when
#'   the preprocess_transform() and reduce_dim_transform()
#'   functions are used to transform a matrix.
#' 
#'   Notes:
#'   * Certain file and directory paths may be stored in the
#'   cell_data_set as identifiers.
#'
#'   * Checksums are calculated using the digest function in
#'   the digest package. The matrix dimensions are stored
#'   with the checksum.
#'
#'   * Matrix transformations such as subsetting and row and
#'   or column reordering do not affect the matrix identity.
#'
#'   * The matrix identity string is stored in the internal
#'   metadata slot of the cell_data_set and the model
#'   identity string is stored in the model object in the
#'   cds@reduce_dim_aux slot of the cell_data_set.
#'
#' @param cds the cell_data_set to use.
#'
#' @return Write identity information to stdout.
#'
#' @importFrom methods is
#' @export
identity_table <- function(cds) {
  assertthat::assert_that(is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  reduction_methods <- c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')
  counts_identity <- get_counts_identity(cds)
  write(sprintf('Count matrix identity'), stdout())
  write(sprintf('  matrix_id: %s', identity_as_string(counts_identity[['matrix_id']])), stdout())
  write(sprintf('  matrix_type: %s', counts_identity[['matrix_type']]), stdout())
  write('', stdout())

  write(sprintf('Reduced dimension matrix identity'), stdout())
  for(reduction_method in reduction_methods) {
    matrix_identity <- get_reduce_dim_matrix_identity(cds, reduction_method)
    if(!is.null(reducedDims(cds)[[reduction_method]]) && matrix_identity[['identity_exists']]) {
      write(sprintf('  %s', reduction_method), stdout())
      write(sprintf('    %s\t%s', 'matrix_type', matrix_identity[['matrix_type']]), stdout())
      write(sprintf('    %s\t%s', 'matrix_id', identity_as_string(matrix_identity[['matrix_id']])), stdout())
      write(sprintf('    %s\t%s', 'prev_matrix_type', matrix_identity[['prev_matrix_type']]), stdout())
      write(sprintf('    %s\t%s', 'prev_matrix_id', identity_as_string(matrix_identity[['prev_matrix_id']])), stdout())
      write(sprintf('    %s\t%s', 'model_type', matrix_identity[['model_type']]), stdout())
      write(sprintf('    %s\t%s', 'model_id', identity_as_string(matrix_identity[['model_id']])), stdout())
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
      write(sprintf('    %s\t%s', 'model_id', identity_as_string(model_identity[['model_id']])), stdout())
      write(sprintf('    %s\t%s', 'prev_model_type', model_identity[['prev_model_type']]), stdout())
      write(sprintf('    %s\t%s', 'prev_model_id', identity_as_string(model_identity[['prev_model_id']])), stdout())
      write(sprintf('    %s\t%s', 'model_path', model_identity[['model_path']]), stdout())
      write(sprintf('    %s\t%s', 'model_version', model_identity[['model_version']]), stdout())
    } else {
       write(sprintf('  %s\n    no identity information', reduction_method), stdout())
    }
    write('', stdout())
  }
}


