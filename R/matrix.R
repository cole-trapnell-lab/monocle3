################################################################################
# BPCells
################################################################################

# Notes:
#   o  the SummarizedExperiment Assay slot information states
#        o  the developer needs to be careful to implement endomorphisms with copy-on-modify semantics
#        o  an Assays concrete subclass needs to implement lossless back and forth coercion from/to
#           SimpleList
# 
# Questions:
#   o  are MatrixDir directories relocatable? The directory is not relocatable while the R MatrixDir
#      object exists. The saved directory that persists after the R session ends can be relocated before
#      being 'reloaded' with the open_matrix_dir() function.
#   o  do BPCells matrices 'respect' copy-on-modify semantics?
#   o  how does one deal with copy-on-semantics and BPCells queued (delayed) operations?
# 
# Ben Parks notes:
# 
#   o  on copy-on-modify semantics
#       First, for the copy-on-modify semantics -- I think BPCells is fine in this respect. As a
#        general rule, BPCells does not keep open file handles sitting around, and files are
#        re-opened each time a data-reading operation is taking place then closed at the end of
#        the operation. Additionally, the only way to modify the contents of an existing matrix
#        on disk or in memory is to write a new matrix while explicitly asking to overwrite the
#        old one. The end result is that the on-disk data will (intentionally) persist across R
#        sessions, and multiple R objects can happily read from the same data source without
#        worrying that one of them will accidentally modify data on disk. Since BPCells objects
#        are basically plain R objects with one field listing a file path, modifications to the
#        BPCells objects automatically support copy-on-modify semantics. The only way to break
#        that is by intentionally changing the contents of the underlying files on disk, which
#        BPCells will never do unless explicitly asked to overwrite data.
# 
#        The only note I'd give is BPCells doesn't support the [<- operator, which would
#        otherwise be the most problematic. Row and column names can be modified, but the
#        modifications take place only in R and only go to disk when explicitly writing the
#        matrix again.
# 
#   o on getting BPCells matrix path/type
# 
#        As for getting the directory path from queued operations, your tactic should mostly
#        work with the exception of sparse matrix multiplies or rbind/cbind operations.
#        There's currently an unexported method called matrix_inputs which does basically the
#        same thing -- it unwraps a single layer of queued operations returning a list of inputs
#        (usually but not always length 1). I might make a helper function that can perform this
#        recursively to get a list of the raw matrix inputs across arbitrarily many operations.
#        For now, I'd recommend trying to directly call the unexported method
#        BPCells:::matrix_inputs within your bpcells_find_base_matrix function, and we can swap
#        that for a more robust and exported method down the line.
# 
# BPCells methods (not all)
#      length(mat)
#      length(mat[,1])
#      length(mat[1,])
#      names(mat)
#      names(mat[,1])
#      names(mat[1,])
#      dim(mat)
#      rownames(mat)
#      colnames(mat)
#      dimnames(mat)
#      storage_order(mat)
#      rowSums(mat)
#      colSums(mat)
#      rowMeans(mat)
#      colMeans(mat)
#      selection_index()
#      rbind2()
#      cbind2()
#      transpose_storage_order(mat)
#      write_matrix_memory()
#      write_matrix_dir()
#      open_matrix_dir()
#      write_matrix_hdf5()
#      open_matrix_hdf5()
#      open_matrix_10x_hdf5()
#      write_matrix_10x_hdf5()
#      open_matrix_anndata_hdf5()
#      convert_matrix_type()
#      matrix_stats()


select_matrix_parameter_value <- function(parameter, matrix_control, matrix_control_default, default_value) {
  if(!is.null(matrix_control[[parameter]])) {
    return(matrix_control[[parameter]])
  }
  else
  if(!is.null(matrix_control_default[[parameter]])) {
    return(matrix_control_default[[parameter]])
  }
  return(default_value)
}


# Usage
#   matrix_class: default: 'dgCMatrix'
#   matrix_class: 'BPCells'
#     matrix_mode: 'mem'  default: 'dir'
#       matrix_type: 'uint32_t', 'float', 'double'  NOTE: 'uint32_t' is valid for assay matrix only
#       matrix_compress: TRUE, FALSE default: FALSE
#     matrix_mode: 'dir'
#       matrix_type: 'uint32_t', 'float', 'double' default: 'uint32_t'
#       matrix_path: <path to directory or file> default: 'NULL' -> temporary directory in pwd
#       matrix_compress: TRUE, FALSE default: FALSE
#       matrix_buffer_size: <integer> default: 8192L
check_matrix_control <- function(matrix_control=list(), control_type=c('any', 'pca'), check_conditional=FALSE) {
  control_type <- match.arg(control_type)
  assertthat::assert_that(is.list(matrix_control))
  assertthat::assert_that(is.logical(check_conditional))

  # If matrix_control list is zero length then return,
  # otherwise, matrix_control[['matrix_class']] must
  # be defined.
  if(length(matrix_control) == 0) {
    return()
  }
  else
  if(is.null(matrix_control[['matrix_class']])) {
    stop('matrix_control[[\'matrix_class\']] is not defined')
  }
  
  allowed_control_parameters <- c('matrix_class',
                                  'matrix_mode',
                                  'matrix_type',
                                  'matrix_compress',
                                  'matrix_path',
                                  'matrix_buffer_size',
                                  'show_values')

  allowed_matrix_class <- c('dgCMatrix', 'BPCells')
  allowed_matrix_type_any <- c('uint32_t', 'float', 'double')
  allowed_matrix_type_pca <- c('float', 'double')
  allowed_matrix_mode <- c('mem', 'dir') 

  control_name <- if(control_type == 'any') 'any_control' else 'pca_control'

  error_string <- ''

  if(!all(names(matrix_control) %in% allowed_control_parameters)) {
    error_string <- paste0('invalid ', control_name, ' control parameter')
  }

  if(check_conditional == FALSE) {
    if(!(is.null(matrix_control[['matrix_class']])) &&
       !(matrix_control[['matrix_class']] %in% allowed_matrix_class)) {
      error_string <- paste0('\ninvalid matrix_class "', matrix_control[['matrix_class']], '"')
    }

    if(!(is.null(matrix_control[['matrix_mode']])) &&
       !(matrix_control[['matrix_mode']] %in% allowed_matrix_mode)) {
      error_string <- paste0('\ninvalid matrix_mode "', matrix_control[['matrix_mode']], '"')
    }

    allowed_values <- if(control_type == 'any') allowed_matrix_type_any else allowed_matrix_type_pca
    if(!(is.null(matrix_control[['matrix_type']])) &&
       !(matrix_control[['matrix_type']] %in% allowed_values)) {
      error_string <- paste0('\ninvalid ', control_name, ' matrix_type "', matrix_control[['matrix_type']], '"')
    }

    if(!(is.null(matrix_control[['matrix_compress']])) &&
       !(is.logical(matrix_control[['matrix_compress']]))) {
      error_string <- paste0('\nmatrix_compress value must be a logical type')
    }

    if(!(is.null(matrix_control[['matrix_path']])) &&
       !(is.character(matrix_control[['matrix_path']]))) {
      error_string <- paste0('\nmatrix_path value must be a character type')
    }

    if(!(is.null(matrix_control[['matrix_buffer_size']])) &&
       !(is.integer(matrix_control[['matrix_buffer_size']]))) {
      error_string <- paste0('\nmatrix_buffer_size value must be an integer type')
    }

  }
  else {
    # Check matrix_class value.
    if(is.null(matrix_control[['matrix_class']])) {
      error_string <- '\nmatrix_class not set'
    }
    else
    if(!(matrix_control[['matrix_class']] %in% allowed_matrix_class)) {
      error_string <- paste0('\ninvalid matrix_class "', matrix_control[['matrix_class']], '\n')
    }
  
    if(matrix_control[['matrix_class']] == 'BPCells') {
      # Check matrix_type value.
      allowed_values <- if(control_type == 'any') allowed_matrix_type_any else allowed_matrix_type_pca
      if(!(matrix_control[['matrix_type']] %in% allowed_values)) {
        error_string <- paste0('\nbad ', control_name, ' matrix_type "', matrix_control[['matrix_type']], '"\n')
      }
  
      # Check matrix_compress value.
      if(!is.logical(matrix_control[['matrix_compress']])) {
        error_string <- '\nmatrix_compress must be as logical type'
      }
  
      # Check matrix_mode value.
      if(!(matrix_control[['matrix_mode']] %in% allowed_matrix_mode)) {
        error_string <- paste0('\ninvalid matrix_mode "', matrix_control[['matrix_mode']], '"')
      }
  
      if(matrix_control[['matrix_mode']] == 'dir') {
        # Check matrix_path value.
        if(!(is.character(matrix_control[['matrix_path']]))) {
          error_string <- paste0('\nbad matrix_path "', matrix_control[['matrix_path']], '"')
        }
  
        # Check matrix_buffer_size.
        if(!(is.integer(matrix_control[['matrix_buffer_size']]))) {
          error_string <- paste0('\nmatrix_buffer_size must be an integer')
        }
      }
    }
  }

  if(error_string != '') {
    stop(stringr::str_trim(error_string))
  }
}


# Usage
#   matrix_class: default: 'dgCMatrix'
#   matrix_class: 'BPCells'
#     matrix_mode: 'mem'  default: 'dir'
#       matrix_type: 'uint32_t', 'float', 'double'
#       matrix_compress: TRUE, FALSE default: FALSE
#     matrix_mode: 'dir'
#       matrix_type: 'uint32_t', 'float', 'double' default: 'uint32_t'
#       matrix_path: <path to directory or file> default: 'NULL' -> temporary directory in pwd
#       matrix_compress: TRUE, FALSE default: FALSE
#       matrix_buffer_size: <integer> default: 8192L
# Notes:
#   o  modification to any of set_assay_control, set_pca_control,
#      or set_pca_control_default may necessitate modifications
#      to all of them.

#' Verify and set the matrix_control parameter list.
#'
#' @description Verifies and sets the list of parameter values
#'   that is used to make the counts matrix that is stored in the
#'   cell_data_set. To see the default values,
#'   call "set_assay_control(assay_control=list(show_values=TRUE))".
#'   "show_values=TRUE" can be used in functions that have the
#'   assay_control list parameter, in which case the function will
#'   show the assay_control values to be used and then stop.
#' @param assay_control Input control list.
#' @return assay_control Output control list.
#'
#' @section assay_control parameters:
#' \describe{
#'   \item{matrix_class}{Specifies the matrix class to use for
#'      matrix storage. The acceptable values are "dgCMatrix"
#'      and "BPCells".}
#'   \item{matrix_type}{Specifies whether to store the matrix
#'      values as unsigned 32-bit integers
#'      (matrix_type="uint32_t"), single precision "floats"
#'      (matrix_type="float"), or double precision "doubles"
#'      (matrix_type="double"). Unsigned 32-bit integers are
#'      suitable for storing non-negative values such as counts.
#'      "matrix_type" is used only for BPCells class matrices.}
#'   \item{matrix_mode}{Specifies whether to store the BPCells
#'      class matrix in memory (matrix_mode="mem") or on
#'      disk (matrix_mode="dir"). "matrix_mode" is used only
#'      for BPCells class matrices.}
#'   \item{matrix_path}{Specifies the directory where the
#'      BPCells on-disk matrix data are stored in a
#'      sub-directory with a randomized name. The default is
#'      in the directory where R is running. "matrix_path" is
#'      used only for BPCells class matrices with
#'      matrix_mode="dir". For example, if matrix_path is set
#'      to "/tmp" Monocle3 will create a directory with a
#'      name that has the form "monocle.bpcells.*.tmp" in
#'      "/tmp". The asterisk represents a random string that
#'      makes the names unique.}
#'   \item{matrix_compress}{Specifies whether to use bit-packing
#'      compression to store BPCells matrix values. This is
#'      most effective for values stored as matrix_type="uint32_t"
#'      but has some benefit for values stored as matrix_type="float"
#'      and matrix_type="double". Compression reduces storage
#'      requirements but increases run time. "matrix_compress" is
#'      used only for BPCells class matrices.}
#'   \item{matrix_buffer_size}{Specifies how many items of
#'      data to buffer in memory before flushing to disk. This
#'      is used for matrix_class="BPCells" with matrix_mode="dir".}
#' @export
set_matrix_control <- function(matrix_control=list(), matrix_control_default=list(), control_type=c('any', 'pca')) {
  control_type <- match.arg(control_type)

  check_matrix_control(matrix_control=matrix_control, control_type=control_type, check_conditional=FALSE)
  check_matrix_control(matrix_control=matrix_control_default, control_type=control_type, check_conditional=FALSE)

  #
  # Last resort fall-back parameter values.
  #
  default_matrix_class <- 'dgCMatrix'
  default_matrix_mode <- NULL
  default_matrix_type <- NULL
  default_matrix_compress <- NULL
  default_matrix_path <- NULL
  default_matrix_buffer_size <- NULL

  matrix_control_out = list()

  matrix_control_out[['matrix_class']] <- select_matrix_parameter_value('matrix_class', matrix_control, matrix_control_default, default_matrix_class)

  if(matrix_control_out[['matrix_class']] == 'BPCells') {
     matrix_control_out[['matrix_mode']] <- select_matrix_parameter_value('matrix_mode', matrix_control, matrix_control_default, default_matrix_mode)

    if(matrix_control_out[['matrix_mode']] == 'mem') {
       matrix_control_out[['matrix_type']] <- select_matrix_parameter_value('matrix_type', matrix_control, matrix_control_default, default_matrix_type)
       matrix_control_out[['matrix_compress']] <- select_matrix_parameter_value('matrix_compress', matrix_control, matrix_control_default, default_matrix_compress)
    }
    else
    if(matrix_control_out[['matrix_mode']] == 'dir') {
       matrix_control_out[['matrix_type']] <- select_matrix_parameter_value('matrix_type', matrix_control, matrix_control_default, default_matrix_type)
       matrix_control_out[['matrix_path']] <- select_matrix_parameter_value('matrix_path', matrix_control, matrix_control_default, default_matrix_path)
       matrix_control_out[['matrix_compress']] <- select_matrix_parameter_value('matrix_compress', matrix_control, matrix_control_default, default_matrix_compress)
       matrix_control_out[['matrix_buffer_size']] <- select_matrix_parameter_value('matrix_buffer_size', matrix_control, matrix_control_default, default_matrix_buffer_size)
    }
  }

  check_matrix_control(matrix_control=matrix_control_out, control_type=control_type, check_conditional=TRUE)

  #
  # Set BPCells out-of-core file/directory name.
  #
  if(matrix_control_out[['matrix_class']] == 'BPCells' &&
     matrix_control_out[['matrix_mode']] == 'dir' &&
     is.null(matrix_control_out[['matrix_path']])) {
        matrix_control_out[['matrix_path']] <- '.'
  }

  #
  # Display matrix_control list if matrix_control[['show_values']] <- TRUE.
  #
  if(!is.null(matrix_control[['show_values']]) && matrix_control[['show_values']] == TRUE)
  {
    show_matrix_control(matrix_control=matrix_control_out, '  matrix_control: ')
    stop_no_noise()
  }

  return(matrix_control_out)
}


# Report matrix_control list values.
show_matrix_control <- function(matrix_control, label=NULL) {
  message('matrix_control:')
  indent <- ''
  if(!is.null(label)) {
    indent <- '  '
  }
  message(indent, '  matrix_class: ', ifelse(!is.null(matrix_control[['matrix_class']]), matrix_control[['matrix_class']], as.character(NA)))

  if(matrix_control[['matrix_class']] == 'dgCMatrix') {
    # nop
    invisible(NULL)
  }
  else
  if(matrix_control[['matrix_class']] == 'BPCells') {
    if(matrix_control[['matrix_mode']] == 'mem') {
      message(indent, '  matrix_mode: ', ifelse(!is.null(matrix_control[['matrix_mode']]), matrix_control[['matrix_mode']], as.character(NA)))
      message(indent, '  matrix_type: ', ifelse(!is.null(matrix_control[['matrix_type']]), matrix_control[['matrix_type']], as.character(NA)))
      message(indent, '  matrix_compress: ', ifelse(!is.null(matrix_control[['matrix_compress']]), matrix_control[['matrix_compress']], as.character(NA)))
    }
    else
    if(matrix_control[['matrix_mode']] == 'dir') {
      message(indent, '  matrix_mode: ', ifelse(!is.null(matrix_control[['matrix_mode']]), matrix_control[['matrix_mode']], as.character(NA)))
      message(indent, '  matrix_type: ', ifelse(!is.null(matrix_control[['matrix_type']]), matrix_control[['matrix_type']], as.character(NA)))
      message(indent, '  matrix_path: ', ifelse(!is.null(matrix_control[['matrix_path']]), matrix_control[['matrix_path']], as.character(NA)))
      message(indent, '  matrix_compress: ', ifelse(!is.null(matrix_control[['matrix_compress']]), matrix_control[['matrix_compress']], as.character(NA)))
      message(indent, '  matrix_buffer_size: ', ifelse(!is.null(matrix_control[['matrix_buffer_size']]), matrix_control[['matrix_buffer_size']], as.character(NA)))
    }
    else {
      stop('show_matrix_control: unsupported matrix class/mode/...\'', matrix_control[['method']], '\'')
    }
  }
  else {
    stop('show_matrix_control: unsupported matrix class/mode/...\'', matrix_control[['method']], '\'')
  }
}


push_matrix_path <- function(mat) {
  matrix_info <- get_matrix_info(mat=mat)
  if(matrix_info[['matrix_class']] == 'BPCells' &&
     matrix_info[['matrix_mode']] == 'dir') {
    matrix_path <- matrix_info[['matrix_path']]
    x <- get_global_variable('monocle_gc_matrix_path')
    x <- append(x, matrix_path)
    set_global_variable('monocle_gc_matrix_path', x)
  }
}


#
# Return the 'base' BPCells matrix.
# Notes:
#   o uses
#      o get base matrix class
#      o get base matrix slots and their information, e.g., @type and @dir..
#
bpcells_find_base_matrix <- function(mat) {
  if(length(BPCells:::matrix_inputs(mat)) > 0) {
    return(bpcells_find_base_matrix(BPCells:::matrix_inputs(mat)[[1]]))
  }
  return(mat)
}


#
# Known matrix types:
#   o  dense matrix: matrix
#   o  sparse matrix: CsparseMatrix, dgCMatrix, ...
#   o  BPCells matrix:
#        o  matrix_mode: 'mem', 'dir'
#        o  matrix_type: uint32_t, float, double
#        o  matrix_compress: TRUE, FALSE
#        o  matrix_path
#        o  matrix_buffer_size
#
get_matrix_class <- function(mat) {
message('get_matrix_class: start')
message('str(mat): ')
message(str(mat))

  nmatch <- 0
  matrix_info <- list()
  if(is(mat, 'matrix')) {
    matrix_info[['matrix_class']] <- 'r_dense_matrix'
    nmatch <- nmatch + 1
  }
  if(any(class(mat) %in% c('dgCMatrix'))) {
    matrix_info[['matrix_class']] <- 'dgCMatrix'
    nmatch <- nmatch + 1
  }
  if(any(class(mat) %in% c('lgCMatrix'))) {
    matrix_info[['matrix_class']] <- 'lgCMatrix'
    nmatch <- nmatch + 1
  }
  if(is(mat, 'dgTMatrix')) {
    matrix_info[['matrix_class']] <- 'dgTMatrix'
    nmatch <- nmatch + 1
  }
  if(is(mat, 'IterableMatrix')) {
    matrix_info[['matrix_class']] <- 'BPCells'
    nmatch <- nmatch + 1
  }

  if(nmatch == 0) {
    stop('get_matrix_class: unrecognized matrix class')
  }
  if(nmatch > 1) {
    stop('get_matrix_class: ambiguous matrix class')
  }
message('get_matrix_class: end')
  return(matrix_info)
}


get_matrix_info <- function(mat) {
  matrix_info <- get_matrix_class(mat)

  if(is.null(matrix_info[['matrix_class']])) {
    message('bad matrix info -- dropping into browser')
    browser()
  }

  if(matrix_info[['matrix_class']] != 'BPCells') {
    return(matrix_info)
  }

  bmat <- bpcells_find_base_matrix(mat)

  # In-memory iterable matrix object classes.
  #   UnpackedMatrixMem_uint32_t
  #   UnpackedMatrixMem_float
  #   UnpackedMatrixMem_double
  #   PackedMatrixMem_uint32_t
  #   PackedMatrixMem_float
  #   PackedMatrixMem_double
  #
  # On-disk iterable matrix object classes and relevant slots.
  #   MatrixDir: slots: dir, compressed, buffer_size, type
  #     compressed: logical
  #     type: character: 'uint32_t', 'float', 'double'
  #     buffer_size: int
  #     dir: character: path
  if(!(class(bmat) %in% c('UnpackedMatrixMem_uint32_t', 'UnpackedMatrixMem_float',
                          'UnpackedMatrixMem_double', 'PackedMatrixMem_uint32_t',
                          'PackedMatrixMem_float', 'PackedMatrixMem_double',
                          'MatrixDir'))) {
    stop('get_matrix_info: unrecognized BPCells matrix class \"', class(mat), '\"')
    return(NULL)
  }

  matrix_info[['matrix_class']] <- 'BPCells'

  if(class(bmat) == 'UnpackedMatrixMem_uint32_t') {
    matrix_info[['matrix_mode']] <- 'mem'
    matrix_info[['matrix_type']] <- 'uint32_t'
    matrix_info[['matrix_compress']] <- FALSE
  }
  else
  if(class(bmat) == 'UnpackedMatrixMem_float') {
    matrix_info[['matrix_mode']] <- 'mem'
    matrix_info[['matrix_type']] <- 'float'
    matrix_info[['matrix_compress']] <- FALSE
  } 
  else
  if(class(bmat) == 'UnpackedMatrixMem_double') {
    matrix_info[['matrix_mode']] <- 'mem'
    matrix_info[['matrix_type']] <- 'double'
    matrix_info[['matrix_compress']] <- FALSE
  } 
  else
  if(class(bmat) == 'PackedMatrixMem_uint32_t') {
    matrix_info[['matrix_mode']] <- 'mem'
    matrix_info[['matrix_type']] <- 'uint32_t'
    matrix_info[['matrix_compress']] <- TRUE
  } 
  else
  if(class(bmat) == 'PackedMatrixMem_float') {
    matrix_info[['matrix_mode']] <- 'mem'
    matrix_info[['matrix_type']] <- 'float'
    matrix_info[['matrix_compress']] <- TRUE
  } 
  else
  if(class(bmat) == 'PackedMatrixMem_double') {
    matrix_info[['matrix_mode']] <- 'mem'
    matrix_info[['matrix_type']] <- 'double'
    matrix_info[['matrix_compress']] <- TRUE
  } 
  else
  if(class(bmat) == 'MatrixDir') {
    matrix_info[['matrix_mode']] <- 'dir'
    matrix_info[['matrix_type']] <- bmat@type
    matrix_info[['matrix_compress']] <- bmat@compressed
    matrix_info[['matrix_buffer_size']] <- bmat@buffer_size
    matrix_info[['matrix_path']] <- bmat@dir
  }

  return(matrix_info)
}


show_matrix_info <- function(matrix_info, indent='') {
  message('matrix_info:')
  if(!is.null(matrix_info[['matrix_class']])) {
    message(paste0(indent, 'class:       ', matrix_info[['matrix_class']]))
  }
  if(!is.null(matrix_info[['matrix_mode']])) { 
    message(paste0(indent, 'mode:        ', matrix_info[['matrix_mode']]))
  }
  if(!is.null(matrix_info[['matrix_type']])) { 
    message(paste0(indent, 'type:        ', matrix_info[['matrix_type']]))
  }
  if(!is.null(matrix_info[['matrix_compress']])) { 
    message(paste0(indent, 'compress:    ', matrix_info[['matrix_compress']]))
  }
  if(!is.null(matrix_info[['matrix_buffer_size']])) { 
    message(paste0(indent, 'buffer_size: ', matrix_info[['matrix_buffer_size']]))
  }
  if(!is.null(matrix_info[['matrix_path']])) { 
    message(paste0(indent, 'path:        ', matrix_info[['matrix_path']]))
  }
}


# set_matrix_class
#  Notes:
#    o  Cast an input matrix into a class given or inferred from the
#       matrix_control, which must be complete in the sense that all
#       values required for the class are given explicitly.
#    o  Always make a new matrix, even if the matrix is a BPCells
#       class and the parameters are the same because we may need to
#       commit the queued operations.
#   
set_matrix_class <- function(mat, matrix_control=list()) {

  # Check matrix_control list.
  check_matrix_control(matrix_control=matrix_control, control_type='any', check_conditional=TRUE)

  # Get input matrix info.
  matrix_info <- get_matrix_info(mat)


message('set_matrix_class: matrix_info: in:')
show_matrix_info(matrix_info, indent='  ')
message('')
message('set_matrix_class: matrix_control in:')
show_matrix_control(matrix_control)
message('')

  if(matrix_info[['matrix_class']] == 'dgTMatrix') {
    mat <- as(mat, 'dgCMatrix')
    matrix_info[['matrix_class']] <- 'dgCMatrix'
  }

  mat_out <- mat
  if(matrix_control[['matrix_class']] == 'dgCMatrix') {
    if(matrix_info[['matrix_class']] != 'dgCMatrix') {
      if(matrix_info[['matrix_class']] == 'BPCells') {
        mat_out <- as(mat, 'dgCMatrix')
      }
      else {
        mat_out <- as(mat, 'dgCMatrix')
      }
    }
  }
  else
  if(matrix_control[['matrix_class']] == 'BPCells') {
    matrix_class_d <- matrix_control[['matrix_class']]
    matrix_mode_d <- matrix_control[['matrix_mode']]
    matrix_type_d <- matrix_control[['matrix_type']]
    matrix_compress_d <- matrix_control[['matrix_compress']]
    matrix_path_d <- matrix_control[['matrix_path']]
    matrix_buffer_size_d <- matrix_control[['matrix_buffer_size']]
    if(matrix_info[['matrix_class']] == 'dgCMatrix') {
      if(matrix_mode_d == 'mem') {
        mat_out <- BPCells::write_matrix_memory(mat=BPCells::convert_matrix_type(mat, type=matrix_type_d),
                                                compress=matrix_compress_d)
      }
      else
      if(matrix_mode_d == 'dir') {
        mat_dir <- tempfile(pattern=paste0('monocle.bpcells.',
                                           format(Sys.Date(), format='%Y%m%d'), '.'),
                            tmpdir=matrix_path_d,
                            fileext='.tmp')[[1]]
        mat_out <- BPCells::write_matrix_dir(mat=BPCells::convert_matrix_type(mat, type=matrix_type_d),
                                             dir=mat_dir,
                                             compress=matrix_compress_d,
                                             buffer_size=matrix_buffer_size_d,
                                             overwrite=FALSE)
        push_matrix_path(mat_out)
      }
      else {
        stop('set_matrix_class: unrecognized matrix_info[[\'matrix_mode\']]')
      }
    }
    else
    if(matrix_info[['matrix_class']] == 'BPCells') {
      if(matrix_mode_d == 'mem') {
        mat_out <- BPCells::write_matrix_memory(mat=BPCells::convert_matrix_type(mat, type=matrix_type_d),
                                                compress=matrix_compress_d)
      }
      else
      if(matrix_mode_d == 'dir') {
        mat_dir <- tempfile(pattern=paste0('monocle.bpcells.', 
                                           format(Sys.Date(), 
                                                  format='%Y%m%d'), '.'), 
                            tmpdir=matrix_path_d,
                            fileext='.tmp')[[1]]
        mat_out <- BPCells::write_matrix_dir(mat=BPCells::convert_matrix_type(mat, type=matrix_type_d),
                                             dir=mat_dir, 
                                             compress=matrix_compress_d, 
                                             buffer_size=matrix_buffer_size_d, 
                                             overwrite=FALSE)
        push_matrix_path(mat_out)
      }
    }
  }
  else {
    stop('set_matrix_class: unrecognized matrix_info[[\'matrix_class\']]: ', matrix_info[['matrix_class']])
  }

message('set_matrix_class: matrix_info: out:')
show_matrix_info(get_matrix_info(mat_out), indent='  ')

  return(mat_out)
}


rm_bpcells_dir <- function(mat) {
  mat_info <- get_matrix_info(mat)
  if(mat_info[['matrix_class']] == 'BPCells' &&
     mat_info[['matrix_mode']] == 'dir') {
    unlink(mat_info[['matrix_path']], recursive=TRUE)
    rm(mat)
  }
}

