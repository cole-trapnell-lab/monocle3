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


# select_matrix_parameter_value() is used by set_matrix_control to select matrix_control values.
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
#     matrix_mode: 'mem'  default: 'dir' NOTE: disallow matrix_mode 'mem'
#       matrix_type: 'float', 'double'  default: 'double'
#       matrix_compress: TRUE, FALSE default: FALSE
#       matrix_bpcells_copy: TRUE, FALSE default: TRUE
#     matrix_mode: 'dir'
#       matrix_type: 'float', 'double' default: 'double'
#       matrix_path: <path to directory or file> default: 'NULL' -> temporary directory in pwd
#       matrix_compress: TRUE, FALSE default: FALSE
#       matrix_buffer_size: <integer> default: 8192L
#       matrix_bpcells_copy: TRUE, FALSE default: TRUE
check_matrix_control <- function(matrix_control=list(), control_type=c('unrestricted', 'pca'), check_conditional=FALSE) {
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
                                  'matrix_bpcells_copy',
                                  'show_values')

  allowed_matrix_class <- c('dgCMatrix', 'BPCells')

  allowed_matrix_mode <- list()
  allowed_matrix_mode[['unrestricted']] <- c('dir')
  allowed_matrix_mode[['pca']] <- c('dir')

  allowed_matrix_type <- list()
  allowed_matrix_type[['unrestricted']] <- c('float', 'double')
  allowed_matrix_type[['pca']] <- c('float', 'double')

  error_string <- ''

  if(!all(names(matrix_control) %in% allowed_control_parameters)) {
    error_string <- paste0(error_string, '\n', 'invalid control parameter')
  }

  if(check_conditional == FALSE) {
    if(!(is.null(matrix_control[['matrix_class']])) &&
       !(matrix_control[['matrix_class']] %in% allowed_matrix_class)) {
      error_string <- paste0(error_string, '\n', 'invalid matrix_class "', matrix_control[['matrix_class']], '"')
    }

    if(matrix_control[['matrix_class']] == 'BPCells') {
      if(control_type == 'unrestricted')
        allowed_values <- allowed_matrix_mode[['unrestricted']]
      else
      if(control_type == 'pca')
        allowed_values <- allowed_matrix_mode[['pca']]
      else
        stop('check_matrix_control: unknown control type \'', control_type, '\'')
      if(!(is.null(matrix_control[['matrix_mode']])) &&
         !(matrix_control[['matrix_mode']] %in% allowed_values)) {
        error_string <- paste0(error_string, '\n', 'invalid matrix_mode "', matrix_control[['matrix_mode']], '"')
      }
    
      if(control_type == 'unrestricted')
        allowed_values <- allowed_matrix_type[['unrestricted']]
      else
      if(control_type == 'pca')
        allowed_values <- allowed_matrix_type[['pca']]
      else
        stop('check_matrix_control: unknown control type \'', control_type, '\'')
      if(!(is.null(matrix_control[['matrix_type']])) &&
         !(matrix_control[['matrix_type']] %in% allowed_values)) {
        error_string <- paste0(error_string, '\n', 'invalid matrix_type "', matrix_control[['matrix_type']], '"')
      }

      if(!(is.null(matrix_control[['matrix_compress']])) &&
         !(is.logical(matrix_control[['matrix_compress']]))) {
        error_string <- paste0(error_string, '\n', 'matrix_compress value must be a logical type')
      }

      if(!(is.null(matrix_control[['matrix_path']])) &&
         !(is.character(matrix_control[['matrix_path']]))) {
        error_string <- paste0(error_string, '\n', 'matrix_path value must be a character type')
      }

      if(!(is.null(matrix_control[['matrix_buffer_size']])) &&
         !(is.integer(matrix_control[['matrix_buffer_size']]))) {
        error_string <- paste0(error_string, '\n', 'matrix_buffer_size value must be an integer type')
      }

      if(!(is.null(matrix_control[['matrix_bpcells_copy']])) &&
         !(is.logical(matrix_control[['matrix_bpcells_copy']]))) {
        error_string <- paste0(error_string, '\n', 'matrix_bpcells_copy value must be a logical type')
      }
    }
  }
  else {
    # Check matrix_class value.
    if(is.null(matrix_control[['matrix_class']])) {
      error_string <- paste0(error_string, '\n', 'matrix_class not set')
    }
    else
    if(!(matrix_control[['matrix_class']] %in% allowed_matrix_class)) {
      error_string <- paste0(error_string, '\n', 'invalid matrix_class "', matrix_control[['matrix_class']], '"')
    }

    if(matrix_control[['matrix_class']] == 'BPCells') {
      # Check matrix_type value.
      if(control_type == 'unrestricted')
        allowed_values <- allowed_matrix_type[['unrestricted']]
      else
      if(control_type == 'pca')
        allowed_values <- allowed_matrix_type[['pca']]
      else
        stop('check_matrix_control: unknown control type \'', control_type, '\'')
      if(!(matrix_control[['matrix_type']] %in% allowed_values)) {
        error_string <- paste0(error_string, '\n', 'bad  matrix_type "', matrix_control[['matrix_type']], '"')
      }
  
      # Check matrix_compress value.
      if(!is.logical(matrix_control[['matrix_compress']])) {
        error_string <- paste0(error_string, '\n', 'matrix_compress must be as logical type')
      }
  
      # Check matrix_mode value.
      if(control_type == 'unrestricted')
        allowed_values <- allowed_matrix_mode[['unrestricted']]
      else
      if(control_type == 'pca')
        allowed_values <- allowed_matrix_mode[['pca']]
      else
        stop('check_matrix_control: unknown control type \'', control_type, '\'')
      if(!(matrix_control[['matrix_mode']] %in% allowed_values)) {
        error_string <- paste0(error_string, '\n', 'invalid matrix_mode "', matrix_control[['matrix_mode']], '"')
      }
  
      if(matrix_control[['matrix_mode']] == 'dir') {
        # Check matrix_path value.
        if(!(is.character(matrix_control[['matrix_path']]))) {
          error_string <- paste0(error_string, '\n', 'bad matrix_path "', matrix_control[['matrix_path']], '"')
        }
  
        # Check matrix_buffer_size.
        if(!(is.integer(matrix_control[['matrix_buffer_size']]))) {
          error_string <- paste0(error_string, '\n', 'matrix_buffer_size must be an integer')
        }
      }
      if(!is.logical(matrix_control[['matrix_bpcells_copy']])) {
        error_string <- paste0(error_string, '\n', 'matrix_bpcells_copy value must be a logical type')
      }
    }
  }

  if(error_string != '') {
    stop(paste0(stringr::str_trim(error_string), '\n'))
  }

  return(TRUE)
}


# Usage
#   matrix_class: default: 'dgCMatrix'
#   matrix_class: 'BPCells'
#     matrix_mode: 'mem'  default: 'dir'
#       matrix_type: 'float', 'double' default: 'double'
#       matrix_compress: TRUE, FALSE default: FALSE
#       matrix_bpcells_copy: TRUE, FALSE default: TRUE
#     matrix_mode: 'dir'
#       matrix_type: 'float', 'double' default: 'double'
#       matrix_path: <path to directory or file> default: 'NULL' -> temporary directory in pwd
#       matrix_compress: TRUE, FALSE default: FALSE
#       matrix_buffer_size: <integer> default: 8192L
#       matrix_bpcells_copy: TRUE, FALSE default: TRUE
# Notes:
#   o  modification to any of set_assay_control, set_pca_control,
#      or set_pca_control_default may necessitate modifications
#      to all of them.

#' Verify and set the matrix_control parameter list.
#'
#' @description Verifies and sets the list of parameter values
#'   that is used to make the count matrix that is stored in the
#'   cell_data_set or certain other matrices that are used
#'   during the Monocle3 run. To see the
#'   default values,
#'   call "set_matrix_control(matrix_control=list(matrix_class='BPCells', show_values=TRUE))".
#'   "show_values=TRUE" can be used in functions that have the
#'   matrix_control list parameter, in which case the function will
#'   show the matrix_control values to be used and then stop.
#' @param matrix_control Input control list.
#' @param matrix_control_default Input default control list.
#' @param control_type A string of either "unrestricted" or "pca". A control_type of "pca"
#'   restricts certain list parameters.
#' @return matrix_control Output control list.
#'
#' @section matrix_control list values:
#' \describe{
#'   \item{matrix_class}{A string that specifies the matrix
#'      class to use for matrix storage. The acceptable
#'      values are "dgCMatrix" and "BPCells".}
#'   \item{matrix_type}{A string that specifies whether to
#'      store the matrix values as single precision "floats"
#'      (matrix_type="float") or double precision "doubles"
#'      (matrix_type="double"). The default is "double".
#'      "matrix_type" is used only for BPCells class matrices.}
#'   \item{matrix_mode}{A string that specifies whether to
#'      store the BPCells class matrix in memory
#'      (matrix_mode="mem") or on disk (matrix_mode="dir").
#'      "matrix_mode" is used only for BPCells class
#'      matrices. At this time, only "dir" is allowed.}
#'   \item{matrix_path}{A string that specifies the directory
#'      where the BPCells on-disk matrix data are stored in a
#'      sub-directory with a randomized name. The default is
#'      in the directory where R is running. "matrix_path" is
#'      used only for BPCells class matrices with
#'      matrix_mode="dir". For example, if matrix_path is set
#'      to "my_dir" Monocle3 will create a directory with a
#'      name that has the form "monocle.bpcells.*.tmp" in
#'      "my_dir". The asterisk represents a random string that
#'      makes the name unique.}
#'   \item{matrix_compress}{A logical that specifies whether
#'      to use bit-packing compression to store BPCells matrix
#'      values. Only the matrix indices are compressed for
#'      matrix_types "float" and "double". "matrix_compress"
#'      is used only for BPCells class matrices. The default is
#'      FALSE, which improves processing speed.}
#'   \item{matrix_buffer_size}{An integer that specifies how
#'      many items of data to buffer in memory before flushing
#'      to disk. This is used for matrix_class="BPCells" with
#'      matrix_mode="dir". The default is 8192L.}
#'   \item{matrix_bpcells_copy}{A logical that specifies
#'      whether the input BPCells matrix is to be copied. This
#'      is relevant only when the input matrix and the desired
#'      output matrix are the same; that is, have the same
#'      matrix_mode, matrix_path, matrix_compress, and
#'      matrix_buffer_size values. If matrix_bpcells_copy is
#'      TRUE, the queued operations are applied to a new
#'      on-disk copy of of the input matrix and the operation
#'      queue is emptied. If FALSE, the queued operations are
#'      not applied and the on-disk storage is unaltered.}
#' }
#' @export
set_matrix_control <- function(matrix_control=list(), matrix_control_default=list(), control_type=c('unrestricted', 'pca')) {
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
  default_matrix_bpcells_copy <- NULL

  matrix_control_out = list()

  matrix_control_out[['matrix_class']] <- select_matrix_parameter_value(parameter='matrix_class', matrix_control=matrix_control, matrix_control_default=matrix_control_default, default_value=default_matrix_class)

  if(matrix_control_out[['matrix_class']] == 'BPCells') {
     matrix_control_out[['matrix_mode']] <- select_matrix_parameter_value(parameter='matrix_mode', matrix_control=matrix_control, matrix_control_default=matrix_control_default, default_value=default_matrix_mode)

    if(matrix_control_out[['matrix_mode']] == 'mem') {
       matrix_control_out[['matrix_type']] <- select_matrix_parameter_value(parameter='matrix_type', matrix_control=matrix_control, matrix_control_default=matrix_control_default, default_value=default_matrix_type)
       matrix_control_out[['matrix_compress']] <- select_matrix_parameter_value(parameter='matrix_compress', matrix_control=matrix_control, matrix_control_default=matrix_control_default, default_value=default_matrix_compress)
       matrix_control_out[['matrix_bpcells_copy']] <- select_matrix_parameter_value(parameter='matrix_bpcells_copy', matrix_control=matrix_control, matrix_control_default=matrix_control_default, default_value=default_matrix_bpcells_copy)
    }
    else
    if(matrix_control_out[['matrix_mode']] == 'dir') {
       matrix_control_out[['matrix_type']] <- select_matrix_parameter_value(parameter='matrix_type', matrix_control=matrix_control, matrix_control_default=matrix_control_default, default_value=default_matrix_type)
       matrix_control_out[['matrix_path']] <- select_matrix_parameter_value(parameter='matrix_path', matrix_control=matrix_control, matrix_control_default=matrix_control_default, default_value=default_matrix_path)
       matrix_control_out[['matrix_compress']] <- select_matrix_parameter_value(parameter='matrix_compress', matrix_control=matrix_control, matrix_control_default=matrix_control_default, default_value=default_matrix_compress)
       matrix_control_out[['matrix_buffer_size']] <- select_matrix_parameter_value(parameter='matrix_buffer_size', matrix_control=matrix_control, matrix_control_default=matrix_control_default, default_value=default_matrix_buffer_size)
       matrix_control_out[['matrix_bpcells_copy']] <- select_matrix_parameter_value(parameter='matrix_bpcells_copy', matrix_control=matrix_control, matrix_control_default=matrix_control_default, default_value=default_matrix_bpcells_copy)
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
    show_matrix_control(matrix_control=matrix_control_out, label='  matrix_control: ')
    stop_no_noise()
  }

  return(matrix_control_out)
}


# Report matrix_control list values.
show_matrix_control <- function(matrix_control=list(), label=NULL) {
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
      message(indent, '  matrix_bpcells_copy: ', ifelse(!is.null(matrix_control[['matrix_bpcells_copy']]), matrix_control[['matrix_bpcells_copy']], as.character(NA)))
    }
    else
    if(matrix_control[['matrix_mode']] == 'dir') {
      message(indent, '  matrix_mode: ', ifelse(!is.null(matrix_control[['matrix_mode']]), matrix_control[['matrix_mode']], as.character(NA)))
      message(indent, '  matrix_type: ', ifelse(!is.null(matrix_control[['matrix_type']]), matrix_control[['matrix_type']], as.character(NA)))
      message(indent, '  matrix_path: ', ifelse(!is.null(matrix_control[['matrix_path']]), matrix_control[['matrix_path']], as.character(NA)))
      message(indent, '  matrix_compress: ', ifelse(!is.null(matrix_control[['matrix_compress']]), matrix_control[['matrix_compress']], as.character(NA)))
      message(indent, '  matrix_buffer_size: ', ifelse(!is.null(matrix_control[['matrix_buffer_size']]), matrix_control[['matrix_buffer_size']], as.character(NA)))
      message(indent, '  matrix_bpcells_copy: ', ifelse(!is.null(matrix_control[['matrix_bpcells_copy']]), matrix_control[['matrix_bpcells_copy']], as.character(NA)))
    }
    else {
      stop('show_matrix_control: unsupported matrix class/mode/...\'', matrix_control[['method']], '\'')
    }
  }
  else {
    stop('show_matrix_control: unsupported matrix class/mode/...\'', matrix_control[['method']], '\'')
  }
}


# Push the directory path of the BPCells on-disk matrix onto the path stack. This is used
# to remove the directories when R is exited.
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
    return(bpcells_find_base_matrix(mat=BPCells:::matrix_inputs(mat)[[1]]))
  }
  return(mat)
}


#
# Known matrix types:
#   o  dense matrix: matrix
#   o  sparse matrix: CsparseMatrix, dgCMatrix, ...
#   o  BPCells matrix:
#        o  matrix_mode: 'dir'
#        o  matrix_type: float, double
#        o  matrix_compress: TRUE, FALSE
#        o  matrix_path
#        o  matrix_buffer_size
#
get_matrix_class <- function(mat) {

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

  return(matrix_info)
}


# Get/infer the matrix information.
get_matrix_info <- function(mat) {
  matrix_info <- get_matrix_class(mat=mat)

  if(is.null(matrix_info[['matrix_class']])) {
    message('bad matrix info -- dropping into browser')
    browser()
  }

  if(matrix_info[['matrix_class']] != 'BPCells') {
    return(matrix_info)
  }

  bmat <- bpcells_find_base_matrix(mat=mat)

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
                          'MatrixDir', 'Iterable_dgCMatrix_wrapper'))) {
    stop('get_matrix_info: unrecognized BPCells matrix class \"', class(mat), '\"')
    return(NULL)
  }

  matrix_info[['matrix_class']] <- 'BPCells'

  if(class(bmat) == 'Iterable_dgCMatrix_wrapper') {
    matrix_info[['matrix_mode']] <- 'dgCMatrix'
  }
  else
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


# Show the matrix information.
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
  if(!is.null(matrix_info[['matrix_path']])) { 
    message(paste0(indent, 'path:        ', matrix_info[['matrix_path']]))
  }
  if(!is.null(matrix_info[['matrix_compress']])) { 
    message(paste0(indent, 'compress:    ', matrix_info[['matrix_compress']]))
  }
  if(!is.null(matrix_info[['matrix_buffer_size']])) { 
    message(paste0(indent, 'buffer_size: ', matrix_info[['matrix_buffer_size']]))
  }
}


# Compare two matrix control lists. Return TRUE if they are the same. There is
# a minor complication. If the matrix_class in both lists is BPCells, the
# matrix storage mode is 'dir, and compare_matrix_path_flag is FALSE, the
# matrix_path values may differ and compare_matrix_control still returns TRUE.
# Otherwise, compare_matrix_control returns FALSE.
compare_matrix_control <- function(matrix_control=list(), matrix_info=list(), compare_matrix_path_flag=FALSE) {
  if(matrix_control[['matrix_class']] == 'dgCMatrix' && matrix_info[['matrix_class']] == 'dgCMatrix') {
    return(TRUE)
  }
  else
  if(matrix_control[['matrix_class']] == 'BPCells' && matrix_info[['matrix_class']] == 'BPCells') {
    if(matrix_control[['matrix_mode']] == 'mem' && matrix_info[['matrix_mode']] == 'mem') {
      if(matrix_control[['matrix_type']] == matrix_info[['matrix_type']] &&
         matrix_control[['matrix_compress']] == matrix_info[['matrix_compress']]) {
        return(TRUE)
      }
    }
    else
    if(matrix_control[['matrix_mode']] == 'dir' && matrix_info[['matrix_mode']] == 'dir') {
      if(matrix_control[['matrix_type']] == matrix_info[['matrix_type']] &&
         matrix_control[['matrix_compress']] == matrix_info[['matrix_compress']] &&
         matrix_control[['matrix_buffer_size']] == matrix_info[['matrix_buffer_size']]) {
        if(compare_matrix_path_flag == FALSE ||
            matrix_control[['matrix_path']] == matrix_info[['matrix_path']]) {
          return(TRUE)
        }
        else {
          return(FALSE)
        }
      }
    }
  }

  return(FALSE)
}


# set_matrix_class
#  Notes:
#    o  set_matrix_class is meant to be the function used
#       nearly exclusively for making BPCells matrices. The
#       exceptions are load_monocle_objects and save_monocle_objects
#       which use BPCells::write_matrix_dir() 'directly'.
#    o  Cast an input matrix into a class given or inferred from the
#       matrix_control, which must be complete in the sense that all
#       values required for the class are given explicitly.
#    o  When the matrix_info values of the input matrix are the same
#       as the matrix_control values AND the matrix_class is BPCells,
#       the queued operations are not applied (without copying the
#       on-disk matrix) when matrix_control[['matrix_bpcells_copy']] is
#       FALSE. Otherwise, the queued BPCells operations are applied
#       and the result is stored in a new on-disk directory. Applying
#       the queued BPCells operations improves performance.
#   
set_matrix_class <- function(mat, matrix_control=list()) {
  # Check matrix_control list.
  check_matrix_control(matrix_control=matrix_control, control_type='unrestricted', check_conditional=TRUE)

  # Get input matrix info.
  matrix_info <- get_matrix_info(mat=mat)

# message('set_matrix_class: matrix_info: in:')
# show_matrix_info(matrix_info, indent='  ')
# message('')
# message('set_matrix_class: matrix_control in:')
# show_matrix_control(matrix_control)
# message('')

  if(matrix_info[['matrix_class']] == 'dgTMatrix') {
    mat <- as(mat, 'CsparseMatrix')
    matrix_info[['matrix_class']] <- 'dgCMatrix'
  }

  mat_out <- mat

  # Can we return the input matrix object immediately? We
  # can when the matrix_info and matrix_control are the same
  # (except for the matrix_path) AND the matrix_class is
  # dgCMatrix OR (the matrix_class is 'BPCells' AND
  # matrix_bpcells_copy is FALSE.
  if(compare_matrix_control(matrix_control, matrix_info, compare_matrix_path_flag=FALSE)) {
    if(matrix_control[['matrix_class']] == 'dgCMatrix' ||
       (matrix_control[['matrix_class']] == 'BPCells' &&
        matrix_control[['matrix_bpcells_copy']] == FALSE)) {
      return(mat_out)
    }
  }

  if(matrix_control[['matrix_class']] == 'dgCMatrix') {
    if(matrix_info[['matrix_class']] != 'dgCMatrix') {
      if(matrix_info[['matrix_class']] == 'BPCells') {
        mat_out <- as(mat, 'dgCMatrix')
      }
      else {
        mat_out <- as(mat, 'CsparseMatrix')
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
        push_matrix_path(mat=mat_out)
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
        push_matrix_path(mat=mat_out)
      }
    }
  }
  else {
    stop('set_matrix_class: unrecognized matrix_info[[\'matrix_class\']]: ', matrix_info[['matrix_class']])
  }

#message('set_matrix_class: matrix_info: out:')
#show_matrix_info(get_matrix_info(mat=mat_out), indent='  ')

  return(mat_out)
}


# Remove the BPCells matrix directory and the matrix R object.
rm_bpcells_dir <- function(mat) {
  mat_info <- get_matrix_info(mat=mat)
  if(mat_info[['matrix_class']] == 'BPCells' &&
     mat_info[['matrix_mode']] == 'dir') {
    unlink(mat_info[['matrix_path']], recursive=TRUE)
    rm(mat)
  }
}


#
# Set row-major order counts matrix in cds with BPCells counts matrix.
#

#' Set the row-major order counts matrix in the assays slot of the CDS when the
#'   CDS has a BPCells counts matrix.
#'
#' @description
#' By default, BPCells stores matrices as one-dimensional vectors in column-major order, as does R.
#' As a result, column access is fast and row access is slow. We use BPCell's ability to also store and
#' access matrices in row-major order, which gives fast row access. The function set_cds_row_order_matrix
#' creates the row-major order matrix and stores it in the assays slot of the CDS with the name
#' "counts_row_order". The two copies of the counts matrix must have the same count values so if you
#' replace or change the CDS's counts matrix, you must also update the `counts_row_order` matrix, which
#' you can do using this function set_cds_row_order_matrix.
#' @param cds cell_data_set A cell_data_set.
#' @return cell_data_set The cell_data_set with the additional row-major order
#'    counts matrix.
#' @examples
#'    cds <- load_a549(matrix_control=list(matrix_class='BPCells'))
#'    cds <- set_cds_row_order_matrix(cds)
#'    str(cds)
#'
#' @export
set_cds_row_order_matrix <- function(cds) {

  mat_c <- counts(cds)
  if(!is(mat_c, 'IterableMatrix')) {
    return(cds)
  }

  matrix_info <- get_matrix_info(mat=mat_c)
  if(matrix_info[['matrix_mode']] == 'dir') {
    bmat <- bpcells_find_base_matrix(mat=mat_c)
    matrix_path <- dirname(bmat@dir)
  }
  else {
    matrix_path <- '.'
  }

  # Remove existing counts_row_order matrix and directory.
  if(!is.null(assays(cds)[['counts_row_order']])) {
    rm_bpcells_dir(mat=assays(cds)[['counts_row_order']])
    assays(cds)[['counts_row_order']] <- NULL
  }

  # Make a BPCells count matrix in row major order.
  outdir <- tempfile(pattern=paste0('monocle.bpcells.',
                                    format(Sys.Date(), format='%Y%m%d'), '.'),
                     tmpdir=matrix_path,
                     fileext='_r.tmp')[[1]]
  tmpdir <- tempfile('monocle.transpose_bpc.', '.', '.tmp')

  # I see no option for choosing compressed matrix and transpose_storage_order appears to
  # compress. This is not a big deal because only the indices are compressed.
  mat_r <- BPCells::transpose_storage_order(matrix=mat_c, outdir=outdir, tmpdir=tmpdir, load_bytes=4194304L, sort_bytes=1073741824L)
  unlink(tmpdir, recursive=TRUE)
  push_matrix_path(mat=mat_r)

  assay(cds, 'counts_row_order') <- mat_r

  return(cds)
}


#' Convert the counts matrix class in the given CDS.
#'
#' @description Converts the counts matrix that is in the
#'   given CDS to the matrix_class specified in the
#'   matrix_control list.
#' @param cds cell_data_set The cell_data_set that has the
#'   counts matrix to be converted.
#' @param matrix_control list A list of matrix control
#'   values used to convert the counts matrix. If the
#'   counts matrix in the cds is the same as the desired
#'   counts matrix, it is not altered.
#' @return cell_data_set The cell_data_set with the converted
#'   counts matrix.
#' @examples
#'    cds <- load_a549()
#'    str(counts(cds))
#'    cds <- convert_counts_matrix(cds, matrix_control=list(matrix_class='BPCells'))
#'    str(counts(cds))
#'
#' @export
convert_counts_matrix <- function(cds, matrix_control=list(matrix_class='BPCells')) {
  if(!is.null(matrix_control[['matrix_class']]) && matrix_control[['matrix_class']] == 'BPCells') {
    matrix_control_default <- get_global_variable('matrix_control_bpcells_unrestricted')
  }
  else {
    matrix_control_default <- get_global_variable('matrix_control_csparsematrix_unrestricted')
  }
  matrix_control_res <- set_matrix_control(matrix_control=matrix_control, matrix_control_default=matrix_control_default, control_type='unrestricted')

  # Do not make a BPCells matrix on-disk copy if the
  # matrix_info values of the counts matrix are the
  # same as the matrix_control values.
  matrix_control_res[['matrix_bpcells_copy']] <- FALSE

  mat <- counts(cds)

  counts(cds, bpcells_warn=FALSE) <- set_matrix_class(mat=mat, matrix_control=matrix_control_res)

  if(matrix_control_res[['matrix_class']] == 'BPCells') {
    push_matrix_path(mat=mat)
    cds <- set_cds_row_order_matrix(cds=cds)
  }

  return(cds)
}

