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


test_bpcells_matrix_dir_assays <- function(cds) {
  if(is(counts(cds), 'IterableMatrix'))
    return(TRUE)
  return(FALSE)
}


select_assay_parameter_value <- function(parameter, assay_control, assay_control_default, default_value) {
  if(!is.null(assay_control[[parameter]])) {
    return(assay_control[[parameter]])
  }
  else
  if(!is.null(assay_control_default[[parameter]])) {
    return(assay_control_default[[parameter]])
  }
  return(default_value)
}


# Usage
#   matrix_assay: default: 'counts'
#   matrix_class: default: 'CsparseMatrix'
#   matrix_class: 'BPCells'
#     matrix_mode: 'mem'  default: 'dir'
#       matrix_type: 'uint32_t', 'float', 'double'
#       matrix_compress: TRUE, FALSE default: FALSE
#     matrix_mode: 'dir'
#       matrix_type: 'uint32_t', 'float', 'double' default: 'uint32_t'
#       matrix_path: <path to directory or file> default: 'NULL' -> temporary directory in pwd
#       matrix_compress: TRUE, FALSE default: FALSE
#       matrix_buffer_size: <integer> default: 8192L
#

#' Verify and set the assay_control parameter list.
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
#'      matrix storage. The acceptable values are "CsparseMatrix"
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
#'   \item{matrix_path}{Specifies the disk directory where the
#'      BPCells on-disk matrix data are to be stored. The
#'      default is a directory, with a randomized name, in
#'      the directory where R is running. "matrix_path" is
#'      used only for BPCells class matrices with
#'      matrix_mode="dir".}
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
set_assay_control <- function(assay_control=list()) {

  assertthat::assert_that(methods::is(assay_control, "list"))

  allowed_control_parameters <- c('matrix_assay',
                                  'matrix_class',
                                  'matrix_mode',
                                  'matrix_type',
                                  'matrix_compress',
                                  'matrix_path',
                                  'matrix_buffer_size',
                                  'show_values')

  allowed_matrix_class <- c('CsparseMatrix', 'BPCells')

  assertthat::assert_that(methods::is(assay_control, "list"))

  if(!is.null(assay_control[['matrix_assay']]) && assay_control[['matrix_assay']] != 'counts') {
    stop('  matrix_assay value must be "counts"')
  }

  if(is.null(assay_control[['matrix_class']])) {
    assay_control[['matrix_class']] <- 'CsparseMatrix'
  }
  else
  if(!(assay_control[['matrix_class']] %in% allowed_matrix_class)) {
    stop('  matrix_class value must be "CsparseMatrix" or "BPCells"')
  }

  if(assay_control[['matrix_class']] == 'CsparseMatrix') {
    assay_control_default <- get_global_variable('assay_control_csparsematrix')
  }
  else
  if(assay_control[['matrix_class']] == 'BPCells') {
    assay_control_default <- get_global_variable('assay_control_bpcells')
  }

  assertthat::assert_that(all(names(assay_control) %in% allowed_control_parameters),
                          msg = "set_assay_control: unknown variable in assay_control")

  assertthat::assert_that(all(names(assay_control_default) %in% allowed_control_parameters),
                          msg = "set_assay_control: unknown variable in assay_control_default")

  #
  # Last resort fall-back parameter values.
  #
  default_matrix_assay <- 'counts'
  default_matrix_class <- 'CsparseMatrix'
  default_matrix_mode <- 'mem'
  default_matrix_type <- 'uint32_t'
  default_matrix_compress <- TRUE
  default_matrix_path <- NULL
  default_matrix_buffer_size <- 8192L


  error_string = list()

  assay_control_out = list()

  assay_control_out[['matrix_assay']] <- select_assay_parameter_value('matrix_assay', assay_control, assay_control_default, default_matrix_assay)
  if(!(assay_control_out[['matrix_assay']] %in% c('counts'))) {
    error_string <- list(error_string, paste0('  ',
                         assay_control_out[['matrix_assay']],
                         ' is not a valid matrix_assay.'))
    stop(error_string)
  }

  assay_control_out[['matrix_class']] <- select_assay_parameter_value('matrix_class', assay_control, assay_control_default, default_matrix_class)

  if(assay_control_out[['matrix_class']] == 'CsparseMatrix') {
    assay_control_out[['matrix_mode']] <- select_assay_parameter_value('matrix_mode', assay_control, assay_control_default, default_matrix_mode)
    if(!(assay_control_out[['matrix_mode']] %in% c('mem'))) {
      error_string <- list(error_string, paste0('  ',
                           assay_control_out[['matrix_mode']],
                           ' is not valid for matrix_class CsparseMatrix.'))
      stop(error_string)
    }
  }
  else
  if(assay_control_out[['matrix_class']] == 'BPCells') {
    assay_control_out[['matrix_mode']] <- select_assay_parameter_value('matrix_mode', assay_control, assay_control_default, default_matrix_mode)
    if(!(assay_control_out[['matrix_mode']] %in% c('mem', 'dir'))) {
      error_string <- list(error_string, paste0('  ',
                           assay_control_out[['matrix_mode']],
                           ' is not a valid matrix_mode for matrix_class "BPCells".'))
      stop(error_string)
    }

    if(assay_control_out[['matrix_mode']] == 'mem') {
      assay_control_out[['matrix_type']] <- select_assay_parameter_value('matrix_type', assay_control, assay_control_default, default_matrix_type)
      if(!(assay_control_out[['matrix_type']] %in% c('uint32_t', 'float', 'double'))) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out[['matrix_type']],
                             ' is not a valid matrix_type for matrix_class "BPCells".'))
      }

      assay_control_out[['matrix_compress']] <- select_assay_parameter_value('matrix_compress', assay_control, assay_control_default, default_matrix_compress)
      if(!is.logical(assay_control_out[['matrix_compress']])) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out[['matrix_compress']],
                             ' is not TRUE or FALSE.'))
      }
      if(length(error_string) > 0) {
        stop(error_string)
      }
    }
    else
    if(assay_control_out[['matrix_mode']] == 'dir') {
      assay_control_out[['matrix_type']] <- select_assay_parameter_value('matrix_type', assay_control, assay_control_default, default_matrix_type)
      if(!(assay_control_out[['matrix_type']] %in% c('uint32_t', 'float', 'double'))) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out[['matrix_type']],
                             ' is not a valid matrix_type for matrix_class "BPCells".'))
      }

      assay_control_out[['matrix_path']] <- select_assay_parameter_value('matrix_path', assay_control, assay_control_default, default_matrix_path)
      if(!is.null(assay_control_out[['matrix_path']]) &&
         !is.character(assay_control_out[['matrix_path']])) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out[['matrix_path']],
                             ' is not a valid matrix_path.'))

      }

      assay_control_out[['matrix_compress']] <- select_assay_parameter_value('matrix_compress', assay_control, assay_control_default, default_matrix_compress)
      if(!is.logical(assay_control_out[['matrix_compress']])) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out[['matrix_compress']],
                             ' is not TRUE or FALSE.'))
      }

      assay_control_out[['matrix_buffer_size']] <- select_assay_parameter_value('matrix_buffer_size', assay_control, assay_control_default, default_matrix_buffer_size)
      if(!is.integer(assay_control_out[['matrix_buffer_size']])) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out[['matrix_buffer_size']],
                             ' is not an integer.'))
      }

      if(length(error_string) > 0) {
        stop(error_string)
      }
    }
  }

  #
  # Set BPCells out-of-core file/directory name, if necessary.
  #
  if(assay_control_out[['matrix_class']] == 'BPCells' && is.null(assay_control_out[['matrix_path']])) {
    if(assay_control_out[['matrix_mode']] == 'dir') {
      assay_control_out[['matrix_path']] <- tempfile(pattern='monocle.bpcells.', tmpdir='.', fileext='.tmp')[[1]]
    }
  }

  #
  # Display assay_control list if assay_control[['show_values']] <- TRUE.
  #
  if(!is.null(assay_control[['show_values']]) && assay_control[['show_values']] == TRUE)
  {
    report_assay_control('  assay_control: ', assay_control=assay_control_out)
    stop_no_noise()
  }

  return(assay_control_out)
}


# Report assay_control list values.
report_assay_control <- function(label=NULL, assay_control) {
  indent <- ''
  if(!is.null(label)) {
    indent <- '  '
  }

  message(ifelse(!is.null(label), label, ''))

  message(indent, '  matrix_class: ', ifelse(!is.null(assay_control[['matrix_class']]), assay_control[['matrix_class']], as.character(NA)))

  if(assay_control[['matrix_class']] == 'CsparseMatrix') {
    message(indent, '  matrix_mode: ', ifelse(!is.null(assay_control[['matrix_mode']]), assay_control[['matrix_mode']], as.character(NA)))
  }
  else
  if(assay_control[['matrix_class']] == 'BPCells') {
    if(assay_control[['matrix_mode']] == 'mem') {
      message(indent, '  matrix_type: ', ifelse(!is.null(assay_control[['matrix_type']]), assay_control[['matrix_type']], as.character(NA)))
      message(indent, '  matrix_compress: ', ifelse(!is.null(assay_control[['matrix_compress']]), assay_control[['matrix_compress']], as.character(NA)))
    }
    else
    if(assay_control[['matrix_mode']] == 'dir') {
      message(indent, '  matrix_type: ', ifelse(!is.null(assay_control[['matrix_type']]), assay_control[['matrix_type']], as.character(NA)))
      message(indent, '  matrix_path: ', ifelse(!is.null(assay_control[['matrix_path']]), assay_control[['matrix_path']], as.character(NA)))
      message(indent, '  matrix_compress: ', ifelse(!is.null(assay_control[['matrix_compress']]), assay_control[['matrix_compress']], as.character(NA)))
      message(indent, '  matrix_buffer_size: ', ifelse(!is.null(assay_control[['matrix_buffer_size']]), assay_control[['matrix_buffer_size']], as.character(NA)))
    }
    else {
      stop('report_assay_control: unsupported assay class/mode/...\'', assay_control[['method']], '\'')
    }
  }
  else {
    stop('report_assay_control: unsupported assay class/mode/...\'', assay_control[['method']], '\'')
  }
}


push_matrix_path <- function(mat, matrix_type) {
  x <- get_global_variable('monocle_gc_matrix_path')
  x <- append(x, mat)
  set_global_variable('monocle_gc_matrix_path', x)
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


bpcells_base_matrix_info <- function(mat, indent='') {
  if(!is(mat, 'IterableMatrix'))
    return(0)
  bmat <- bpcells_find_base_matrix(mat)
  slot_names <- slotNames(bmat)
  message(paste0(indent, 'class: ', indent, class(bmat)))
  if('type' %in% slot_names)
    message(paste0(indent, 'type: ', bmat@type))
  if('dir' %in% slot_names)
    message(paste0(indent, 'dir: ', bmat@dir))
  if('compressed' %in% slot_names)
    message(paste0(indent, 'compressed: ', bmat@compressed))
  if('buffer_size' %in% slot_names)
    message(paste0(indent, 'buffer_size: ', bmat@buffer_size))
  if('version' %in% slot_names)
    message(paste0(indent, 'version: ', bmat@version))
}

