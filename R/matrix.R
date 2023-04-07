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
check_matrix_control <- function(matrix_control, control_source=c('assay', 'pca'), check_conditional=FALSE) {
  assertthat::assert_that(is.list(matrix_control))
  assertthat::assert_that(is.logical(check_conditional))

  control_source <- match.arg(control_source)

  allowed_control_parameters <- c('matrix_assay',
                                  'matrix_class',
                                  'matrix_mode',
                                  'matrix_type',
                                  'matrix_compress',
                                  'matrix_path',
                                  'matrix_buffer_size',
                                  'show_values')

  allowed_matrix_class <- c('dgCMatrix', 'BPCells')
  allowed_matrix_type_assay <- c('uint32_t', 'float', 'double')
  allowed_matrix_type_pca <- c('float', 'double')
  allowed_matrix_mode <- c('mem', 'dir') 

  control_name <- if(control_source == 'assay') 'assay_control' else 'pca_control'

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

    allowed_values <- if(control_source == 'assay') allowed_matrix_type_assay else allowed_matrix_type_pca
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
      allowed_values <- if(control_source == 'assay') allowed_matrix_type_assay else allowed_matrix_type_pca
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
#   matrix_assay: default: 'counts'
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
#

# Notes:
#   o  modification to any of set_assay_control, set_pca_control,
#      or set_pca_control_default may necessitate modifications
#      to all of them.

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
set_assay_control <- function(assay_control=list()) {

  check_matrix_control(matrix_control=assay_control, control_source='assay', check_conditional=FALSE)

  if(is.null(assay_control[['matrix_class']]) ||
     assay_control[['matrix_class']] == 'dgCMatrix') {
    assay_control_default <- get_global_variable('assay_control_csparsematrix')
  }
  else
  if(assay_control[['matrix_class']] == 'BPCells') {
    assay_control_default <- get_global_variable('assay_control_bpcells')
  }

  check_matrix_control(matrix_control=assay_control_default, control_source='assay', check_conditional=FALSE)

  #
  # Last resort fall-back parameter values.
  #
  default_matrix_assay <- 'counts'
  default_matrix_class <- 'dgCMatrix'
  default_matrix_mode <- 'mem'
  default_matrix_type <- 'uint32_t'
  default_matrix_compress <- TRUE
  default_matrix_path <- '.'
  default_matrix_buffer_size <- 8192L

  assay_control_out = list()

  assay_control_out[['matrix_assay']] <- select_assay_parameter_value('matrix_assay', assay_control, assay_control_default, default_matrix_assay)
  if(!(assay_control_out[['matrix_assay']] %in% c('counts'))) {
    stop(paste0(assay_control_out[['matrix_assay']], ' is not a valid matrix_assay.'))
  }

  assay_control_out[['matrix_class']] <- select_assay_parameter_value('matrix_class', assay_control, assay_control_default, default_matrix_class)

  if(assay_control_out[['matrix_class']] == 'BPCells') {
     assay_control_out[['matrix_mode']] <- select_assay_parameter_value('matrix_mode', assay_control, assay_control_default, default_matrix_mode)

    if(assay_control_out[['matrix_mode']] == 'mem') {
       assay_control_out[['matrix_type']] <- select_assay_parameter_value('matrix_type', assay_control, assay_control_default, default_matrix_type)
       assay_control_out[['matrix_compress']] <- select_assay_parameter_value('matrix_compress', assay_control, assay_control_default, default_matrix_compress)
    }
    else
    if(assay_control_out[['matrix_mode']] == 'dir') {
       assay_control_out[['matrix_type']] <- select_assay_parameter_value('matrix_type', assay_control, assay_control_default, default_matrix_type)
       assay_control_out[['matrix_path']] <- select_assay_parameter_value('matrix_path', assay_control, assay_control_default, default_matrix_path)
       assay_control_out[['matrix_compress']] <- select_assay_parameter_value('matrix_compress', assay_control, assay_control_default, default_matrix_compress)
       assay_control_out[['matrix_buffer_size']] <- select_assay_parameter_value('matrix_buffer_size', assay_control, assay_control_default, default_matrix_buffer_size)
    }
  }

  check_matrix_control(matrix_control=assay_control_out, control_source='assay', check_conditional=TRUE)

  #
  # Set BPCells out-of-core file/directory name.
  #
  if(assay_control_out[['matrix_class']] == 'BPCells') {
    if(assay_control_out[['matrix_mode']] == 'dir') {
      tmp_dir <- if(is.null(assay_control_out[['matrix_path']])) '.' else assay_control_out[['matrix_path']]
      assay_control_out[['matrix_path']] <- tempfile(pattern=paste0('monocle.bpcells.', format(Sys.Date(), format='%Y%m%d'), '.'), tmpdir=tmp_dir, fileext='.tmp')[[1]]
    }
  }

  #
  # Display assay_control list if assay_control[['show_values']] <- TRUE.
  #
  if(!is.null(assay_control[['show_values']]) && assay_control[['show_values']] == TRUE)
  {
    report_assay_control(assay_control=assay_control_out, '  assay_control: ')
    stop_no_noise()
  }

  return(assay_control_out)
}


# Report assay_control list values.
report_assay_control <- function(assay_control, label=NULL) {
  indent <- ''
  if(!is.null(label)) {
    indent <- '  '
  }

  message(ifelse(!is.null(label), label, ''))

  message(indent, '  matrix_class: ', ifelse(!is.null(assay_control[['matrix_class']]), assay_control[['matrix_class']], as.character(NA)))

  if(assay_control[['matrix_class']] == 'dgCMatrix') {
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

