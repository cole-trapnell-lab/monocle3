################################################################################
# BPCells
################################################################################


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


#' Check and set the assay_control list values.
#' @param assay_control Input control list.
#' @return assay_control Output control list.
# Usage
#   matrix_assay: default: 'counts'
#   matrix_class: default: 'CsparseMatrix'
#     matrix_mode: default: 'mem'
#   matrix_class: 'BPCells'
#     matrix_mode: 'mem'  default: 'dir'
#       matrix_type: 'uint32_t', 'float', 'double'
#       matrix_compress: TRUE, FALSE default: FALSE
#     matrix_mode: 'dir'
#       matrix_type: 'uint32_t', 'float', 'double' default: 'uint32_t'
#       matrix_path: <path to directory or file> default: 'NULL' -> temporary directory in pwd
#       matrix_compress: TRUE, FALSE default: FALSE
#       matrix_buffer_size: <integer> default: 8192L
#       matrix_overwrite: TRUE, FALSE default: FALSE
set_assay_control <- function(assay_control=list()) {

  assertthat::assert_that(methods::is(assay_control, "list"))

  allowed_control_parameters <- c('matrix_assay',
                                  'matrix_class',
                                  'matrix_mode',
                                  'matrix_type',
                                  'matrix_compress',
                                  'matrix_path',
                                  'matrix_group',
                                  'matrix_buffer_size',
                                  'matrix_chunk_size',
                                  'matrix_overwrite',
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
  default_matrix_group <- '/'
  default_matrix_buffer_size <- 8192L
  default_matrix_chunk_size <- 1024L
  default_matrix_overwrite <- FALSE


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

      assay_control_out[['matrix_overwrite']] <- select_assay_parameter_value('matrix_overwrite', assay_control, assay_control_default, default_matrix_overwrite)
      if(!is.logical(assay_control_out[['matrix_overwrite']])) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out[['matrix_overwrite']],
                             ' is not TRUE or FALSE.'))
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
      message(indent, '  matrix_overwrite: ', ifelse(!is.null(assay_control[['matrix_overwrite']]), assay_control[['matrix_overwrite']], as.character(NA)))
    }
    else {
      stop('report_nn_control: unsupported pca class/mode/...\'', assay_control[['method']], '\'')
    }
  }
  else {
    stop('report_nn_control: unsupported pca class/mode/...\'', assay_control[['method']], '\'')
  }
}


push_matrix_path <- function(mat, matrix_type) {
  message('**** push_matrix_path: ', mat, ' -- ', matrix_type)
  x <- get_global_variable('monocle_gc_matrix_path')
  message('monocle_gc_matrix_path: start:\n', paste0(x, collapse='\n'))
  x <- append(x, mat)
  set_global_variable('monocle_gc_matrix_path', x)

  x <- get_global_variable('monocle_gc_matrix_path')
  message('monocle_gc_matrix_path: end:\n', paste0(x, collapse='\n'))
}


