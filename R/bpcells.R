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
#'@noRD
# Usage
#   storage_class: 'CsparseMatrix'
#     storage_mode: 'mem'
#   storage_class: 'BPCells'
#     storage_mode: 'mem'
#       storage_type_values: 'uint32_t', 'float', 'double'
#       storage_compress: TRUE, FALSE
#     storage_mode: 'dir'
#       storage_type_values: 'uint32_t', 'float', 'double' default: 'uint32_t'
#       storage_path: <path to directory or file> default: temporary directory in pwd
#       storage_compress: TRUE, FALSE
#       storage_buffer_size: <integer> default: 8192L
#       storage_overwrite: TRUE, FALSE default: FALSE
#     storage_mode: '10xhdf5'
#       storage_path: <path to directory or file> default: temporary file in pwd
#       storage_group: default '/'
#       storage_compress: TRUE, FALSE
#       storage_buffer_size: <integer> default: 8192L
#       storage_chunk_size: <integer> default: 1024L
#       storage_overwrite: TRUE, FALSE default: FALSE
set_assay_control <- function(assay_control=list(), assay_control_default=list()) {

  assertthat::assert_that(methods::is(assay_control, "list"))
  assertthat::assert_that(methods::is(assay_control_default, "list"))

  allowed_control_parameters <- c('storage_class',
                                  'storage_mode',
                                  'storage_type',
                                  'storage_compress',
                                  'storage_path',
                                  'storage_group',
                                  'storage_buffer_size',
                                  'storage_chunk_size',
                                  'storage_overwrite')

  assertthat::assert_that(all(names(assay_control) %in% allowed_control_parameters),
                          msg = "set_assay_control: unknown variable in assay_control")
  assertthat::assert_that(all(names(assay_control_default) %in% allowed_control_parameters),
                          msg = "set_assay_control: unknown variable in assay_control_default")


  default_storage_class <- 'CsparseMatrix'
  default_storage_mode <- 'mem'
  default_storage_type <- 'uint32_t'
  default_storage_compress <- TRUE
  default_storage_path <- NULL
  default_storage_group <- '/'
  default_storage_buffer_size <- 8192L
  default_storage_chunk_size <- 1024L
  default_storage_overwrite <- FALSE

  assay_control_out = list()

  error_string = list()

  assay_control_out$storage_class <- select_assay_parameter_value('storage_class', assay_control, assay_control_default, default_storage_class)
  if(!(assay_control_out$storage_class %in% c('CsparseMatrix', 'BPCells'))) {
    error_string <- list(error_string, paste0('  ',
                         assay_control_out$storage_class,
                         ' is not a valid storage_class.'))
    stop(error_string)
  }

  if(assay_control_out$storage_class == 'CsparseMatrix') {
    assay_control_out$storage_mode <- select_assay_parameter_value('storage_mode', assay_control, assay_control_default, default_storage_mode)
    if(!(assay_control_out$storage_mode %in% c('mem'))) {
      error_string <- list(error_string, paste0('  ',
                           assay_control_out$storage_mode,
                           ' is not valid for storage_class CsparseMatrix.'))
      stop(error_string)
    }
  }

  if(assay_control_out$storage_class== 'BPCells') {
    assay_control_out$storage_mode <- select_assay_parameter_value('storage_mode', assay_control, assay_control_default, default_storage_mode)
    if(!(assay_control_out$storage_mode %in% c('mem', 'dir', '10xhdf5'))) {
      error_string <- list(error_string, paste0('  ',
                           assay_control_out$storage_mode,
                           ' is not a valid storage_mode for storage_class "BPCells".'))
      stop(error_string)
    }

    if(assay_control_out$storage_mode == 'mem') {
      assay_control_out$storage_type <- select_assay_parameter_value('storage_type', assay_control, assay_control_default, default_storage_type)
      if(!(assay_control_out$storage_type %in% c('uint32_t', 'float', 'double'))) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out$storage_type,
                             ' is not a valid storage_type for storage_class "BPCells".'))
      }

      assay_control_out$storage_compress <- select_assay_parameter_value('storage_compress', assay_control, assay_control_default, default_storage_compress)
      if(!is.logical(assay_control_out$storage_compress)) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out$storage_compress,
                             ' is not TRUE or FALSE.'))
      }
      if(length(error_string) > 0) {
        stop(error_string)
      }
    } else
    if(assay_control_out$storage_mode == 'dir') {
      assay_control_out$storage_type <- select_assay_parameter_value('storage_type', assay_control, assay_control_default, default_storage_type)
      if(!(assay_control_out$storage_type %in% c('uint32_t', 'float', 'double'))) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out$storage_type,
                             ' is not a valid storage_type for storage_class "BPCells".'))
      }

      assay_control_out$storage_path <- select_assay_parameter_value('storage_path', assay_control, assay_control_default, default_storage_path)
      if(!is.null(assay_control_out$storage_path) &&
         !is.character(assay_control_out$storage_path)) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out$storage_path,
                             ' is not a valid storage_path.'))

      }

      assay_control_out$storage_compress <- select_assay_parameter_value('storage_compress', assay_control, assay_control_default, default_storage_compress)
      if(!is.logical(assay_control_out$storage_compress)) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out$storage_compress,
                             ' is not TRUE or FALSE.'))
      }

      assay_control_out$storage_buffer_size <- select_assay_parameter_value('storage_buffer_size', assay_control, assay_control_default, default_storage_buffer_size)
      if(!is.integer(assay_control_out$storage_buffer_size)) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out$storage_buffer_size,
                             ' is not an integer.'))
      }

      assay_control_out$storage_overwrite <- select_assay_parameter_value('storage_overwrite', assay_control, assay_control_default, default_storage_overwrite)
      if(!is.logical(assay_control_out$storage_overwrite)) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out$storage_overwrite,
                             ' is not TRUE or FALSE.'))
      }

      if(length(error_string) > 0) {
        stop(error_string)
      }
    } else
      if(assay_control_out$storage_mode == '10xhdf5') {
      assay_control_out$storage_path <- select_assay_parameter_value('storage_path', assay_control, assay_control_default, default_storage_path)
      if(!is.null(assay_control_out$storage_path) &&
         !is.character(assay_control_out$storage_path)) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out$storage_path,
                             ' is not a valid storage_path.'))

      }

      assay_control_out$storage_group <- select_assay_parameter_value('storage_group', assay_control, assay_control_default, default_storage_group)
      if(!is.character(assay_control_out$storage_group)) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out$storage_group,
                             ' is not a valid storage_path.'))

      }

      assay_control_out$storage_compress <- select_assay_parameter_value('storage_compress', assay_control, assay_control_default, default_storage_compress)
      if(!is.logical(assay_control_out$storage_compress)) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out$storage_compress,
                             ' is not TRUE or FALSE.'))
      }

      assay_control_out$storage_buffer_size <- select_assay_parameter_value('storage_buffer_size', assay_control, assay_control_default, default_storage_buffer_size)
      if(!is.integer(assay_control_out$storage_buffer_size)) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out$storage_buffer_size,
                             ' is not an integer.'))
      }

      assay_control_out$storage_chunk_size <- select_assay_parameter_value('storage_chunk_size', assay_control, assay_control_default, default_storage_chunk_size)
      if(!is.integer(assay_control_out$storage_chunk_size)) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out$storage_chunk_size,
                             ' is not an integer.'))
      }

      assay_control_out$storage_overwrite <- select_assay_parameter_value('storage_overwrite', assay_control, assay_control_default, default_storage_overwrite)
      if(!is.logical(assay_control_out$storage_overwrite)) {
        error_string <- list(error_string, paste0('  ',
                             assay_control_out$storage_overwrite,
                             ' is not TRUE or FALSE.'))
      }

      if(length(error_string) > 0){
        stop(error_string)
      }
    }
  }

  return(assay_control_out)
}


