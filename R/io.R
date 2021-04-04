# Save Annoy model files.
# Complexities
#   o  uwot annoy model object is declared to be subject to change
#   o  uwot annoy object changes some time between releases
#      0.1.8 and 0.1.10: the annoy index moved from
#      object$nn_index to object$nn_index$ann
#   o  the uwot annoy object may have more than one metric
#   o  see uwot/R/uwot.R functions save_uwot and load_uwot
# Notes
#   o  uwot 0.1.10 has nn_index$type='annoyv1', which may be a
#      umap annoy object version number. It does not exist for
#      uwot 0.1.8. nn_index$ann does not exist either.
#   o  the nn_index parameter is the nn_index object in the
#      uwot umap model and returned by uwot::annoy_build()
#         uwot v0.1.8::annoy_build() returns
#           ann <- uwot::create_ann()
#             ...
#           return ann
#         uwot v0.1.10::annoy_build() returns
#           annoy <- annoy_create(metric, nc)
#                 ...
#               ann <- annoy$ann
#                 ...
#               for (i in 1:nr) {
#                 ann$addItem(i - 1, X[i, ])
#               }
#                 ...
#           return annoy
#        where annoy is assigned to nn_index in either the
#        UMAP model or by uwot::annoy_build().
# untested code when is.null(nn_index[['type']])
save_annoy_index <- function(nn_index, file_name) {
  if( !is.null(nn_index[['type']])) {
    if(nn_index[['type']] == 'annoyv1') {
      nn_index[['ann']]$save(file_name)
    } else {
      stop('Unrecognized uwot annoy index type')
    }
  } else {
    nn_index$save(file_name)
  }
}

# see comments for save_annoy_index
# untested code when is.null(nn_index[['type']])
load_annoy_index <- function(nn_index, file_name, metric, ndim) {
  if( !is.null(nn_index[['type']])) {
    if(nn_index[['type']] == 'annoyv1') {
      nn_index[['ann']] <- uwot:::create_ann(metric, ndim)
      nn_index[['ann']]$load(file_name)
    } else {
      stop('Unrecognized uwot annoy index type')
    }
  } else {
    nn_index <- uwot:::create_ann(metric, ndim)
    nn_index$load(file_name)
  }
  nn_index
}


# Save umap annoy indexes to files and return md5sum
# value(s) as either a character string, in case of
# one metric, or a list, in case of more than one matric.
save_umap_nn_indexes <- function(umap_model, file_name) {
  metrics <- names(umap_model[['metric']])
  n_metrics <- length(metrics)
  if(n_metrics == 1) {
    save_annoy_index(umap_model[['nn_index']], file_name)
    md5sum_umap_index <- tools::md5sum(file_name)
  } else {
    warn('save_umap_nn_indexes is untested with more than one umap metric')
    md5sum_vec <- character()
    for(i in 1:n_metrics) {
      file_name_expand <- paste0(file_name, i)
      save_annoy_index(umap_model[['nn_index']][[i]], file_name_expand)
      md5sum <- tools::md5sum(file_name_expand)
      append(md5sum_vec, md5sum)
    }
    md5sum_umap_index <- paste(md5sum_vec, collapse='_')
  }
  md5sum_umap_index
}


# Load umap annoy indexes into umap_model and return umap_model.
load_umap_nn_indexes <- function(umap_model, file_name, md5sum_umap_index) {
  metrics <- names(umap_model[['metric']])
  n_metrics <- length(metrics)
  if(n_metrics == 1) {
    md5sum <- tools::md5sum(file_name)
    if(!is.null(md5sum_umap_index) && md5sum != md5sum_umap_index) {
      stop('The UMAP annoy index file, \'', file_name, '\', differs from the file made using the save_reduce_dimension_model() function.')
    }
    metric <- metrics[[1]]
    annoy_metric <- if(metric == 'correlation') 'cosine' else metric
    annoy_ndim <- umap_model[['metric']][[1]][['ndim']]
    umap_model[['nn_index']] <- load_annoy_index(umap_model[['nn_index']], file_name, annoy_metric, annoy_ndim)
  } else {
    warn('load_umap_nn_indexes is untested with more than one umap metric')
    if(!is.null(md5sum_umap_index)) {
      md5sum_vec <- unlist(strsplit(md5sum_umap_index, '_', fixed=TRUE))
    }
    for(i in 1:n_metrics) {
      file_name_expand <- paste0(file_name, i)
      md5sum <- tools::md5sum(file_name_expand)
      if(!is.null(md5sum_umap_index) && md5sum != md5sum_vec[[i]]) {
        stop('The UMAP annoy index file, \'', file_name_expand, '\', differs from the file made using the save_reduce_dimension_model() function.')
      }
      metric <- metrics[[i]]
      annoy_metric <- if(metric == 'correlation') 'cosine' else metric
      annoy_ndim <- length(umap_model[['metric']][[i]])
      umap_model[['nn_index']][[i]] <- load_annoy_index(umap_model[['nn_index']][[i]], file_name_expand, annoy_metric, annoy_ndim)
    }
  }
  umap_model
}


#
# object_name_to_string() is used to save the object
# name as a string in file_index.rds.
#
object_name_to_string <- function( object ) {
  str <- deparse(substitute(object))
  return( str )
}


#
# Report files saved.
#
report_files_saved <- function(file_index) {
  appendLF <- TRUE
  processes <- list()
  files <- file_index[['files']]
  for( i in seq_along(files[['cds_object']])) {
    cds_object <- files[['cds_object']][[i]]
    method <- files[['method']][[i]]
    if(cds_object == 'cds') {
      process <- 'cell_data_set'
      method <- 'full_cds'
    } else
    if(cds_object == 'preprocess_aux') {
      if(method == 'Aligned') {
        process <- 'align_cds'
      } else
      if(method == 'PCA' || method == 'LSI') {
        process <- 'preprocess_cds'
      } else {
        stop('Unrecognized preprocess method \'', method, '\'')
      }
    } else
    if(cds_object == 'reduce_dim_aux') {
      process <- 'reduce_dimension'
    } else {
      stop('Unrecognized cds_object value \'', files[['cds_object']][[i]], '\'')
    }
    file_format <- files[['file_format']][[i]]
    if(file_format == 'rds') {
      file_type <- 'RDS'
    } else
    if(file_format == 'hdf5') {
      file_type <- 'RDS_HDF5'
    } else
    if(file_format == 'annoy_index') {
      file_type <- 'NN_index'
    } else
    if(file_format == 'umap_annoy_index') {
      file_type <- 'UMAP_NN_index'
    } else {
      stop('Unrecognized file_format value \'', file_format, '\'')
    }

    file_name <- basename(files[['file_path']][[i]])
    message('  ', file_name, '  (', method, '  ', file_type, '  from  ', process, ')', appendLF=appendLF)
  }
}


#
#' Save cell_data_set transform models.
#'
#' Save the transform models to the specified directory
#' by writing the R objects to RDS files and
#' the Annoy nearest neighbor indexes to individual index files. 
#' save_transform_models saves transform models made by running
#' the preprocess_cds, align_cds, and reduce_dimension functions
#' on an initial cell_data_set. Subsequent cell_data_sets are
#' transformed into the reduced dimension space of the initial
#' cds by loading the new data into a new cds, loading the
#' initial data set transform models into the new cds using
#' the load_transform_models function, and applying those transform models
#' to the new data set using the preprocess_transform,
#' align_transform, and reduce_dimension_transform functions.
#' In this case, do not run the preprocess_cds, align_cds, or
#' reduce_dimension functions on the new cds. Additionally,
#' save_transform_models saves Annoy nearest neighbor indexes
#' when the preprocess_cds, align_cds, and reduce_dimension
#' functions are run with the build_nn_index=TRUE parameter. These
#' indexes are used to find matches between cells in the new
#' processed cds and the initial cds using the xxx functions.
#' save_transform_models scans the initial cell_data_set for
#' models and Annoy nearest neighbor indexes and saves those
#' that it finds. It does not save transform models
#' and indexes that were not made using preprocess_cds,
#' align_cds, and reduce_dimension.
#'
#' @param cds A cell_data_set with existing models.
#' @param directory_path A string giving the name of the directory
#'   in which to write the model files.
#' @param comment A string with optional notes that is saved with 
#'   the objects.
#' @param verbose A boolean determining whether to print information
#'   about the saved files.
#'
#' @return None.
#'
#' @export
save_transform_models <- function( cds, directory_path, comment="", verbose=TRUE) {
  appendLF <- TRUE
  # file information is written to an RDS file
  # in directory_path
  #   cds_object: preprocess_aux |r reduce_dim_aux
  #   method: PCA | LSI | Aligned | UMAP ...
  #   object_spec: ex: cds@reduce_dim_aux[[method]]
  #   file_format: rds | annoy_index | umap_annoy_index
  #   file_path: path within directory_path (need file name)
  #   file_md5sum: md5sum of file(s)
  file_index <- list( 'save_function' = 'save_transform_models',
                      'archive_date' = Sys.time(),
                      'r_version' = R.Version()$version.string,
                      'uwot_version' = packageVersion('uwot'),
                      'monocle_version' = packageVersion('monocle3'),
                      'cds_version' = metadata(cds)$cds_version,
                      'archive_version' = '1.0.0',
                      'directory' = directory_path,
                      'comment' = comment,
                      'files' = data.frame(cds_object = character(0),
                                           method = character(0),
                                           object_spec = character(0),
                                           file_format = character(0),
                                           file_path = character(0),
                                           file_md5sum = character(0),
                                           stringsAsFactors = FALSE))

  # Gather preprocess methods and whether each has an annoy index.
  methods_preprocess <- list()
  for(method in names(cds@preprocess_aux)) {
    methods_preprocess[[method]] <- list()
    methods_preprocess[[method]][['rds_path']] <- file.path(directory_path, paste0('ppc_', tolower(method), '_transform_model.rds'))
    methods_preprocess[[method]][['nn_index_path']] <- file.path(directory_path, paste0('ppc_', tolower(method), '_transform_model_nn.idx'))
    if(!is.null(cds@preprocess_aux[[method]][['nn_index']][['annoy_index']]))
      methods_preprocess[[method]][['has_nn_index']] <- TRUE
    else
      methods_preprocess[[method]][['has_nn_index']] <- FALSE
  }

  # Gather reduce_dimension methods and whether each has an annoy index.
  methods_reduce_dim <- list()
  for( method in names(cds@reduce_dim_aux)) {
    methods_reduce_dim[[method]] <- list()
    methods_reduce_dim[[method]][['rds_path']] <- file.path(directory_path, paste0('rdd_', tolower(method), '_transform_model.rds'))
    methods_reduce_dim[[method]][['umap_index_path']] <- file.path(directory_path, paste0('rdd_', tolower(method), '_transform_model_umap.idx'))
    methods_reduce_dim[[method]][['nn_index_path']] <- file.path(directory_path, paste0('rdd_', tolower(method), '_transform_model_nn.idx'))
    if(!is.null(cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']]))
      methods_reduce_dim[[method]][['has_nn_index']] <- TRUE
    else
      methods_reduce_dim[[method]][['has_nn_index']] <- FALSE
  }

  # Make directory if necessary.
  dir.create(path = directory_path, showWarnings=FALSE, recursive=TRUE, mode='0700')

  # Remove files, if they exist.
  for(method in names(methods_preprocess)) {
    if(file.exists(methods_preprocess[[method]][['rds_path']]))
      file.remove(methods_preprocess[[method]][['rds_path']])
    if(file.exists(methods_preprocess[[method]][['nn_index_path']]))
      file.remove(methods_preprocess[[method]][['nn_index_path']])
  }
  for(method in names(methods_reduce_dim)) {
    if(file.exists(methods_reduce_dim[[method]][['rds_path']]))
       file.remove(methods_reduce_dim[[method]][['rds_path']])
    if(file.exists(methods_reduce_dim[[method]][['umap_index_path']]))
       file.remove(methods_reduce_dim[[method]][['umap_index_path']])
    if(file.exists(methods_reduce_dim[[method]][['nn_index_path']]))
       file.remove(methods_reduce_dim[[method]][['nn_index_path']])
  }

  # Save preprocess_cds annoy indexes.
  # Notes:
  #   o  save RDS files before the corresponding index files in
  #      order to enable loading.
  #
  for(method in names(methods_preprocess)) {
    tryCatch(
      {
        saveRDS(cds@preprocess_aux[[method]], file=methods_preprocess[[method]][['rds_path']])
      },
      error = function(cnd) {
                     message('Error writing file \'', methods_preprocess[[method]][['rds_path']], '\': ', cnd, appendLF=appendLF)
                     return(NULL)
      },
      finally = {
        md5sum <- tools::md5sum(methods_preprocess[[method]][['rds_path']])
        file_index[['files']] <- rbind(file_index[['files']],
                                       data.frame(cds_object = 'preprocess_aux',
                                                  method = method,
                                                  object_spec = object_name_to_string(cds@preprocess_aux[[method]]),
                                                  file_format = 'rds',
                                                  file_path = methods_preprocess[[method]][['rds_path']],
                                                  file_md5sum = md5sum,
                                                  stringsAsFactors = FALSE))
      })
    if(methods_preprocess[[method]][['has_nn_index']]) {
      tryCatch(
        {
          save_annoy_index(cds@preprocess_aux[[method]][['nn_index']][['annoy_index']], methods_preprocess[[method]][['nn_index_path']])
        },
        error = function(cnd) {
                       message('Error writing file \'', methods_preprocess[[method]][['nn_index_path']], '\': ', cnd, appendLF=appendLF)
                       return(NULL)
        },
        finally = {
          md5sum <- tools::md5sum(methods_preprocess[[method]][['nn_index_path']])
          file_index[['files']] <- rbind(file_index[['files']],
                                         data.frame(cds_object = 'preprocess_aux',
                                                    method = method,
                                                    object_spec = object_name_to_string(cds@preprocess_aux[[method]][['nn_index']][['annoy_index']]),
                                                    file_format = 'annoy_index',
                                                    file_path = methods_preprocess[[method]][['nn_index_path']],
                                                    file_md5sum = md5sum,
                                                    stringsAsFactors = FALSE))
        })
    }
  }

  # Save reduce_dimension annoy indexes.
  for(method in names(methods_reduce_dim)) {
    tryCatch(
      {
        saveRDS(cds@reduce_dim_aux[[method]], file=methods_reduce_dim[[method]][['rds_path']])
      },
      error = function(cnd) {
                     message('Error writing file \'', methods_reduce_dim[[method]][['rds_path']], '\': ', cnd, appendLF=appendLF)
                     return(NULL)
      },
      finally = {
        md5sum <- tools::md5sum(methods_reduce_dim[[method]][['rds_path']])
        file_index[['files']] <- rbind(file_index[['files']],
                                       data.frame(cds_object = 'reduce_dim_aux',
                                                  method = method,
                                                  object_spec = object_name_to_string(cds@reduce_dim_aux[[method]]),
                                                  file_format = 'rds',
                                                  file_path = methods_reduce_dim[[method]][['rds_path']],
                                                  file_md5sum = md5sum,
                                                  stringsAsFactors = FALSE))
      })
    if(method == 'UMAP') {
      tryCatch(
        {
          md5sum <- save_umap_nn_indexes(cds@reduce_dim_aux[[method]][['model']][['umap_model']], methods_reduce_dim[[method]][['umap_index_path']])
        },
        error = function(cnd) {
                       message('Error writing file \'', methods_reduce_dim[[method]][['umap_index_path']], '\': ', cnd, appendLF=appendLF)
                       return(NULL)
        },
        finally = {
          file_index[['files']] <- rbind(file_index[['files']],
                                         data.frame(cds_object = 'reduce_dim_aux',
                                                    method = method,
                                                    object_spec = object_name_to_string(cds@reduce_dim_aux[[method]][['model']][['umap_model']]),
                                                    file_format = 'umap_annoy_index',
                                                    file_path = methods_reduce_dim[[method]][['umap_index_path']],
                                                    file_md5sum = md5sum,
                                                    stringsAsFactors = FALSE))
        })
    }
    if(methods_reduce_dim[[method]][['has_nn_index']]) {
      tryCatch(
        {
          save_annoy_index(cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']], methods_reduce_dim[[method]][['nn_index_path']])
        },
        error = function(cnd) {
                       message('Error writing file \'', methods_reduce_dim[[method]][['nn_index_path']], '\': ', cnd, appendLF=appendLF)
                       return(NULL)
        },
        finally = {
          md5sum <- tools::md5sum(methods_reduce_dim[[method]][['nn_index_path']])
          file_index[['files']] <- rbind(file_index[['files']],
                                         data.frame(cds_object = 'reduce_dim_aux',
                                                    method = method,
                                                    object_spec = object_name_to_string(cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']]),
                                                    file_format = 'annoy_index',
                                                    file_path = methods_reduce_dim[[method]][['nn_index_path']],
                                                    file_md5sum = md5sum,
                                                    stringsAsFactors = FALSE))
        })
    }
  }

  # Save file_index.rds.
  saveRDS(file_index, file=file.path(directory_path, 'file_index.rds'))

  if(verbose) {
    report_files_saved(file_index)
  }
}


#' Load transform models into a cell_data_set.
#'
#' Load transform models, which were saved using save_transform_models,
#' into a cell_data_set. This function over-writes existing models in
#' the cell_data_set. For more information read the help information
#' for save_transform_models.
#'
#' @param cds A cell_data_set to be transformed using the models.
#' @param directory_path A string giving the name of the directory
#'   from which to read the model files.
#'
#' @return A cell_data_set with the transform models loaded by
#'   load_transform_models.
#'
#' @export
load_transform_models <- function(cds, directory_path) {
  appendLF <- TRUE
  # Check for directory.
  if(!file.exists(directory_path))
    stop('Directory \'', directory_path, '\' does not exist.')

  # Check for file_index.rds.
  file_index_path <- file.path(directory_path, 'file_index.rds')
  if(!file.exists(file_index_path))
    stop('Missing file index file \'', file_index_path, '\'')

  # Read file index.
  file_index <- tryCatch(
    {
      readRDS(file_index_path)
    },
    error = function(cnd) {
              message('Error reading file \'', file_index_path, '\': ', cnd, appendLF=appendLF);
              return(NULL)
    }
  )

  # Check that this is a save_transform_models archive.
  if(file_index[['save_function']] != 'save_transform_models') {
    stop('The files in ', directory_path, ' are not from save_transform_models.')
  }

  # Write stored comment field.
  if(length(file_index[['comment']]) > 1) {
    message('File comment: ', file_index[['comment']], appendLF=appendLF)
  }

  # Loop through the files in file_index.rds in order
  # to restore objects.
  for(ifile in seq_along(file_index[['files']][['cds_object']])) {
    file_path <- file_index[['files']][['file_path']][[ifile]]
    file_format <- file_index[['files']][['file_format']][[ifile]]
    cds_object <- file_index[['files']][['cds_object']][[ifile]]
    method <- file_index[['files']][['method']][[ifile]]
    md5sum <- file_index[['files']][['file_md5sum']][[ifile]]

    #
    # For UWOT UMAP annoy index, the function load_umap_nn_indexes
    # checks md5sums internally so don't check here.
    #
    md5sum_file <- tools::md5sum(file_path)
    if(!(cds_object == 'reduce_dim_aux' &&
         method == 'UMAP' &&
         file_format == 'umap_nn_index' &&
         nchar(md5sum) > 32)) {
      if(md5sum_file != md5sum) {
        stop('md5sum mismatch for file \'', file_path, '\'')
      }
    }

    #
    # Note:
    #   o  expect that the RDS file for a method
    #      appears before index files for the method
    #
    if(cds_object == 'preprocess_aux') {
      if(file_format == 'rds') {
        cds@preprocess_aux[[method]] <- tryCatch(
          {
            readRDS(file_path)
          },
          error = function(cnd) {
            message('Error reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
          })
     
      } else
      if(file_format == 'annoy_index') {
        metric <- cds@preprocess_aux[[method]][['nn_index']][['annoy_metric']]
        ndim <- cds@preprocess_aux[[method]][['nn_index']][['annoy_ndim']]
        cds@preprocess_aux[[method]][['nn_index']][['annoy_index']] <- tryCatch(
         {
           load_annoy_index(cds@preprocess_aux[[method]][['nn_index']][['annoy_index']], file_path, metric, ndim)
         },
         error = function(cds) {
           message('Error reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
         })
      } else {
        stop('Unrecognized file format value \'', file_format, '\'')
      }
    } else
    if(cds_object == 'reduce_dim_aux') {
      if(file_format == 'rds') {
        cds@reduce_dim_aux[[method]] <- tryCatch(
          { 
            readRDS(file_path)
          },
          error = function(cnd) {
            message('Error reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
          })
      } else
      if(file_format == 'annoy_index') {
        metric <- cds@reduce_dim_aux[[method]][['nn_index']][['annoy_metric']]
        ndim <- cds@reduce_dim_aux[[method]][['nn_index']][['annoy_ndim']]
        cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']] <- tryCatch(
          {
            load_annoy_index(cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']], file_path, metric, ndim)
          },
          error = function(cnd) {
            message('Error reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
          })
      } else
      if(method == 'UMAP' && file_format == 'umap_annoy_index') {
        cds@reduce_dim_aux[[method]][['model']][['umap_model']] <- tryCatch(
          {
            load_umap_nn_indexes(cds@reduce_dim_aux[[method]][['model']][['umap_model']], file_path, md5sum)
          },
          error = function(cnd) {
            message('Error reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
         })
      } else {
        stop('Unrecognized file format value \'', file_format, '\'')
      }
    } else {
      stop('Unrecognized cds_object value \'', cds_object, '\'')
    }
  }

  cds
}


#
# Check cds assays for HDF5Array objects.
#
test_hdf5_assays <- function(cds) {
  assays <- assays(cds)
  for( idx in seq_along(assays)) {
    asyl <- getListElement(assays, idx)
    hdf5_test<- unlist(DelayedArray:::seedApply(asyl, is, "HDF5ArraySeed"))
    if(any(unlist(hdf5_test))) return(TRUE)
  } 
  FALSE
}


#
#' Save a Monocle3 full cell_data_set.
#'
#' Save a Monocle3 full cell_data_set to a specified directory
#' by writing the R objects to an RDS file and
#' the Annoy nearest neighbor indexes to individual index files.
#' The assays objects are saved as HDF5Array files when
#' hdf5_assays=TRUE or when the cell_data_set assays are
#' HDF5Array objects. If any assay in the cell_data set is an
#' HDF5 object, all assays must be. When save_monocle_objects is
#' run with hdf5_assays=TRUE, the load_monocle_objects function
#' loads the saved assays into HDF5Array objects in the resulting
#' cell_data_set. Note: operations such as preprocess_cds that
#' are run on assays stored as HDF5Arrays are much, much slower
#' than the same operations run on assays stored as in-memory
#' matrices. You may want to investigate parameters related to
#' the Bioconductor DelayArray and BiocParallel packages in this
#' case.
#'
#' @param cds A cell_data_set to save.
#' @param directory_path A string giving the name of the directory
#'   in which to write the object files.
#' @param hdf5_assays A boolean determining whether the
#'   non-HDF5Array assays objects are saved as HDF5 files. At this
#'   time cell_data_set HDF5Array assay objects are stored as
#'   HDF5Assay files regardless of the hdf5_assays parameter value.
#' @param comment A string with optional notes that is saved with
#'   the objects.
#' @param verbose A boolean determining whether to print information
#'   about the saved files.
#'
#' @return None.
#'
#' @export
save_monocle_objects <- function(cds, directory_path, hdf5_assays=FALSE, comment="", verbose=TRUE) {
  appendLF <- TRUE
  # file information is written to an RDS file
  # in directory_path
  #   cds_object: cds | preprocess_aux | reduce_dim_aux
  #   method: PCA | LSI | Aligned | UMAP  ...
  #   object_spec: ex: cds@reduce_dim_aux[[method]]
  #   file_format: rds | annoy_index | umap_annoy_index
  #   file_path: path within directory_path (need file name)
  #   file_md5sum: md5sum of file(s)
  file_index <- list( 'save_function' = 'save_monocle_objects',
                      'archive_date' = Sys.time(),
                      'r_version' = R.Version()$version.string,
                      'uwot_version' = packageVersion('uwot'),
                      'hdf5array_version' = packageVersion('HDF5Array'),
                      'monocle_version' = packageVersion('monocle3'),
                      'cds_version' = metadata(cds)$cds_version,
                      'archive_version' = '1.0.0',
                      'directory' = directory_path,
                      'comment' = comment,
                      'files' = data.frame(cds_object = character(0),
                                           method = character(0),
                                           object_spec = character(0),
                                           file_format = character(0),
                                           file_path = character(0),
                                           file_md5sum = character(0),
                                           stringsAsFactors = FALSE))

  # Save assays as HDF5Array objects?
  hdf5_assay_flag <- hdf5_assays || test_hdf5_assays(cds)

  # Path of cds object file.
  rds_path <- file.path(directory_path, 'cds_object.rds')
  hdf5_path <- file.path(directory_path, 'hdf5_object')

  # Gather preprocess methods and whether each has an annoy index.
  methods_preprocess <- list()
  for(method in names(cds@preprocess_aux)) {
    methods_preprocess[[method]] <- list()
    methods_preprocess[[method]][['nn_index_path']] <- file.path(directory_path, paste0('ppc_', tolower(method), '_transform_model_nn.idx'))
    if(!is.null(cds@preprocess_aux[[method]][['nn_index']][['annoy_index']]))
      methods_preprocess[[method]][['has_nn_index']] <- TRUE
    else
      methods_preprocess[[method]][['has_nn_index']] <- FALSE
  }

  # Gather reduce_dimension method names for which indexes exist.
  methods_reduce_dim <- list()
  for(method in names(cds@reduce_dim_aux)) {
    methods_reduce_dim[[method]] <- list()
    methods_reduce_dim[[method]][['umap_index_path']] <- file.path(directory_path, paste0('rdd_', tolower(method), '_transform_model_umap.idx'))
    methods_reduce_dim[[method]][['nn_index_path']] <- file.path(directory_path, paste0('rdd_', tolower(method), '_transform_model_nn.idx'))
    if(!is.null(cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']]))
      methods_reduce_dim[[method]][['has_nn_index']] <- TRUE
    else
      methods_reduce_dim[[method]][['has_nn_index']] <- FALSE
  }

  # Make directory if necessary.
  dir.create(path = directory_path, showWarnings=FALSE, recursive=TRUE, mode='0700')

  # Remove files, if they exist.
  for(method in names(methods_preprocess)) {
    if(file.exists(methods_preprocess[[method]][['nn_index_path']]))
      file.remove(methods_preprocess[[method]][['nn_index_path']])
  }
  for(method in names(methods_reduce_dim)) {
    if(file.exists(methods_reduce_dim[[method]][['umap_index_path']]))
       file.remove(methods_reduce_dim[[method]][['umap_index_path']])
    if(file.exists(methods_reduce_dim[[method]][['nn_index_path']]))
       file.remove(methods_reduce_dim[[method]][['nn_index_path']])
  }

  #
  # Save cds object.
  # Notes:
  #   o  allow for HDF5Array assay objects.
  #
  if(!hdf5_assay_flag) {
    tryCatch(
      {
        saveRDS(cds, rds_path)
      },
      error = function(cnd) {
                       message('Error writing file \'', rds_path, '\': ', cnd, appendLF=appendLF)
                       return(NULL)
      },
      finally = {
        md5sum <- tools::md5sum(rds_path)
        file_index[['files']] <- rbind(file_index[['files']],
                                       data.frame(cds_object = 'cds',
                                                  method = NA,
                                                  object_spec = object_name_to_string(cds),
                                                  file_format = 'rds',
                                                  file_path = rds_path,
                                                  file_md5sum = md5sum,
                                                  stringsAsFactors = FALSE))
      })
  } else {
    tryCatch(
      {
        HDF5Array::saveHDF5SummarizedExperiment(cds, hdf5_path, replace=TRUE)
      },
      error = function(cnd) {
                       message('Error writing file \'', hdf5_path, '\': ', cnd, appendLF=appendLF)
                       return(NULL)
      },
      finally = {
        md5sum <- tools::md5sum(file.path(hdf5_path, 'se.rds'))
        file_index[['files']] <- rbind(file_index[['files']],
                                       data.frame(cds_object = 'cds',
                                                  method = NA,
                                                  object_spec = object_name_to_string(cds),
                                                  file_format = 'hdf5',
                                                  file_path = hdf5_path,
                                                  file_md5sum = md5sum,
                                                  stringsAsFactors = FALSE))
      })
  }

  # Save preprocess_cds annoy indexes.
  # Notes:
  #   o  save RDS files before the corresponding index files in
  #      order to enable loading.
  #
  for(method in names(methods_preprocess)) {
    if(methods_preprocess[[method]][['has_nn_index']]) {
      tryCatch(
        {
          save_annoy_index(cds@preprocess_aux[[method]][['nn_index']][['annoy_index']], methods_preprocess[[method]][['nn_index_path']])
        },
        error = function(cnd) {
                       message('Error writing file \'', methods_preprocess[[method]][['nn_index_path']], '\': ', cnd, appendLF=appendLF)
                       return(NULL)
        },
        finally = {
          md5sum <- tools::md5sum(methods_preprocess[[method]][['nn_index_path']])
          file_index[['files']] <- rbind(file_index[['files']],
                                         data.frame(cds_object = 'preprocess_aux',
                                                    method = method,
                                                    object_spec = object_name_to_string(cds@preprocess_aux[[method]][['nn_index']][['annoy_index']]),
                                                    file_format = 'annoy_index',
                                                    file_path = methods_preprocess[[method]][['nn_index_path']],
                                                    file_md5sum = md5sum,
                                                    stringsAsFactors = FALSE))
        })
    }
  }

  # Save reduce_dimension annoy indexes.
  for(method in names(methods_reduce_dim)) {
    if(method == 'UMAP') {
      tryCatch(
        {
          md5sum <- save_umap_nn_indexes(cds@reduce_dim_aux[[method]][['model']][['umap_model']], methods_reduce_dim[[method]][['umap_index_path']])
        },
        error = function(cnd) {
                       message('Error writing file \'', methods_reduce_dim[[method]][['umap_index_path']], '\': ', cnd, appendLF=appendLF)
                       return(NULL)
        },
        finally = {
          file_index[['files']] <- rbind(file_index[['files']],
                                         data.frame(cds_object = 'reduce_dim_aux',
                                                    method = method,
                                                    object_spec = object_name_to_string(cds@reduce_dim_aux[[method]][['model']][['umap_model']]),
                                                    file_format = 'umap_annoy_index',
                                                    file_path = methods_reduce_dim[[method]][['umap_index_path']],
                                                    file_md5sum = md5sum,
                                                    stringsAsFactors = FALSE))
        })
    }
    if(methods_reduce_dim[[method]][['has_nn_index']]) {
      tryCatch(
        {
          save_annoy_index(cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']], methods_reduce_dim[[method]][['nn_index_path']])
        },
        error = function(cnd) {
                       message('Error writing file \'', methods_reduce_dim[[method]][['nn_index_path']], '\': ', cnd, appendLF=appendLF)
                       return(NULL)
        },
        finally = {
          md5sum <- tools::md5sum(methods_reduce_dim[[method]][['nn_index_path']])
          file_index[['files']] <- rbind(file_index[['files']],
                                         data.frame(cds_object = 'reduce_dim_aux',
                                                    method = method,
                                                    object_spec = object_name_to_string(cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']]),
                                                    file_format = 'annoy_index',
                                                    file_path = methods_reduce_dim[[method]][['nn_index_path']],
                                                    file_md5sum = md5sum,
                                                    stringsAsFactors = FALSE))
        })
    }
  }

  # Save file_index.rds.
  saveRDS(file_index, file=file.path(directory_path, 'file_index.rds'))

  if(verbose) {
    report_files_saved(file_index)
  }
}


#
#' Load a full Monocle3 cell_data_set.
#'
#' Load a full Monocle3 cell_data_set, which was saved using
#' save_monocle_objects. For more information read the help
#' information for save_monocle_objects.
#'
#' @param directory_path A string giving the name of the directory
#'   from which to read the saved cell_data_set files.
#'
#' @return A cell_data_set.
#'
#' @export
load_monocle_objects <- function(directory_path) {
  appendLF <- TRUE
  # Check for directory.
  if(!file.exists(directory_path))
    stop('Directory \'', directory_path, '\' does not exist.')

  # Check for file_index.rds.
  file_index_path <- file.path(directory_path, 'file_index.rds')
  if(!file.exists(file_index_path))
    stop('Missing file index file \'', file_index_path, '\'')

  # Read file index.
  file_index <- tryCatch(
    {
      readRDS(file_index_path)
    },
    error = function(cnd) {
              message('Error reading file \'', file_index_path, '\': ', cnd, appendLF=appendLF);
              return(NULL)
    }
  )

  # Check that this is a save_monocle_objects archive.
  if(file_index[['save_function']] != 'save_monocle_objects') {
    stop('The files in ', directory_path, ' are not from save_monocle_objects.')
  }

  # Write stored comment field.
  if(length(file_index[['comment']]) > 1) {
    message('File comment: ', file_index[['comment']], appendLF=appendLF)
  }

  # Loop through the files in file_index.rds in order
  # to restore objects.
  for(ifile in seq_along(file_index[['files']][['cds_object']])) {
    file_path <- file_index[['files']][['file_path']][[ifile]]
    file_format <- file_index[['files']][['file_format']][[ifile]]
    cds_object <- file_index[['files']][['cds_object']][[ifile]]
    method <- file_index[['files']][['method']][[ifile]]
    md5sum <- file_index[['files']][['file_md5sum']][[ifile]]

    #
    # The functions load_umap_nn_indexes and
    # loadHDF5SummarizedExperiment check md5sums
    # internally so don't check here.
    #
    if(!(cds_object == 'reduce_dim_aux' &&
         method == 'UMAP' &&
         file_format == 'umap_nn_index' &&
         nchar(md5sum) > 32) &&
       file_format != 'hdf5') {
      md5sum_file <- tools::md5sum(file_path)
      if(md5sum_file != md5sum) {
        stop('md5sum mismatch for file \'', file_path, '\'')
      }
    }

    #
    # Note:
    #   o  expect that the RDS file for a method
    #      appears before index files for the method
    #
    if(cds_object == 'cds') {
      if(file_format == 'rds') {
        cds <- tryCatch(
          {
            readRDS(file_path)
          },
          error = function(cnd) {
            message('Error reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
          })
      } else
      if(file_format == 'hdf5') {
        cds <- tryCatch(
          {
            HDF5Array::loadHDF5SummarizedExperiment(file_path)
          },
          error = function(cnd) {
            message('Error reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
          })
      } else {
        stop('Unrecognized cds format value \'', file_format, '\'')
      }
    } else
    if(cds_object == 'preprocess_aux') {
      if(file_format == 'annoy_index') {
        metric <- cds@preprocess_aux[[method]][['nn_index']][['annoy_metric']]
        ndim <- cds@preprocess_aux[[method]][['nn_index']][['annoy_ndim']]
        cds@preprocess_aux[[method]][['nn_index']][['annoy_index']] <- tryCatch(
         {
           load_annoy_index(cds@preprocess_aux[[method]][['nn_index']][['annoy_index']], file_path, metric, ndim)
         },
         error = function(cds) {
           message('Error reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
         })
      } else {
        stop('Unrecognized file format value \'', file_format, '\'')
      }
    } else
    if(cds_object == 'reduce_dim_aux') {
      if(file_format == 'annoy_index') {
        metric <- cds@reduce_dim_aux[[method]][['nn_index']][['annoy_metric']]
        ndim <- cds@reduce_dim_aux[[method]][['nn_index']][['annoy_ndim']]
        cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']] <- tryCatch(
          {
            load_annoy_index(cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']], file_path, metric, ndim)
          },
          error = function(cnd) {
            message('Error reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
          })
      } else
      if(method == 'UMAP' && file_format == 'umap_annoy_index') {
        cds@reduce_dim_aux[[method]][['model']][['umap_model']] <- tryCatch(
          {
            load_umap_nn_indexes(cds@reduce_dim_aux[[method]][['model']][['umap_model']], file_path, md5sum)
          },
          error = function(cnd) {
            message('Error reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
         })
      } else {
        stop('Unrecognized file format value \'', file_format, '\'')
      }
    } else {
      stop('Unrecognized cds_object value \'', cds_object, '\'')
    }
  }

  cds
}

