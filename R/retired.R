#
#' Save a Monocle3 full cell_data_set.
#'
#' Save a Monocle3 full cell_data_set to a specified directory
#' by writing the R objects to an RDS file and
#' the Annoy nearest neighbor indexes to individual index files.
#'
#' @param cds cell_data_set to save.
#' @param directory_path A string giving the directory in
#'   which to write the files.
#'
#' @return None.
#'
#' @export
save_monocle_objects_retired <- function(cds, directory_path) {
  md5sum_list = list()
  # Gather preprocess method names for which indexes exist.
  methods_preprocess <- c()
  for(method in names(cds@preprocess_aux)) {
    if(!is.null(cds@preprocess_aux[[method]][['nn_index']][['annoy_index']]))
      methods_preprocess <- c(methods_preprocess, method)
  }

  # Gather reduce_dimension method names for which indexes exist.
  methods_reduce_dim <- c()
  for( method in names(cds@reduce_dim_aux)) {
    if(!is.null(cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']]))
      methods_reduce_dim <- c(methods_reduce_dim, method)
  }

  # Make file paths.
  path_cds_rds <- file.path(directory_path, 'cell_data_set.rds')
  path_umap_ann_index <- file.path(directory_path, 'rdd_umap_rd_ann.idx')

  paths_preprocess = list()
  for(method in methods_preprocess)
    paths_preprocess[[method]] <- file.path(directory_path, paste0('ppc_', tolower(method), '_nn_ann.idx'))

  paths_reduce_dim = list()
  for(method in methods_reduce_dim)
    paths_reduce_dim[[method]] <- file.path(directory_path, paste0('rdd_', tolower(method), '_nn_ann.idx'))

  # Make directory if necessary.
  dir.create(path = directory_path, showWarnings=FALSE, recursive=TRUE, mode='0700')

  # Remove files, if they exist.
  for(method in methods_preprocess) {
    if(file.exists(paths_preprocess[[method]]))
      file.remove( paths_preprocess[[method]])
  }
  for(method in methods_reduce_dim) {
    if(file.exists(paths_reduce_dim[[method]]))
       file.remove(paths_reduce_dim[[method]])
  }
  if(file.exists(path_cds_rds))
    file.remove(path_cds_rds)
  if(file.exists(path_umap_ann_index))
    file.remove(path_umap_ann_index)

  # save CDS in RDS file
  saveRDS(cds, path_cds_rds)
  md5sum_file <- tools::md5sum(path_cds_rds)
  md5sum_list <- c(md5sum_list, md5sum_file)

  # save preprocess_cds annoy indexes.
  for(method in methods_preprocess) {
    save_annoy_index(cds@preprocess_aux[[method]][['nn_index']][['annoy_index']], paths_preprocess[[method]])
    md5sum_file <- tools::md5sum(paths_preprocess[[method]])
    md5sum_list <- c(md5sum_list, md5sum_file)
  }

  # save reduce_dimension annoy indexes.
  for(method in methods_reduce_dim) {
    save_annoy_index(cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']], paths_reduce_dim[[method]])
    md5sum_file <- tools::md5sum(paths_reduce_dim[[method]])
    md5sum_list <- c(md5sum_list, md5sum_file)
  }

  # save reduce_dimension UMAP annoy index
  if(!is.null(cds@reduce_dim_aux[[method]][['model']][['umap_model']])) {
    md5sum_file <- save_umap_nn_indexes(cds@reduce_dim_aux[[method]][['model']][['umap_model']], path_umap_ann_index)
    md5sum_list <- c(md5sum_list, md5sum_file)
  }

  names(md5sum_list) <- basename(names(md5sum_list))
  write.table(md5sum_list, file=file.path(directory_path,'md5sums.dat'), row.names=FALSE)
}


#
#' Load Monocle3 cell_data_set objects.
#'
#' Load a full Monocle3 cell_data_set from the directory
#' made by the previously run save_monocle_objects function.
#'
#' @param directory_path A string giving the directory from
#'   which to read the previously saved cell_data_set files.
#'
#' @return A cell_data_set.
#'
#' @export
load_monocle_objects_retired <- function(directory_path) {
  appendLF <- TRUE

  # Make file paths.
  path_cds_rds <- file.path(directory_path, 'cell_data_set.rds')
  path_umap_ann_index <- file.path(directory_path, 'rdd_umap_rd_ann.idx')

  # Check whether directory exists.
  if(!file.exists(directory_path))
    stop('Directory \'', directory_path, '\' does not exist.')

  # Check whether cds rds file exists.
  if(!file.exists(path_cds_rds))
    stop('Missing RDS file \'', path_cds_rds, '\'')

  # Read cds rds file.
  cds <- tryCatch(
           {
             readRDS(path_cds_rds)
           },
           error = function(cnd) {
                     message('Error reading file \'', path_cds_rds, '\': ', cnd, appendLF=appendLF);
                     return(NULL)
           },
           finally = {
                       message('loaded RDS file.', appendLF=appendLF)
           }
  )

  # Gather preprocess method names for which indexes exist.
  methods_preprocess <- c()
  for(method in names(cds@preprocess_aux)) {
    if(!is.null(cds@preprocess_aux[[method]][['nn_index']][['annoy_index']]))
      methods_preprocess <- c(methods_preprocess, method)
  }

  # Gather reduce_dimension method names for which indexes exist.
  methods_reduce_dim <- c()
  for( method in names(cds@reduce_dim_aux)) {
    if(!is.null(cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']]))
      methods_reduce_dim <- c(methods_reduce_dim, method)
  }

  paths_preprocess = list()
  for(method in methods_preprocess)
    paths_preprocess[[method]] <- file.path(directory_path, paste0('ppc_', tolower(method), '_nn_ann.idx'))

  paths_reduce_dim = list()
  for(method in methods_reduce_dim)
    paths_reduce_dim[[method]] <- file.path(directory_path, paste0('rdd_', tolower(method), '_nn_ann.idx'))

  # Load preprocess annoy indexes.
  for(method in methods_preprocess) {
    metric <- cds@preprocess_aux[[method]][['nn_index']][['annoy_metric']]
    ndim <- cds@preprocess_aux[[method]][['nn_index']][['annoy_ndim']]
    cds@preprocess_aux[[method]][['nn_index']][['annoy_index']] <- tryCatch(
             {
               load_annoy_index(cds@preprocess_aux[[method]][['nn_index']][['annoy_index']], paths_preprocess[[method]], metric, ndim)
             },
             error = function(cnd) {
                       message('Error reading file \'', paths_preprocess[[method]], '\': ', cnd, appendLF=appendLF)
                       return(NULL)
             },
             finally = {
                         message('Loaded ', method, ' Annoy nearest neighbor index file.', appendLF=appendLF)
             })
  }

  # Load reduce_dimension indexes.
  for(method in methods_reduce_dim) {
    metric <- cds@reduce_dim_aux[[method]][['nn_index']][['annoy_metric']]
    ndim <- cds@reduce_dim_aux[[method]][['nn_index']][['annoy_ndim']]
    cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']] <- tryCatch(
             {
               load_annoy_index(cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']], paths_reduce_dim[[method]], metric, ndim)
             },
             error = function(cnd) {
                       message('Error reading file \'', paths_preprocess[[method]], '\': ', cnd, appendLF=appendLF)
                       return(NULL)
             },
             finally = {
                         message('Loaded ', method, ' Annoy nearest neighbor index file.', appendLF=appendLF)
             })
  }

  # Load UMAP annoy index.
  cds@reduce_dim_aux[['UMAP']][['model']][['umap_model']] <- tryCatch(
           {
             load_umap_nn_indexes(cds@reduce_dim_aux[['UMAP']][['model']][['umap_model']], path_umap_ann_index, NULL)
           },
           error = function(cnd) {
                     message('Error reading file \'', path_umap_ann_index, '\': ', cnd, appendLF=appendLF)
                     return(NULL)
           },
           finally = {
                       message('Loaded UMAP internal Annoy nearest neighbor index file.', appendLF=appendLF)
           })
  cds
}


#' Save preprocess model.
#'
#' Save a preprocess model consisting of dimensionality reduction
#' model and preprocess nearest neighbor index objects, which
#' are used to transform new data sets.
#'
#' @param cds cell_data_set with an existing model, which was created using
#'   cds <- preprocess_cds(cds, build_nn_index=TRUE).
#' @param method A string indicating the method used to build the model
#'   that you want to save. Default is "PCA".
#' @param file_name_root. A string with the root of names given to the two
#'   files to which the model will be saved. The files are named
#'   <file_name_root>.ppc_rds and <file_name_root>.ppc_nn_index. Note: the
#'   <file_name_root>.ppc_rds file contains only the information needed to load
#'   the preprocess model into the cds object given in
#'   load_preprocess_model(cds,...).
#' @param comment An optional string that describes the model.
#'
#' @return None
#' @export
save_preprocess_model_retired <- function(cds, method = c('PCA', 'LSI'), file_name_root = NULL, comment = NULL) {
  model_file_version <- '1.0.0'
  method <- match.arg(method)
  if(is.null(cds@preprocess_aux[[method]])) {
    stop('preprocess_cds was not run on this cds.')
  } else if(is.null(cds@preprocess_aux[[method]][['nn_index']])) {
    stop('There is no nearest neighbor index -- run preprocess_cds() with build_nn_index=TRUE')
  }

  file_name_rds <- paste0(file_name_root, '.ppc_rds')
  file_name_nn_index1 <- paste0(file_name_root, '.ppc_nn_index')


  object <- list()
  object[['model_file_version']] <- model_file_version
  object[['model_file']] <- 'preprocess'
  object[['method']] <- method
  object[['bundle']] <- cds@preprocess_aux[[method]]
  object[['comment']] <- comment

  save_annoy_index(cds@preprocess_aux[[method]][['nn_index']][['annoy_index']], file_name_nn_index1)
  object[['md5sum_nn_index']] <- tools::md5sum(file_name_nn_index1)

  appendLF <- FALSE
  message('Models saved:', appendLF=appendLF)
  message('  preprocess_cds', appendLF=appendLF)
  message('Files written:', appendLF=appendLF)
  message('  R object information:', appendLF=appendLF)
  message('    models: ', file_name_rds, appendLF=appendLF)
  message('  Nearest neighbor indexes:', appendLF=appendLF)
  message('    preprocess_cds: ', file_name_nn_index1, appendLF=appendLF)

  saveRDS(object, file_name_rds)
}


#' Load preprocess model.
#'
#' Load into an existing cell_data_set a preprocess model consisting
#' of a dimensionality reduction model and nearest neighbor nearest neighbor index
#' objects, which are used to transform new data sets.
#'
#' @param cds cell_data_set into which the model is to be loaded.
#' @param file_name_root. A string with the root of names given to the two
#'   files to which the model was saved. The files are named
#'   <file_name_root>.ppc_rds and <file_name_root>.ppc_nn_index.
#'
#' @return a cell_data_set
#' @export
load_preprocess_model_retired <- function(cds, file_name_root = NULL) {
  model_file_version <- '1.0.0'
  file_name_rds <- paste0(file_name_root, '.ppc_rds')
  file_name_nn_index1 <- paste0(file_name_root, '.ppc_nn_index')

  object <- readRDS(file_name_rds)
  if( is.null(object[['model_file']]) || object[['model_file']] != 'preprocess') {
    stop('File \'', file_name_rds, '\' was not made using the save_preprocess_model() function.')
  }

  method <- object[['method']]
  metric_index <- object[['bundle']][['nn_index']][['annoy_metric']]
  ndim_index <- object[['bundle']][['nn_index']][['annoy_ndim']]

  md5sum_nn_index <- tools::md5sum(file_name_nn_index1)
  if(md5sum_nn_index != object[['md5sum_nn_index']]) {
    stop('The annoy index file, \'', file_name_nn_index1, '\', differs from the file made using the save_preprocess_model() function.')
  }

  object[['bundle']][['nn_index']][['annoy_index']] <- load_annoy_index(object[['bundle']][['nn_index']][['annoy_index']], file_name_nn_index1, metric_index, ndim_index)

  if(length(object[['comment']]) > 0 ) {
    message('Comment: ', object[['comment']], appendLF=appendLF)
  }

  cds@preprocess_aux[[method]] <- object[['bundle']]

  cds
}


#' Save align_cds model.
#'
#' Save an align_cds model consisting of dimensionality reduction
#' models and align_cds nearest neighbor index objects, which are
#' used to transform new data sets. Additionally, it saves the
#' preprocess_cds nearest neighbor index, if preprocess_cds() was run with
#' build_nn_index=TRUE.
#'
#' @param cds cell_data_set with an existing model, which was created using
#'   cds <- preprocess_cds(cds, build_nn_index=TRUE).
#' @param file_name_root. A string with the root of names given to the
#'   files to which the model will be saved. The files are named
#'   <file_name_root>.aln_rds and <file_name_root>.aln_nn_index<n>. Note: the
#'   <file_name_root>.aln_rds file contains only the information needed to
#'   load the align_cds model into the cds object given in
#'   save_align_cds_model(cds,...).
#' @param comment An optional string that describes the model.
#'
#' @return None
#' @export
save_align_cds_model_retired <- function(cds, method = c('Aligned'), file_name_root = NULL, comment = NULL) {
  model_file_version <- '1.0.0'
  if(is.null(cds@preprocess_aux[['Aligned']])) {
    stop('align_cds was not run on this cds.')
  } else if(is.null(cds@preprocess_aux[['Aligned']][['nn_index']])) {
    stop('There is no nearest neighbor index -- run align_cds() with build_nn_index=TRUE')
  }
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "method must be one of 'Aligned'")
  method <- match.arg(method)

  file_name_rds <- paste0(file_name_root, '.aln_rds')
  # preprocess nearest neighbor index (PCA or LSI).
  file_name_nn_index1 <- paste0(file_name_root, '.aln_nn_index1')
  # align_cds nearest neighbor index.
  file_name_nn_index2 <- paste0(file_name_root, '.aln_nn_index2')

  # Was preprocess_method 'PCA' or 'LSI', as specified in align_cds(..., preprocess_method=xx,...)?
  preprocess_method <- cds@preprocess_aux[[method]][['model']][['preprocess_method']]

  # Is there a preprocess_cds() nearest neighbor index?
  exists_index1 <- !(is.null(cds@preprocess_aux[[preprocess_method]][['nn_index']]))

  object <- list()
  object[['model_file_version']] <- model_file_version
  object[['model_file']] <- 'align_cds'
  object[['preprocess_method']] <- preprocess_method
  object[['method']] <- method
  object[['bundle']] <- list(preprocess_aux=cds@preprocess_aux[[preprocess_method]], align_aux=cds@preprocess_aux[[method]])
  object[['comment']] <- comment

  if(exists_index1) {
    save_annoy_index(cds@preprocess_aux[[preprocess_method]][['nn_index']][['annoy_index']], file_name_nn_index1)
    object[['md5sum_nn_index1']] <- tools::md5sum(file_name_nn_index1)
  }

  save_annoy_index(cds@preprocess_aux[[method]][['nn_index']][['annoy_index']], file_name_nn_index2)
  object[['md5sum_nn_index2']] <- tools::md5sum(file_name_nn_index2)

  appendLF <- FALSE
  message('Models saved:', appendLF=appendLF)
  message('  preprocess_cds\n  align_cds', appendLF=appendLF)
  message('Files written:', appendLF=appendLF)
  message('  R object information:', appendLF=appendLF)
  message('    models: ', file_name_rds, appendLF=appendLF)
  message('  Nearest neighbor indexes:', appendLF=appendLF)
  message('    preprocess_cds: ', ifelse(exists_index1, file_name_nn_index1, 'No preprocess_cds index'), appendLF=appendLF)
  message('    align_cds: ', file_name_nn_index2, appendLF=appendLF)

  saveRDS(object, file_name_rds)
}


#' Load align_cds model.
#' 
#' Load into an existing cell_data_set an align_cds model consisting
#' of dimensionality reduction models and nearest neighbor index
#' objects, which are used to transform new data sets.
#' 
#' @param cds cell_data_set into which the model is to be loaded.
#' @param file_name_root. A string with the root of names given to the two
#'   files to which the model was saved. The files are named
#'   <file_name_root>.aln_rds
#'   and <file_name_root>.aln_nn_index<n>.
#'   
#' @return a cell_data_set
#' @export
load_align_cds_model_retired <- function(cds, file_name_root = NULL) {
  model_file_version <- '1.0.0'
  file_name_rds <- paste0(file_name_root, '.aln_rds')
  file_name_nn_index1 <- paste0(file_name_root, '.aln_nn_index1')
  file_name_nn_index2 <- paste0(file_name_root, '.aln_nn_index2')

  object <- readRDS(file_name_rds)
  if( is.null(object[['model_file']]) || object[['model_file']] != 'align_cds') {
    stop('File \'', file_name_rds, '\' was not made using the save_align_cds_model() function.')
  } 
    
  preprocess_method <- object[['preprocess_method']]
  method <- object[['method']]
    
  preprocess_aux <- object[['bundle']][['preprocess_aux']]
  align_aux <- object[['bundle']][['align_aux']]
  
  exists_index1 <- !is.null(object[['md5sum_nn_index1']])
  
  if(exists_index1) {
    metric_index1 <- preprocess_aux[['nn_index']][['annoy_metric']]
    ndim_index1 <- preprocess_aux[['nn_index']][['annoy_ndim']]
  }
  metric_index2 <- align_aux[['nn_index']][['annoy_metric']]
  ndim_index2 <- align_aux[['nn_index']][['annoy_ndim']]
  
  if(exists_index1) {
    md5sum_nn_index1 <- tools::md5sum(file_name_nn_index1)
    if(md5sum_nn_index1 != object[['md5sum_nn_index1']]) {
      stop('The annoy index file, \'', file_name_nn_index1, '\', differs from the file made using the save_align_cds_model() function.')
    }
  }
  
  md5sum_nn_index2 <- tools::md5sum(file_name_nn_index2)
  if(md5sum_nn_index2 != object[['md5sum_nn_index2']]) {
    stop('The annoy index file, \'', file_name_nn_index2, '\', differs from the file made using the save_align_cds_model() function.')
  }

  if(exists_index1) {
    preprocess_aux[['nn_index']][['annoy_index']] <- load_annoy_index(preprocess_aux[['nn_index']][['annoy_index']], file_name_nn_index1, metric_index1, ndim_index1)
  } 
  
  align_aux[['nn_index']][['annoy_index']] <- load_annoy_index(align_aux[['nn_index']][['annoy_index']], file_name_nn_index2, metric_index2, ndim_index2)
  
  cds@preprocess_aux[[preprocess_method]] <- preprocess_aux
  cds@preprocess_aux[[method]] <- align_aux
  
  if(length(object[['comment']]) > 0 ) {
    message('Comment: ', object[['comment']], appendLF=appendLF)
  }
  
  cds
} 
  
  
#' Save reduce_dimension model.
#'
#' Save a reduce_dimension model consisting of dimensionality
#' reduction models and reduce_dimension nearest neighbor
#' index objects, which are used to transform new data
#' sets. Additionally, it saves preprocess and align_cds nearest neighbor
#' indexes, when preprocess_cds() and align_cds() were run with
#' build_nn_index=TRUE.
#'
#' @param cds cell_data_set with an existing model, which was created using
#'   cds <- reduce_dimension(cds, build_nn_index=TRUE).
#' @param method A string indicating the method used to build the model
#'   that you want saved. This must be "UMAP". Default = "UMAP".
#' @param file_name_root. A string with the root of names given to the
#'   files to which the model will be saved. The files are named
#'   <file_name_root>.rdd_rds, <file_name_root>.rdd_umap_nn_index, and
#'   <file_name_root>._rdd_nn_index<n>. Note: the <file_name_root>.rdd_rds file
#'   contains only the information needed to load the reduce_dimension
#'   model into the cds object given in load_reduce_dimension_model(cds,...).
#' @param comment An optional string that describes the model.
#'
#' @return None
#' @export
save_reduce_dimension_model_retired <- function(cds, method = c('UMAP'), file_name_root = NULL, comment = NULL) {
  model_file_version <- '1.0.0'
  if(is.null(cds@reduce_dim_aux[['UMAP']])) {
    stop('Reduce_dimension was not run on this cds.')
  } else if(is.null(cds@reduce_dim_aux[['UMAP']][['nn_index']])) {
    stop('There is no nearest neighbor index -- run reduce_dimension() with build_nn_index=TRUE')
  }
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "method must be 'UMAP'")
  method <- match.arg(method)

  # Set up filenames.
  file_name_rds <- paste0(file_name_root, '.rdd_rds')
  # uwot::umap annoy index.
  file_name_umap_index <- paste0(file_name_root, '.rdd_umap_nn_index')
  # preprocess nearest neighbor index (PCA or LSI).
  file_name_nn_index1 <- paste0(file_name_root, '.rdd_nn_index1')
  # align_cds Aligned nearest neighbor index.
  file_name_nn_index2 <- paste0(file_name_root, '.rdd_nn_index2')
  # reduce_dimension nearest neighbor index.
  file_name_nn_index3 <- paste0(file_name_root, '.rdd_nn_index3')

  # Get preprocess information.
  umap_preprocess_method <- cds@reduce_dim_aux[['UMAP']][['model']][['umap_preprocess_method']]
  # preprocess_method2 is given in reduce_dimension(...,preprocess_method=xx,...)
  # preprocess_method1 is given in align_cds(...,preprocess_method=xx,...) if align_cds()
  #                    was used; otherwise, it is given in preprocess_cds(...,method=xx,...)
  if(umap_preprocess_method != 'Aligned') {
    preprocess_method1 <- umap_preprocess_method
    preprocess_method2 <- NULL
  } else {
    preprocess_method1 <- cds@preprocess_aux[['Aligned']][['model']][['preprocess_method']]
    preprocess_method2 <- 'Aligned'
  }

  # Annoy index flags.
  exists_index1 <- !(is.null(cds@preprocess_aux[[preprocess_method1]][['nn_index']]))
  exists_index2 <- !(is.null(preprocess_method2)) && !(is.null(cds@preprocess_aux[[preprocess_method2]][['nn_index']]))

  # Initialize output object.
  object <- list()
  object[['model_file_version']] <- model_file_version
  object[['model_file']] <- 'reduce_dimension'
  object[['preprocess_method1']] <- preprocess_method1
  object[['preprocess_method2']] <- preprocess_method2
  object[['method']] <- method

  # Write annoy indexes to files and set up R objects bundle.
  if(exists_index1) {
    save_annoy_index(cds@preprocess_aux[[preprocess_method1]][['nn_index']][['annoy_index']], file_name_nn_index1)
    object[['md5sum_nn_index1']] <- tools::md5sum(file_name_nn_index1)
  }

  if(is.null(preprocess_method2)) {
    object[['bundle']] <- list(preprocess_aux=cds@preprocess_aux[[preprocess_method1]],
                               reduce_dim_aux=cds@reduce_dim_aux[[method]])
  } else {
    if(exists_index2) {
      save_annoy_index(cds@preprocess_aux[[preprocess_method2]][['nn_index']][['annoy_index']], file_name_nn_index2)
      object[['md5sum_nn_index2']] <- tools::md5sum(file_name_nn_index2)
    }
    object[['bundle']] <- list(preprocess_aux=cds@preprocess_aux[[preprocess_method1]],
                               align_aux=cds@preprocess_aux[['Aligned']],
                               reduce_dim_aux=cds@reduce_dim_aux[[method]])
  }

  save_annoy_index(cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']], file_name_nn_index3)
  object[['md5sum_nn_index3']] <- tools::md5sum(file_name_nn_index3)

  object[['md5sum_umap_index']] <- save_umap_nn_indexes(cds@reduce_dim_aux[[method]][['model']][['umap_model']], file_name_umap_index)

  # Final set up and information message.
  object[['comment']] <- comment

  appendLF <- FALSE
  message('Models saved:', appendLF=appendLF)
  if(is.null(preprocess_method2)) {
    message('  preprocess_cds\n  reduce_dimension', appendLF=appendLF)
  } else {
    message('  preprocess_cds\n  align_cds\n  reduce_dimension', appendLF=appendLF)
  }
  message('Files written:', appendLF=appendLF)
  message('  R object information:', appendLF=appendLF)
  message('    models: ', file_name_rds, appendLF=appendLF)
  message('  Nearest neighbor indexes:', appendLF=appendLF)
  message('    preprocess_cds: ', ifelse(exists_index1, file_name_nn_index1, 'No preprocess_cds index'), appendLF=appendLF)
  message('    align_cds: ', ifelse(exists_index2, file_name_nn_index2, 'No align_cds index'), appendLF=appendLF)
  message('    reduce_dimension: ', file_name_nn_index3, appendLF=appendLF)
  message('    UMAP: ', file_name_umap_index, appendLF=appendLF)

  saveRDS(object, file_name_rds)
}


#' Load reduce_dimension model.
#'
#' Load into an existing cell_data_set a reduce_dimension model consisting
#' of dimensionality reduction models and nearest neighbor index objects,
#' which are used to transform new data sets. 
#' 
#' @param cds cell_data_set into which the model is to be loaded.
#' @param file_name_root. A string with the root of names given to the two
#'   files to which the model was saved. The files are named
#'   <file_name_root>.rdd_rds, <file_name_root>.rdd_umap_nn_index, and
#'   <file_name_root>._rdd_nn_index<n>.
#' 
#' @return A cell_data_set.
#' 
#' @export 
load_reduce_dimension_model_retired <- function(cds, file_name_root = NULL) {
  model_file_version <- '1.0.0'
  file_name_rds <- paste0(file_name_root, '.rdd_rds') 
  file_name_umap_index <- paste0(file_name_root, '.rdd_umap_nn_index')
  file_name_nn_index1 <- paste0(file_name_root, '.rdd_nn_index1')
  file_name_nn_index2 <- paste0(file_name_root, '.rdd_nn_index2')
  file_name_nn_index3 <- paste0(file_name_root, '.rdd_nn_index3')

  # Read R object from RDS file.
  object <- readRDS(file_name_rds)
  if( is.null(object[['model_file']]) || object[['model_file']] != 'reduce_dimension') {
    stop('File \'', file_name_rds, '\' was not made using the save_reduce_dimension_model() function.')
  } 
    
  # Get preprocess method information.
  preprocess_method1 <- object[['preprocess_method1']]
  preprocess_method2 <- object[['preprocess_method2']]
  method <- object[['method']]
    
  # Did the user run align_cds?
  preprocess_aux <- object[['bundle']][['preprocess_aux']]
  align_aux <- if(!is.null(preprocess_method2)) object[['bundle']][['align_aux']] else NULL
  reduce_dim_aux <- object[['bundle']][['reduce_dim_aux']]
  
  exists_index1 <- !is.null(object[['md5sum_nn_index1']])
  exists_index2 <- !is.null(object[['md5sum_nn_index2']])
  
  # Get annoy metric and ndim values.
  metric_umap_index <- reduce_dim_aux[['model']][['umap_model']][['nn_index']][['metric']]
  ndim_umap_index <- reduce_dim_aux[['model']][['umap_model']][['metric']][[metric_umap_index]][['ndim']]
  
  if(exists_index1) {
    metric_nn_index1 <- preprocess_aux[['nn_index']][['annoy_metric']]
    ndim_nn_index1 <- preprocess_aux[['nn_index']][['annoy_ndim']]
  } 
  
  if(exists_index2) {
    metric_nn_index2 <- align_aux[['nn_index']][['annoy_metric']]
    ndim_nn_index2 <- align_aux[['nn_index']][['annoy_ndim']]
  } 

  metric_nn_index3 <- reduce_dim_aux[['nn_index']][['annoy_metric']]
  ndim_nn_index3 <- reduce_dim_aux[['nn_index']][['annoy_ndim']]

  # Check the annoy index file md5sums.
  if(exists_index1) {
    md5sum_nn_index1 <- tools::md5sum(file_name_nn_index1)
    if(md5sum_nn_index1 != object[['md5sum_nn_index1']] ) { 
      stop('The annoy index file, \'', file_name_nn_index1, '\', differs from the file made using the save_reduce_dimension_model() function.')
    }
  }

  if(exists_index2) { 
    md5sum_nn_index2 <- tools::md5sum(file_name_nn_index2)
    if(md5sum_nn_index2 != object[['md5sum_nn_index2']] ) {
      stop('The annoy index file, \'', file_name_nn_index2, '\', differs from the file made using the save_reduce_dimension_model() function.')
    }
  } 

  md5sum_nn_index3 <- tools::md5sum(file_name_nn_index3)
  if(md5sum_nn_index3 != object[['md5sum_nn_index3']] ) {
    stop('The annoy index file, \'', file_name_nn_index3, '\', differs from the file made using the save_reduce_dimension_model() function.')
  }

  # Load annoy indexes.
  reduce_dim_aux[['model']][['umap_model']] <- load_umap_nn_indexes(reduce_dim_aux[['model']][['umap_model']], file_name_umap_index, object[['md5sum_umap_index']])

  if(exists_index1) {
    preprocess_aux[['nn_index']][['annoy_index']] <- load_annoy_index(preprocess_aux[['nn_index']][['annoy_index']], file_name_nn_index1, metric_nn_index1, ndim_nn_index1)
  }

  if(exists_index2) {
    align_aux[['nn_index']][['annoy_index']] <- load_annoy_index(align_aux[['nn_index']][['annoy_index']], file_name_nn_index2, metric_nn_index2, ndim_nn_index2)
  }

  reduce_dim_aux[['nn_index']][['annoy_index']] <- load_annoy_index(reduce_dim_aux[['nn_index']][['annoy_index']], file_name_nn_index3, metric_nn_index3, ndim_nn_index3)

  # Final assignments and comment message, if it exists.
  cds@preprocess_aux[[preprocess_method1]] <- preprocess_aux
  if(!is.null(preprocess_method2)) {
    cds@preprocess_aux[['Aligned']] <- align_aux
  }
  cds@reduce_dim_aux[[method]] <- reduce_dim_aux

  if(length(object[['comment']]) > 0 ) {
    message('Comment: ', object[['comment']], appendLF=appendLF)
  }

  cds
}


