#' Build a cell_data_set from the data stored in inst/extdata directory.
#' @export
load_a549 <- function(){
  small_a549_colData_df <- readRDS(system.file("extdata",
          "small_a549_dex_pdata.rda",
          package = "monocle3"))
  small_a549_rowData_df <- readRDS(system.file("extdata",
          "small_a549_dex_fdata.rda",
          package = "monocle3"))
  small_a549_exprs <- readRDS(system.file("extdata",
          "small_a549_dex_exprs.rda",
          package = "monocle3"))
  small_a549_exprs <- small_a549_exprs[,row.names(small_a549_colData_df)]
  
  cds <- new_cell_data_set(expression_data = small_a549_exprs,
      cell_metadata = small_a549_colData_df,
      gene_metadata = small_a549_rowData_df)
  cds
}


#' Build a cell_data_set from the data stored in inst/extdata directory.
#' @keywords internal
load_worm_embryo <- function(){
  expression_matrix <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_expression.rds"))
  cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_colData.rds"))
  gene_annotation <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_rowData.rds"))
  gene_annotation$use_for_ordering <- NULL
  
  cds <- new_cell_data_set(expression_matrix,
      cell_metadata = cell_metadata,
      gene_metadata = gene_annotation)
  cds <- estimate_size_factors(cds)
  cds
}


#' Test if a file has a Matrix Market header.
#' @param matpath Path to test file.
#' @return TRUE if matpath file has Matrix Market header.
#' @noRd
is_matrix_market_file <- function( matpath )
{
  first_line <- read.table( matpath, nrows=1 )
  grepl( "%%MatrixMarket", first_line$V1 )
}


#' Load matrix dimension metadata from file.
#' @param anno_path Path to matdim annotation file.
#' @return list consisting of matdim_names, which label matrix dimension,
#' and metadata, if present in the file, which are
#' additional dimension metadata.
#' @noRd
load_annotations_data <- function( anno_path, metadata_column_names=NULL, header=FALSE, sep="", quote="\"'", annotation_type=NULL )
{
  assertthat::assert_that( ! is.null( annotation_type ) )
  tryCatch(
      {
        annotations <- read.table( anno_path, header=header, sep=sep, quote=quote, stringsAsFactors=FALSE )
      }, error = function( emsg )
      {
        stop( 'load_mm_data: bad status reading ', annotation_type, ' file \'', anno_path, '\'\n  ',
            emsg, '\n',
            '  note: possible problems include the wrong filename, a missing file,\n',
            '  and incorrect file format parameters, for example \'header\', \'sep\', and \'quote\'' )
      })
  
  metadata = NULL
  if( .row_names_info( annotations ) < 0 )
  {
    names <- annotations[,1]
    if( ncol( annotations ) > 1 )
      metadata <- annotations[,-1,drop=FALSE]
  }
  else
  {
    names <- rownames( annotations )
    metadata <- annotations
  }
  
  if( ! ( is.null( metadata_column_names ) || is.null( metadata ) ) )
  {
    assertthat::assert_that( length( metadata_column_names ) == ncol( metadata ),
        msg=paste0( annotation_type,
            ' metadata column name count (',
            length( metadata_column_names ),
            ') != ',
            annotation_type,
            'annotation column count (',
            ncol( metadata ),
            ')' ) )
    colnames( metadata ) <- metadata_column_names
  }
  
  if( ! is.null( metadata ) )
    rownames(metadata)<-names
  
  list( names=names, metadata=metadata )
}


#' Load data from matrix market format files.
#'
#' @param mat_path Path to the Matrix Market .mtx matrix file. The
#' values are read and stored  as a sparse matrix with nrows and ncols,
#' as inferred from the file. Required.
#' @param feature_anno_path Path to a feature annotation file. The
#' feature_anno_path file must have nrows lines and at least one column.
#' The values in the first column label the matrix rows and each must be
#' distinct in the column. Values in additional columns are stored in
#' the cell_data_set 'gene' metadata. For gene features, we urge use of
#' official gene IDs for labels, such as Ensembl or Wormbase IDs. In this
#' case, the second column has typically a 'short' gene name. Additional
#' information such as gene_biotype may be stored in additional columns
#' starting with column 3. Required.
#' @param cell_anno_path Path to a cell annotation file. The cell_anno_path
#' file must have ncols lines and at least one column. The values in the
#' first column label the matrix columns and each must be distinct in the
#' column. Values in additional columns are stored in the cell_data_set
#' cells metadata. Required.
#' @param header Logical set to TRUE if both feature_anno_path and
#' cell_anno_path files have column headers, or set to FALSE if both
#' files do not have column headers (only these cases are supported).
#' The files may have either ncols or ncols-1 header fields. In both
#' cases, the first column is used as the matrix dimension names. The
#' default is FALSE.
#' @param feature_metadata_column_names A character vector of feature
#' metadata column names. The number of names must be one less than the
#' number of columns in the feature_anno_path file. These values
#' will replace those read from the feature_anno_path file header,
#' if present. The default is NULL.
#' @param cell_metadata_column_names A character vector of cell
#' metadata column names. The number of names must be one less than the
#' number of columns in the cell_anno_path file. These values will
#' replace those read from the cell_anno_path file header, if present.
#' The default is NULL.
#' @param quote A character string specifying the quoting characters
#' used in the feature_anno_path and cell_anno_path files. The default
#' is "\"'".
#' @param umi_cutoff UMI per cell cutoff. Columns (cells) with less
#' than umi_cutoff total counts are removed from the matrix. The
#' default is 100.
#' @param sep field separator character in the annotation files. If
#' sep = "", the separator is white space, that is, one or more spaces,
#' tabs, newlines, or carriage returns. The default is the tab
#' character for tab-separated-value files.
#'
#' @return cds object
#'
#' @section Comments:
#' * load_mm_data estimates size factors.
#'
#' @examples
#' \dontrun{
#' library( monocle3 )
#' pmat<-system.file("extdata", "matrix.mtx.gz", package = "monocle3")
#' prow<-system.file("extdata", "features_c3h0.txt", package = "monocle3")
#' pcol<-system.file("extdata", "barcodes_c2h0.txt", package = "monocle3")
#' cds <- load_mm_data( pmat, prow, pcol,
#'                      feature_metadata_column_names = c('gene_short_name',
#'                      'gene_biotype'), sep='' )
#'
#' In this example, the features_c3h0.txt file has three columns,
#' separated by spaces. The first column has official gene names, the
#' second has short gene names, and the third has gene biotypes.
#' }
#'
#' @export
#'
load_mm_data <- function( mat_path,
    feature_anno_path,
    cell_anno_path,
    header = FALSE,
    feature_metadata_column_names = NULL,
    cell_metadata_column_names = NULL,
    umi_cutoff = 100,
    quote="\"'",
    sep="\t") {
  assertthat::assert_that(assertthat::is.readable(mat_path), msg='unable to read matrix file')
  assertthat::assert_that(assertthat::is.readable(feature_anno_path), msg='unable to read feature annotation file')
  assertthat::assert_that(assertthat::is.readable(cell_anno_path), msg='unable to read cell annotation file')
  assertthat::assert_that(is.numeric(umi_cutoff))
  
  feature_annotations <- load_annotations_data( feature_anno_path, feature_metadata_column_names, header, sep, quote=quote, annotation_type='features' )
  cell_annotations <- load_annotations_data( cell_anno_path, cell_metadata_column_names, header, sep, quote=quote, annotation_type='cells' )
  
  assertthat::assert_that( ! any( duplicated( feature_annotations$names ) ), msg='duplicate feature names in feature annotation file' )
  assertthat::assert_that( ! any( duplicated( cell_annotations$names ) ), msg='duplicate cell names in cell annotation file' )
  
  mat <- Matrix::readMM( mat_path )
  
  assertthat::assert_that( length( feature_annotations$names ) == nrow( mat ),
      msg=paste0( 'feature name count (',
          length( feature_annotations$names ),
          ') != matrix row count (',
          nrow( mat ),
          ')' ) )
  assertthat::assert_that( length( cell_annotations$names ) == ncol( mat ),
      msg=paste0( 'cell name count (',
          length( cell_annotations$names ),
          ') != matrix column count (',
          ncol( mat ),
          ')' ) )
  
  rownames( mat ) <- feature_annotations$names
  colnames( mat ) <- cell_annotations$names
  
  cds <- new_cell_data_set( mat,
      cell_metadata = cell_annotations$metadata,
      gene_metadata = feature_annotations$metadata )
  
  colData(cds)$n.umi <- Matrix::colSums(exprs(cds))
  cds <- cds[,colData(cds)$n.umi >= umi_cutoff]
  cds <- estimate_size_factors(cds)

  cds <- initialize_counts_metadata(cds)
  matrix_id <- get_unique_id()
  cds <- set_counts_identity(cds, mat_path, matrix_id)
 
  return( cds )
}


#' Load data from matrix market format
#'
#' @param mat_path Path to the .mtx matrix market file.
#' @param gene_anno_path Path to gene annotation file.
#' @param cell_anno_path Path to cell annotation file.
#' @param umi_cutoff UMI per cell cutoff, default is 100.
#'
#' @return cds object
#' @export
#'
load_mtx_data <- function( mat_path,
    gene_anno_path,
    cell_anno_path,
    umi_cutoff = 100) {
  assertthat::assert_that(assertthat::is.readable(mat_path))
  assertthat::assert_that(assertthat::is.readable(gene_anno_path))
  assertthat::assert_that(assertthat::is.readable(cell_anno_path))
  assertthat::assert_that(is.numeric(umi_cutoff))
  
  if( is_matrix_market_file( mat_path ) )
  {
    #
    # Read an feature annotation file with two tab-separated
    # columns where the second column has short gene names.
    # Read a cell annotation file with one column that has the
    # matrix row names.
    #
    cds <- load_mm_data( mat_path,
        gene_anno_path,
        cell_anno_path,
        feature_metadata_column_names=c('gene_short_name'),
        umi_cutoff=umi_cutoff,
        sep="\t" )
    return( cds )
  }
  
  df <- utils::read.table(mat_path, col.names = c("gene.idx", "cell.idx", "count"),
      colClasses = c("integer", "integer", "integer"))
  
  gene.annotations <- utils::read.table(gene_anno_path,
      col.names = c("id", "gene_short_name"),
      colClasses = c("character", "character"))
  
  cell.annotations <- utils::read.table(cell_anno_path, col.names = c("cell"),
      colClasses = c("character"))
  
  rownames(gene.annotations) <- gene.annotations$id
  rownames(cell.annotations) <- cell.annotations$cell
  
  # add a dummy cell to ensure that all genes are included in the matrix
  # even if a gene isn't expressed in any cell
  df <- rbind(df, data.frame(gene.idx = c(1, nrow(gene.annotations)),
          cell.idx = rep(nrow(cell.annotations) + 1, 2),
          count = c(1, 1)))
  
  mat <- Matrix::sparseMatrix(i = df$gene.idx, j = df$cell.idx, x = df$count)
  
  if(ncol(mat) == 1) {
    mat <- mat[,0, drop=FALSE]
  } else {
    mat <- mat[, 1:(ncol(mat)-1), drop=FALSE]
  }
  
  rownames(mat) <- gene.annotations$id
  colnames(mat) <- cell.annotations$cell
  
  cds <- new_cell_data_set(mat, cell_metadata = cell.annotations,
      gene_metadata = gene.annotations)
  colData(cds)$n.umi <- Matrix::colSums(exprs(cds))
  cds <- cds[,colData(cds)$n.umi >= umi_cutoff]
  cds <- estimate_size_factors(cds)
  return(cds)
}


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
  return(nn_index)
}


# Save umap annoy indexes to files and return md5sum
# value(s) as either a character string, in case of
# one metric, or a list, in case of more than one matric.
save_umap_nn_indexes <- function(umap_model, file_name) {
  if(is.null(umap_model[['metric']])) {
    stop('No UMAP model in this CDS.')
  }
  metrics <- names(umap_model[['metric']])
  n_metrics <- length(metrics)
  if(n_metrics == 1) {
    save_annoy_index(umap_model[['nn_index']], file_name)
    md5sum_umap_index <- tools::md5sum(file_name)
  } else {
    warning('save_umap_nn_indexes is untested with more than one umap metric')
    md5sum_vec <- character()
    for(i in 1:n_metrics) {
      file_name_expand <- paste0(file_name, i)
      save_annoy_index(umap_model[['nn_index']][[i]], file_name_expand)
      md5sum <- tools::md5sum(file_name_expand)
      append(md5sum_vec, md5sum)
    }
    md5sum_umap_index <- paste(md5sum_vec, collapse='_')
  }
  return(md5sum_umap_index)
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
    warning('load_umap_nn_indexes is untested with more than one umap metric')
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
  return(umap_model)
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
    if(cds_object == 'reduce_dim_aux') {
      if(method == 'Aligned') {
        process <- 'align_cds'
      } else
      if(method == 'PCA' || method == 'LSI') {
        process <- 'preprocess_cds'
      } else
      if(method == 'tSNE' || method == 'UMAP') {
        process <- 'reduce_dimension'
      } else {
        stop('Unrecognized preprocess method \'', method, '\'')
      }
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
  #   cds_object: reduce_dim_aux
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

  # Gather reduce_dimension methods and whether each has an annoy index.
  methods_reduce_dim <- list()
  for( method in names(cds@reduce_dim_aux)) {
    methods_reduce_dim[[method]] <- list()
    methods_reduce_dim[[method]][['rds_path']] <- paste0('rdd_', tolower(method), '_transform_model.rds')
    methods_reduce_dim[[method]][['umap_index_path']] <- paste0('rdd_', tolower(method), '_transform_model_umap.idx')
    methods_reduce_dim[[method]][['nn_index_path']] <- paste0('rdd_', tolower(method), '_transform_model_nn.idx')
    if(!is.null(cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']]))
      methods_reduce_dim[[method]][['has_nn_index']] <- TRUE
    else
      methods_reduce_dim[[method]][['has_nn_index']] <- FALSE
  }

  # Make directory if necessary.
  dir.create(path = directory_path, showWarnings=FALSE, recursive=TRUE, mode='0700')

  # Remove files, if they exist.
  for(method in names(methods_reduce_dim)) {
    if(file.exists(file.path(directory_path, methods_reduce_dim[[method]][['rds_path']])))
      file.remove(file.path(directory_path, methods_reduce_dim[[method]][['rds_path']]))
    if(file.exists(file.path(directory_path, methods_reduce_dim[[method]][['nn_index_path']])))
      file.remove(file.path(directory_path, methods_reduce_dim[[method]][['nn_index_path']]))
    if(method == 'UMAP') {
      if(file.exists(file.path(directory_path, methods_reduce_dim[[method]][['umap_index_path']])))
         file.remove(file.path(directory_path, methods_reduce_dim[[method]][['umap_index_path']]))
    }
  }

  # Save reduce_dimension annoy indexes.
  # Notes:
  #   o  save RDS files before the corresponding index files in
  #      order to enable loading.
  #
  for(method in names(methods_reduce_dim)) {
    tryCatch(
      {
        saveRDS(cds@reduce_dim_aux[[method]], file=file.path(directory_path, methods_reduce_dim[[method]][['rds_path']]))
      },
      error = function(cnd) {
                     message('Error writing file \'', file.path(directory_path, methods_reduce_dim[[method]][['rds_path']]), '\': ', cnd, appendLF=appendLF)
                     return(NULL)
      },
      finally = {
        md5sum <- tools::md5sum(file.path(directory_path, methods_reduce_dim[[method]][['rds_path']]))
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
          md5sum <- save_umap_nn_indexes(cds@reduce_dim_aux[[method]][['model']][['umap_model']], file.path(directory_path, methods_reduce_dim[[method]][['umap_index_path']]))
        },
        error = function(cnd) {
                       message('Error writing file \'', file.path(directory_path, methods_reduce_dim[[method]][['umap_index_path']]), '\': ', cnd, appendLF=appendLF)
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
          save_annoy_index(cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']], file.path(directory_path, methods_reduce_dim[[method]][['nn_index_path']]))
        },
        error = function(cnd) {
                       message('Error writing file \'', file.path(directory_path, methods_reduce_dim[[method]][['nn_index_path']]), '\': ', cnd, appendLF=appendLF)
                       return(NULL)
        },
        finally = {
          md5sum <- tools::md5sum(file.path(directory_path, methods_reduce_dim[[method]][['nn_index_path']]))
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
    file_path <- file.path(directory_path, file_index[['files']][['file_path']][[ifile]])
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
      cds <- set_model_identity_path(cds, method, directory_path) 
    } else {
      stop('Unrecognized cds_object value \'', cds_object, '\'')
    }
  }

  return(cds)
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
  #   cds_object: cds | reduce_dim_aux
  #   method: PCA | LSI | Aligned | tSNE | UMAP  ...
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
  rds_path <- 'cds_object.rds'
  hdf5_path <- 'hdf5_object'

  # Gather reduce_dimension method names for which indexes exist.
  methods_reduce_dim <- list()
  for(method in names(cds@reduce_dim_aux)) {
    methods_reduce_dim[[method]] <- list()
    if(method == 'UMAP') {
      methods_reduce_dim[[method]][['umap_index_path']] <- paste0('rdd_', tolower(method), '_transform_model_umap.idx')
    }
    methods_reduce_dim[[method]][['nn_index_path']] <- paste0('rdd_', tolower(method), '_transform_model_nn.idx')
    if(!is.null(cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']]))
      methods_reduce_dim[[method]][['has_nn_index']] <- TRUE
    else
      methods_reduce_dim[[method]][['has_nn_index']] <- FALSE
  }

  # Make directory if necessary.
  dir.create(path = directory_path, showWarnings=FALSE, recursive=TRUE, mode='0700')

  # Remove files, if they exist.
  for(method in names(methods_reduce_dim)) {
    if(method == 'UMAP') {
      if(file.exists(file.path(directory_path, methods_reduce_dim[[method]][['umap_index_path']])))
         file.remove(file.path(directory_path, methods_reduce_dim[[method]][['umap_index_path']]))
    }
    if(file.exists(file.path(directory_path, methods_reduce_dim[[method]][['nn_index_path']])))
       file.remove(file.path(directory_path, methods_reduce_dim[[method]][['nn_index_path']]))
  }

  #
  # Save cds object.
  # Notes:
  #   o  allow for HDF5Array assay objects.
  #
  if(!hdf5_assay_flag) {
    tryCatch(
      {
        saveRDS(cds, file.path(directory_path, rds_path))
      },
      error = function(cnd) {
                       message('Error writing file \'', file.path(directory_path, rds_path), '\': ', cnd, appendLF=appendLF)
                       return(NULL)
      },
      finally = {
        md5sum <- tools::md5sum(file.path(directory_path, rds_path))
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
        HDF5Array::saveHDF5SummarizedExperiment(cds, file.path(directory_path, hdf5_path), replace=TRUE)
      },
      error = function(cnd) {
                       message('Error writing file \'', file.path(directory_path, hdf5_path), '\': ', cnd, appendLF=appendLF)
                       return(NULL)
      },
      finally = {
        md5sum <- tools::md5sum(file.path(directory_path, hdf5_path, 'se.rds'))
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

  # Save reduce_dimension annoy indexes.
  # Notes:
  #   o  save RDS files before the corresponding index files in
  #      order to enable loading.
  #
  for(method in names(methods_reduce_dim)) {
    if(method == 'UMAP') {
      tryCatch(
        {
          md5sum <- save_umap_nn_indexes(cds@reduce_dim_aux[[method]][['model']][['umap_model']], file.path(directory_path, methods_reduce_dim[[method]][['umap_index_path']]))
        },
        error = function(cnd) {
                       message('Error writing file \'', file.path(directory_path, methods_reduce_dim[[method]][['umap_index_path']]), '\': ', cnd, appendLF=appendLF)
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
          save_annoy_index(cds@reduce_dim_aux[[method]][['nn_index']][['annoy_index']], file.path(directory_path, methods_reduce_dim[[method]][['nn_index_path']]))
        },
        error = function(cnd) {
                       message('Error writing file \'', file.path(directory_path, methods_reduce_dim[[method]][['nn_index_path']]), '\': ', cnd, appendLF=appendLF)
                       return(NULL)
        },
        finally = {
          md5sum <- tools::md5sum(file.path(directory_path, methods_reduce_dim[[method]][['nn_index_path']]))
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
    file_path <- file.path(directory_path, file_index[['files']][['file_path']][[ifile]])
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
      cds <- set_model_identity_path(cds, method, directory_path)
    } else {
      stop('Unrecognized cds_object value \'', cds_object, '\'')
    }
  }

  return(cds)
}


load_monocle_objects_1_0 <- function(file) {
  cds_tmp <- readRDS(file)
  cds_preprocess_aux <- cds_tmp[['preprocess_aux']]
  cds_reduce_dim_aux <- cds_tmp@reduce_dim_aux
  cds_tmp[['preprocess_aux']] <- NULL
  cds_tmp@reduce_dim_aux <- NULL
  reduced_dims <-list()
  if(!is.null(reducedDim(cds_tmp)[['PCA']])) {
    reduced_dims[['PCA']] <- TRUE
#need to fuss with gene_loadings -> svd_v and svd_sdev
#transfer prop_var_expl
  }
  if(!is.null(reducedDim(cds_tmp)[['LSI']])) {
    reduced_dims[['LSI']] <- TRUE
#need to fuss with gene_loadings -> svd_v and svd_sdev
  }
  if(!is.null(reducedDim(cds_tmp)[['Aligned']])) {
    reduced_dims[['Aligned']] <- TRUE
#if there is a beta, then used linear regression, if not used only MNN
  }
  if(!is.null(reducedDim(cds_tmp)[['tSNE']])) {
    reduced_dims[['tSNE']] <- TRUE
  }
  if(!is.null(reducedDim(cds_tmp)[['UMAP']])) {
    reduced_dims[['UMAP']] <- TRUE
  }
}
