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
    stop( 'Bad status reading ', annotation_type, ' file \'', anno_path, '\'\n  ',
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
  nn_index
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
save_preprocess_model <- function(cds, method = c('PCA', 'LSI'), file_name_root = NULL, comment = NULL) {
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
load_preprocess_model <- function(cds, file_name_root = NULL) {
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
    message('Comment: ', object[['comment']])
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
save_align_cds_model <- function(cds, method = c('Aligned'), file_name_root = NULL, comment = NULL) {
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
  message('    align_cds: ', file_name_nn_index2)

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
load_align_cds_model <- function(cds, file_name_root = NULL) {
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
    message('Comment: ', object[['comment']])
  }

  cds
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
save_reduce_dimension_model <- function(cds, method = c('UMAP'), file_name_root = NULL, comment = NULL) {
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
#' @return a cell_data_set
#' @export
load_reduce_dimension_model <- function(cds, file_name_root = NULL) {
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
    message('Comment: ', object[['comment']])
  }

  cds
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
# Make model string from file_index[['file']] entry.
#
file_index_to_model_string <- function(file_index, i) {
  if(file_index[['cds_object']][[i]] == 'preprocess_aux') {
    model <- paste0('preprocess_cds ', file_index[['method']][[i]])
  } else
  if(file_index[['cds_object']][[i]] == 'reduce_dim_aux') {
    model <- paste0('reduce_dimension ', file_index[['method']][[i]])
  } else {
    stop('Unrecognized cds_object value \'', file_index[['cds_object']][[i]], '\'')
  }
  model
}


#
# Report files saved.
#
report_files_saved <- function(file_index) {
  appendLF <- FALSE
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
    message('  ', file_name, '  (', method, '  ', file_type, '  from  ', process, ')')
  }
}


#
#' Save the set of cell_data_set transform models.
#'
#' Save the transform models to a specified directory
#' by writing the R objects to RDS files and
#' the Annoy nearest neighbor indexes to individual index files. 
#' save_transform_models saves transforms made by the
#' preprocess_cds, align_cds, and reduce_dimension functions.
#' From the preprocess_cds, align_cds, and reduce_dimension
#' transforms, save_transform_models saves the objects required
#' to transform new count matrices into their reduced dimension
#' spaces. If the preprocess_cds, align_cds, or reduce_dimension
#' functions are run with the build_nn_index=TRUE parameter,
#' those nearest neighbor indexes are saved as well. The model
#' transforms are made by running preprocess_cds, align_cds,
#' and reduce_dimension on a primary data set and saving the
#' the transform models. In order to transform a new count
#' matrix, load it into a new cell_data_set, load the saved
#' transform models, and apply the models to the new count
#' matrix without running preprocess_cds, align_cds, or
#' reduce_dimension.
#'
#' @param cds A cell_data_set with existing models.
#' @param directory_path A string giving the name of the directory
#'   in which to write the model files.
#' @param comment A string with optional notes that is saved with 
#'   the objects.
#' @param verbose A boolean determining whether to print information
#'   about the saved files.
#' @return
#'
#' @export
save_transform_models <- function( cds, directory_path, comment="", verbose=TRUE) {
  appendLF <- FALSE
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
#' the cell_data_set.
#'
#' @param cds A cell_data_set to be transformed using the models.
#' @param directory_path A string giving the name of the directory
#'   from which to read the model files.
#' @return
#'
#' @export
load_transform_models <- function(cds, directory_path) {
  appendLF <- FALSE
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
    message('File comment: ', file_index[['comment']])
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
            message('Error reading file \'', file_path, '\'')
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
           message('Error reading file \'', file_path, '\'')
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
            message('Error reading file \'', file_path, '\'')
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
            message('Error reading file \'', file_path, '\'')
            return(NULL)
          })
      } else
      if(method == 'UMAP' && file_format == 'umap_annoy_index') {
        cds@reduce_dim_aux[[method]][['model']][['umap_model']] <- tryCatch(
          {
            load_umap_nn_indexes(cds@reduce_dim_aux[[method]][['model']][['umap_model']], file_path, md5sum)
          },
          error = function(cnd) {
            message('Error reading file \'', file_path, '\'')
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
#' @return NA
#'
# # @export
save_monocle_objects_old <- function(cds, directory_path) {
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
#' @return cell_data_set
#'
# # @export
load_monocle_objects_old <- function(directory_path) {
  appendLF <- FALSE

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
#' HDF5 object, all assays must be.
#'
#' @param cds A cell_data_set to save.
#' @param directory_path A string giving the name of the directory
#'   in which to write the object files.
#' @param hdf5_assays A boolean determining whether the
#'   non-HDF5Array assays objects are saved as HDF5 files. At this
#'   time HDF5Array assay objects are stored as HDF5Assay files
#'   regardless of the hdf5_assay parameter value.
#' @param comment A string with optional notes that is saved with
#'   the objects.
#' @param verbose A boolean determining whether to print information
#'   about the saved files.
#' @return 
#'
#' @export
save_monocle_objects <- function(cds, directory_path, hdf5_assays=FALSE, comment="", verbose=TRUE) {
  appendLF <- FALSE
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
#' save_monocle_objects.
#'
#' @param directory_path A string giving the name of the directory
#'   from which to read the saved cell_data_set files.
#' @return
#'
#' @export
load_monocle_objects <- function(directory_path) {
  appendLF <- FALSE
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
    message('File comment: ', file_index[['comment']])
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
            message('Error reading file \'', file_path, '\'')
            return(NULL)
          })
      } else
      if(file_format == 'hdf5') {
        cds <- tryCatch(
          {
            HDF5Array::loadHDF5SummarizedExperiment(file_path)
          },
          error = function(cnd) {
            message('Error reading file \'', file_path, '\'')
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
           message('Error reading file \'', file_path, '\'')
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
            message('Error reading file \'', file_path, '\'')
            return(NULL)
          })
      } else
      if(method == 'UMAP' && file_format == 'umap_annoy_index') {
        cds@reduce_dim_aux[[method]][['model']][['umap_model']] <- tryCatch(
          {
            load_umap_nn_indexes(cds@reduce_dim_aux[[method]][['model']][['umap_model']], file_path, md5sum)
          },
          error = function(cnd) {
            message('Error reading file \'', file_path, '\'')
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

