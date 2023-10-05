#' Build a small cell_data_set.
#' @param matrix_control A list used to control how the counts matrix is stored
#'    in the CDS. By default, Monocle3 stores the counts matrix in-memory as a
#'    sparse matrix. Setting 'matrix_control=list(matrix_class="BPCells")',
#'    stores the matrix on-disk as a sparse matrix.
#' @return cds object
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'   }
#'
#' @export
load_a549 <- function(matrix_control=list()){

  if(!is.null(matrix_control[['matrix_class']]) && matrix_control[['matrix_class']] == 'BPCells') {
    matrix_control_default <- get_global_variable('matrix_control_bpcells_unrestricted')
  }
  else {
    matrix_control_default <- get_global_variable('matrix_control_csparsematrix_unrestricted')
  }
  matrix_control_res <- set_matrix_control(matrix_control=matrix_control, matrix_control_default=matrix_control_default, control_type='unrestricted')

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

  expression_matrix <- set_matrix_class(mat=small_a549_exprs, matrix_control=matrix_control_res)

  cds <- new_cell_data_set(expression_data = expression_matrix,
                           cell_metadata = small_a549_colData_df,
                           gene_metadata = small_a549_rowData_df)
  cds
}


#' Build a cell_data_set from C. elegans embryo data.
#' @param matrix_control A list used to control how the counts matrix is stored
#'    in the CDS. By default, Monocle3 stores the counts matrix in-memory as a
#'    sparse matrix. Setting 'matrix_control=list(matrix_class="BPCells")',
#'    stores the matrix on-disk as a sparse matrix.
#' @return cds object
#' @importFrom SingleCellExperiment counts
#' @export
load_worm_embryo <- function(matrix_control=list()) {

  if(!is.null(matrix_control[['matrix_class']]) && matrix_control[['matrix_class']] == 'BPCells') {
    matrix_control_default <- get_global_variable('matrix_control_bpcells_unrestricted')
  }
  else {
    matrix_control_default <- get_global_variable('matrix_control_csparsematrix_unrestricted')
  }
  matrix_control_res <- set_matrix_control(matrix_control=matrix_control, matrix_control_default=matrix_control_default, control_type='unrestricted')

  expression_matrix <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_expression.rds"))
  cell_metadata <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_colData.rds"))
  gene_annotation <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_rowData.rds"))
  gene_annotation$use_for_ordering <- NULL

  expression_matrix <- set_matrix_class(mat=expression_matrix, matrix_control=matrix_control_res)

  cds <- new_cell_data_set(expression_matrix,
      cell_metadata = cell_metadata,
      gene_metadata = gene_annotation)
  cds <- estimate_size_factors(cds)

  cds <- initialize_counts_metadata(cds)
  matrix_id <- get_unique_id(counts(cds))
  cds <- set_counts_identity(cds, 'URL: https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_expression.rds', matrix_id)

  cds
}


#' Build a cell_data_set from C. elegans L2 data.
#' @param matrix_control A list used to control how the counts matrix is stored
#'    in the CDS. By default, Monocle3 stores the counts matrix in-memory as a
#'    sparse matrix. Setting 'matrix_control=list(matrix_class="BPCells")',
#'    stores the matrix on-disk as a sparse matrix.
#' @return cds object
#' @importFrom SingleCellExperiment counts
#' @export
load_worm_l2 <- function(matrix_control=list()) {

  if(!is.null(matrix_control[['matrix_class']]) && matrix_control[['matrix_class']] == 'BPCells') {
    matrix_control_default <- get_global_variable('matrix_control_bpcells_unrestricted')
  }
  else {
    matrix_control_default <- get_global_variable('matrix_control_csparsematrix_unrestricted')
  }
  matrix_control_res <- set_matrix_control(matrix_control=matrix_control, matrix_control_default=matrix_control_default, control_type='unrestricted')

  expression_matrix <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_expression.rds"))
  cell_metadata <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_colData.rds"))
  gene_annotation <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_rowData.rds"))

  expression_matrix <- set_matrix_class(mat=expression_matrix, matrix_control=matrix_control_res)

  cds <- new_cell_data_set(expression_matrix,
      cell_metadata = cell_metadata,
      gene_metadata = gene_annotation)
  cds <- estimate_size_factors(cds)

  cds <- initialize_counts_metadata(cds)
  matrix_id <- get_unique_id(counts(cds))
  cds <- set_counts_identity(cds, 'URL: https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_expression.rds', matrix_id)

  cds
}

#' Test if a file has a Matrix Market header.
#' @param matpath Path to test file.
#' @return TRUE if matpath file has Matrix Market header.
#' @noRd
is_matrix_market_file <- function( matpath )
{
  first_line <- utils::read.table( matpath, nrows=1 )
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
        annotations <- utils::read.table( anno_path, header=header, sep=sep, quote=quote, stringsAsFactors=FALSE )
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
    rownames(metadata) <- names

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
#' @param verbose a logical value that determines whether or not the
#' function writes diagnostic information.
#' @param matrix_control an optional list of values that control how
#' matrices are stored in the cell_data_set assays slot. Typically,
#' matrices are stored in-memory as dgCMatrix class (compressed sparse
#' matrix) objects using matrix_class="dgCMatrix". This is the
#' default. A very large matrix can be stored in a file and accessed
#' by Monocle3 as if it were in-memory. For this, Monocle3 uses the
#' BPCells R package. Here the matrix_control list values are set to
#' matrix_class="BPCells" and matrix_mode="dir". Then the counts matrix
#' is stored in a directory, on-disk, which is created by Monocle3 in
#' the directory where you run Monocle3. This directory has a name
#' with the form "monocle.bpcells.*.tmp" where the asterisk is a
#' string of random characters that makes the name unique. Do not
#' remove this directory while Monocle3 is running! If you choose to
#' store the counts matrix as an on-disk BPCells object, you must use
#' the "save_monocle_objects" and "load_monocle_objects" functions
#' to save and restore the cell_data_set. Monocle3 tries to remove
#' the BPCells matrix directory when your R session ends; however,
#' sometimes a matrix directory may persist after the session ends.
#' In this case, the user must remove the directory after the
#' session ends. For additional information about the matrix_control
#' list, see the examples below and the set_matrix_control help.
#' Note that for the load_mm_data function the BPCells matrix_mode
#' is "dir", the matrix_type is "double", and the matrix_compress is
#' FALSE.
#' @return cds object
#'
#' @section Comments:
#' * load_mm_data estimates size factors.
#'
#' @examples
#'   \donttest{
#'     pmat<-system.file("extdata", "matrix.mtx.gz", package = "monocle3")
#'     prow<-system.file("extdata", "features_c3h0.txt", package = "monocle3")
#'     pcol<-system.file("extdata", "barcodes_c2h0.txt", package = "monocle3")
#'     cds <- load_mm_data( pmat, prow, pcol,
#'                          feature_metadata_column_names =
#'                          c('gene_short_name', 'gene_biotype'), sep='' )
#'
#'     # In this example, the features_c3h0.txt file has three columns,
#'     # separated by spaces. The first column has official gene names, the
#'     # second has short gene names, and the third has gene biotypes.
#'     #
#'     # For typical count matrices with a small to medium number of cells,
#'     # we suggest that you use the default matrix_control list by not
#'     # not setting the matrix_control parameter. In this case, the
#'     # counts matrix is stored in-memory as a sparse matrix in the
#'     # dgCMatrix format, as it has in the past. It is also possible to
#'     # set the matrix_control list explicitly to use this in-memory
#'     # dgCMatrix format by setting the matrix_control parameter to
#'     #
#'       load_mm_data(..., matrix_control=list(matrix_class='dgCMatrix'))
#'     #
#'     # For large matrices, we suggest that you try storing the count
#'     # matrix as a BPCells object on-disk by setting the matrix_control
#'     # parameter list as follows
#'     #
#'       load_mm_data(..., matrix_control=list(matrix_class='BPCells'))
#'     #
#'   }
#' 
#' @importFrom SingleCellExperiment counts
#' @export
load_mm_data <- function( mat_path,
    feature_anno_path,
    cell_anno_path,
    header = FALSE,
    feature_metadata_column_names = NULL,
    cell_metadata_column_names = NULL,
    umi_cutoff = 100,
    quote="\"'",
    sep="\t",
    verbose=FALSE,
    matrix_control=list()) {
  assertthat::assert_that(assertthat::is.readable(mat_path), msg='unable to read matrix file')
  assertthat::assert_that(assertthat::is.readable(feature_anno_path), msg='unable to read feature annotation file')
  assertthat::assert_that(assertthat::is.readable(cell_anno_path), msg='unable to read cell annotation file')
  assertthat::assert_that(is.numeric(umi_cutoff))

  if(!is.null(matrix_control[['matrix_class']]) && matrix_control[['matrix_class']] == 'BPCells') {
    matrix_control_default <- get_global_variable('matrix_control_bpcells_unrestricted')
  }
  else {
    matrix_control_default <- get_global_variable('matrix_control_csparsematrix_unrestricted')
  }

  matrix_control_res <- set_matrix_control(matrix_control=matrix_control, matrix_control_default=matrix_control_default, control_type='unrestricted')

  feature_annotations <- load_annotations_data( feature_anno_path, feature_metadata_column_names, header, sep, quote=quote, annotation_type='features' )
  cell_annotations <- load_annotations_data( cell_anno_path, cell_metadata_column_names, header, sep, quote=quote, annotation_type='cells' )

  assertthat::assert_that( ! any( duplicated( feature_annotations$names ) ), msg='duplicate feature names in feature annotation file' )
  assertthat::assert_that( ! any( duplicated( cell_annotations$names ) ), msg='duplicate cell names in cell annotation file' )

  if(matrix_control_res[['matrix_class']] != 'BPCells') {
    # Read MatrixMarket file and convert to dgCMatrix format.
    mat <- Matrix::readMM(mat_path)
    mat <- set_matrix_class(mat=mat, matrix_control=matrix_control_res)
  
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
  }
  else {
    toutdir <- tempfile(pattern=paste0('monocle.bpcells.',
                                       format(Sys.Date(), format='%Y%m%d'), '.'),
                        tmpdir=matrix_control_res[['matrix_path']],
                        fileext='.tmp')[[1]]
    tmpdir <- tempfile('monocle.import_mm.', '.', '.tmp')
    tmat <- BPCells::import_matrix_market(mtx_path=mat_path,
                                          outdir=toutdir,
                                          row_names=feature_annotations$names,
                                          col_names=cell_annotations$names,
                                          row_major=FALSE,
                                          tmpdir=tmpdir,
                                          load_bytes=4194304L,
                                          sort_bytes=1073741824L)
    unlink(tmpdir, recursive=TRUE)
    outdir_c <- tempfile(pattern=paste0('monocle.bpcells.',
                                        format(Sys.Date(), format='%Y%m%d'), '.'),
                         tmpdir=matrix_control_res[['matrix_path']],
                         fileext='.tmp')[[1]]
    mat <- BPCells::write_matrix_dir(BPCells::convert_matrix_type(tmat, 'double'), outdir_c, compress=FALSE, buffer_size=8192L, overwrite=FALSE)
    unlink(toutdir, recursive=TRUE)
    push_matrix_path(mat=mat)
  }

  cds <- new_cell_data_set(mat,
                           cell_metadata = cell_annotations$metadata,
                           gene_metadata = feature_annotations$metadata,
                           verbose = verbose)

  if(is(counts(cds), 'CsparseMatrix')) {
    colData(cds)$n.umi <- Matrix::colSums(counts(cds))
  }
  else
  if(is(counts(cds), 'IterableMatrix')) {
    colData(cds)$n.umi <- BPCells::colSums(counts(cds))
  }

  cds <- cds[,colData(cds)$n.umi >= umi_cutoff]
  cds <- estimate_size_factors(cds)

  cds <- initialize_counts_metadata(cds)
  matrix_id <- get_unique_id(counts(cds))
  cds <- set_counts_identity(cds, mat_path, matrix_id)

  if(verbose) {
    message('load_mm_data: matrix_info: out:')
    show_matrix_info(matrix_info=get_matrix_info(mat=mat), '  ')
  }

  return(cds)
}


#' Load data from matrix market format
#'
#' @param mat_path Path to the .mtx matrix market file.
#' @param gene_anno_path Path to gene annotation file.
#' @param cell_anno_path Path to cell annotation file.
#' @param umi_cutoff UMI per cell cutoff, default is 100.
#' @param matrix_control A list used to control how the counts matrix is stored
#'    in the CDS. By default, Monocle3 stores the counts matrix in-memory as a
#'    sparse matrix. Setting 'matrix_control=list(matrix_class="BPCells")',
#'    stores the matrix BPCells on-disk as a sparse matrix.
#'
#' @return cds object
#' @importFrom SingleCellExperiment counts
#'
#' @examples
#'   \donttest{
#'     pmat<-system.file("extdata", "matrix.mtx.gz", package = "monocle3")
#'     prow<-system.file("extdata", "features_c3h0.txt", package = "monocle3")
#'     pcol<-system.file("extdata", "barcodes_c2h0.txt", package = "monocle3")
#'     cds <- load_mtx_data( pmat, prow, pcol)
#'   }
#'
#' @export
load_mtx_data <- function( mat_path,
    gene_anno_path,
    cell_anno_path,
    umi_cutoff = 100,
    matrix_control=list()) {

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
        sep="\t",
        verbose=FALSE,
        matrix_control=matrix_control)
    return( cds )
  }

  if(!is.null(matrix_control[['matrix_class']]) && matrix_control[['matrix_class']] == 'BPCells') {
    matrix_control_default <- get_global_variable('matrix_control_bpcells_unrestricted')
  }
  else {
    matrix_control_default <- get_global_variable('matrix_control_csparsematrix_unrestricted')
  }
  matrix_control_res <- set_matrix_control(matrix_control=matrix_control, matrix_control_default=matrix_control_default, control_type='unrestricted')

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
  }
  else {
    if(ncol(mat) < 2) warning('bad loop: ncol(mat) < 2')
    mat <- mat[, 1:(ncol(mat)-1), drop=FALSE]
  }

  mat <- set_matrix_class(mat=mat, matrix_control=matrix_control_res)

  rownames(mat) <- gene.annotations$id
  colnames(mat) <- cell.annotations$cell

  cds <- new_cell_data_set(mat, cell_metadata = cell.annotations,
      gene_metadata = gene.annotations)
  colData(cds)$n.umi <- Matrix::colSums(counts(cds))
  cds <- cds[,colData(cds)$n.umi >= umi_cutoff]
  cds <- estimate_size_factors(cds)

  cds <- initialize_counts_metadata(cds)
  matrix_id <- get_unique_id(counts(cds))
  cds <- set_counts_identity(cds, mat_path, matrix_id)

  return(cds)
}


update_annoy_index <- function(annoy) {
  if(!is.null(annoy[['nn_index']][['version']]) && annoy[['nn_index']][['version']] == 2) {
    annoy_out <- annoy
  }
  else
  if(!is.null(annoy[['nn_index']][['version']]) && annoy[['nn_index']][['version']] == 1) {
    annoy_out <- S4Vectors::SimpleList()
    annoy_out[['nn_index']] <- list()
    annoy_out[['nn_index']][['method']] <- 'annoy'
    annoy_out[['nn_index']][['annoy_index']] <- annoy[['nn_index']][['annoy_index']]
    annoy_out[['nn_index']][['version']] <- 2
    annoy_out[['nn_index']][['annoy_index_version']] <- annoy[['index_version']]
    annoy_out[['nn_index']][['metric']] <- annoy[['metric']]
    annoy_out[['nn_index']][['n_trees']] <- annoy[['n_trees']]
    annoy_out[['nn_index']][['nrow']] <- annoy[['nrow']]
    annoy_out[['nn_index']][['ncol']] <- annoy[['ncol']]
    annoy_out[['nn_index']][['checksum_rownames']] <- NA_character_
    annoy_out[['matrix_id']] <- annoy[['matrix_id']]
    annoy_out[['updated']] <- TRUE
  }
  else
  if(!is.null(annoy[['nn_index']][['type']]) && annoy[['nn_index']][['type']] == 'annoyv1') {
    annoy_out <- S4Vectors::SimpleList()
    annoy_out[['nn_index']] <- list()
    annoy_out[['nn_index']][['method']] <- 'annoy'
    annoy_out[['nn_index']][['annoy_index']] <- annoy[['nn_index']][['ann']]
    annoy_out[['nn_index']][['version']] <- 2
    annoy_out[['nn_index']][['annoy_index_version']] <- annoy[['index_version']]
    annoy_out[['nn_index']][['metric']] <- annoy[['metric']]
    annoy_out[['nn_index']][['n_trees']] <- annoy[['n_trees']]
    annoy_out[['nn_index']][['nrow']] <- annoy[['nrow']]
    annoy_out[['nn_index']][['ncol']] <- annoy[['ncol']]
    annoy_out[['nn_index']][['checksum_rownames']] <- NA_character_
    annoy_out[['matrix_id']] <- annoy[['matrix_id']]
    annoy_out[['updated']] <- TRUE
  }
  else
  if(is.null(annoy[['nn_index']][['type']])) {
    annoy_out <- S4Vectors::SimpleList()
    annoy_out[['nn_index']] <- list()
    annoy_out[['nn_index']][['method']] <- 'annoy'
    annoy_out[['nn_index']][['annoy_index']] <- annoy[['nn_index']]
    annoy_out[['nn_index']][['version']] <- 2
    annoy_out[['nn_index']][['annoy_index_version']] <- NA_character_
    annoy_out[['nn_index']][['metric']] <- annoy[['metric']]
    annoy_out[['nn_index']][['metric']] <- ifelse(!is.null(annoy[['metric']]), annoy[['metric']], NA_character_)
    annoy_out[['nn_index']][['n_trees']] <- ifelse(!is.null(annoy[['n_trees']]), annoy[['n_trees']], NA_integer_)
    annoy_out[['nn_index']][['nrow']] <- ifelse(!is.null(annoy[['n_trees']]), annoy[['n_trees']], NA_integer_)
    annoy_out[['nn_index']][['ncol']] <- ifelse(!is.null(annoy[['n_trees']]), annoy[['n_trees']], NA_integer_)
    annoy_out[['nn_index']][['checksum_rownames']] <- NA_character_
    annoy_out[['matrix_id']] <- NA_character_
    annoy_out[['updated']] <- TRUE
  }

  return(annoy_out)
}


update_hnsw_index <- function(hnsw) {
  if(!is.null(hnsw[['nn_index']][['version']]) && hnsw[['nn_index']][['version']] == 1) {
    hnsw_out <- hnsw
  }
  else
  if(!is.null(hnsw[['nn_index']][['version']]) && hnsw[['nn_index']][['version']] == 1) {
    hnsw_out <- S4Vectors::SimpleList()
    annoy_out[['nn_index']] <- list()
    hnsw_out[['nn_index']][['method']] <- 'hnsw'
    hnsw_out[['nn_index']][['hnsw_index']] <- hnsw[['nn_index']]
    hnsw_out[['nn_index']][['version']] <- 1
    hnsw_out[['nn_index']][['hnsw_index_version']] <- hnsw[['index_version']]
    hnsw_out[['nn_index']][['metric']] <- hnsw[['metric']]
    hnsw_out[['nn_index']][['M']] <- hnsw[['M']]
    hnsw_out[['nn_index']][['ef_construction']] <- hnsw[['ef_construction']]
    hnsw_out[['nn_index']][['nrow']] <- hnsw[['nrow']]
    hnsw_out[['nn_index']][['ncol']] <- hnsw[['ncol']]
    hnsw_out[['nn_index']][['checksum_rownames']] <- NA_character_
    hnsw_out[['matrix_id']] <- hnsw[['matrix_id']]
    annoy_out[['updated']] <- TRUE
  }

  return(hnsw_out)
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
  if(!is.null(nn_index[['version']])) {
    if(nn_index[['version']] == 1 || nn_index[['version']] == 2) {
      tryCatch( nn_index[['annoy_index']]$save(file_name),
                error = function(e) {message('Unable to save annoy index: it may not exist in this cds: error message is ', e)})
    }
    else {
      stop('Unrecognized Monocle3 annoy index type')
    }
  }
  else
  if(!is.null(nn_index[['type']])) {
    if(nn_index[['type']] == 'annoyv1') {
      tryCatch( nn_index[['ann']]$save(file_name),
                error = function(e) {message('Unable to save annoy index: it may not exist in this cds: error message is ', e)})
    }
    else {
      stop('Unrecognized uwot annoy index type')
    }
  }
  else {
    nn_index$save(file_name)
  }
}


# see comments for save_annoy_index
# modified for RcppAnnoy
load_annoy_index <- function(nn_index, file_name, metric, ndim) {
  if(!is.null(nn_index[['version']])) {
    if(nn_index[['version']] == 1 || nn_index[['version']] == 2) {
      annoy_index <- new_annoy_index(metric, ndim)
      tryCatch(
        {
          annoy_index$load(file_name)
        }, error = function(emsg)
        {
          stop('load_annoy_index: bad status reading annoy index file')
        }
      )
      nn_index[['annoy_index']] <- annoy_index
    }
    else {
      stop('Unrecognized annoy index type')
    }
  }
  else
  if(!is.null(nn_index[['type']])) {
    if(nn_index[['type']] == 'annoyv1') {
      annoy_index <- new_annoy_index(metric, ndim)
      tryCatch(
        {
          annoy_index$load(file_name)
        }, error = function(emsg)
        {
          stop('load_annoy_index: bad status reading annoy index file')
        }
      )
      nn_index[['ann']] <- annoy_index
    }
    else {
      stop('Unrecognized annoy index type')
    }
  }
  else {
    # Assume to be an older uwot annoy index version.
    nn_index <- new_annoy_index(metric, ndim)
    nn_index$load(file_name)
  }
  return(nn_index)
}


save_umap_annoy_index <- function(nn_index, file_name) {
  if(!is.null(nn_index[['type']])) {
    if(nn_index[['type']] == 'annoyv1') {
      tryCatch( nn_index[['ann']]$save(file_name),
                error = function(e) {message('Unable to save annoy index: it may not exist in this cds: error message is ', e)})
    }
    else {
      stop('Unrecognized umap annoy index type')
    }
  }
  else {
    nn_index$save(file_name)
  }
}


load_umap_annoy_index <- function(nn_index, file_name, metric, ndim) {
  if(!is.null(nn_index[['type']])) {
    if(nn_index[['type']] == 'annoyv1') {
      annoy_index <- new_annoy_index(metric, ndim)
      tryCatch(
        {
          annoy_index$load(file_name)
        }, error = function(emsg)
        {
          stop('load_annoy_index: bad status reading annoy index file')
        }
      )
      nn_index[['ann']] <- annoy_index
    }
    else {
      stop('Unrecognized umap annoy index type')
    }
  }
  else {
    # Assume to be an older uwot annoy index version.
    nn_index <- new_annoy_index(metric, ndim)
    tryCatch(
      {
        nn_index$load(file_name)
      }, error = function(emsg)
      {
        stop('load_annoy_index: bad status reading annoy index file')
      }
    )
  }
  return(nn_index)
}


save_hnsw_index <- function(nn_index, file_name) {
  if(is.null(nn_index)) return()

  if(!is.null(nn_index[['version']])) {
    out_index <- nn_index[['hnsw_index']]
  }
  else {
    out_index <- nn_index
  }
  tryCatch(out_index$save(file_name),
           error = function(e) {message('Unable to save hnsw index: it may not exist in this cds: error message is ', e)})
}


load_hnsw_index <- function(nn_index, file_name, metric, ndim) {
  if(metric == 'l2') {
    tryCatch(
      {
        new_index <- methods::new(RcppHNSW::HnswL2, ndim, file_name)
      }, error = function(emsg)
      {
        stop('load_hnsw_index: bad status reading hnsw index file')
      }
    )
    unlink(file_name)
  }
  else
  if(metric == 'euclidean') {
    tryCatch(
      {
        new_index <- methods::new(RcppHNSW::HnswL2, ndim, file_name)
      }, error = function(emsg)
      {
        stop('load_hnsw_index: bad status reading hnsw index file')
      }
    )
    attr(new_index, "distance") <- "euclidean"
    unlink(file_name)
  }
  else
    if(metric == 'cosine') {
    tryCatch(
      {
        new_index <- methods::new(RcppHNSW::HnswCosine, ndim, file_name)
      }, error = function(emsg)
      {
        stop('load_hnsw_index: bad status reading hnsw index file')
      }
    )
    unlink(file_name)
  }
  else
  if(metric == 'ip') {
    tryCatch(
      {
        new_index <- methods::new(RcppHNSW::HnswIp, ndim, file_name)
      }, error = function(emsg)
      {
        stop('load_hnsw_index: bad status reading hnsw index file')
      }
    )
    unlink(file_name)
  }
  else
    stop('Unrecognized HNSW metric ', metric)

  if(!is.null(nn_index[['version']]))
    nn_index[['hnsw_index']] <- new_index
  else
    nn_index <- new_index

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
    save_umap_annoy_index(umap_model[['nn_index']], file_name)
    md5sum_umap_index <- tools::md5sum(file_name)
  }
  else {
    warning('save_umap_nn_indexes is untested with more than one umap metric')
    md5sum_vec <- character()
    for(i in seq(1, n_metrics, 1)) {
      file_name_expand <- paste0(file_name, i)
      save_umap_annoy_index(umap_model[['nn_index']][[i]], file_name_expand)
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
    umap_model[['nn_index']] <- load_umap_annoy_index(umap_model[['nn_index']], file_name, annoy_metric, annoy_ndim)
  }
  else {
    warning('load_umap_nn_indexes is untested with more than one umap metric')
    if(!is.null(md5sum_umap_index)) {
      md5sum_vec <- unlist(strsplit(md5sum_umap_index, '_', fixed=TRUE))
    }
    for(i in seq(1, n_metrics, 1)) {
      file_name_expand <- paste0(file_name, i)
      md5sum <- tools::md5sum(file_name_expand)
      if(!is.null(md5sum_umap_index) && md5sum != md5sum_vec[[i]]) {
        stop('The UMAP annoy index file, \'', file_name_expand, '\', differs from the file made using the save_reduce_dimension_model() function.')
      }
      metric <- metrics[[i]]
      annoy_metric <- if(metric == 'correlation') 'cosine' else metric
      annoy_ndim <- length(umap_model[['metric']][[i]])
      umap_model[['nn_index']][[i]] <- load_umap_annoy_index(umap_model[['nn_index']][[i]], file_name_expand, annoy_metric, annoy_ndim)
    }
  }
  return(umap_model)
}

# This is a specialized function for use in load_monocle_objects. There are
# no matrix_control checks and it requires the path to an existing
# BPCells matrix stored in a directory. The matrix control is used only
# to set the resulting the matrix_path.
load_bpcells_matrix_dir <- function(file_name, md5sum, matrix_control=list()) {
  matrixDirTmp <- BPCells::open_matrix_dir(dir=file_name, buffer_size=8192L)
  matrix_info <- get_matrix_info(mat=matrixDirTmp)
  matrix_control_res <- list(matrix_class=matrix_info[['matrix_class']],
                             matrix_mode=matrix_info[['matrix_mode']],
                             matrix_type=matrix_info[['matrix_type']],
                             matrix_compress=matrix_info[['matrix_compress']],
                             matrix_path=matrix_control[['matrix_path']],
                             matrix_buffer_size=matrix_info[['matrix_buffer_size']],
                             matrix_bpcells_copy=TRUE)
  matrixDir <- set_matrix_class(mat=matrixDirTmp, matrix_control=matrix_control_res)
  return(matrixDir)
}


#
# Report files saved.
#
report_files_saved <- function(file_index) {
  appendLF <- TRUE
  files <- file_index[['files']]
  for( i in seq_along(files[['cds_object']])) {
    cds_object <- files[['cds_object']][[i]]
    reduction_method <- files[['reduction_method']][[i]]
    if(cds_object == 'cds') {
      process <- 'cell_data_set'
      reduction_method <- 'full_cds'
    }
    else
    if(cds_object == 'reduce_dim_aux') {
      if(reduction_method == 'Aligned') {
        process <- 'align_cds'
      }
      else
      if(reduction_method == 'PCA' || reduction_method == 'LSI') {
        process <- 'preprocess_cds'
      }
      else
      if(reduction_method == 'tSNE' || reduction_method == 'UMAP') {
        process <- 'reduce_dimension'
      }
      else {
        stop('Unrecognized preprocess reduction_method \'', reduction_method, '\'')
      }
    }
    else
    if(cds_object == 'bpcells_matrix_dir') {
      process <- 'BPCells MatrixDir'
      reduction_method <- 'full_counts_matrix'
    }
    else {
      stop('Unrecognized cds_object value \'', files[['cds_object']][[i]], '\'')
    }
    file_format <- files[['file_format']][[i]]
    if(file_format == 'rds') {
      file_type <- 'RDS'
    }
    else
    if(file_format == 'hdf5') {
      file_type <- 'RDS_HDF5'
    }
    else
    if(file_format == 'annoy_index' || file_format == 'hnsw_index') {
      file_type <- 'NN_index'
    }
    else
    if(file_format == 'umap_annoy_index') {
      file_type <- 'UMAP_NN_index'
    }
    else
    if(file_format == 'BPCells:MatrixDir') {
      file_type <- 'BPCells:MatrixDir'
    }
    else {
      stop('Unrecognized file_format value \'', file_format, '\'')
    }

    file_name <- basename(files[['file_path']][[i]])
    message('  ', file_name, '  (', reduction_method, '  ', file_type, '  from  ', process, ')', appendLF=appendLF)
  }
}


#
#' Save cell_data_set transform models.
#'
#' Save the transform models in the cell_data_set to the
#' specified directory by writing the R objects to RDS
#' files and the nearest neighbor indexes to
#' index files. save_transform_models saves transform
#' models made by running the preprocess_cds and
#' reduce_dimension functions on an initial cell_data_set.
#' Subsequent cell_data_sets are transformed into the
#' reduced dimension space of the initial cell_data_set by
#' loading the new data into a new cell_data_set, loading
#' the initial data set transform models into the new
#' cell_data_set using the load_transform_models function,
#' and applying those transform models to the new data set
#' using the preprocess_transform and
#' reduce_dimension_transform functions. In this case, do
#' not run the preprocess_cds or reduce_dimension
#' functions on the new cell_data_set. Additionally,
#' save_transform_models saves nearest neighbor indexes
#' when the preprocess_cds and reduce_dimension
#' functions are run with the make_nn_index=TRUE parameter.
#' These indexes are used to find matches between cells in
#' the new processed cell_data_set and the initial
#' cell_data_set using index search functions. For more
#' information see the help for transfer_cell_labels.
#' save_transform_models saves the models to a directory
#' given by directory_path.
#'
#' @param cds a cell_data_set with existing models.
#' @param directory_path a string giving the name of the directory
#'   in which to write the model files.
#' @param comment a string with optional notes that is saved with
#'   the objects.
#' @param verbose a boolean determining whether to print information
#'   about the saved files.
#'
#' @return none.
#'
#' @examples
#'   \dontrun{
#'     cds <- load_a549()
#'     cds <- preprocess_cds(cds)
#'     cds <- reduce_dimension(cds)
#'     save_transform_models(cds, 'tm')
#'   }
#'
#' @export
# Bioconductor forbids writing to user directories so examples
# is not run.
save_transform_models <- function( cds, directory_path, comment="", verbose=TRUE) {
  appendLF <- TRUE
  # file information is written to an RDS file
  # in directory_path
  #   cds_object: reduce_dim_aux
  #   reduction_method: PCA | LSI | Aligned | UMAP ...
  #   nn_method: annoy | hnsw
  #   object_spec: ex: cds@reduce_dim_aux[[reduction_method]]
  #   file_format: rds | annoy_index | umap_annoy_index | hnsw_index
  #   file_path: path within directory_path (need file name)
  #   file_md5sum: md5sum of file(s)
  file_index <- list( 'save_function' = 'save_transform_models',
                      'archive_date' = Sys.time(),
                      'r_version' = R.Version()$version.string,
                      'uwot_version' = utils::packageVersion('uwot'),
                      'hnsw_version' = utils::packageVersion('RcppHNSW'),
                      'monocle_version' = utils::packageVersion('monocle3'),
                      'cds_version' = S4Vectors::metadata(cds)$cds_version,
                      'archive_version' = get_global_variable('transform_models_version'),
                      'directory' = directory_path,
                      'comment' = comment,
                      'files' = data.frame(cds_object = character(0),
                                           reduction_method = character(0),
                                           object_spec = character(0),
                                           file_format = character(0),
                                           file_path = character(0),
                                           file_md5sum = character(0),
                                           stringsAsFactors = FALSE))

  # Gather reduce_dimension reduction_methods and whether each has an annoy index.
  methods_reduce_dim <- list()
  for( reduction_method in names(cds@reduce_dim_aux)) {
    methods_reduce_dim[[reduction_method]] <- list()
    methods_reduce_dim[[reduction_method]][['rds_path']] <- paste0('rdd_', tolower(reduction_method), '_transform_model.rds')

    methods_reduce_dim[[reduction_method]][['annoy_index_path']] <- paste0('rdd_', tolower(reduction_method), '_transform_model_annoy.idx')
    methods_reduce_dim[[reduction_method]][['hnsw_index_path']] <- paste0('rdd_', tolower(reduction_method), '_transform_model_hnsw.idx')

    if(reduction_method == 'UMAP') {
      if(has_nn_index(cds, 'umap_model_annoy')) {
        methods_reduce_dim[[reduction_method]][['has_model_index']] <- TRUE
        methods_reduce_dim[[reduction_method]][['umap_index_path']] <- paste0('rdd_', tolower(reduction_method), '_transform_model_umap.idx')
      }
      else
        methods_reduce_dim[[reduction_method]][['has_model_index']] <- FALSE
    }

    if(has_nn_index(cds, paste0(tolower(reduction_method), '_search_annoy')))
      methods_reduce_dim[[reduction_method]][['has_annoy_index']] <- TRUE
    else
      methods_reduce_dim[[reduction_method]][['has_annoy_index']] <- FALSE

    if(has_nn_index(cds, paste0(tolower(reduction_method), '_search_hnsw')))
      methods_reduce_dim[[reduction_method]][['has_hnsw_index']] <- TRUE
    else
      methods_reduce_dim[[reduction_method]][['has_hnsw_index']] <- FALSE

  }

  # Make directory if necessary.
  dir.create(path = directory_path, showWarnings=FALSE, recursive=TRUE, mode='0700')

  # Remove files, if they exist.
  for(reduction_method in names(methods_reduce_dim)) {
    if(file.exists(file.path(directory_path, methods_reduce_dim[[reduction_method]][['rds_path']])))
      file.remove(file.path(directory_path, methods_reduce_dim[[reduction_method]][['rds_path']]))

    if(file.exists(file.path(directory_path, methods_reduce_dim[[reduction_method]][['annoy_index_path']])))
      file.remove(file.path(directory_path, methods_reduce_dim[[reduction_method]][['annoy_index_path']]))

    if(file.exists(file.path(directory_path, methods_reduce_dim[[reduction_method]][['hnsw_index_path']])))
      file.remove(file.path(directory_path, methods_reduce_dim[[reduction_method]][['hnsw_index_path']]))

    if(reduction_method == 'UMAP') {
      if(methods_reduce_dim[[reduction_method]][['has_model_index']]) {
        if(file.exists(file.path(directory_path, methods_reduce_dim[[reduction_method]][['umap_index_path']])))
           file.remove(file.path(directory_path, methods_reduce_dim[[reduction_method]][['umap_index_path']]))
      }
    }
  }

  # Save reduce_dimension annoy indexes.
  # Notes:
  #   o  save RDS files before the corresponding index files in
  #      order to enable loading.
  #
  for(reduction_method in names(methods_reduce_dim)) {
    tryCatch(
      {
        base::saveRDS(cds@reduce_dim_aux[[reduction_method]], file=file.path(directory_path, methods_reduce_dim[[reduction_method]][['rds_path']]))
      },
      error = function(cond) {
                     message('problem writing file \'', file.path(directory_path, methods_reduce_dim[[reduction_method]][['rds_path']]), '\': ', cond, appendLF=appendLF)
                     return(NULL)
      },
      finally = {
        md5sum <- tools::md5sum(file.path(directory_path, methods_reduce_dim[[reduction_method]][['rds_path']]))
        file_index[['files']] <- rbind(file_index[['files']],
                                       data.frame(cds_object = 'reduce_dim_aux',
                                                  reduction_method = reduction_method,
                                                  object_spec = object_name_to_string(cds@reduce_dim_aux[[reduction_method]]),
                                                  file_format = 'rds',
                                                  file_path = methods_reduce_dim[[reduction_method]][['rds_path']],
                                                  file_md5sum = md5sum,
                                                  stringsAsFactors = FALSE))
    })
    if(methods_reduce_dim[[reduction_method]][['has_annoy_index']]) {
      tryCatch(
        {
          save_annoy_index(cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy']][['nn_index']], file.path(directory_path, methods_reduce_dim[[reduction_method]][['annoy_index_path']]))
        },
        error = function(cond) {
                       message('problem writing file \'', file.path(directory_path, methods_reduce_dim[[reduction_method]][['annoy_index_path']]), '\': ', cond, appendLF=appendLF)
                       return(NULL)
        },
        finally = {
          md5sum <- tools::md5sum(file.path(directory_path, methods_reduce_dim[[reduction_method]][['annoy_index_path']]))
          file_index[['files']] <- rbind(file_index[['files']],
                                         data.frame(cds_object = 'reduce_dim_aux',
                                                    reduction_method = reduction_method,
                                                    object_spec = object_name_to_string(cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy']][['nn_index']]),
                                                    file_format = 'annoy_index',
                                                    file_path = methods_reduce_dim[[reduction_method]][['annoy_index_path']],
                                                    file_md5sum = md5sum,
                                                    stringsAsFactors = FALSE))
        })
    }

    if(methods_reduce_dim[[reduction_method]][['has_hnsw_index']]) {
      tryCatch(
        {
          save_hnsw_index(cds@reduce_dim_aux[[reduction_method]][['nn_index']][['hnsw']][['nn_index']], file.path(directory_path, methods_reduce_dim[[reduction_method]][['hnsw_index_path']]))
        },
        error = function(cond) {
                       message('problem writing file \'', file.path(directory_path, methods_reduce_dim[[reduction_method]][['hnsw_index_path']]), '\': ', cond, appendLF=appendLF)
                       return(NULL)
        },
        finally = {
          md5sum <- tools::md5sum(file.path(directory_path, methods_reduce_dim[[reduction_method]][['hnsw_index_path']]))
          file_index[['files']] <- rbind(file_index[['files']],
                                         data.frame(cds_object = 'reduce_dim_aux',
                                                    reduction_method = reduction_method,
                                                    object_spec = object_name_to_string(cds@reduce_dim_aux[[reduction_method]][['nn_index']][['hnsw']][['nn_index']]),
                                                    file_format = 'hnsw_index',
                                                    file_path = methods_reduce_dim[[reduction_method]][['hnsw_index_path']],
                                                    file_md5sum = md5sum,
                                                    stringsAsFactors = FALSE))
        })
    }
    if(reduction_method == 'UMAP' && methods_reduce_dim[[reduction_method]][['has_model_index']]) {
      tryCatch(
        {
          md5sum <- save_umap_nn_indexes(cds@reduce_dim_aux[[reduction_method]][['model']][['umap_model']], file.path(directory_path, methods_reduce_dim[[reduction_method]][['umap_index_path']]))
        },
        error = function(cond) {
                       message('problem writing file \'', file.path(directory_path, methods_reduce_dim[[reduction_method]][['umap_index_path']]), '\': ', cond, appendLF=appendLF)
                       return(NULL)
        },
        finally = {
          file_index[['files']] <- rbind(file_index[['files']],
                                         data.frame(cds_object = 'reduce_dim_aux',
                                                    reduction_method = reduction_method,
                                                    object_spec = object_name_to_string(cds@reduce_dim_aux[[reduction_method]][['model']][['umap_model']]),
                                                    file_format = 'umap_annoy_index',
                                                    file_path = methods_reduce_dim[[reduction_method]][['umap_index_path']],
                                                    file_md5sum = md5sum,
                                                    stringsAsFactors = FALSE))
        })
    }
  }

  # Save file_index.rds.
  base::saveRDS(file_index, file=file.path(directory_path, 'file_index.rds'))

  if(verbose) {
    report_files_saved(file_index)
  }
}


#
# Copy reduce_dim_aux slot objects from cds_src to cds_dst. The copy
# runs on the R objects only. It does not try to copy objects about
# which R knows nothing, such as annoy and hnsw indices. At this
# time, this function is used only by load_transform_models().
#
copy_reduce_dim_aux <- function(cds_dst, cds_src) {
  for( reduction_method in names(cds_src@reduce_dim_aux@listData)) {
    cds_dst@reduce_dim_aux[[reduction_method]] <- cds_src@reduce_dim_aux[[reduction_method]]
  }
  return(cds_dst)
}


#' Load transform models into a cell_data_set.
#'
#' Load transform models, which were saved using save_transform_models,
#' into a cell_data_set. This function over-writes existing models in
#' the cell_data_set. For more information read the help information
#' for save_transform_models.
#'
#' @param cds a cell_data_set to be transformed using the models.
#' @param directory_path a string giving the name of the directory
#'   from which to read the model files. The model file directory
#'   is made by either save_transform_models() or
#'   save_monocle_objects().
#'
#' @return a cell_data_set with the transform models loaded by
#'   load_transform_models.
#'
#' @examples
#'   \dontrun{
#'     cds <- load_a549()
#'     cds <- preprocess_cds(cds)
#'     cds <- reduce_dimension(cds)
#'     save_transform_models(cds, 'tm')
#'     cds1 <- load_a549()
#'     cds1 <- load_transform_models(cds1, 'tm')
#'   }
#' @export
# Bioconductor forbids writing to user directories so examples
# is not run.
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
    error = function(cond) {
              message('problem reading file \'', file_index_path, '\': ', cond, appendLF=appendLF);
              return(NULL)
    }
  )

  # Check that this is a save_transform_models archive.
  if(file_index[['save_function']] != 'save_transform_models' && file_index[['save_function']] != 'save_monocle_objects') {
    stop('The files in ', directory_path, ' are not from save_transform_models or save_monocle_objects.')
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
    reduction_method <- file_index[['files']][['reduction_method']][[ifile]]
    md5sum <- file_index[['files']][['file_md5sum']][[ifile]]

    if(!file.exists(file_path))
      stop('load_transform_models: missing file \'', file_path, '\'')

    #
    # For UWOT UMAP annoy index, the function load_umap_nn_indexes
    # checks md5sums internally so don't check here.
    #
    md5sum_file <- tools::md5sum(file_path)
    if(!(cds_object == 'reduce_dim_aux' &&
         reduction_method == 'UMAP' &&
         file_format == 'umap_nn_index' &&
         nchar(md5sum) > 32)) {
      if(is.na(md5sum_file) || (md5sum_file != md5sum)) {
        stop('md5sum mismatch for file \'', file_path, '\'')
      }
    }

    #
    # Note:
    #   o  expect that the RDS file for a reduction_method
    #      (or cds) appears before index files for the
    #      reduction_method
    #

    if(cds_object == 'cds') {
      if(file_format == 'rds') {
        cds_tmp <- tryCatch(
          {
            readRDS(file_path)
          },
          error = function(cond) {
            message('problem reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
          })
        cds <- copy_reduce_dim_aux(cds, cds_tmp)
        rm(cds_tmp)
      }
    }
    else

    if(cds_object == 'reduce_dim_aux') {
      if(file_format == 'rds') {
        cds@reduce_dim_aux[[reduction_method]] <- tryCatch(
          {
            readRDS(file_path)
          },
          error = function(cond) {
            message('problem reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
          })
      }
      else
      if(file_format == 'annoy_index') {
        cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy']] <- update_annoy_index(cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy']])

        metric <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy']][['nn_index']][['metric']]
        ncolumn <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy']][['nn_index']][['ncol']]
        cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy']][['nn_index']] <- tryCatch(
          {
            load_annoy_index(cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy']][['nn_index']], file_path, metric, ncolumn)
          },
          error = function(cond) {
            message('problem reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
          })
      }
      else
      if(file_format == 'hnsw_index') {
        cds@reduce_dim_aux[[reduction_method]][['nn_index']][['hnsw']] <- update_hnsw_index(cds@reduce_dim_aux[[reduction_method]][['nn_index']][['hnsw']])

        metric <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][['hnsw']][['nn_index']][['metric']]
        ncolumn <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][['hnsw']][['nn_index']][['ncol']]
        cds@reduce_dim_aux[[reduction_method]][['nn_index']][['hnsw']][['nn_index']] <- tryCatch(
          {
            load_hnsw_index(cds@reduce_dim_aux[[reduction_method]][['nn_index']][['hnsw']][['nn_index']], file_path, metric, ncolumn)
          },
          error = function(cond) {
            message('problem reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
          })
      }
      else
      if(reduction_method == 'UMAP' && file_format == 'umap_annoy_index') {
        cds@reduce_dim_aux[[reduction_method]][['model']][['umap_model']] <- tryCatch(
          {
            load_umap_nn_indexes(cds@reduce_dim_aux[[reduction_method]][['model']][['umap_model']], file_path, md5sum)
          },
          error = function(cond) {
            message('problem reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
         })
      }
      else {
        stop('Unrecognized file format value \'', file_format, '\'')
      }
      cds <- set_model_identity_path(cds, reduction_method, directory_path)
    }
    else {
      stop('Unrecognized cds_object value \'', cds_object, '\'')
    }
  }

  return(cds)
}


#
# Check cds assays for HDF5Array objects.
#
#' @importFrom S4Vectors getListElement
test_hdf5_assays <- function(cds) {
  assays <- assays(cds)
  for( idx in seq_along(assays)) {
    asyl <- getListElement(assays, idx)
    hdf5_test <- unlist(DelayedArray::seedApply(asyl, methods::is, "HDF5ArraySeed"))
    if(any(unlist(hdf5_test))) return(TRUE)
  }
  FALSE
}


#
#' Save a Monocle3 full cell_data_set.
#'
#' Save a Monocle3 full cell_data_set to a specified directory
#' by writing the R objects to RDS files and the nearest
#' neighbor indexes to index files. The assays
#' objects are saved as HDF5Array files when hdf5_assays=TRUE
#' or when the cell_data_set assays are HDF5Array objects. If
#' any assay in the cell_data set is an HDF5 object, all assays
#' must be. When save_monocle_objects is run with hdf5_assays=TRUE,
#' the load_monocle_objects function loads the saved assays into
#' HDF5Array objects in the resulting cell_data_set. Note:
#' operations such as preprocess_cds that are run on assays stored
#' as HDF5Arrays are much, much slower than the same operations
#' run on assays stored as in-memory matrices. You may want to
#' investigate parameters related to the Bioconductor DelayedArray
#' and BiocParallel packages in this case.
#'
#' @param cds a cell_data_set to save.
#' @param directory_path a string giving the name of the directory
#'   in which to write the object files.
#' @param hdf5_assays a boolean determining whether the
#'   non-HDF5Array assay objects are saved as HDF5 files. At this
#'   time cell_data_set HDF5Array assay objects are stored as
#'   HDF5Assay files regardless of the hdf5_assays parameter value.
#' @param comment a string with optional notes that is saved with
#'   the objects.
#' @param verbose a boolean determining whether to print information
#'   about the saved files.
#' @param archive_control a list that is used to control archiving
#'   the output directory. The archive_control parameters are
#'   \describe{
#'     \item{archive_type}{a string giving the method used to
#'        archive the directory. The acceptable values are
#'        "tar" and "none". The directory is not archived when
#'        archive_type is "none". The default is "tar".}
#'     \item{archive_compression}{a string giving the type of
#'        compression applied to the archive file. The acceptable
#'        values are "none", "gzip", "bzip2", and "xz". The
#'        default is "none".}
#'   }
#'   Note: the output directory is not removed after it is
#'         archived.
#'
#' @return none.
#'
#' @examples
#'   \dontrun{
#'     cds <- load_a549()
#'     save_monocle_objects(cds, 'mo')
#'   }
#'
#' @export
# Bioconductor forbids writing to user directories so examples
# is not run.
save_monocle_objects <- function(cds, directory_path, hdf5_assays=FALSE, comment="", verbose=TRUE, archive_control=list(archive_type="tar", archive_compression="none")) {

  if(is.null(archive_control[['archive_type']])) archive_control[['archive_type']] <- 'tar'
  if(is.null(archive_control[['archive_compression']])) archive_control[['archive_compression']] <- 'none'

  assertthat::assert_that(archive_control[['archive_type']] %in% c('tar', 'none'),
    msg=paste0("archive_type must be either \'none\' or \'tar\'"))
  assertthat::assert_that(archive_control[['archive_compression']] %in% c('gzip', 'bzip2', 'xz', 'none'),
    msg=paste0("archive_compression must be \'none\', \'gzip\', \'bzip2\', or \'xz\'."))
  assertthat::assert_that(archive_control[['archive_compression']] %in% c('gzip', 'bzip2', 'xz', 'none'))

  appendLF <- TRUE
  # file information is written to an RDS file
  # in directory_path
  #   cds_object: cds | reduce_dim_aux
  #   reduction_method: PCA | LSI | Aligned | tSNE | UMAP  ...
  #   nn_method: annoy | hnsw
  #   object_spec: ex: cds@reduce_dim_aux[[reduction_method]]
  #   file_format: rds | annoy_index | umap_annoy_index | hnsw_index
  #   file_path: path within directory_path (need file name)
  #   file_md5sum: md5sum of file(s)
  file_index <- list( 'save_function' = 'save_monocle_objects',
                      'archive_date' = Sys.time(),
                      'r_version' = R.Version()$version.string,
                      'uwot_version' = utils::packageVersion('uwot'),
                      'hnsw_version' = utils::packageVersion('RcppHNSW'),
                      'hdf5array_version' = utils::packageVersion('HDF5Array'),
                      'bpcells_version' = utils::packageVersion('BPCells'),
                      'monocle_version' = utils::packageVersion('monocle3'),
                      'cds_version' = S4Vectors::metadata(cds)$cds_version,
                      'archive_version' = get_global_variable('monocle_objects_version'),
                      'directory' = directory_path,
                      'comment' = comment,
                      'files' = data.frame(cds_object = character(0),
                                           reduction_method = character(0),
                                           object_spec = character(0),
                                           file_format = character(0),
                                           file_path = character(0),
                                           file_md5sum = character(0),
                                           stringsAsFactors = FALSE))

  # Save assays as HDF5Array objects?
  hdf5_assay_flag <- hdf5_assays || test_hdf5_assays(cds)

  # Save assays as BPCells Matrix_dir or BPCells 10xHDF5 file?
  bpcells_matrix_dir_flag <- FALSE
  matrix_info <- get_matrix_info(mat=counts(cds))
  if(matrix_info[['matrix_class']] == 'BPCells' &&
     matrix_info[['matrix_mode']] == 'dir') {
    #
    # Check that BPCells matrix directory exists.
    # We use the bpcells_matrix_dir_flag to save the
    # directory (later) and so we don't try to save
    # the BPCells matrix files if the directory
    # doesn't exist.
    #
    if(dir.exists(matrix_info[['matrix_path']])) {
      bpcells_matrix_dir_flag <- TRUE
    }
    else {
      message('save_monocle_objects: warning: the CDS has a BPCells count matrix but\nbut the BPCells count matrix directory is missing, which will likely\ncause problems in the future.\nI\'m continuing without it.')
    }
  }

  # Path of cds object file.
  rds_path <- 'cds_object.rds'
  hdf5_path <- 'hdf5_object'
  bpcells_matrix_dir <- 'bpcells_matrix_dir'
  
  # Gather reduce_dimension reduction_method names for which indexes exist.
  methods_reduce_dim <- list()
  for(reduction_method in names(cds@reduce_dim_aux)) {
    methods_reduce_dim[[reduction_method]] <- list()

    # These are the nn search indices make when build_nn_index=TRUE.
    # These names are confusing...
    methods_reduce_dim[[reduction_method]][['annoy_index_path']] <- paste0('rdd_', tolower(reduction_method), '_transform_model_annoy.idx')
    methods_reduce_dim[[reduction_method]][['hnsw_index_path']] <- paste0('rdd_', tolower(reduction_method), '_transform_model_hnsw.idx')

    # This is the annoy index used internally by UMAP.
    # These names are confusing...
    if(reduction_method == 'UMAP') {
      if(has_nn_index(cds, 'umap_model_annoy')) {
        methods_reduce_dim[[reduction_method]][['has_model_index']] <- TRUE
        methods_reduce_dim[[reduction_method]][['umap_index_path']] <- paste0('rdd_', tolower(reduction_method), '_transform_model_umap.idx')
      }
      else
        methods_reduce_dim[[reduction_method]][['has_model_index']] <- FALSE
    }

    if(has_nn_index(cds, paste0(tolower(reduction_method), '_search_annoy')))
      methods_reduce_dim[[reduction_method]][['has_annoy_index']] <- TRUE
    else
      methods_reduce_dim[[reduction_method]][['has_annoy_index']] <- FALSE

    if(has_nn_index(cds, paste0(tolower(reduction_method), '_search_hnsw')))
      methods_reduce_dim[[reduction_method]][['has_hnsw_index']] <- TRUE
    else
      methods_reduce_dim[[reduction_method]][['has_hnsw_index']] <- FALSE
  }

  # Make directory if necessary.
  dir.create(path = directory_path, showWarnings=FALSE, recursive=TRUE, mode='0700')

  # Remove files, if they exist.

  # BPCells MatrixDir directory.
  if(file.exists(file.path(directory_path, bpcells_matrix_dir)))
    unlink(file.path(directory_path, bpcells_matrix_dir), recursive=TRUE)

  # Reduction method related files.
  for(reduction_method in names(methods_reduce_dim)) {
    if(file.exists(file.path(directory_path, rds_path)))
      file.remove(file.path(directory_path, rds_path))

    if(file.exists(file.path(directory_path, methods_reduce_dim[[reduction_method]][['annoy_index_path']])))
       file.remove(file.path(directory_path, methods_reduce_dim[[reduction_method]][['annoy_index_path']]))

    if(file.exists(file.path(directory_path, methods_reduce_dim[[reduction_method]][['hnsw_index_path']])))
       file.remove(file.path(directory_path, methods_reduce_dim[[reduction_method]][['hnsw_index_path']]))

    if(reduction_method == 'UMAP') {
      if(methods_reduce_dim[[reduction_method]][['has_model_index']]) {
        if(file.exists(file.path(directory_path, methods_reduce_dim[[reduction_method]][['umap_index_path']])))
           file.remove(file.path(directory_path, methods_reduce_dim[[reduction_method]][['umap_index_path']]))
      }
    }
  }

  #
  # Save cds object.
  # Notes:
  #   o  allow for HDF5Array assay objects.
  #
  if(!hdf5_assay_flag) {
    tryCatch(
      {
        base::saveRDS(cds, file.path(directory_path, rds_path))
      },
      error = function(cond) {
                       message('problem writing file \'', file.path(directory_path, rds_path), '\': ', cond, appendLF=appendLF)
                       return(NULL)
      },
      finally = {
        md5sum <- tools::md5sum(file.path(directory_path, rds_path))
        file_index[['files']] <- rbind(file_index[['files']],
                                       data.frame(cds_object = 'cds',
                                                  reduction_method = NA,
                                                  object_spec = object_name_to_string(cds),
                                                  file_format = 'rds',
                                                  file_path = rds_path,
                                                  file_md5sum = md5sum,
                                                  stringsAsFactors = FALSE))
        # Save BCells MatrixDir, if required.
        if(bpcells_matrix_dir_flag) {
          bpcells_matrix_path <- file.path(directory_path, bpcells_matrix_dir)
          mat <- counts(cds)
          BPCells::write_matrix_dir(mat=mat, dir=bpcells_matrix_path, compress=FALSE, buffer_size=8192L, overwrite=FALSE)

          file_index[['files']] <- rbind(file_index[['files']],
                                         data.frame(cds_object = 'bpcells_matrix_dir',
                                                    reduction_method = NA,
                                                    object_spec = object_name_to_string(mat),
                                                    file_format = 'BPCells:MatrixDir',
                                                    file_path = bpcells_matrix_dir,
                                                    file_md5sum = NA,
                                                    stringsAsFactors = FALSE))
        }
      })
  }
  else {
    tryCatch(
      {
        HDF5Array::saveHDF5SummarizedExperiment(cds, file.path(directory_path, hdf5_path), replace=TRUE)
      },
      error = function(cond) {
                       message('problem writing file \'', file.path(directory_path, hdf5_path), '\': ', cond, appendLF=appendLF)
                       return(NULL)
      },
      finally = {
        md5sum <- tools::md5sum(file.path(directory_path, hdf5_path, 'se.rds'))
        file_index[['files']] <- rbind(file_index[['files']],
                                       data.frame(cds_object = 'cds',
                                                  reduction_method = NA,
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
  for(reduction_method in names(methods_reduce_dim)) {
    if(methods_reduce_dim[[reduction_method]][['has_annoy_index']]) {
      tryCatch(
        {
          save_annoy_index(cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy']][['nn_index']], file.path(directory_path, methods_reduce_dim[[reduction_method]][['annoy_index_path']]))
        },
        error = function(cond) {
                       message('problem writing file \'', file.path(directory_path, methods_reduce_dim[[reduction_method]][['annoy_index_path']]), '\': ', cond, appendLF=appendLF)
                       return(NULL)
        },
        finally = {
          md5sum <- tools::md5sum(file.path(directory_path, methods_reduce_dim[[reduction_method]][['annoy_index_path']]))
          file_index[['files']] <- rbind(file_index[['files']],
                                         data.frame(cds_object = 'reduce_dim_aux',
                                                    reduction_method = reduction_method,
                                                    object_spec = object_name_to_string(cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy']][['nn_index']]),
                                                    file_format = 'annoy_index',
                                                    file_path = methods_reduce_dim[[reduction_method]][['annoy_index_path']],
                                                    file_md5sum = md5sum,
                                                    stringsAsFactors = FALSE))
        })
    }
    if(methods_reduce_dim[[reduction_method]][['has_hnsw_index']]) {
      tryCatch(
        {
          save_hnsw_index(cds@reduce_dim_aux[[reduction_method]][['nn_index']][['hnsw']][['nn_index']], file.path(directory_path, methods_reduce_dim[[reduction_method]][['hnsw_index_path']]))
        },
        error = function(cond) {
                       message('problem writing file \'', file.path(directory_path, methods_reduce_dim[[reduction_method]][['hswn_index_path']]), '\': ', cond, appendLF=appendLF)
                       return(NULL)
        },
        finally = {
          md5sum <- tools::md5sum(file.path(directory_path, methods_reduce_dim[[reduction_method]][['hnsw_index_path']]))
          file_index[['files']] <- rbind(file_index[['files']],
                                         data.frame(cds_object = 'reduce_dim_aux',
                                                    reduction_method = reduction_method,
                                                    object_spec = object_name_to_string(cds@reduce_dim_aux[[reduction_method]][['nn_index']][['hnsw']][['nn_index']]),
                                                    file_format = 'hnsw_index',
                                                    file_path = methods_reduce_dim[[reduction_method]][['hnsw_index_path']],
                                                    file_md5sum = md5sum,
                                                    stringsAsFactors = FALSE))
        })
    }
    if(reduction_method == 'UMAP' && methods_reduce_dim[[reduction_method]][['has_model_index']]) {
      tryCatch(
        {
          md5sum <- save_umap_nn_indexes(cds@reduce_dim_aux[[reduction_method]][['model']][['umap_model']], file.path(directory_path, methods_reduce_dim[[reduction_method]][['umap_index_path']]))
        },
        error = function(cond) {
                       message('problem writing file \'', file.path(directory_path, methods_reduce_dim[[reduction_method]][['umap_index_path']]), '\': ', cond, appendLF=appendLF)
                       return(NULL)
        },
        finally = {
          file_index[['files']] <- rbind(file_index[['files']],
                                         data.frame(cds_object = 'reduce_dim_aux',
                                                    reduction_method = reduction_method,
                                                    object_spec = object_name_to_string(cds@reduce_dim_aux[[reduction_method]][['model']][['umap_model']]),
                                                    file_format = 'umap_annoy_index',
                                                    file_path = methods_reduce_dim[[reduction_method]][['umap_index_path']],
                                                    file_md5sum = md5sum,
                                                    stringsAsFactors = FALSE))
        })
    }
  }

  # Save file_index.rds.
  base::saveRDS(file_index, file=file.path(directory_path, 'file_index.rds'))

  if(verbose) {
    report_files_saved(file_index)
  }

  if(archive_control[['archive_type']] == 'tar') {
    if(archive_control[['archive_compression']] == 'gzip') {
      archive_name <- paste0(directory_path, '.tar.gz')
    }
    else
    if(archive_control[['archive_compression']] == 'bzip2') {
      archive_name <- paste0(directory_path, '.tar.bz2')
    }
    else
    if(archive_control[['archive_compression']] == 'xz') {
      archive_name <- paste0(directory_path, '.tar.xz')
    }
    else {
      archive_name <- paste0(directory_path, '.tar')
    }
    tryCatch({
      tar(tarfile=archive_name,
          files=directory_path,
          compression=archive_control[['archive_compression']])
      },
      error=function(cond) {
              message('problem writing the archive file \'', archive_name, '\': ', cond, appendLF=appendLR)
              return(NULL)
      },
      finally={
        message(paste0('save_monocle_objects made an archive file called \"', archive_name, '\"'))
      }
    ) # tryCatch
  } # if(archive_control...
}


#
#' Load a full Monocle3 cell_data_set.
#'
#' Load a full Monocle3 cell_data_set, which was saved using
#' save_monocle_objects. For more information read the help
#' information for save_monocle_objects.
#'
#' @param directory_path a string giving the name of the directory
#'   from which to read the saved cell_data_set files.
#' @param matrix_control a list that is used only to set the
#'   matrix path when the saved monocle objects has the counts matrix
#'   stored as a BPCells on-disk matrix. By default, the BPCells matrix
#'   directory path is set to the current working directory.
#' @return a cell_data_set.
#'
#' @examples
#'   \dontrun{
#'     cds <- load_a549()
#'     save_monocle_objects(cds, 'mo')
#'     cds1 <- load_monocle_objects('mo')
#'   }
#'
#' @export
# Bioconductor forbids writing to user directories so examples
# is not run.
load_monocle_objects <- function(directory_path, matrix_control=list(matrix_path='.')) {
  appendLF <- FALSE
  # Check for directory.
  if(!file.exists(directory_path))
    stop('Directory or file \'', directory_path, '\' does not exist.')

  # Check for file_index.rds. If none, assume that 'directory_path' is
  # an RDS file.
  file_index_path <- file.path(directory_path, 'file_index.rds')
  if(!file.exists(file_index_path)) {
    cds <- load_monocle_rds(directory_path)
    return(cds)
  }

  # Read file index.
  catch_error <- FALSE
  file_index <- tryCatch(
    {
      readRDS(file_index_path)
    },
    error = function(cond) {
              message('problem reading file \'', file_index_path, '\': ', cond, appendLF=appendLF);
              catch_error <<- TRUE
    },
    warning = function(cond) {
              message('problem reading file \'', file_index_path, '\': ', cond, appendLF=appendLF);
              catch_error <<- TRUE
    }
  )

  if(catch_error) {
    return(NULL)
  }

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
    reduction_method <- file_index[['files']][['reduction_method']][[ifile]]
    md5sum <- file_index[['files']][['file_md5sum']][[ifile]]

    if(!file.exists(file_path)) 
      stop('load_monocle_objects: missing file \'', file_path, '\'')

    #
    # The functions load_umap_nn_indexes and
    # loadHDF5SummarizedExperiment check md5sums
    # internally so don't check here.
    #
    if(!(cds_object == 'reduce_dim_aux' &&
         reduction_method == 'UMAP' &&
         file_format == 'umap_nn_index' &&
         nchar(md5sum) > 32) &&
       file_format != 'hdf5' &&
       cds_object != 'bpcells_matrix_dir') {
      md5sum_file <- tools::md5sum(file_path)
      if(is.na(md5sum_file) || (md5sum_file != md5sum)) {
        stop('md5sum mismatch for file \'', file_path, '\'')
      }
    }

    #
    # Note:
    #   o  expect that the RDS file for a reduction_method
    #      appears before index files for the reduction_method
    #
    if(cds_object == 'cds') {
      if(file_format == 'rds') {
        cds <- tryCatch(
          {
            readRDS(file_path)
          },
          error = function(cond) {
            message('problem reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
          })
      }
      else
      if(file_format == 'hdf5') {
        cds <- tryCatch(
          {
            HDF5Array::loadHDF5SummarizedExperiment(file_path)
          },
          error = function(cond) {
            message('problem reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
          })
      }
      else {
        stop('Unrecognized cds format value \'', file_format, '\'')
      }
    }
    else
    if(cds_object == 'reduce_dim_aux') {
      if(file_format == 'annoy_index') {
        cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy']] <- update_annoy_index(cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy']])

        metric <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy']][['nn_index']][['metric']]
        ncolumn <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy']][['nn_index']][['ncol']]
        cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy']][['nn_index']] <- tryCatch(
          {
            load_annoy_index(cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy']][['nn_index']], file_path, metric, ncolumn)
          },
          error = function(cond) {
            message('problem reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
          })
      }
      else
      if(file_format == 'hnsw_index') {
        cds@reduce_dim_aux[[reduction_method]][['nn_index']][['hnsw']] <- update_hnsw_index(cds@reduce_dim_aux[[reduction_method]][['nn_index']][['hnsw']])

        metric <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][['hnsw']][['nn_index']][['metric']]
        ncolumn <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][['hnsw']][['nn_index']][['ncol']]
        cds@reduce_dim_aux[[reduction_method]][['nn_index']][['hnsw']][['nn_index']] <- tryCatch(
          {
            load_hnsw_index(cds@reduce_dim_aux[[reduction_method]][['nn_index']][['hnsw']][['nn_index']], file_path, metric, ncolumn)
          },
          error = function(cond) {
            message('problem reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
          })
      }
      else
      if(reduction_method == 'UMAP' && file_format == 'umap_annoy_index') {
        cds@reduce_dim_aux[[reduction_method]][['model']][['umap_model']] <- tryCatch(
          {
            load_umap_nn_indexes(cds@reduce_dim_aux[[reduction_method]][['model']][['umap_model']], file_path, md5sum)
          },
          error = function(cond) {
            message('problem reading file \'', file_path, '\'', appendLF=appendLF)
            return(NULL)
         })
      }
      else {
        stop('Unrecognized file format value \'', file_format, '\'')
      }
      cds <- set_model_identity_path(cds, reduction_method, directory_path)
    }
    else
    if(cds_object == 'bpcells_matrix_dir') {
      if(!is.null(assay(cds, 'counts_row_order'))) {
        assay(cds, 'counts_row_order') <- NULL
      }
      counts(cds, bpcells_warn=FALSE ) <- tryCatch(
        {
          load_bpcells_matrix_dir(file_path, md5sum, matrix_control=matrix_control)
        },
        error = function(cond) {
          message('problem reading file \'', file_path, '\'', appendLF=appendLF)
          return(NULL)
        })
        # Rebuild the BPCells row-major order counts matrix.
        cds <- set_cds_row_order_matrix(cds=cds)
    }
    else {
      stop('Unrecognized cds_object value \'', cds_object, '\'')
    }
  }

  return(cds)
}


#
#  Operations from Monocle3 documentation
#   o  expression matrix
#   o  PCA: svd_v svd_sdev and gene_loadings and prop_var_expl
#   o  LSI: ?
#   o  Aligned: beta array and ?
#   o  cell clustering
#   o  marker genes
#   o  cell annotations
#   o  find_gene_modules
#   o  Garnett annotations?
#   o  learn_graph()
#   o  order_cells()
#   o  fit_models() DE regression analysis
#   o  graph_test() DE graph-autocorrelation analysis
#   o  Cicero annotations?
#
# Objects on loading (initial)
#   o  @ reduce_dim_aux
#   o  @ principal_graph_aux
#   o  @ principal_graph
#   o  @ clusters
#   o  @ int_elementMetadata
#   o  @ int_colData
#   o  @ int_metadata
#   o  @ rowRanges
#   o  @ colData
#   o  @ assays
#   o  @ NAMES
#   o  @ elementMetadata
#   o  @ metadata
#   o  $ preprocess_aux (vestigial now)
#
# Objects after analysis consisting of
#   cds <- preprocess_cds(cds, num_dim = 100)
#   cds <- align_cds(cds, num_dim = 100, alignment_group = "plate")
#   cds <- reduce_dimension(cds)
#   cds <- cluster_cells(cds, resolution=1e-5)
#   marker_test_res <- top_markers(cds, group_cells_by="partition", reference_cells=1000, cores=8)
#   colData(cds)$assigned_cell_type <- as.character(partitions(cds))
#   cds_subset <- choose_cells(cds)
#   pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="knn", cores=8)
#   pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))
#   gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=1e-3)
#   cds_subset <- cluster_cells(cds_subset, resolution=1e-2)
#
#   o  @ reduce_dim_aux
#   o  @ principal_graph_aux
#   o  @ principal_graph
#   o  @ clusters
#   o  @ int_elementMetadata
#   o  @ int_colData
#   o  @ int_metadata
#   o  @ rowRanges
#   o  @ colData
#   o  @ assays
#   o  @ NAMES
#   o  @ elementMetadata
#   o  @ metadata
#   o  $ preprocess_aux (empty)
#
# Objects after analysis of
#   cds <- load_monocle_rds('packer_embryo.load.rds')
#   cds <- preprocess_cds(cds, num_dim = 50)
#   cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")
#   cds <- reduce_dimension(cds)
#   cds <- cluster_cells(cds)
#   cds <- learn_graph(cds)
#   cds <- order_cells(cds)
#   ciliated_genes <- c("che-1", "hlh-17", "nhr-6", "dmd-6", "ceh-36", "ham-1")
#   cds_subset <- cds[rowData(cds)$gene_short_name %in% ciliated_genes,]
#   gene_fits <- fit_models(cds_subset, model_formula_str = "~embryo.time")
#   subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
#   pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
#
#   o  @ reduce_dim_aux
#   o  @ principal_graph_aux
#   o  @ principal_graph
#   o  @ clusters
#   o  @ int_elementMetadata
#   o  @ int_colData
#   o  @ int_metadata
#   o  @ rowRanges
#   o  @ colData
#   o  @ assays
#   o  @ NAMES
#   o  @ elementMetadata
#   o  @ metadata
#
# Strategy
#   o  copy slots
#        o  @ assays
#        o  @ colData
#        o  @ rowRanges
#        o  @ int_metadata ?
#        o  @ int_colData
#        o  @ int_elementMetadata
#        o  @ clusters
#        o  @ principal_graph
#        o  @ principal_graph_aux
#        o  @ reduce_dim_aux
#   o  transfer
#        o  @ preprocess_aux -> reduce_dim_aux (done)
#
# Issues
#   o  deal with missing values/objects; e.g., projection related and nearest neighbors
#        o  no annoy indices
#        o  no models initially; partial models after loading
#        o  no model identities
#        o  affected functions
#             o  projection functions
#             o  save/load models
#             o  save/load monocle objects?
#   o  gene_loadings to svd_v and svd_sdev conversion? (done)
#   o  the user may save a cds made using the recent Monocle3 version
#      and load it using load_monocle_objects. In that case, load the
#      cds using readRDS and return it unmodified. Use a model_version
#      check?
#   o  can there be an inconsistency between the model and matrix identity
#      version information?

load_monocle_rds <- function(file_path) {
  appendLF <- TRUE
  catch_error <- FALSE
  cds_tmp <- tryCatch(
                       {
                         readRDS(file_path)
                       },
                       error=function(cond) {
                         message('problem reading file \'', file_path, '\': ', cond, appendLF=appendLF);
                         catch_error <<- TRUE
                       },
                       warning=function(cond) {
                         message('problem reading file \'', file_path, '\': ', cond, appendLF=appendLF);
                         catch_error <<- TRUE
                       }
                     )

  if(catch_error) {
    return(NULL)
  }

  cds <- cds_tmp
  if(!is.null(SingleCellExperiment::reducedDims(cds_tmp)[['PCA']])) {
    if(is.null(cds_tmp@reduce_dim_aux[['PCA']][['model']][['identity']])) {
      cds <- initialize_reduce_dim_model_identity(cds, 'PCA')
    }
    if(!is.null(attr(cds, which='preprocess_aux')) && !is.null(cds_tmp@preprocess_aux[['gene_loadings']])) {
      gene_loadings <- cds_tmp@preprocess_aux[['gene_loadings']]
      svd_sdev <- apply(gene_loadings, 2, function(x) {sqrt(sum(x^2))})
      svd_v <- gene_loadings %*% diag(1.0/svd_sdev)
      colnames(svd_v) <- paste("PC", seq(1, ncol(svd_v)), sep="")
      cds@reduce_dim_aux[['PCA']][['model']][['svd_sdev']] <- svd_sdev
      cds@reduce_dim_aux[['PCA']][['model']][['svd_v']] <- svd_v
      cds@reduce_dim_aux[['PCA']][['model']][['prop_var_expl']] <- cds_tmp@preprocess_aux[['prop_var_expl']]
    }
  }

  if(!is.null(SingleCellExperiment::reducedDims(cds_tmp)[['LSI']])) {
    if(is.null(cds_tmp@reduce_dim_aux[['LSI']][['model']][['identity']])) {
      cds <- initialize_reduce_dim_model_identity(cds, 'LSI')
    }
    if(!is.null(attr(cds, which='preprocess_aux')) && !is.null(cds_tmp@preprocess_aux[['gene_loadings']])) {
      gene_loadings <- cds_tmp@preprocess_aux[['gene_loadings']]
      svd_sdev <- apply(gene_loadings, 2, function(x) {sqrt(sum(x^2))})
      cds@reduce_dim_aux[['LSI']][['model']][['svd_v']] <- gene_loadings %*% diag(1.0/svd_sdev)
      cds@reduce_dim_aux[['LSI']][['model']][['svd_sdev']] <- svd_sdev
    }
  }

  if(!is.null(SingleCellExperiment::reducedDims(cds_tmp)[['Aligned']])) {
    if(is.null(cds_tmp@reduce_dim_aux[['Aligned']][['model']][['identity']])) {
      cds <- initialize_reduce_dim_model_identity(cds, 'Aligned')
    }
    if(!is.null(attr(cds, which='preprocess_aux')) && !is.null(cds_tmp@preprocess_aux$beta)) {
      cds@reduce_dim_aux[['Aligned']][['model']][['beta']] <- cds_tmp@preprocess_aux$beta
    }
  }

  if(!is.null(SingleCellExperiment::reducedDims(cds_tmp)[['tSNE']])) {
    if(is.null(cds_tmp@reduce_dim_aux[['tSNE']][['model']][['identity']])) {
      cds <- initialize_reduce_dim_model_identity(cds, 'tSNE')
    }
  }

  if(!is.null(SingleCellExperiment::reducedDims(cds_tmp)[['UMAP']])) {
    if(is.null(cds_tmp@reduce_dim_aux[['UMAP']][['model']][['identity']])) {
      cds <- initialize_reduce_dim_model_identity(cds, 'UMAP')
    }
  }

  if(is.null(SingleCellExperiment::int_metadata(cds_tmp)[['counts_metadata']])) {
    cds <- initialize_counts_metadata(cds)
  }

  # The RDS file may have a preprocess_aux slot, which doesn't
  # exist in the current cell_data_set class. The attr command
  # appears to remove it whereas setting cds@preprocess_cds
  # and cds$preprocess_cds to NULL do not.
  if(!is.null(attr(cds, which='preprocess_aux'))) {
    attr(cds, which='preprocess_aux') <- NULL
  }

  rm(cds_tmp)

  return(cds)
}


