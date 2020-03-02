setOldClass(c("igraph"), prototype=structure(list(), class="igraph"))

#' The cell_data_set class
#'
#' The main class used by Monocle3 to hold single-cell expression data.
#' cell_data_set extends the Bioconductor SingleCellExperiment class.
#'
#' This class is initialized from a matrix of expression values along with cell
#' and feature metadata.
#'
#' @field preprocess_aux SimpleList, auxiliary information from preprocessing.
#' @field reduce_dim_aux SimpleList, auxiliary information from reduced
#'   dimension.
#' @field principal_graph_aux SimpleList, auxiliary information from principal
#'   graph construction
#' @field principal_graph SimpleList of igraph objects containing principal
#'   graphs for different dimensionality reduction.
#' @field clusters SimpleList of cluster information for different
#'   dimensionality reduction.
#' @name cell_data_set
#' @rdname cell_data_set
#' @aliases cell_data_set-class
#' @exportClass cell_data_set
#' @importFrom SingleCellExperiment SingleCellExperiment colData rowData
#' @importFrom SingleCellExperiment reducedDim<- reducedDim reducedDims<-
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom SummarizedExperiment Assays colData<- rowData<- assays assays<-
#' @importFrom S4Vectors metadata metadata<- SimpleList
setClass("cell_data_set",
          contains = c("SingleCellExperiment"),
          slots = c(preprocess_aux = "SimpleList",
                    reduce_dim_aux = "SimpleList",
                    principal_graph_aux="SimpleList",
                    principal_graph = "SimpleList",
                    clusters = "SimpleList")
)


#' Create a new cell_data_set object.
#'
#' @param expression_data expression data matrix for an experiment, can be a
#'   sparseMatrix.
#' @param cell_metadata data frame containing attributes of individual cells,
#'  where \code{row.names(cell_metadata) = colnames(expression_data)}.
#' @param gene_metadata data frame containing attributes of features
#'   (e.g. genes), where
#'   \code{row.names(gene_metadata) = row.names(expression_data)}.
#' @return a new cell_data_set object
#' @export
#' @examples
#' small_a549_colData_df <- readRDS(system.file("extdata",
#'                                              "small_a549_dex_pdata.rda",
#'                                              package = "monocle3"))
#' small_a549_rowData_df <- readRDS(system.file("extdata",
#'                                              "small_a549_dex_fdata.rda",
#'                                              package = "monocle3"))
#' small_a549_exprs <- readRDS(system.file("extdata",
#'                                         "small_a549_dex_exprs.rda",
#'                                         package = "monocle3"))
#' small_a549_exprs <- small_a549_exprs[,row.names(small_a549_colData_df)]
#'
#' cds <- new_cell_data_set(expression_data = small_a549_exprs,
#'                          cell_metadata = small_a549_colData_df,
#'                          gene_metadata = small_a549_rowData_df)
#'
new_cell_data_set <- function(expression_data,
                              cell_metadata = NULL,
                              gene_metadata = NULL) {

  assertthat::assert_that(class(expression_data) == "matrix" ||
                            is_sparse_matrix(expression_data),
                          msg = paste("Argument expression_data must be a",
                                      "matrix - either sparse from the",
                                      "Matrix package or dense"))
  if (!is.null(cell_metadata)) {
    assertthat::assert_that(nrow(cell_metadata) == ncol(expression_data),
                            msg = paste("cell_metadata must be NULL or have",
                                        "the same number of rows as columns",
                                        "in expression_data"))
    assertthat::assert_that(!is.null(row.names(cell_metadata)) &
      all(row.names(cell_metadata) == colnames(expression_data)),
      msg = paste("row.names of cell_metadata must be equal to colnames of",
                  "expression_data"))
  }

  if (!is.null(gene_metadata)) {
    assertthat::assert_that(nrow(gene_metadata) == nrow(expression_data),
                            msg = paste("gene_metadata must be NULL or have",
                                        "the same number of rows as rows",
                                        "in expression_data"))
    assertthat::assert_that(!is.null(row.names(gene_metadata)) & all(
      row.names(gene_metadata) == row.names(expression_data)),
      msg = paste("row.names of gene_metadata must be equal to row.names of",
                  "expression_data"))
  }

  if (is.null(cell_metadata)) {
    cell_metadata <- data.frame(cell = colnames(expression_data),
                                row.names = colnames(expression_data))
  }

  if(!('gene_short_name' %in% colnames(gene_metadata))) {
    warning(paste("Warning: gene_metadata must contain a column verbatim",
                  "named 'gene_short_name' for certain functions."))
  }

  sce <- SingleCellExperiment(list(counts=methods::as(expression_data, "dgCMatrix")),
                              rowData = gene_metadata,
                              colData = cell_metadata)

  cds <- methods::new("cell_data_set",
             assays = SummarizedExperiment::Assays(
               list(counts=methods::as(expression_data, "dgCMatrix"))),
             colData = colData(sce),
             int_elementMetadata =int_elementMetadata(sce),
             int_colData = int_colData(sce),
             int_metadata = int_metadata(sce),
             metadata = metadata(sce),
             NAMES = NULL,
             elementMetadata = elementMetadata(sce)[,0],
             rowRanges = rowRanges(sce))

  metadata(cds)$cds_version <- Biobase::package.version("monocle3")
  clusters <- stats::setNames(SimpleList(), character(0))
  cds <- estimate_size_factors(cds)
  cds
}
