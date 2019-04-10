setOldClass(c("igraph"), prototype=structure(list(), class="igraph"))

#' The cell_data_set class
#'
#' The main class used by Monocle3 to hold single cell expression data.
#' cell_data_set extends the basic Bioconductor SingleCellExperiment class.
#'
#' This class is initialized from a matrix of expression values Methods that
#' operate on cell_data_set objects constitute the basic Monocle workflow.
#'
#' @field principal_graph_aux
#' @field principal_graph
#' @field partitions
#' @name cell_data_set
#' @rdname cell_data_set
#' @aliases cell_data_set-class
#' @exportClass cell_data_set
#' @importFrom SingleCellExperiment SingleCellExperiment colData rowData
#' @importFrom SingleCellExperiment reducedDim<- reducedDim reducedDims<-
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom SummarizedExperiment Assays colData<- rowData<- assays assays<-
#' @importFrom S4Vectors metadata metadata<- SimpleList
setClass( "cell_data_set",
          contains = c("SingleCellExperiment"),
          slots = c(principal_graph_aux="SimpleList",
                    principal_graph = "SimpleList",
                    partitions = "SimpleList")
)


#' Creates a new cell_data_set object.
#'
#' @param expression_data expression data matrix for an experiment
#' @param cell_metadata data frame containing attributes of individual cells
#' @param gene_metadata data frame containing attributes of features (e.g. genes)
#' @param lower_detection_limit the minimum expression level that consistitutes
#'   true expression
#' @param expression_family character, the expression family function to be
#'   used for expression response variables
#' @return a new cell_data_set object
#' @export
#' @examples
#' \dontrun{
#' sample_sheet_small <- read.delim("../data/sample_sheet_small.txt",
#'                                  row.names=1)
#' sample_sheet_small$Time <- as.factor(sample_sheet_small$Time)
#' gene_annotations_small <- read.delim("../data/gene_annotations_small.txt",
#'                                       row.names=1)
#' fpkm_matrix_small <- read.delim("../data/fpkm_matrix_small.txt")
#' pd <- new("AnnotatedDataFrame", data = sample_sheet_small)
#' fd <- new("AnnotatedDataFrame", data = gene_annotations_small)
#' HSMM <- new("cell_data_set", exprs = as.matrix(fpkm_matrix_small),
#'             phenoData = pd, featureData = fd)
#' }
new_cell_data_set <- function(expression_data,
                              cell_metadata = NULL,
                              gene_metadata = NULL,
                              lower_detection_limit = 0.1,
                              expression_family="negbinomial.size") {

  if(!('gene_short_name' %in% colnames(gene_metadata))) {
    warning(paste("Warning: gene_metadata must contain a column verbatim named",
                  "'gene_short_name' for certain functions"))
  }

  if (class(expression_data) != "matrix" && is_sparse_matrix(expression_data) == FALSE){
    stop(paste("Error: argument expression_data must be a matrix (either sparse from",
               "the Matrix package or dense)"))
  }

  sce <- SingleCellExperiment(list(exprs=expression_data),
                              rowData = gene_metadata,
                              colData = cell_metadata)

  cds <- new("cell_data_set",
             assays = SummarizedExperiment::Assays(list(exprs=expression_data)),
             colData = colData(sce),
             int_elementMetadata =sce@int_elementMetadata,
             int_colData = sce@int_colData,
             int_metadata = sce@int_metadata,
             metadata = sce@metadata,
             NAMES = sce@NAMES,
             elementMetadata = sce@elementMetadata,
             rowRanges = sce@rowRanges)

  metadata(cds)$lower_detection_limit <- lower_detection_limit
  metadata(cds)$cds_version <- Biobase::package.version("monocle3")
  partitions <- setNames(SimpleList(), character(0))
  cds
}
