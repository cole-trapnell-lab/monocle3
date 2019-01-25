setOldClass(c("igraph"), prototype=structure(list(), class="igraph"))

#' The cell_data_set class
#'
#' The main class used by Monocle3 to hold single cell expression data.
#' cell_data_set extends the basic Bioconductor ExpressionSet class.
#'
#' This class is initialized from a matrix of expression values Methods that
#' operate on cell_data_set objects constitute the basic Monocle workflow.
#'
#'
#' @field reducedDimS Matrix of class numeric, containing the source values computed by Independent Components Analysis.
#' @field reducedDimW Matrix of class numeric, containing the whitened expression values computed during Independent Components Analysis.
#' @field reducedDimA Matrix of class numeric, containing the weight values computed by Independent Components Analysis.
#' @field reducedDimK A Matrix of class numeric, containing the pre-whitening matrix computed by Independent Components Analysis.
#' @field principal_graph An Object of class igraph, containing the minimum spanning tree used by Monocle to order cells according to progress through a biological process.
#' @field cellPairwiseDistances A Matrix of class numeric, containing the pairwise distances between cells in the reduced dimension space.
#' @field expression_family An Object of class character, specifying the expression family function used for expression responses.
#' @field lower_detection_limit A numeric value specifying the minimum expression level considered to be true expression.
#' @field disp_fit_info An environment containing lists, one for each set of estimated dispersion values. See estimateDispersions.
#' @field dim_reduce_type A string encoding how this cell_data_set has been reduced in dimensionality
#' @field rge_method A string encoding how this cell_data_set has been fitted with a principal graph
#' @field aux_ordering_data An environment of auxilliary data structures used by various steps in Monocle. Not to be accessed by users directly.
#' @name cell_data_set
#' @rdname cell_data_set
#' @aliases cell_data_set-class
#' @exportClass cell_data_set
#' @importFrom Biobase ExpressionSet pData fData exprs pData<- fData<-
setClass( "cell_data_set",
          contains = "ExpressionSet",
          slots = c(reducedDimS = "matrix",
                    reducedDimW = "matrix",
                    reducedDimA = "matrix",
                    reducedDimK = "matrix",
                    principal_graph="igraph",
                    #cellPairwiseDistances="matrix",
                    expression_family="character",
                    lower_detection_limit="numeric",
                    disp_fit_info = "environment",
                    dim_reduce_type="character",
                    #rge_method="character",
                    aux_ordering_data = "environment",
                    aux_clustering_data = "environment",
                    normalized_data_projection = "matrix"
          ),
          prototype = prototype( new( "VersionedBiobase",
                                      versions = c( Biobase::classVersion("ExpressionSet"), cell_data_set = "1.2.0" ) ))
)


#' Creates a new cell_data_set object.
#'
#' @param cellData expression data matrix for an experiment
#' @param phenoData data frame containing attributes of individual cells
#' @param featureData data frame containing attributes of features (e.g. genes)
#' @param lower_detection_limit the minimum expression level that consistitutes true expression
#' @param expression_family character, the expression family function to be used for expression response variables
#' @return a new cell_data_set object
#' @export
#' @examples
#' \dontrun{
#' sample_sheet_small <- read.delim("../data/sample_sheet_small.txt", row.names=1)
#' sample_sheet_small$Time <- as.factor(sample_sheet_small$Time)
#' gene_annotations_small <- read.delim("../data/gene_annotations_small.txt", row.names=1)
#' fpkm_matrix_small <- read.delim("../data/fpkm_matrix_small.txt")
#' pd <- new("AnnotatedDataFrame", data = sample_sheet_small)
#' fd <- new("AnnotatedDataFrame", data = gene_annotations_small)
#' HSMM <- new("cell_data_set", exprs = as.matrix(fpkm_matrix_small), phenoData = pd, featureData = fd)
#' }
new_cell_data_set <- function(cellData,
                            phenoData = NULL,
                            featureData = NULL,
                            lower_detection_limit = 0.1,
                            expression_family="negbinomial.size")
{

  if(!('gene_short_name' %in% colnames(featureData))) {
    warning("Warning: featureData must contain a column verbatim named 'gene_short_name' for certain functions")
  }

  if (class(cellData) != "matrix" && is_sparse_matrix(cellData) == FALSE){
    stop("Error: argument cellData must be a matrix (either sparse from the Matrix package or dense)")
  }

  if(!('gene_short_name' %in% colnames(featureData))) {
    warning("Warning: featureData must contain a column verbatim named 'gene_short_name' for certain functions")
  }

  sizeFactors <- rep( NA_real_, ncol(cellData) )


  if( is.null( phenoData ) )
    phenoData <- annotatedDataFrameFrom( cellData, byrow=FALSE )
  if( is.null( featureData ) )
    featureData <- annotatedDataFrameFrom(cellData, byrow=TRUE)

  if(!('gene_short_name' %in% colnames(featureData))) {
    warning("Warning: featureData must contain a column verbatim named 'gene_short_name' for certain functions")
  }

  phenoData$`Size_Factor` <- sizeFactors

  cds <- new( "cell_data_set",
              assayData = assayDataNew( "environment", exprs=cellData ),
              phenoData=phenoData,
              featureData=featureData,
              lower_detection_limit=lower_detection_limit,
              expression_family=expression_family,
              disp_fit_info = new.env( hash=TRUE ))

  validObject( cds )
  cds
}
