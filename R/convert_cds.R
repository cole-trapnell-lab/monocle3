#' Import a seurat object and convert it to a monocle cds.
#'
#' This function takes a Seurat object and converts it to a monocle3 cds.
#'
#' @param seurat_cds the Seurat object you would like to convert into a monocle3 cds
#' @param import_all logical, whether to import all the slots in seurat.
#'   Default is FALSE (or only keep minimal dataset).
#' @return a new monocle3 cell_data_set
#' @export
#' @examples
seurat_to_cds <- function(seurat_cds, import_all = FALSE) {
  if(!class(seurat_cds)[1] == 'seurat')
    stop("Input must be a Seurat object")

  #requireNamespace("Seurat")
  data <- seurat_cds@raw.data

  if(class(data) == "data.frame") {
    data <- as(as.matrix(data), "sparseMatrix")
  }

  pd <- tryCatch( {
    pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
    pd
  },
  error = function(e) {
    pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
    pd <- new("AnnotatedDataFrame", data = pData)

    message("This Seurat object doesn't provide any meta data");
    pd
  })

  # remove filtered cells from Seurat
  if(length(setdiff(colnames(data), rownames(pd))) > 0) {
    data <- data[, rownames(pd)]
  }

  fData <- data.frame(gene_short_name = row.names(data),
                      row.names = row.names(data))
  fd <- new("AnnotatedDataFrame", data = fData)
  lowerDetectionLimit <- otherCDS@is.expr

  if(all(data == floor(data))) {
    expression_family <- "negbinomial.size"
  } else if(any(data < 0)){
    expression_family <- "uninormal"
  } else {
    expression_family <- "tobit"
  }

  valid_data <- data[, row.names(pd)]

  monocle_cds <- new_cell_data_set(data,
                                phenoData = pd,
                                featureData = fd,
                                lower_detection_limit=lowerDetectionLimit,
                                expression_family=expression_family)

  if(import_all) {
    if("Monocle" %in% names(otherCDS@misc)) {
      # if(slotNames(lung) == )
      # monocle_cds@reducedDimS = otherCDS@misc$Monocle@reducedDimS
      # monocle_cds@reducedDimW = otherCDS@misc$Monocle@reducedDimW
      # monocle_cds@reducedDimA = otherCDS@misc$Monocle@reducedDimA
      # monocle_cds@reducedDimK = otherCDS@misc$Monocle@reducedDimK
      # monocle_cds@minSpanningTree = otherCDS@misc$Monocle@minSpanningTree
      # monocle_cds@cellPairwiseDistances = otherCDS@misc$Monocle@cellPairwiseDistances
      # monocle_cds@expression_family = otherCDS@misc$Monocle@expression_family
      # monocle_cds@disp_fit_info = otherCDS@misc$Monocle@disp_fit_info
      # monocle_cds@dim_reduce_type = otherCDS@misc$Monocle@dim_reduce_type
      # monocle_cds@auxOrderingData = otherCDS@misc$Monocle@auxOrderingData
      # monocle_cds@auxClusteringData = otherCDS@misc$Monocle@auxClusteringData
      # monocle_cds@experimentData = otherCDS@misc$Monocle@experimentData
      # monocle_cds@classVersion = otherCDS@misc$Monocle@.__classVersion__
      # monocle_cds@annotation = otherCDS@misc$Monocle@annotation
      # monocle_cds@protocolData = otherCDS@misc$Monocle@protocolData
      # monocle_cds@featureData = otherCDS@misc$Monocle@featureData

      # clean all conversion related slots
      otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
      otherCDS@misc$Monocle@auxClusteringData$scran <- NULL

      monocle_cds <- otherCDS@misc$Monocle
      mist_list <- otherCDS

    } else {
      # mist_list <- list(ident = ident,
      #                   project.name = project.name,
      #                   dr = otherCDS@dr,
      #                   assay = otherCDS@assay,
      #                   hvg.info = otherCDS@hvg.info,
      #                   imputed = otherCDS@imputed,
      #                   cell.names = otherCDS@cell.names,
      #                   cluster.tree = otherCDS@cluster.tree,
      #                   snn = otherCDS@snn,
      #                   kmeans = otherCDS@kmeans,
      #                   spatial = otherCDS@spatial,
      #                   misc = otherCDS@misc
      # )
      mist_list <- otherCDS
    }
  } else {
    mist_list <- list()
  }

  if("var.genes" %in% slotNames(otherCDS)) {
    monocle_cds <- set_ordering_filter(monocle_cds, otherCDS@var.genes)

  }
  monocle_cds@aux_clustering_data$seurat <- mist_list

return(monocle_cds)
}
