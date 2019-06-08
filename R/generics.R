#' Generic to extract pseudotime from CDS object
#'
#' @param x A cell_data_set object.
#' @param reduction_method Reduced dimension to extract pseudotime for.
#'
#' @export
setGeneric("pseudotime", function(x, reduction_method = "UMAP")
  standardGeneric("pseudotime"))

#' Method to extract pseudotime from CDS object
#' @param cell_data_set A cell_data_set object.
#'
#' @export
setMethod("pseudotime", "cell_data_set",
          function(x, reduction_method = "UMAP") {
            value <- x@principal_graph_aux[[
              reduction_method]]$pseudotime[colnames(exprs(x))]
            if (is.null(value)) {
              stop(paste0("No pseudotime calculated for reduction_method = ",
                          reduction_method, ". Please first run ",
                          "order_cells with reduction_method = ",
                          reduction_method, "."))
            }
            return(value)
          })

#' Generic to extract clusters from CDS object
#'
#' @param x A cell_data_set object.
#' @param reduction_method Reduced dimension to extract clusters for.
#' @export
setGeneric("clusters", function(x, reduction_method = "UMAP")
  standardGeneric("clusters"))

#' Method to extract clusters from CDS object
#' @param cell_data_set A cell_data_set object.
#'
#' @export
setMethod("clusters", "cell_data_set",
          function(x, reduction_method = "UMAP") {
            value <- x@clusters[[
              reduction_method]]$clusters[colnames(exprs(x))]
            if (is.null(value)) {
              stop(paste0("No clusters calculated for reduction_method = ",
                          reduction_method, ". Please first run ",
                          "cluster_cells with reduction_method = ",
                          reduction_method, "."))
            }
            return(value)
          })

#' Generic to extract partitions from CDS object
#'
#' @param x A cell_data_set object.
#' @param reduction_method Reduced dimension to partitions clusters for.
#'
#' @export
setGeneric("partitions", function(x, reduction_method = "UMAP")
  standardGeneric("partitions"))

#' Method to extract partitions from CDS object
#' @param cell_data_set A cell_data_set object.
#'
#' @export
setMethod("partitions", "cell_data_set",
          function(x, reduction_method = "UMAP") {
            value <- x@clusters[[
              reduction_method]]$partitions[colnames(exprs(x))]
            if (is.null(value)) {
              stop(paste0("No partitions calculated for reduction_method = ",
                          reduction_method, ". Please first run ",
                          "cluster_cells with reduction_method = ",
                          reduction_method, "."))
            }
            return(value)
          })

#' Generic to extract principal graph from CDS
#' @param x A cell_data_set object.
#'
#' @export
setGeneric("principal_graph", function(x) standardGeneric("principal_graph"))

#' Generic to set principal graph to CDS
#' @param x A cell_data_set object.
#' @param value A principal graph object.
#'
#' @export
setGeneric("principal_graph<-", function(x, value)
  standardGeneric("principal_graph<-"))

#' Method to extract principal graph from CDS
#' @param x A cell_data_set object.
#'
#' @export
setMethod("principal_graph", "cell_data_set", function(x) {
  value <- x@principal_graph
  return(value)
})

#' Generic to set principal graph to CDS
#' @param x A cell_data_set object.
#' @param value A principal graph object.
#'
#' @export
#' @importClassesFrom S4Vectors List
setReplaceMethod("principal_graph", "cell_data_set", function(x, value) {
  value <- methods::as(value, "List")
  if (is.null(names(value))) {
    names(value) <- character(length(value))
  }
  x@principal_graph <- value
  validObject(x)
  return(x)
})

#' Generic to extract principal graph auxilliary information from CDS
#' @param x A cell_data_set object.
#' @export
setGeneric("principal_graph_aux", function(x)
  standardGeneric("principal_graph_aux"))

#' Generic to set principal graph auxilliary information into CDS
#' @param x A cell_data_set object.
#' @param value A SimpleList of principal graph auxilliary information.
#' @export
setGeneric("principal_graph_aux<-", function(x, value)
  standardGeneric("principal_graph_aux<-"))

#' Method to extract principal graph auxilliary information from CDS
#' @param x A cell_data_set object.
#' @export
#' @export
setMethod("principal_graph_aux", "cell_data_set", function(x) {
  value <- x@principal_graph_aux
  return(value)
})

#' Method to set principal graph auxilliary information into CDS
#' @param x A cell_data_set object.
#' @param value A SimpleList of principal graph auxilliary information.
#' @export
#' @importClassesFrom S4Vectors List
setReplaceMethod("principal_graph_aux", "cell_data_set", function(x, value) {
  value <- methods::as(value, "List")
  if (is.null(names(value))) {
    names(value) <- character(length(value))
  }
  x@principal_graph_aux <- value
  validObject(x)
  return(x)
})


# Set of wrappers for easy transition from monocle.

#' Generic to access cds count matrix
#' @param x A cell_data_set object.
#'
#' @export
setGeneric("exprs", function(x) standardGeneric("exprs"))

#' Method to access cds count matrix
#' @param x A cell_data_set object.
#'
#' @export
setMethod("exprs", "cell_data_set", function(x) {
  value <- assays(x)$counts
  return(value)
})

#' Generic to access cds colData table
#' @param x A cell_data_set object.
#'
#' @export
setGeneric("pData", function(x) standardGeneric("pData"))

#' Generic to set cds colData table
#' @param x A cell_data_set object.
#' @param value A data frame to set to colData table.
#'
#' @export
setGeneric("pData<-", function(x, value) standardGeneric("pData<-"))

#' Method to access cds colData table
#' @param x A cell_data_set object.
#'
#' @export
setMethod("pData", "cell_data_set", function(x) {
  value <- colData(x)
  return(value)
})

#' Method to set cds colData table
#' @param x A cell_data_set object.
#' @param value A data frame to set to colData table.
#'
#' @export
#' @importClassesFrom S4Vectors List
setReplaceMethod("pData", "cell_data_set", function(x, value) {
  colData(x) <- value
  validObject(x)
  return(x)
})

#' Generic to access cds rowData table
#' @param x A cell_data_set object.
#'
#' @export
setGeneric("fData", function(x) standardGeneric("fData"))

#' Generic to set cds rowData table
#' @param x A cell_data_set object.
#' @param value A data frame to set to colData table.
#'
#' @export
setGeneric("fData<-", function(x, value) standardGeneric("fData<-"))

#' Generic to access cds rowData table
#' @param x A cell_data_set object.
#'
#' @export
setMethod("fData", "cell_data_set", function(x) {
  value <- rowData(x)
  return(value)
})

#' Method to set cds rowData table
#' @param x A cell_data_set object.
#' @param value A data frame to set to colData table.
#'
#' @export
#' @importClassesFrom S4Vectors List
setReplaceMethod("fData", "cell_data_set", function(x, value) {
  rowData(x) <- value
  validObject(x)
  return(x)
})
