% Generated manually for S4 method documentation
\docType{methods}
\name{monocle3-s4-methods}
\alias{counts<-,SingleCellExperiment-method}
\alias{saveRDS,IterableMatrix-method}
\alias{saveRDS,cell_data_set-method}
\title{S4 methods in monocle3}
\description{
Monocle3 provides S4 methods for several generics to support its
core classes and interoperability with Bioconductor classes.
}
\section{Methods}{
\describe{
\item{\code{counts<-(object = "SingleCellExperiment")}}{
Sets the \code{counts} assay and emits a warning when assigning a BPCells
matrix, because the row-order matrix must also be updated.
}
\item{\code{saveRDS(object = "IterableMatrix")}}{
Calls \code{base::saveRDS} and warns that BPCells out-of-core matrices are
not preserved in the serialized object.
}
\item{\code{saveRDS(object = "cell_data_set")}}{
Calls \code{base::saveRDS} and warns that nearest-neighbor indices and
out-of-core matrices are not preserved; use \code{save_monocle_objects}
to persist all assets.
}
}
}
