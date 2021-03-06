% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_cellranger_data.R
\name{load_cellranger_data}
\alias{load_cellranger_data}
\title{Load data from the 10x Genomics Cell Ranger pipeline}
\usage{
load_cellranger_data(
  pipestance_path = NULL,
  genome = NULL,
  barcode_filtered = TRUE,
  umi_cutoff = 100
)
}
\arguments{
\item{pipestance_path}{Path to the output directory produced by Cell Ranger}

\item{genome}{The desired genome (e.g., 'hg19' or 'mm10')}

\item{barcode_filtered}{Load only the cell-containing barcodes}

\item{umi_cutoff}{Numeric, desired cutoff to include a cell. Default is 100.}
}
\value{
a new cell_data_set object
}
\description{
Loads cellranger data into a cell_data_set object.  Note that if your dataset
is from version 3.0 and contains non-Gene-Expression data (e.g. Antibodies
or CRISPR features), only the Gene Expression data is returned.
}
\details{
\itemize{
\item the \emph{pipestance_path} argument takes the name of a Cell Ranger
output directory, in which it looks for the required data files,
for example, \emph{pipestance_path=10x_data}
\item for Cell Ranger version 2 data, \emph{load_cellranger_data} expects to
find the required files \emph{barcodes.tsv}, \emph{genes.tsv}, and
\emph{matrix.mtx}
in the directories as
\itemize{
\item \emph{10x_data/outs/filtered_gene_bc_matrices/<genome>/barcodes.tsv}
\item \emph{10x_data/outs/filtered_gene_bc_matrices/<genome>/genes.tsv}
\item \emph{10x_data/outs/filtered_gene_bc_matrices/<genome>/matrix.mtx}
}

where <genome> is the name of a genome. \emph{load_cellranger_data}
expects to find either a single \emph{genome} directory in
\emph{10x_data/outs/filtered_gene_bc_matrices} or a \emph{genome}
directory with the name given with the \emph{genome} argument.
\item for Cell Ranger version 3 data, \emph{load_cellranger_data} expects to
find the required files \emph{barcodes.tsv.gz}, \emph{features.tsv.gz},
and \emph{matrix.mtx.gz} in the directories as
\itemize{
\item \emph{10x_data/outs/filtered_feature_bc_matrix/barcodes.tsv.gz}
\item \emph{10x_data/outs/filtered_feature_bc_matrix/features.tsv.gz}
\item \emph{10x_data/outs/filtered_feature_bc_matrix/matrix.mtx.gz}
}
\item if any of the files is not in the expected directory,
\emph{load_cellranger_data} will terminate with an error
}
}
\examples{
\dontrun{
# Load from a Cell Ranger output directory
gene_bc_matrix <- load_cellranger_data("/home/user/cellranger_output")
}
}
