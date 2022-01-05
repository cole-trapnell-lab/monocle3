# Contributed by Nigel Delaney with adaptations

#' Get a genome from Cell Ranger output
#'
#' @param matrix_path Path to a matrices directory produced by the Cell Ranger
#'   pipeline
#' @param genome Genome to specifically check for, otherwise will check for
#'   whatever genome(s) exist there
#' @return A string representing the genome found
get_genome_in_matrix_path <- function(matrix_path, genome=NULL) {
  genomes <- dir(matrix_path)
  if (is.null(genome)) {
    if (length(genomes) == 1) {
      genome <- genomes[1]
    } else {
      stop(sprintf(paste("Multiple genomes found; please specify one. \n",
                         "Genomes present: %s"),
                   paste(genomes, collapse=", ")))
    }
  } else if (!(genome %in% genomes)) {
    stop(sprintf("Could not find specified genome: '%s'. Genomes present: %s",
                 genome,paste(genomes, collapse=", ")))
  }
  return(genome)
}

#' Load data from the 10x Genomics Cell Ranger pipeline
#'
#' Loads cellranger data into a cell_data_set object.  Note that if your dataset
#' is from version 3.0 and contains non-Gene-Expression data (e.g. Antibodies
#' or CRISPR features), only the Gene Expression data is returned.
#'
#' @details
#' *  the \emph{pipestance_path} argument takes the name of a Cell Ranger
#'     output directory, in which it looks for the required data files,
#'     for example, \emph{pipestance_path=10x_data}
#' *  for Cell Ranger version 2 data, \emph{load_cellranger_data} expects to
#'    find the required files \emph{barcodes.tsv}, \emph{genes.tsv}, and
#'    \emph{matrix.mtx}
#'    in the directories as
#'      -  \emph{10x_data/outs/filtered_gene_bc_matrices/<genome>/barcodes.tsv}
#'      -  \emph{10x_data/outs/filtered_gene_bc_matrices/<genome>/genes.tsv}
#'      -  \emph{10x_data/outs/filtered_gene_bc_matrices/<genome>/matrix.mtx}
#'
#'    where <genome> is the name of a genome. \emph{load_cellranger_data}
#'    expects to find either a single \emph{genome} directory in
#'     \emph{10x_data/outs/filtered_gene_bc_matrices} or a \emph{genome}
#'    directory with the name given with the \emph{genome} argument.
#' *  for Cell Ranger version 3 data, \emph{load_cellranger_data} expects to
#'    find the required files \emph{barcodes.tsv.gz}, \emph{features.tsv.gz},
#'    and \emph{matrix.mtx.gz} in the directories as
#'      -  \emph{10x_data/outs/filtered_feature_bc_matrix/barcodes.tsv.gz}
#'      -  \emph{10x_data/outs/filtered_feature_bc_matrix/features.tsv.gz}
#'      -  \emph{10x_data/outs/filtered_feature_bc_matrix/matrix.mtx.gz}
#'
#' * if any of the files is not in the expected directory,
#'   \emph{load_cellranger_data} will terminate with an error
#'
#' @param pipestance_path Path to the output directory produced by Cell Ranger
#' @param genome The desired genome (e.g., 'hg19' or 'mm10')
#' @param barcode_filtered Load only the cell-containing barcodes
#' @param umi_cutoff Numeric, desired cutoff to include a cell. Default is 100.
#' @return a new cell_data_set object
#'
#' @examples
#'   \donttest{
#'     cell_ranger_data <- system.file("extdata", "cell_ranger_3", package = "monocle3")
#'     gene_bc_matrix <- load_cellranger_data(cell_ranger_data)
#'   }
#'
#' @export
load_cellranger_data <- function(pipestance_path=NULL, genome=NULL,
                                 barcode_filtered=TRUE, umi_cutoff = 100) {
  # check for correct directory structure
  if (!dir.exists(pipestance_path))
    stop("Could not find the pipestance path: '", pipestance_path,"'.
         Please double-check if the directory exists.\n")
  od = file.path(pipestance_path, "outs")
  if (!dir.exists(od))
    stop("Could not find the pipestance output directory: '",
         file.path(pipestance_path,'outs'),
         "'. Please double-check if the directory exists.\n")

  v3p = file.path(od, "filtered_feature_bc_matrix")
  v2p = file.path(od, "filtered_gene_bc_matrices")
  v3d = dir.exists(v3p)
  if(barcode_filtered) {
    matrix_dir = ifelse(v3d, v3p, v2p)
  } else {
    matrix_dir = ifelse(v3d, file.path(od, "raw_feature_bc_matrix"),
                        file.path(od, "raw_gene_bc_matrices"))
  }
  if(!dir.exists(matrix_dir))
    stop("Could not find directory: ", matrix_dir)

  if(v3d) {
    features.loc <- file.path(matrix_dir, "features.tsv.gz")
    barcode.loc <- file.path(matrix_dir, "barcodes.tsv.gz")
    matrix.loc <- file.path(matrix_dir, "matrix.mtx.gz")
    summary.loc <- file.path(od, "metrics_summary_csv.csv")
  } else {
    genome = get_genome_in_matrix_path(matrix_dir, genome)
    barcode.loc <- file.path(matrix_dir, genome, "barcodes.tsv")
    features.loc <- file.path(matrix_dir, genome, "genes.tsv")
    matrix.loc <- file.path(matrix_dir, genome, "matrix.mtx")
    summary.loc <- file.path(od, "metrics_summary.csv")
  }
  if (!file.exists(barcode.loc)){
    stop("Barcode file missing")
  }
  if (!file.exists(features.loc)){
    stop("Gene name or features file missing")
  }
  if (!file.exists(matrix.loc)){
    stop("Expression matrix file missing")
  }
  # Not importing for now.
  #if(!file.exists(summary.loc)) {
  #  stop("Metrics summary file missing")
  #}
  data <- Matrix::readMM(matrix.loc)

  feature.names = utils::read.delim(features.loc,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  # Duplicate row names not allowed
  feature.names$V1 = make.unique(feature.names$V1)
  if(dim(data)[1] != length(feature.names[,1])) {
    stop(sprintf(paste("Mismatch dimension between gene file: \n\t %s\n and",
                       "matrix file: \n\t %s\n"), features.loc, matrix.loc))
  }
  if(v3d) {
    # We will only load GEX data for the relevant genome
    data_types = factor(feature.names$V3)
    allowed = data_types == "Gene Expression"
    if(!is.null(genome)) {
      # If not multigenome, no prefix will be added and we won't filter out
      # the one genome
      gfilter = grepl(genome, feature.names$V1)
      if(any(gfilter)) {
        allowed = allowed & grepl(genome, feature.names$V1)
      } else {
        message(paste("Data does not appear to be from a multi-genome sample,",
                      "simply returning all gene feature data without",
                      "filtering by genome."))
      }

    }
    data = data[allowed, ]
    feature.names = feature.names[allowed,1:2]
  }
  colnames(feature.names) = c("id", "gene_short_name")
  rownames(data) = feature.names[,"id"]
  rownames(feature.names) = feature.names[,"id"]

  barcodes <- utils::read.delim(barcode.loc, stringsAsFactors=FALSE, header=FALSE)
  if (dim(data)[2] != length(barcodes[,1])) {
    stop(sprintf(paste("Mismatch dimension between barcode file: \n\t %s\n",
                       "and matrix file: \n\t %s\n"), barcode.loc,matrix.loc))
  }
  barcodes$V1 = make.unique(barcodes$V1)
  colnames(data) = barcodes[,1]
  pd = data.frame(barcode=barcodes[,1], row.names=barcodes[,1])
  data <- data[,Matrix::colSums(data) > umi_cutoff]
  pd <- pd[colnames(data),, drop=FALSE]
  gbm <- new_cell_data_set(data,
                        cell_metadata = pd,
                        gene_metadata =  feature.names)

  return(gbm)
}
