#' Identify the genes most specifically expressed in groups of cells
#'
#' @param cds A cell_data_set object to calculate top markers for.
#' @param group_cells_by String indicating what to group cells by for
#'   comparison. Default is "cluster".
#' @param genes_to_test_per_group Numeric, how many genes of the top ranked
#'   specific genes by Jenson-Shannon to do the more expensive regression test
#'   on.
#' @param reduction_method String indicating the method used for dimensionality
#'   reduction. Currently only "UMAP" is supported.
#' @param marker_sig_test A flag indicating whether to assess the discriminative
#' power of each marker through logistic regression. Can be slow, consider
#' disabling to speed up top_markers().
#' @param reference_cells If provided, top_markers will perform the marker
#' significance test against a "reference set" of cells. Must be either a list
#' of cell ids from colnames(cds), or a positive integer. If the latter, top_markers()
#' will randomly select the specified number of reference cells. Accelerates
#' the marker significance test at some cost in sensitivity.
#' @param speedglm.maxiter Maximum number of iterations allowed for fitting GLM
#' models when testing markers for cell group.
#' @param cores Number of cores to use.
#' @param verbose Whether to print verbose progress output.
#'
#' @return a data.frame where the rows are genes and the columns are
#' * gene_id vector of gene names
#' * gene_short_name vector of gene short names
#' * cell_group character vector of the cell group to which the cell belongs
#' * marker_score numeric vector of marker scores as the fraction expressing scaled by the specificity. The value ranges from 0 to 1.
#' * mean_expression numeric vector of mean normalized expression of the gene in the cell group
#' * fraction_expressing numeric vector of fraction of cells expressing the gene within the cell group
#' * specificity numeric vector of a measure of how specific the gene's expression is to the cell group based on the Jensen-Shannon divergence. The value ranges from 0 to 1.
#' * pseudo_R2 numeric vector of pseudo R-squared values, a measure of how well the gene expression model fits the categorical data relative to the null model. The value ranges from 0 to 1.
#' * marker_test_p_value numeric vector of likelihood ratio p-values
#' * marker_test_q_value numeric vector of likelihood ratio q-values
#'
#' @export
top_markers <- function(cds,
                        group_cells_by="cluster",
                        genes_to_test_per_group=25,
                        reduction_method="UMAP",
                        marker_sig_test=TRUE,
                        reference_cells=NULL,
                        speedglm.maxiter=25,
                        cores=1,
                        verbose=FALSE
){

  # Yes, it's stupid we have cell ids both as a column and as the rownames.
  cell_group_df <- data.frame(row.names=row.names(colData(cds)),
                              cell_id=row.names(colData(cds)))

  # Set up the table that partitions the cells into groups.
  # Must be either a column in colData or one of "cluster" or "partition.
  # FIXME: Should check its not a column you can't really use for grouping
  # (i.e. a floating point value)
  if (group_cells_by == "cluster"){
    cell_group_df$cell_group <-
      tryCatch({clusters(cds, reduction_method = reduction_method)},
               error = function(e) {NULL})
  } else if (group_cells_by == "partition") {
    cell_group_df$cell_group <-
      tryCatch({partitions(cds, reduction_method = reduction_method)},
               error = function(e) {NULL})
  } else{
    cell_group_df$cell_group <- colData(cds)[,group_cells_by]
  }
  cell_group_df$cell_group <- as.character(cell_group_df$cell_group)

  if (verbose)
    message("Aggregating gene expression values for groups")
  # For each gene compute the fraction of cells expressing it within each group
  # in a matrix thats genes x cell groups

  cluster_binary_exprs = as.matrix(aggregate_gene_expression(cds,
                                                             cell_group_df=cell_group_df,
                                                             norm_method="binary"))

  cluster_mean_exprs = as.matrix(aggregate_gene_expression(cds,
                                                           cell_group_df=cell_group_df,
                                                           norm_method="size_only"))

  if (verbose)
    message("Computing Jensen-Shannon specificities")

  # Now compute a Jensen Shannon specificity score for each gene w.r.t each group
  cluster_spec_mat = specificity_matrix(cluster_mean_exprs, cores=cores)
  cluster_marker_score_mat = as.matrix(cluster_binary_exprs * cluster_spec_mat)

  if (verbose)
    message("Gathering score tables")
  cluster_marker_score_table = tibble::rownames_to_column(as.data.frame(cluster_marker_score_mat))
  cluster_marker_score_table = tidyr::gather(cluster_marker_score_table, "cell_group", "marker_score", -rowname)

  cluster_spec_table = tibble::rownames_to_column(as.data.frame(cluster_spec_mat))
  cluster_spec_table = tidyr::gather(cluster_spec_table, "cell_group", "specificity", -rowname)

  cluster_expr_table = tibble::rownames_to_column(as.data.frame(cluster_mean_exprs))
  cluster_expr_table = tidyr::gather(cluster_expr_table, "cell_group", "mean_expression", -rowname)

  cluster_fraction_expressing_table = tibble::rownames_to_column(as.data.frame(cluster_binary_exprs))
  cluster_fraction_expressing_table = tidyr::gather(cluster_fraction_expressing_table, "cell_group", "fraction_expressing", -rowname)

  # spec_model_df = data.frame(rowname=row.names(cluster_spec_score_mat),
  #                            num_expressing=Matrix::rowSums(SingleCellExperiment::counts(cds) > 0),
  #                            mean_exprs=Matrix::rowMeans(cluster_agg_exprs),
  #                            max_spec=rowMaxs(cluster_spec_score_mat))

  # # Now compute the expected max specificity as a function of how many cells express a given gene
  # # Genes that are expressed in few cells tend to have very high specificity, so we want to
  # # control for this trend when ranking genes by specificity later on
  # spec_model_df = spec_model_df %>% dplyr::mutate(quantile = dplyr::ntile(num_expressing, expression_bins))
  # spec_summary = spec_model_df %>% dplyr::group_by(quantile) %>% dplyr::summarize(log_spec_mean = mean(log(max_spec)), log_spec_sd = sd(log(max_spec)))
  # spec_model_df = dplyr::left_join(spec_model_df, spec_summary)
  #
  # # Compute the "specifity above expectation" for each gene w.r.t. each group:
  # cluster_spec_table = dplyr::left_join(cluster_spec_table, spec_model_df)
  # cluster_spec_table = cluster_spec_table %>% dplyr::mutate(log_spec=log(specificity),
  #                                                    pval_excess_spec = pnorm(log(specificity),log_spec_mean, log_spec_sd, lower.tail=FALSE))

  cluster_marker_score_table$specificity = cluster_spec_table$specificity
  cluster_marker_score_table$mean_expression = cluster_expr_table$mean_expression
  cluster_marker_score_table$fraction_expressing = cluster_fraction_expressing_table$fraction_expressing

  cluster_marker_score_table = cluster_marker_score_table %>%
    #filter(num_expressing > 10) %>%
    dplyr::group_by(cell_group) %>%
    dplyr::top_n(genes_to_test_per_group, marker_score)

  cell_group_df$cell_id <- as.character(cell_group_df$cell_id)
  cell_group_df$cell_group <- as.character(cell_group_df$cell_group)


  if (marker_sig_test){
    if (verbose)
      message("Running marker significance tests")
    # Temporarily disable OpenMP threading in functions to be run in parallel
    old_omp_num_threads = as.numeric(Sys.getenv("OMP_NUM_THREADS"))
    if (is.na(old_omp_num_threads)){
      old_omp_num_threads = 1
    }
    RhpcBLASctl::omp_set_num_threads(1)

    # Temporarily set the number of threads the BLAS library can use to be 1
    old_blas_num_threads = as.numeric(Sys.getenv("OPENBLAS_NUM_THREADS"))
    if (is.na(old_omp_num_threads)){
      old_blas_num_threads = 1
    }
    RhpcBLASctl::blas_set_num_threads(1)

    # Set up a balanced "reference" panel of cells from each group
    if (is.null(reference_cells) == FALSE){
      if(is.numeric(reference_cells)){
        num_ref_cells_per_group = reference_cells / length(unique(cell_group_df$cell_group))
        reference_cells = cell_group_df %>% dplyr::group_by(cell_group) %>%
          dplyr::sample_n(min(num_ref_cells_per_group, dplyr::n())) %>%
          dplyr::pull(cell_id)
        #reference_cells = sample(colnames(cds), reference_cells)
      } else {
        # TODO: check that reference cells is a list of valid cell ids.
      }
    }

    marker_test_res = tryCatch({pbmcapply::pbmcmapply(test_marker_for_cell_group,
                                                      cluster_marker_score_table$rowname,
                                                      cluster_marker_score_table$cell_group,
                                                      MoreArgs=list(cell_group_df, cds, reference_cells,
                                                      speedglm.maxiter),
                                                      ignore.interactive = TRUE,
                                                      mc.cores=cores)},
                               finally = {
                                 RhpcBLASctl::omp_set_num_threads(old_omp_num_threads)
                                 RhpcBLASctl::blas_set_num_threads(old_blas_num_threads)
                               })

    #
    # Check for possible convergence failure or other problems. Issue: #383
    #
    if('warning' %in% marker_test_res) {
        warning('test_marker_for_cell_group() caught warning: possible convergence failure.')
    }

    marker_test_res = t(marker_test_res)
    marker_test_res = as.matrix(marker_test_res)
    colnames(marker_test_res) = c("pseudo_R2", "lrtest_p_value")

    #marker_test_res = as.data.frame(marker_test_res)

    marker_test_res = dplyr::bind_cols(cluster_marker_score_table, as.data.frame(marker_test_res))
    marker_test_res$lrtest_q_value = stats::p.adjust(marker_test_res$lrtest_p_value,
                                              method="bonferroni",
                                              n=length(cluster_spec_mat))
    marker_test_res = marker_test_res %>% dplyr::select(rowname,
                                                        cell_group,
                                                        marker_score,
                                                        mean_expression,
                                                        fraction_expressing,
                                                        specificity,
                                                        pseudo_R2,
                                                        lrtest_p_value,
                                                        lrtest_q_value)
    
    marker_test_res = marker_test_res %>% dplyr::rename(gene_id=rowname, marker_test_p_value=lrtest_p_value,  marker_test_q_value=lrtest_q_value)
    marker_test_res$pseudo_R2 = unlist(marker_test_res$pseudo_R2)
    marker_test_res$marker_test_p_value = unlist(marker_test_res$marker_test_p_value)
    if ("gene_short_name" %in% colnames(rowData(cds))){
      marker_test_res = rowData(cds) %>%
      as.data.frame %>%
      tibble::rownames_to_column() %>%
      dplyr::select(rowname, gene_short_name) %>%
      dplyr::inner_join(marker_test_res, by=c("rowname"="gene_id"))
      marker_test_res = marker_test_res %>% dplyr::rename(gene_id=rowname)
    }
  } else {
    marker_test_res = cluster_marker_score_table
    marker_test_res = marker_test_res %>% dplyr::rename(gene_id=rowname)
  }

  if (verbose)
    message("Done")

  return(marker_test_res)
}


# Calculate the probability vector
makeprobsvec <- function(p) {
  phat <- p/sum(p)
  phat[is.na(phat)] = 0
  phat
}

# Calculate the probability matrix for a relative abundance matrix
makeprobs <- function(a) {
  colSums<-apply(a,2,sum)
  b <- Matrix::t(Matrix::t(a)/colSums)
  b[is.na(b)] = 0
  b
}

# Calculate the Shannon entropy based on the probability vector
shannon.entropy <- function(p) {
  if (min(p) < 0 || sum(p) <=0)
    return(Inf)
  p.norm <- p[p>0]/sum(p)
  -sum(log2(p.norm)*p.norm)
}

# Calculate the Jessen-Shannon distance for two probability distribution
JSdistVec <- function (p, q)
{
  JSdiv <- shannon.entropy((p + q)/2) - (shannon.entropy(p) +
                                           shannon.entropy(q)) * 0.5
  JSdiv[is.infinite(JSdiv)] <- 1
  JSdiv[JSdiv < 0] <- 0
  JSdist <- sqrt(JSdiv)
  JSdist
}

specificity_matrix <- function(agg_expr_matrix, cores=1){
  specificity_mat <-
    pbmcapply::pbmclapply(row.names(agg_expr_matrix),
                          FUN = function(x) {
                            agg_exprs = as.numeric(agg_expr_matrix[x,])
                            agg_exprs = makeprobsvec(agg_exprs)
                            perfect_spec_matrix = diag(ncol(agg_expr_matrix))
                            sapply(1:ncol(agg_expr_matrix), function(col_idx) {
                              1 - JSdistVec(agg_exprs,
                                            perfect_spec_matrix[,col_idx])
                            })
                          }, mc.cores=cores,
                          ignore.interactive = TRUE)
  specificity_mat = do.call(rbind, specificity_mat)
  colnames(specificity_mat) = colnames(agg_expr_matrix)
  row.names(specificity_mat) = row.names(agg_expr_matrix)
  return(specificity_mat)
  #
}

enrichment_matrix <- function(agg_expr_matrix, cores=1){
  specificity_mat = pbmcapply::pbmclapply(row.names(agg_expr_matrix),
                                          FUN = function(x)
                                          {
                                            agg_exprs = as.numeric(agg_expr_matrix[x,])
                                            agg_exprs = makeprobsvec(agg_exprs)
                                            perfect_spec_matrix = diag(ncol(agg_expr_matrix))
                                            sapply(1:ncol(agg_expr_matrix), function(col_idx) {
                                              1 - JSdistVec(agg_exprs, perfect_spec_matrix[,col_idx])
                                            }
                                            )
                                          }, mc.cores=cores,
                                          ignore.interactive = TRUE)
  specificity_mat = do.call(rbind, specificity_mat)
  colnames(specificity_mat) = colnames(agg_expr_matrix)
  row.names(specificity_mat) = row.names(agg_expr_matrix)
  return(specificity_mat)
  #
}

test_marker_for_cell_group = function(gene_id, cell_group, cell_group_df, cds,
                                      reference_cells=NULL, speedglm.maxiter=25){
  #print(gene_id)
  #print(cell_group)
  #print (length(reference_cells))
  results <- tryCatch({
    f_expression <-
      log(as.numeric(SingleCellExperiment::counts(cds)[gene_id,]) / size_factors(cds) + 0.1)
    #print(sum(SingleCellExperiment::counts(cds)[gene_id,] > 0))
    is_member <-
      as.character(cell_group_df[colnames(cds),2]) == as.character(cell_group)
    names(is_member) = names(f_expression) = colnames(cds)
    is_member[is.na(is_member)] = FALSE
    is_member[is.null(is_member)] = FALSE

    if (is.null(reference_cells) == FALSE){
      # Exclude cells that aren't in either the cell_group or the
      # reference_panel
      f_expression <- f_expression[is_member | names(f_expression) %in%
                                     reference_cells]
      is_member <- is_member[is_member | names(is_member) %in% reference_cells]
    }

    if (sum(is.na(f_expression)) > 0 || sum(is.na(is_member)) > 0){
      stop("Expression and group membership can't be NA")
    }
    model <- speedglm::speedglm(is_member ~ f_expression,
                                acc=1e-3, model=FALSE,
                                y=FALSE,
                                verbose=TRUE,
                                family=stats::binomial(),
                                maxit=speedglm.maxiter)
    assertthat::assert_that(model$convergence == TRUE, msg=paste0('speedglm model failed to converge in ',speedglm.maxiter, ' iterations.'))
    null_model <- speedglm::speedglm(is_member ~ 1,
                                     acc=1e-3, model=FALSE,
                                     y=FALSE,
                                     verbose=TRUE,
                                     family=stats::binomial(),
                                     maxit=speedglm.maxiter)
    assertthat::assert_that(model$convergence == TRUE, msg=paste0('speedglm null model failed to converge in ',speedglm.maxiter, ' iterations.'))
    lr.stat <- lmtest::lrtest(null_model, model)
    #print (summary(model))
    # #print(summary(null_model))
    # #print (lr.stat)
    #print (str(lr.stat))
    n=ncol(cds)
    pseudo_R2 <-
      (1-exp(-as.numeric(lr.stat$Chisq[2])/n)) /
      (1-exp(2*as.numeric(stats::logLik(null_model)/n)))
    LR_test_pval = lr.stat$`Pr(>Chisq)`[2]
    # model_summary = summary(model)
    # #print(model_summary)
    # #pval = as.numeric(as.character(model_summary$coefficients[2,4]))
    # pseudo_R2
    #pval
    return (list(pseudo_R2, LR_test_pval))
  }, error = function(e) { return(list(0.0, 1.0)) })

  #print(pval)
  return(results)
}

#' Generate a Garnett marker file from top_markers output.
#'
#' @param marker_test_res Tibble of top markers, output of
#'   \code{\link{top_markers}}.
#' @param file Path to the marker file to be generated. Default is
#'   "./marker_file.txt".
#' @param max_genes_per_group Numeric, the maximum number of genes to output
#'   per cell type entry. Default is 10.
#' @param remove_duplicate_genes Logical indicating whether marker genes that
#'   mark multiple cell groups should be excluded. Default is FALSE. When
#'   FALSE, a message will be emitted when duplicates are present.
#'
#' @return None, marker file is written to \code{file} parameter location.
#' @export
#'
generate_garnett_marker_file <- function(marker_test_res,
                                         file = "./marker_file.txt",
                                         max_genes_per_group = 10,
                                         remove_duplicate_genes = FALSE) {
  marker_test_res <- as.data.frame(marker_test_res)
  if(is.null(marker_test_res$group_name)) {
    marker_test_res$group_name <- paste("Cell type", marker_test_res$cell_group)
  }
  group_list <- unique(marker_test_res$group_name)

  good_markers <- marker_test_res %>% dplyr::group_by(group_name) %>%
    dplyr::top_n(max_genes_per_group, marker_score)

  dups <- good_markers$gene_id[duplicated(good_markers$gene_id)]

  if ("gene_short_name" %in% colnames(good_markers)) {
    dups_gsn <- good_markers$gene_short_name[duplicated(good_markers$gene_short_name)]
  }


  if(remove_duplicate_genes) {
    good_markers <- good_markers[!good_markers$gene_id %in% dups,]
  } else {
    if (length(dups) > 0) {
      if("gene_short_name" %in% colnames(good_markers)) {
        message(paste("The following marker genes mark multiple cell groups.",
                      "Prior to using Garnett, we recommend either excluding",
                      "these genes using remove_duplicate_genes = TRUE, or",
                      "modifying your marker file to make the cell types with",
                      "the shared marker subtypes in a hierarchy.",
                      paste(dups_gsn, collapse = ", ")))
      } else {
        message(paste("The following marker genes mark multiple cell groups.",
                      "Prior to using Garnett, we recommend either excluding",
                      "these genes using remove_duplicate_genes = TRUE, or",
                      "modifying your marker file to make the cell types with",
                      "the shared marker subtypes in a hierarchy.",
                      paste(dups, collapse = ", ")))
      }
    }
  }

  output <- list()

  for (group in group_list) {
    if (sum(good_markers$group_name == group) == 0) {
      message(paste(group, "did not have any markers above the q-value",
                    "threshold. It will be skipped."))
      next
    }

    good_name <- gsub("\\(|\\)|:|>|,|#", ".", group)
    if (good_name != group) {
      warning(paste("Group name contained an illegal character for a Garnett ",
                    "marker file. ", group, " will be substituted for ",
                    good_name, "."))
    }

    sub <- good_markers[good_markers$group_name == group,]
    if (nrow(sub) > max_genes_per_group) {
      sub <- sub[order(sub$marker_test_q_value),][1:max_genes_per_group,]
    }
    if ("gene_short_name" %in% colnames(sub)){
      entry <- paste0("> ", good_name, "\n", "expressed: ",
                      paste(sub$gene_short_name, collapse = ", "), "\n")
    } else {
      entry <- paste0("> ", good_name, "\n", "expressed: ",
                      paste(sub$gene_id, collapse = ", "), "\n")
    }

    output <- append(output, entry)
  }

  all <- paste(output, collapse = "\n")

  write(all, file=file)
  message(paste("Garnett marker file written to", file))
}

