#' Identify the genes most specifically expressed in specified groups of cells
#'
#' @export
top_markers <- function(cds,
                        group_cells_by="cluster",
                        genes_to_test_per_group=50,
                        reduction_method="UMAP",
                        expression_bins=5,
                        cores=1
                        ){

  # Yes, it's stupid we have cell ids both as a column and as the rownames.
  cell_group_df = data.frame(row.names=row.names(colData(cds)),
                             cell_id=row.names(colData(cds)))

  # Set up the table that partitions the cells into groups.
  # Must be either a column in colData or one of "cluster" or "partition.
  # FIXME: Should check its not a column you can't really use for grouping (i.e. a floating point value)
  if (group_cells_by == "cluster"){
    cell_group_df$cell_group = tryCatch({clusters(cds, reduction_method = reduction_method)}, error = function(e) {NULL})
  } else if (group_cells_by == "partition") {
    cell_group_df$cell_group = tryCatch({partitions(cds, reduction_method = reduction_method)}, error = function(e) {NULL})
  } else{
    cell_group_df$cell_group = colData(cds)[,group_cells_by]
  }

  # For each gene ggregate expression values across cells within each group
  # in a matrix thats genes x cell groups
  cluster_agg_exprs = aggregate_gene_expression(cds,
                                                cell_group_df=cell_group_df,
                                                norm_method="size_only")

  # Now compute a Jensen Shannon specificity score for each gene w.r.t each group
  cluster_spec_mat = specificity_matrix(cluster_agg_exprs, cores=cores)
  cluster_spec_table = tibble::rownames_to_column(as.data.frame(cluster_spec_mat))
  cluster_spec_table = gather(cluster_spec_table, "cell_group", "specificity", -rowname)
  spec_model_df = data.frame(rowname=row.names(cluster_spec_mat),
                             num_expressing=Matrix::rowSums(counts(cds) > 0),
                             mean_exprs=Matrix::rowMeans(cluster_agg_exprs),
                             max_spec=rowMaxs(cluster_spec_mat))

  # Now compute the expected max specificity as a function of how many cells express a given gene
  # Genes that are expressed in few cells tend to have very high specificity, so we want to
  # control for this trend when ranking genes by specificity later on
  spec_model_df = spec_model_df %>% mutate(quantile = ntile(num_expressing, expression_bins))
  spec_summary = spec_model_df %>% group_by(quantile) %>% summarize(log_spec_mean = mean(log(max_spec)), log_spec_sd = sd(log(max_spec)))
  spec_model_df = left_join(spec_model_df, spec_summary)

  # Compute the "specifity above expectation" for each gene w.r.t. each group:
  cluster_spec_table = dplyr::left_join(cluster_spec_table, spec_model_df)
  cluster_spec_table = cluster_spec_table %>% mutate(log_spec=log(specificity),
                                                     pval_excess_spec = pnorm(log(specificity),log_spec_mean, log_spec_sd, lower.tail=FALSE))


  cluster_spec_table = cluster_spec_table %>%
    #filter(num_expressing > 10) %>%
    group_by(cell_group) %>%
    top_n(genes_to_test_per_group, -pval_excess_spec)


  old_omp_num_threads = Sys.getenv("OMP_NUM_THREADS")
  Sys.setenv(OMP_NUM_THREADS = 1)
  marker_test_res = pbmcapply::pbmcmapply(test_marker_for_cell_group,
                                          cluster_spec_table$rowname,
                                          cluster_spec_table$cell_group,
                                          MoreArgs=list(cell_group_df, cds),
                                          mc.cores=cores)
  marker_test_res = t(marker_test_res)
  marker_test_res = as.matrix(marker_test_res)
  colnames(marker_test_res) = c("pseudo_R2", "lrtest_p_value")

  #marker_test_res = as.data.frame(marker_test_res)
  Sys.setenv(OMP_NUM_THREADS = old_omp_num_threads)

  marker_test_res = dplyr::bind_cols(cluster_spec_table, as.data.frame(marker_test_res))
  marker_test_res$lrtest_q_value = p.adjust(marker_test_res$lrtest_p_value,
                                            method="bonferroni",
                                            n=length(cluster_spec_mat))
  marker_test_res = marker_test_res %>% dplyr::select(rowname, cell_group, specificity, pseudo_R2, lrtest_p_value, lrtest_q_value)
  marker_test_res = marker_test_res %>% dplyr::rename(gene_id=rowname, marker_test_p_value=lrtest_p_value,  marker_test_q_value=lrtest_q_value)
  marker_test_res$pseudo_R2 = unlist(marker_test_res$pseudo_R2)
  marker_test_res$marker_test_p_value = unlist(marker_test_res$marker_test_p_value)
  return(marker_test_res)
}


# Calculate the probability vector
makeprobsvec<-function(p){
  phat<-p/sum(p)
  phat[is.na(phat)] = 0
  phat
}

# Calculate the probability matrix for a relative abundance matrix
makeprobs<-function(a){
  colSums<-apply(a,2,sum)
  b<-Matrix::t(Matrix::t(a)/colSums)
  b[is.na(b)] = 0
  b
}

# Calculate the Shannon entropy based on the probability vector
shannon.entropy <- function(p) {
  if (min(p) < 0 || sum(p) <=0)
    return(Inf)
  p.norm<-p[p>0]/sum(p)
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
  specificity_mat = pbmcapply::pbmclapply(row.names(agg_expr_matrix),
                                          FUN = function(x)
                                          {
                                            agg_exprs = as.numeric(agg_expr_matrix[x,])
                                            agg_exprs = makeprobsvec(agg_exprs)
                                            perfect_spec_matrix = diag(ncol(agg_expr_matrix))
                                            sapply(1:ncol(agg_expr_matrix), function(col_idx) {
                                              1 - JSdistVec(agg_exprs, perfect_spec_matrix[,col_idx])
                                            })
                                          }, mc.cores=cores,
                                          ignore.interactive = TRUE)
  specificity_mat = do.call(rbind, specificity_mat)
  colnames(specificity_mat) = colnames(agg_expr_matrix)
  row.names(specificity_mat) = row.names(agg_expr_matrix)
  return(specificity_mat)
  #
}


test_marker_for_cell_group = function(gene_id, cell_group, cell_group_df, cds){
  #print(gene_id)
  #print(cell_group)
  results = tryCatch({
    f_expression = log(as.numeric(counts(cds)[gene_id,]) / size_factors(cds) + 0.1)
    #print(sum(counts(cds)[gene_id,] > 0))
    is_member = as.character(cell_group_df[colnames(cds),2]) == as.character(cell_group)
    #print (is_member)

    model = speedglm::speedglm(is_member ~ f_expression,
                               acc=1e-3, model=FALSE,
                               y=FALSE,
                               verbose=TRUE,
                               family=binomial())
    null_model = speedglm::speedglm(is_member ~ 1,
                                    acc=1e-3, model=FALSE,
                                    y=FALSE,
                                    verbose=TRUE,
                                    family=binomial())
    lr.stat <- lmtest::lrtest(null_model, model)
    # #print (summary(model))
    # #print(summary(null_model))
    # #print (lr.stat)
    #print (str(lr.stat))
    n=ncol(cds)
    pseudo_R2 = (1-exp(-as.numeric(lr.stat$Chisq[2])/n))/(1-exp(2*as.numeric(logLik(null_model)/n)))
    LR_test_pval = lr.stat$`Pr(>Chisq)`[2]
    # model_summary = summary(model)
    # #print(model_summary)
    # #pval = as.numeric(as.character(model_summary$coefficients[2,4]))
    # pseudo_R2
    #pval
    return (list(pseudo_R2, LR_test_pval))
  }, error = function(e) { print(e); return(list(0.0, 1.0)) })

  #print(pval)
  return(results)
}

