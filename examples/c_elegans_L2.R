rm(list = ls())  # Clear the environment
#options(warn=-1) # Turn off warning message globally
library(monocle3) # Load Monocle
library(tidymodels)
theme_set(theme_gray(base_size = 6))

expression_matrix = readRDS("L2_expression.rda")
cell_metadata = readRDS("L2_cell_metadata.rda")
gene_annotation = readRDS("L2_gene_metadata.rda")

# expression_matrix = Matrix::readMM(gzcon(url("http://jpacker-data.s3.amazonaws.com/public/worm/L2/Cao_et_al_2017_data_2019_update.exprs.mm.gz")))
# gene_annotation = read.delim(url("http://jpacker-data.s3.amazonaws.com/public/worm/L2/Cao_et_al_2017_data_2019_update.fData.tsv"))
# cell_metadata  = read.delim(url("http://jpacker-data.s3.amazonaws.com/public/worm/L2/Cao_et_al_2017_data_2019_update.pData.tsv"))


cds = new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


## Step 1: Normalize and pre-process the data
cds = preprocess_cds(cds, num_dim = 100)

## Step 2: Reduce the dimensionality of the data
#### Without batch correction:
cds = reduce_dimension(cds)
plot_cells(cds) + ggsave("L2_umap_no_color.png", width=5, height=4, dpi = 600)
plot_cells(cds, color_cells_by="cao_cell_type") + ggsave("L2_umap_color_cells_by_cao_type.png", width=5, height=4, dpi = 600)
plot_cells(cds, color_cells_by="plate", label_cell_groups=FALSE) + ggsave("L2_umap_plate.png", width=5, height=4, dpi = 600)

#### With batch correction:
cds = preprocess_cds(cds, num_dim = 100, residual_model_formula_str = "~ plate")
cds = reduce_dimension(cds)
plot_cells(cds, color_cells_by="plate", label_cell_groups=FALSE) + ggsave("L2_umap_corrected_plate.png", width=5, height=4, dpi = 600)

## Step 3: (Optional) Cluster cells
cds = cluster_cells(cds, resolution=c(10^seq(-6,-1)))
plot_cells(cds) + ggsave("L2_umap_color_cells_by_cluster.png", width=5, height=4, dpi = 600)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition") + ggsave("L2_umap_color_cells_by_partition.png", width=5, height=4, dpi = 600)

pheatmap::pheatmap(log(table(clusters(cds), colData(cds)$cao_cell_type)+1),
                   clustering_method="ward.D2",
                   fontsize=6, width=5, height=8, filename="L2_cell_type_by_cluster.png")

plot_cells(cds, color_cells_by="cao_cell_type")  + ggsave("L2_umap_corrected_cao_type.png", width=5, height=4, dpi = 600)
plot_cells(cds, color_cells_by="cao_cell_type", label_groups_by_cluster=FALSE)  + ggsave("L2_umap_corrected_cao_type_no_cluster_label.png", width=5, height=4, dpi = 600)

## Step 4: Find genes expressed in each cluster or cell type
pr_graph_test_res = graph_test(cds, neighbor_graph="knn", cores=8)
pr_deg_ids = row.names(subset(pr_graph_test_res, q_value < 0.05))

# Make a matrix with average gene expression values for each cluster:
cluster_agg_exprs = aggregate_gene_expression(cds[pr_deg_ids,], 
                                             cell_group_df=data.frame(row.names(colData(cds)), clusters(cds)),
                                             norm_method="size_only")

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

cluster_spec_mat = specificity_matrix(cluster_agg_exprs[pr_deg_ids,], cores=8)
cluster_spec_table = tibble::rownames_to_column(as.data.frame(cluster_spec_mat))
cluster_spec_table = gather(cluster_spec_table, "cell_group", "specificity", -rowname)
spec_model_df = data.frame(num_expressing=Matrix::rowSums(counts(cds[pr_deg_ids,]) > 0), 
                           mean_exprs=rowMeans(cluster_agg_exprs), 
                           max_spec=rowMaxs(cluster_spec_mat))

# spec_model = glm(max_spec ~ splines::ns(log(mean_exprs), df=3), data=spec_model_df)
# spec_model_predictions = predict(spec_model, type = "link", se.fit = TRUE)
# spec_model$family$linkinv(spec_model_predictions$fit)

# mod <- glm(max_spec ~ splines::ns(mean_exprs, df=3),, data=spec_model_df)
# preddata <- with(spec_model_df, data.frame(mean_exprs = seq(min(mean_exprs), max(mean_exprs), length = 100)))
# preds <- predict(mod, newdata = preddata, type = "response", se.fit = TRUE)
# critval <- 1.96 ## approx 95% CI
# upr <- preds$fit + (critval * preds$se.fit)
# lwr <- preds$fit - (critval * preds$se.fit)
# fit <- preds$fit

# fit2 <- mod$family$linkinv(fit)
# upr2 <- mod$family$linkinv(upr)
# lwr2 <- mod$family$linkinv(lwr)

# preddata$lwr <- lwr2 
# preddata$upr <- upr2 
# ggplot(data=spec_model_df, mapping=aes(x=mean_exprs,y=max_spec)) + geom_point() +         
#    #stat_smooth(method="glm") + 
#    geom_line(data=preddata, mapping=aes(x=mean_exprs, y=upr), col="red") + 
#    geom_line(data=preddata, mapping=aes(x=mean_exprs, y=lwr), col="red") 

spec_model = loess(max_spec ~ log(num_expressing), data=spec_model_df)
spec_model_predictions = predict(spec_model, se=TRUE)
spec_model_df$fitted_spec = spec_model_predictions$fit
spec_model_df$spec_residual = residuals(spec_model) 
spec_model_df$spec_std_residual = spec_model_df$spec_residual / sqrt(spec_model_predictions$se.fit)

cluster_spec_table = dplyr::left_join(cluster_spec_table, tibble::rownames_to_column(spec_model_df))
cluster_spec_table = cluster_spec_table %>% mutate(excess_spec = specificity - fitted_spec,
                                                   spec_score = 0.1*log(num_expressing) + excess_spec)
top_specific_markers = unique(cluster_spec_table %>% 
    filter(num_expressing > 10) %>%
    group_by(cell_group) %>% 
    top_n(1, spec_score) %>% pull(rowname))
pheatmap::pheatmap(cluster_spec_mat[top_specific_markers,], cluster_rows=FALSE, cluster_cols=FALSE)



plot_genes_by_group <- function(cds, 
                                markers, 
                                group_cells_by="cluster",
                                reduction_method = "UMAP",
                                norm_method = c("log", "size_only"),
                                  lower_threshold = 0,
                                  max.size = 10, 
                                  ordering_type = c('cluster_row_col', 'maximal_on_diag', 'none'), # maybe be also do the maximum color on the diagonal; the axis change be switched too 
                                  axis_order = c('group_marker', 'marker_group'), 
                                  flip_percentage_mean = FALSE, 
                                  pseudocount = 1, 
                                  scale_max = 3, 
                                  scale_min = -3,
                                  ...) {
  #if(!(group_by %in% colnames(pData(cds)))) 
  #  stop(paste0(group_by, ' is not in the pData, please first perform cell clustering (clusterCells) before running this function!'))
  
  assertthat::assert_that(is(cds, "cell_data_set"))
  
  if(!is.null(group_cells_by)) {
    assertthat::assert_that(group_cells_by %in% c("cluster", "partition") | group_cells_by %in% names(colData(cds)),
                            msg = paste("group_cells_by must be a column in the",
                                        "colData table."))
  }

  norm_method = match.arg(norm_method)
  #group_cells_by=match.arg(group_cells_by)

  gene_ids = as.data.frame(fData(cds)) %>% 
                tibble::rownames_to_column() %>% 
                dplyr::filter(rowname %in% markers | gene_short_name %in% markers) %>% dplyr::pull(rowname)
  if(length(gene_ids) < 1) 
    stop('Please make sure markers are included in the gene_short_name column of the fData!')          
  
  if(flip_percentage_mean == FALSE){
    major_axis <- 1
    minor_axis <- 2
  } else if (flip_percentage_mean == TRUE){
    major_axis <- 2
    minor_axis <- 1
  }
  
  exprs_mat <- t(as.matrix(exprs(cds)[gene_ids, ]))
  exprs_mat <- reshape2::melt(exprs_mat)
  colnames(exprs_mat) <- c('Cell', 'Gene', 'Expression')
  exprs_mat$Gene <- as.character(exprs_mat$Gene)


  if (group_cells_by == "cluster"){
    cell_group <- tryCatch({clusters(cds, reduction_method = reduction_method)}, error = function(e) {NULL})
  } else if (group_cells_by == "partition") {
    cell_group <- tryCatch({partitions(cds, reduction_method = reduction_method)}, error = function(e) {NULL})
  } else{
    cell_group <- colData(cds)[,group_cells_by]
  }
  names(cell_group) = colnames(cds)

  exprs_mat$Group <- cell_group[exprs_mat$Cell]
  exprs_mat = exprs_mat %>% filter(is.na(Group) == FALSE)
  ExpVal <- exprs_mat %>% dplyr::group_by(Group, Gene) %>% dplyr::summarize(mean = mean(log(Expression + pseudocount)), percentage = sum(Expression > lower_threshold) / length(Expression))
  ExpVal$mean <- ifelse(ExpVal$mean < scale_min, scale_min, ExpVal$mean)
  ExpVal$mean <- ifelse(ExpVal$mean > scale_max, scale_max, ExpVal$mean)
  
  ExpVal$Gene <- fData(cds)[ExpVal$Gene, 'gene_short_name']
  
  res <- reshape2::dcast(ExpVal[, 1:4], Group ~ Gene, value.var = colnames(ExpVal)[2 + major_axis])
  group_id <- res[, 1]
  res <- res[, -1]
  row.names(res) <- group_id
  
  if(ordering_type == 'cluster_row_col') {
    row_dist <- as.dist((1 - cor(t(res[, -1])))/2)
    row_dist[is.na(row_dist)] <- 1
    
    col_dist <- as.dist((1 - cor(res[, -1]))/2)
    col_dist[is.na(col_dist)] <- 1
    
    ph <- pheatmap::pheatmap(res[, -1], 
                   useRaster = T,
                   cluster_cols=TRUE, 
                   cluster_rows=TRUE, 
                   show_rownames=F, 
                   show_colnames=F, 
                   clustering_distance_cols=col_dist,
                   clustering_distance_rows=row_dist,
                   clustering_method = 'ward.D2',
                   silent=TRUE,
                   filename=NA)
    
    ExpVal$Gene <- factor(ExpVal$Gene, levels = colnames(res)[-1][ph$tree_col$order])
    ExpVal$Group <- factor(ExpVal$Group, levels = row.names(res)[ph$tree_row$order])
    
  } else if(ordering_type == 'maximal_on_diag'){
    
    order_mat <- t(apply(res, major_axis, order))
    max_ind_vec <- c()
    for(i in 1:nrow(order_mat)) {
      tmp <- max(which(!(order_mat[i, ] %in% max_ind_vec)))
      max_ind_vec <- c(max_ind_vec, order_mat[i, tmp])
    }
    max_ind_vec <- max_ind_vec[!is.na(max_ind_vec)]
    
    if(major_axis == 1){
      # ExpVal$Gene <- factor(ExpVal %>% dplyr::pull(colnames(ExpVal)[minor_axis]) , levels = dimnames(res)[[minor_axis]][max_ind_vec])
      max_ind_vec <- c(max_ind_vec, setdiff(1:length(markers), max_ind_vec))
      ExpVal$Gene <- factor(ExpVal$Gene , levels = dimnames(res)[[2]][max_ind_vec])
    }
    else{
      max_ind_vec <- c(max_ind_vec, setdiff(1:length(unique(exprs_mat$Group)), max_ind_vec))
      ExpVal$Group <- factor(ExpVal$Group , levels = dimnames(res)[[1]][max_ind_vec])
    }
  } else if(ordering_type == 'none'){
    ExpVal$Gene <- factor(ExpVal$Gene, levels = markers)
  }
  
  # + scale_color_gradient(low ="blue",   high = "red", limits=c(1, max(ExpVal$mean) ))
  if(flip_percentage_mean){
    g <- ggplot(ExpVal, aes(y = Gene,  x = Group)) + geom_point(aes(colour = percentage,  size = mean)) + viridis::scale_color_viridis(name = 'percentage') + scale_size(name = 'log(mean + 0.1)', range = c(0, max.size))
  } else {
    g <- ggplot(ExpVal, aes(y = Gene,  x = Group)) + geom_point(aes(colour = mean,  size = percentage)) + viridis::scale_color_viridis(name = 'log(mean + 0.1)') + scale_size(name = 'percentage', range = c(0, max.size))
  }
  
  g <- g + xlab("Cluster") + ylab("Gene") + monocle3:::monocle_theme_opts() + xlab("Cluster") + ylab("Gene") + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  if(axis_order == 'marker_group') {
    g <- g + coord_flip()
  }
  
  g
}
#debug(plot_genes_by_group)
plot_genes_by_group(cds, top_specific_markers, group_cells_by="cao_cell_type")

gene_cluster_df = cluster_genes(cds[pr_deg_ids,], resolution=0.001)

# text_df = gene_cluster_df %>% dplyr::group_by(cluster) %>% dplyr::summarize(text_x = median(x = dim_1),
#                                                                              text_y = median(x = dim_2))
# ggplot2::qplot(dim_1, dim_2, color=cluster, data=gene_cluster_df) +
#   ggplot2::geom_text(data=text_df, mapping = ggplot2::aes_string(x = "text_x", y = "text_y", label = "cluster"), color=I("black"),  size = 4)

gene_cluster_df = cluster_genes(cds[pr_deg_ids,], resolution=c(0,10^seq(-6,-1)))
cell_group_df = tibble::tibble(cell=row.names(colData(cds)), cell_group=clusters(cds))
agg_mat = aggregate_gene_expression(cds, gene_cluster_df, cell_group_df)
#agg_mat
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                    scale="column", clustering_method="ward.D2", fontsize=6)


png("module_graph.png", res=600, width=8, height=8, units="in")
plot_cells(cds, genes=gene_cluster_df, cell_size=0.5, show_backbone=FALSE, label_branch_points=FALSE, label_leaves=FALSE, label_roots=FALSE)
dev.off()




## Step 5: Annotate cells by type:
#ceModel = readRDS(url("https://cole-trapnell-lab.github.io/garnett/classifiers/ceWhole"))
#library(org.Ce.eg.db)
library(garnett)
load(url("https://cole-trapnell-lab.github.io/garnett/classifiers/ceWhole"))

colData(cds)$garnett_cluster = clusters(cds)
cds = classify_cells(cds, ceWhole,
                           db = org.Ce.eg.db::org.Ce.eg.db,
                           cluster_extend = TRUE,
                           cds_gene_id_type = "ENSEMBL")

plot_cells(cds, group_cells_by="partition", color_cells_by="cluster_ext_type") + ggsave("L2_umap_corrected_garnett_ext_type.png", width=5, height=4, dpi = 600)

########
# Look at a specific set of clusters

cds_subset = choose_cells(cds)

pr_graph_test_res = graph_test(cds_subset, neighbor_graph="knn", cores=8)
pr_deg_ids = row.names(subset(pr_graph_test_res, q_value < 0.05))

gene_cluster_df = monocle3:::cluster_genes(cds[pr_deg_ids,], resolution=c(0,10^seq(-6,-1)), verbose=TRUE)
cell_group_df = tibble::tibble(cell=row.names(colData(cds)), cell_group=clusters(cds))
agg_mat = aggregate_gene_expression(cds, gene_cluster_df, cell_group_df)
#agg_mat
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                    scale="column", clustering_method="ward.D2", fontsize=6)

png("branch_modules.png", res=600, width=6, height=6, units="in")
plot_cells(cds_subset, genes=gene_cluster_df, show_trajectory_graph=FALSE, label_cell_groups=FALSE)
dev.off()


mod_df = data.frame(resolution= c(0,10^seq(-6,-1)),
                    modularity = c(0.0463080808623619, 0.0463080808623619, 0.519517905711304, 0.722441186959141, 0.893506647690171, 0.926503643050994, 0.861890178533872),
                    clusters = c(5, 5, 6, 8, 18, 60, 179))

qplot(clusters, modularity, data=mod_df)

qplot(diff(mod_df$modularity), diff(mod_df$clusters))

qplot(10^seq(-6,-1), diff(mod_df$modularity)/diff(mod_df$clusters), log="x")
