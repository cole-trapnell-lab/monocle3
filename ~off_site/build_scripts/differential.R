#
# This script creates the plots and tables used by
# the differential.html page and stores the new
# files in the directories given by the bash
# variables 'manual_images_dir' and 'data_de_dir', as
# defined in globals.R.
#
# Call this script from R using the command
#   source('differential.R', local=TRUE, echo=TRUE)
#
# There is information about using ggsave in the blog
# page at
#
#   URL: https://minimaxir.com/2017/08/ggplot2-web/
#
# A plot template looks like
#   plot_file_name <- 'L2_umap_no_color.png'
#   message('plot ', plot_file_name)
#   plot_cmd <- plot_cells(cds)
#   plot_cmd + theme(text=element_text(size=6))
#   ggsave(file.path(manual_images_dir, plot_file_name), units='in', width=5, height=4, dpi=600, device='png')
#
# Notes:
#   o  the png() driver makes the text and cell points too small
#   o  the ggsave() function works except the text is too large
#      without setting 'theme(text=element_text(size=6))'
#
# A png driver-based template looks like
#   plot_cmd <- plot_cells(cds, color_cells_by="cao_cell_type")
#   plot_cmd
#   plot_file_name <- 'xxxx'
#   message('plot ', plot_file_name)
#   png(file.path(manual_images_dir, plot_file_name), width=<nx>, height=<ny>, res=300)
#   plot_cmd
#   dev.off()
#
# but I abandoned this form except for heatmaps, which are
# not made using ggplot2.


source('globals.R', local=TRUE, echo=TRUE)
devtools::load_all(monocle3_git_dir, quiet=TRUE)

library(dplyr)

# Number of significant digits in real values in tables.
nsigdigits <- 3

# This is copied from trajectories because the differential page appears
# to assume that the embryo data are loaded.
#
expression_matrix <- readRDS(url("https://staff.washington.edu/hpliner/data/packer_embryo_expression.rds"))
cell_metadata <- readRDS(url("https://staff.washington.edu/hpliner/data/packer_embryo_colData.rds"))
gene_annotation <- readRDS(url("https://staff.washington.edu/hpliner/data/packer_embryo_rowData.rds"))

cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

# Note: Your data will not have the loading batch information demonstrated here, you will correct batch using your own batch information.
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")

# However, unlike clustering, which works well with both UMAP and t-SNE, here we strongly urge you to use UMAP, the default method:
cds <- reduce_dimension(cds)

# Starting on differential expression page...
# Let's begin with a small set of genes that we know are important in ciliated neurons to demonstrate Monocle's capabilities:
ciliated_genes <- c("che-1",
                    "hlh-17",
                    "nhr-6",
                    "dmd-6",
                    "ceh-36",
                    "ham-1")
cds_subset <- cds[rowData(cds)$gene_short_name %in% ciliated_genes,]

#  To do so, we first call the fit_models() function:
gene_fits <- fit_models(cds_subset, model_formula_str = "~embryo.time")

# Now let's see which of these genes have time-dependent expression. First, we extract a table of coefficients from each model using the coefficient_table() function:
fit_coefs <- coefficient_table(gene_fits)

# fit_coefs looks like this: ...
fit_coefs
out_emb_terms <- as.data.frame(fit_coefs[,c('id', 'gene_short_name', 'num_cells_expressed', 'status', 'term', 'estimate', 'std_err', 'test_val', 'p_value', 'normalized_effect', 'model_component', 'q_value')])
out_emb_terms <- rapply(out_emb_terms, f=formatC, classes="numeric", how="replace", digits=nsigdigits, format='g')
csv_file_name <- 'emb_terms.csv'
write.csv(out_emb_terms, file.path(data_de_dir, csv_file_name), row.names=TRUE)

# We generally don't care about the intercept term \beta_0, so we can easily just extract the time terms:
emb_time_terms <- fit_coefs %>% filter(term == "embryo.time")

# We can filter the results and control the false discovery rate as follows:
out_emb_time_terms <- as.data.frame(emb_time_terms %>% filter (q_value < 0.05) %>% select(gene_short_name, term, q_value, estimate))
out_emb_time_terms <- rapply(out_emb_time_terms, f=formatC, classes="numeric", how="replace", digits=nsigdigits, format='g')
csv_file_name <- 'emb_time_sig_terms.csv'
write.csv(out_emb_time_terms, file.path(data_de_dir, csv_file_name), row.names=TRUE)

# One type of plot is a "violin" plot.
plot_file_name <- 'embryo_ciliated_markers_violin.png'
plot_cmd <- plot_genes_violin(cds_subset, group_cells_by="embryo.time.bin", ncol=2) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# Controlling for batch effects and other factors
gene_fits <- fit_models(cds_subset, model_formula_str = "~embryo.time + batch")
fit_coefs <- coefficient_table(gene_fits)
out_emb_plus_batch_terms <- as.data.frame(fit_coefs %>% filter(term != "(Intercept)") %>% select(gene_short_name, term, q_value, estimate))
out_emb_plus_batch_terms <- rapply(out_emb_plus_batch_terms, f=formatC, classes="numeric", how="replace", digits=nsigdigits, format='g')
csv_file_name <- 'emb_plus_batch_terms.csv'
write.csv(out_emb_plus_batch_terms, file.path(data_de_dir, csv_file_name), row.names=TRUE)

# We can evaluate the fits of each model using the evaluate_fits() function:
evaluate_gene_fits <- evaluate_fits(gene_fits)
out_emb_time_evals <- as.data.frame(evaluate_gene_fits)
out_emb_time_evals <- rapply(out_emb_time_evals, f=formatC, classes="numeric", how="replace", digits=nsigdigits, format='g')
csv_file_name <- 'emb_time_evals.csv'
write.csv(out_emb_time_evals, file.path(data_de_dir, csv_file_name), row.names=TRUE)

# You run compare_models() like this:
time_batch_models <- fit_models(cds_subset, model_formula_str = "~embryo.time + batch", expression_family="negbinomial")
time_models <- fit_models(cds_subset, model_formula_str = "~embryo.time", expression_family="negbinomial")
out_emb_model_lr_test <- compare_models(time_batch_models, time_models) %>% select(gene_short_name, q_value)
out_emb_model_lr_test <- rapply(out_emb_model_lr_test, f=formatC, classes="numeric", how="replace", digits=nsigdigits, format='g')
csv_file_name <- 'emb_model_lr_test.csv'
write.csv(out_emb_model_lr_test, file.path(data_de_dir, csv_file_name), row.names=TRUE)


# In the L2 worm data, we identified a number of clusters that were very distinct as neurons:
# reload and reprocess the data as described in the 'Clustering and classifying your cells' section
expression_matrix <- readRDS(url("https://staff.washington.edu/hpliner/data/cao_l2_expression.rds"))
cell_metadata <- readRDS(url("https://staff.washington.edu/hpliner/data/cao_l2_colData.rds"))
gene_annotation <- readRDS(url("https://staff.washington.edu/hpliner/data/cao_l2_rowData.rds"))

# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds, resolution=1e-5)

colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$assigned_cell_type,
                                                 "1"="Body wall muscle",
                                                 "2"="Germline",
                                                 "3"="Motor neurons",
                                                 "4"="Seam cells",
                                                 "5"="Sex myoblasts",
                                                 "6"="Socket cells",
                                                 "7"="Marginal_cell",
                                                 "8"="Coelomocyte",
                                                 "9"="Am/PH sheath cells",
                                                 "10"="Ciliated neurons",
                                                 "11"="Intestinal/rectal muscle",
                                                 "12"="Excretory gland",
                                                 "13"="Chemosensory neurons",
                                                 "14"="Interneurons",
                                                 "15"="Unclassified eurons",
                                                 "16"="Ciliated neurons",
                                                 "17"="Pharyngeal gland cells",
                                                 "18"="Unclassified neurons",
                                                 "19"="Chemosensory neurons",
                                                 "20"="Ciliated neurons",
                                                 "21"="Ciliated neurons",
                                                 "22"="Inner labial neuron",
                                                 "23"="Ciliated neurons",
                                                 "24"="Ciliated neurons",
                                                 "25"="Ciliated neurons",
                                                 "26"="Hypodermal cells",
                                                 "27"="Mesodermal cells",
                                                 "28"="Motor neurons",
                                                 "29"="Pharyngeal gland cells",
                                                 "30"="Ciliated neurons",
                                                 "31"="Excretory cells",
                                                 "32"="Amphid neuron",
                                                 "33"="Pharyngeal muscle")

# Subset just the neurons:
neurons_cds <- cds[,grepl("neurons", colData(cds)$assigned_cell_type, ignore.case=TRUE)]
plot_file_name <- 'L2_umap_neurons.png'
plot_cmd <- plot_cells(neurons_cds, color_cells_by="partition")
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# You can run graph_test() like this:
pr_graph_test_res <- graph_test(neurons_cds, neighbor_graph="knn", cores=4)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

# You can call find_gene_modules(), which essentially runs UMAP on the genes (as opposed to the cells) and then groups them into modules using Louvain community analysis:
gene_module_df <- find_gene_modules(neurons_cds[pr_deg_ids,], resolution=1e-2)

# The first is just to make a simple table that shows the aggregate expression of all genes in each module across all the clusters. Monocle provides a simple utility function called aggregate_gene_expression for this purpose:
cell_group_df <- tibble::tibble(cell=row.names(colData(neurons_cds)), cell_group=partitions(cds)[colnames(neurons_cds)])
agg_mat <- aggregate_gene_expression(neurons_cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))
plot_file_name <- 'L2_neuron_module_heatmap.png'
message('plot ', plot_file_name)
plot_cmd <- pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE, scale="column", clustering_method="ward.D2", fontsize=6)
plot_cmd
png(file.path(manual_images_dir, plot_file_name), width=1500, height=2400, res=300)
plot_cmd
dev.off()

# If there are many modules, it can be hard to see where each one is expressed, so we'll just look at a subset of them:
plot_file_name <- 'L2_umap_neurons_selected_modules.png'
plot_cmd <- plot_cells(neurons_cds, genes=gene_module_df %>% filter(module %in% c(8, 28, 33, 37)), group_cells_by="partition", color_cells_by="partition", show_trajectory_graph=FALSE)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# Let's return to the embryo data:

# We will load it as we did with the L2 data:
expression_matrix <- readRDS(url("https://staff.washington.edu/hpliner/data/packer_embryo_expression.rds"))
cell_metadata <- readRDS(url("https://staff.washington.edu/hpliner/data/packer_embryo_colData.rds"))
gene_annotation <- readRDS(url("https://staff.washington.edu/hpliner/data/packer_embryo_rowData.rds"))

cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

# Note: Your data will not have the loading batch information demonstrated here, you will correct batch using your own batch information.
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")

# However, unlike clustering, which works well with both UMAP and t-SNE, here we strongly urge you to use UMAP, the default method:
cds <- reduce_dimension(cds)

# I believe I need to cluster cells and learn_graph here, so...
ciliated_genes <- c("che-1", "hlh-17", "nhr-6", "dmd-6", "ceh-36", "ham-1")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)

plot_file_name <- 'embryo_pr_graph_packer_cell_type.png'
plot_cmd <- plot_cells(cds, color_cells_by = "cell.type", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# Once again, we turn to graph_test(), this time passing it neighbor_graph="principal_graph", which tells it to test whether cells at similar positions on the trajectory have correlated expression:
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))

# Here are a couple of interesting genes that score as highly significant according to graph_test():
plot_file_name <- 'embryo_umap_selected_markers.png'
plot_cmd <- plot_cells(cds, genes=c("hlh-4", "gcy-8", "dac-1", "oig-8"), show_trajectory_graph=FALSE, label_cell_groups=FALSE, label_leaves=FALSE)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# As before, we can collect the trajectory-variable genes into modules:
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6,-1)))

# Here we plot the aggregate module scores within each group of cell types as annotated by Packer & Zhu et al:
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$cell.type)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
plot_file_name <- 'emb_pseudotime_module_heatmap.png'
message('plot ', plot_file_name)
plot_cmd <- pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2", fontsize=6)
plot_cmd
png(file.path(manual_images_dir, plot_file_name), width=1500, height=2400, res=300)
plot_cmd
dev.off()

# We can also pass gene_module_df to plot_cells() as we did when we compared clusters in the L2 data above.
plot_file_name <- 'embryo_umap_selected_pseudotime_modules.png'
plot_cmd <- plot_cells(cds, genes=gene_module_df %>% filter(module %in% c(27, 10, 7, 30)), label_cell_groups=FALSE, show_trajectory_graph=FALSE)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# You can select a path with choose_cells() or by subsetting the cell data set by cluster, cell type, or other annotation that's restricted to the path. Let's pick one such path, the AFD cells:
AFD_genes <- c("gcy-8", "dac-1", "oig-8")
AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes, colData(cds)$cell.type %in% c("AFD")]

# The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:
plot_file_name <- 'embryo_AFD_dynamic_genes.png'
plot_cmd <- plot_genes_in_pseudotime(AFD_lineage_cds, color_cells_by="embryo.time.bin", min_expr=0.5)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)


# Doing so is as simple as selecting the cells (and branch point) of interest with choose_cells():
cds_subset <- choose_cells(cds)

# This will identify genes with interesting patterns of expression that fall only within the region of the trajectory you selected, giving you a more refined and relevant set of genes.
subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))

# Grouping these genes into modules can reveal fate specific genes or those that are activate immediate prior to or following the branch point:
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.001)

# We will organize the modules by their similarity (using hclust) over the trajectory so it's a little easier to see which ones come on before others:
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, levels = row.names(agg_mat)[module_dendro$order])
plot_file_name <- 'embryo_umap_AFD_modules.png'
plot_cmd <- plot_cells(cds_subset, genes=gene_module_df, label_cell_groups=FALSE, show_trajectory_graph=FALSE)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

