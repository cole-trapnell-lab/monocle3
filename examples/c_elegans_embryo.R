rm(list = ls())  # Clear the environment
#options(warn=-1) # Turn off warning message globally
library(monocle3) # Load Monocle
library(tidymodels)

theme_set(theme_gray(base_size = 6))

expression_matrix = readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_expression.rds"))
cell_metadata = readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_colData.rds"))
gene_annotation = readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_rowData.rds"))

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

## Step 1: Normalize and pre-process the data
#cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = 50)

## Step 2: Correct for batch effects (optional)
cds <- align_cds(cds, 
                 residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading",
                 alignment_group="batch")

set.seed(42)
## Step 2: Reduce the dimensionality of the data
cds <- reduce_dimension(cds, umap.fast_sgd = FALSE, cores=1)
plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "cell.type") + ggsave("embryo_umap_packer_cell_type.png", width=5, height=4, dpi = 600)
#plot_cells(cds, color_cells_by = "batch", label_cell_groups=FALSE) + ggsave("embryo_umap_batch.png", width=5, height=4, dpi = 600)

## Step 3: Cluster cells
#cds <- cluster_cells(cds, resolution=c(0, 1e-5, 1e-4, 1e-3, 1e-2))
cds <- cluster_cells(cds, resolution=5e-6)
plot_cells(cds, group_cells_by="partition", color_cells_by = "partition")+ ggsave("embryo_umap_partition.png", width=5, height=4, dpi = 600)

## Step 4: Learn cell trajectories
#
cds <- learn_graph(cds)

# Default is too coarse, need to set ncenter explicitly:
cds <- learn_graph(cds, learn_graph_control=list(ncenter=1000))

plot_cells(cds,
           color_cells_by = "cell.type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE) + ggsave("embryo_pr_graph_packer_cell_type.png", width=5, height=4, dpi = 600)


plot_cells(cds,
           color_cells_by = "embryo.time.bin",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5) + ggsave("embryo_pr_graph_by_time.png", width=5, height=4, dpi = 600)

cds = order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots=FALSE,
           graph_label_size=1.5) + ggsave("embryo_pr_graph_by_pseudotime.png", width=5, height=4, dpi = 600)


plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

root_group = colnames(cds)[clusters(cds) == 1]
cds = order_cells(cds, root_cells = root_group)


# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="130-170"){
  cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)

  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
      (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}

cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5) + ggsave("embryo_pr_graph_by_pseudotime_programmatically_ordered.png", width=5, height=4, dpi = 600)

#plot_cells(cds, genes=c("egl-21", "egl-1"))
#plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name %in% c("egl-21", "egl-1"),])

ciliated_genes = c("che-1",
                   "hlh-17",
                   "nhr-6",
                   "dmd-6",
                   "ceh-36",
                   "ham-1")
plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE) + ggsave("embryo_ciliated_markers.png", width=5, height=4, dpi = 600)

#graph_test(cds_subset)


plot_percent_cells_positive(cds[rowData(cds)$gene_short_name %in% ciliated_genes,], group_cells_by="time.point") +
    guides(fill=FALSE) +
    theme(axis.text.x=element_text(angle=45, hjust=1))

# plot_genes_violin(cds[rowData(cds)$gene_short_name %in% ciliated_genes,], color_cells_by="cell.type") +
#     guides(fill=FALSE) +
#     theme(axis.text.x=element_text(angle=45, hjust=1))

### Basic differential expression:
cds_subset = cds[rowData(cds)$gene_short_name %in% ciliated_genes,]

gene_fits = fit_models(cds_subset, model_formula_str = "~embryo.time")
fit_coefs = coefficient_table(gene_fits)
write.csv(fit_coefs, "emb_terms.csv")

emb_time_terms = fit_coefs %>% filter(term == "embryo.time")
write.csv(emb_time_terms, "emb_time_terms.csv")

sig_emb_time_terms = emb_time_terms %>% filter (q_value < 0.05) %>% select(gene_short_name, term, q_value, estimate)
write.csv(sig_emb_time_terms, "emb_time_sig_terms.csv")

plot_genes_violin(cds[rowData(cds)$gene_short_name %in% c("che-1",
                   "nhr-6",
                   "dmd-6",
                   "ceh-36")], group_cells_by="embryo.time.bin", ncol=2) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) + ggsave("embryo_ciliated_markers_violin.png", width=4, height=4, dpi = 600)

### Subtracting unwanted effects:

gene_fits = fit_models(cds_subset, model_formula_str = "~embryo.time + batch")
fit_coefs = coefficient_table(gene_fits)
emb_terms = fit_coefs %>% filter(term != "(Intercept)")
write.csv(emb_terms, "emb_plus_batch_terms.csv")

### Evaluating fit

fit_evals = evaluate_fits(gene_fits)
write.csv(fit_evals, "emb_time_evals.csv")

gene_fits = fit_models(cds[rowData(cds)$gene_short_name %in% ciliated_genes,], model_formula_str = "~embryo.time + batch")
fit_coefs = coefficient_table(gene_fits)
emb_time_terms = fit_coefs %>% filter(term == "embryo.time")
emb_time_terms = emb_time_terms %>% mutate(q_value = stats::p.adjust(p_value))
emb_time_terms %>% filter (q_value < 0.05) %>% select(gene_short_name, term, q_value, estimate)

### Comparing models

time_models = fit_models(cds_subset, model_formula_str = "~embryo.time", expression_family="negbinomial")
time_batch_models = fit_models(cds_subset, model_formula_str = "~embryo.time + batch", expression_family="negbinomial")

emb_model_lr_test = compare_models(time_batch_models, time_models)
write.csv(emb_model_lr_test, "emb_model_lr_test.csv")


cds_subset = cds[rowData(cds)$gene_short_name %in% ciliated_genes,
                 is.finite(colData(cds)$pseudotime)]

gene_fits = fit_models(cds_subset, model_formula_str = "~splines::ns(pseudotime, df=3)", verbose=TRUE)
fit_coefs = coefficient_table(gene_fits)
emb_time_terms = fit_coefs %>% filter(grepl("pseudotime", term))
emb_time_terms = emb_time_terms %>% mutate(q_value = stats::p.adjust(p_value))
emb_time_terms %>% filter (q_value < 0.05) %>% select(gene_short_name, term, q_value, estimate)

# Principal graph test:

ciliated_cds_pr_test_res = graph_test(cds,  neighbor_graph="principal_graph", cores=4)
pr_deg_ids = row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))

gene_module_df = monocle3:::find_gene_modules(cds[pr_deg_ids,], resolution=1e-3)
cell_group_df = tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$cell.type)
agg_mat = aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                    scale="column", clustering_method="ward.D2",
                    fontsize=6, width=5, height=8,
                    file="emb_pseudotime_module_heatmap.png")

plot_cells(cds,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE) + ggsave("embryo_umap_selected_pseudotime_modules.png", width=5, height=4, dpi = 600)

rowData(cds)[gene_module_df %>% filter(module == 29) %>% pull(id),]

plot_cells(cds, genes=gene_module_df, color_cells_by="cell.type", show_trajectory_graph=FALSE)

plot_cells(cds, genes=c("hlh-4", "gcy-8", "dac-1", "oig-8"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE) + ggsave("embryo_umap_selected_markers.png", width=5, height=4, dpi = 600)


plot_cells(cds, show_trajectory_graph=FALSE) + ggsave("embryo_umap_cluster.png", width=5, height=4, dpi = 600)


AFD_lineage_cds = cds[,clusters(cds) %in% c(22, 28, 35)]
plot_genes_in_pseudotime(AFD_lineage_cds[rowData(AFD_lineage_cds)$gene_short_name %in% c("gcy-8", "dac-1", "oig-8"),],
                         color_cells_by="embryo.time.bin",
                         min_expr=0.5) +
    ggsave("embryo_AFD_dynamic_genes.png", width=3, height=3, dpi = 600)


plot_cells(cds, genes=c("F29C4.2", "grld-1", "csn-5", "vamp-7"), label_branch_points=FALSE, label_roots=FALSE, label_leaves=FALSE)


######## Branch analysis:

cds_subset = choose_cells(cds)

subset_pr_test_res = graph_test(cds_subset, cores=4)
pr_deg_ids = row.names(subset(subset_pr_test_res, q_value < 0.05))

gene_module_df = monocle3:::find_gene_modules(cds_subset[pr_deg_ids,], resolution=1e-2)
cell_group_df = tibble::tibble(cell=row.names(colData(cds_subset)), cell_group=clusters(cds_subset)[colnames(cds_subset)])
agg_mat = aggregate_gene_expression(cds_subset, gene_module_df, cell_group_df)
gene_module_df$module <- factor(gene_module_df$module, levels = row.names(agg_mat)[module_dendro$order])

plot_cells(cds_subset,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE) + ggsave("embryo_umap_AFD_modules.png", width=5, height=4, dpi = 600)

subset_pr_test_res = inner_join(subset_pr_test_res, gene_module_df)
#top_diff_genes = subset_pr_test_res %>% group_by(module) %>% filter(q_value < 0.05) %>% top_n(5, morans_I)
top_AFD_genes = subset_pr_test_res %>% filter(q_value < 0.05 & module == 3) %>% top_n(3, morans_I)
top_AWC_genes = subset_pr_test_res %>% filter(q_value < 0.05 & module %in% c(11, 7)) %>% top_n(3, morans_I)

plot_cells(cds_subset,
           genes=c(top_AFD_genes %>% pull(gene_short_name), top_AWC_genes %>% pull(gene_short_name)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE) + ggsave("embryo_umap_AFD_AWC_branch_genes.png", width=5, height=4, dpi = 600)

###### Making 3D trajectories:

cds_3d = reduce_dimension(cds, max_components = 3)
cds_3d = cluster_cells(cds_3d)
cds_3d = learn_graph(cds_3d)
cds_3d = order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))

cds_3d_plot_obj = plot_cells_3d(cds_3d, color_cells_by="partition")

htmlwidgets::saveWidget(cds_3d_plot_obj, "emb_3d_by_partition.html")
