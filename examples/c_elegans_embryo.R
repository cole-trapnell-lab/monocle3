rm(list = ls())  # Clear the environment
#options(warn=-1) # Turn off warning message globally
library(monocle3) # Load Monocle
library(tidymodels)

theme_set(theme_gray(base_size = 6))

# color_metadata = read.delim("/Users/coletrap/dropbox_lab/Analysis/worm-lineage/animations/ciliated_neurons_3D_umap_metadata.txt")
# color_metadata = color_metadata %>% dplyr::select(plot.cell.type, color.for.movie) %>% distinct()
# color_map = as.character(color_metadata$color.for.movie)
# names(color_map) = as.character(color_metadata$plot.cell.type)

# expression_matrix = readRDS("for-cole.ciliated.amphid.neurons.exprs.rds")
# cell_metadata = readRDS("for-cole.ciliated.amphid.neurons.pData.rds")
# gene_annotation = readRDS("for-cole.ciliated.amphid.neurons.fData.rds")

expression_matrix = readRDS(url("http://jpacker-data.s3.amazonaws.com/for-cole/for-cole.ciliated.amphid.neurons.exprs.rds"))
cell_metadata = readRDS(url("http://jpacker-data.s3.amazonaws.com/for-cole/for-cole.ciliated.amphid.neurons.pData.rds"))
gene_annotation = readRDS(url("http://jpacker-data.s3.amazonaws.com/for-cole/for-cole.ciliated.amphid.neurons.fData.rds"))
gene_annotation$use_for_ordering = NULL

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

## Step 1: Normalize and pre-process the data
#cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = 75, residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")

set.seed(42)
## Step 2: Reduce the dimensionality of the data
cds <- reduce_dimension(cds, umap.fast_sgd = FALSE, cores=1)
plot_cells(cds, color_by = "cell.type")
plot_cells(cds, color_by = "batch", label_cell_groups=FALSE)

## Step 3: Cluster cells
#cds <- cluster_cells(cds, resolution=c(0, 1e-5, 1e-4, 1e-3, 1e-2))
cds <- cluster_cells(cds)
plot_cells(cds, color_by = "Cluster")

## Step 4: Learn cell trajectories
#cds <- learn_graph(cds, learn_graph_control=list(ncenter=1000))
cds <- learn_graph(cds)
plot_cells(cds, color_by = "cell.type")

plot_cells(cds, color_by = "embryo.time.bin", label_cell_groups=FALSE)


cds = order_cells(cds)
plot_cells(cds, color_by = "Pseudotime")

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

plot_cells(cds, color_by = "Pseudotime")

## Step 5: Visualize the trajectory


plot_cells(cds, color_by = "embryo.time.bin", label_cell_groups=FALSE)


#plot_cells(cds, genes=c("egl-21", "egl-1"))
#plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name %in% c("egl-21", "egl-1"),])

ciliated_genes = c("che-1",
                   "hlh-17",
                   "nhr-6",
                   "dmd-6",
                   "ceh-36",
                   "ham-1")
cds_subset = cds[rowData(cds)$gene_short_name %in% ciliated_genes,]

plot_cells(cds, genes=ciliated_genes, show_trajectory_graph=FALSE)

#graph_test(cds_subset)


# plot_percent_cells_positive(cds_subset, color_by="cell.type") +
#     guides(fill=FALSE) +
#     theme(axis.text.x=element_text(angle=45, hjust=1))

# plot_genes_violin(cds[rowData(cds)$gene_short_name %in% ciliated_genes,], color_by="cell.type") +
#     guides(fill=FALSE) +
#     theme(axis.text.x=element_text(angle=45, hjust=1))

### Basic differential expression:

gene_fits = fit_models(cds_subset, model_formula_str = "~embryo.time")
fit_coefs = coefficient_table(gene_fits)
write.csv(fit_coefs, "emb_terms.csv")

emb_time_terms = emb_terms %>% filter(term == "embryo.time")
write.csv(emb_time_terms, "emb_time_terms.csv")

sig_emb_time_terms = emb_time_terms %>% filter (q_value < 0.05) %>% select(gene_short_name, term, q_value, estimate)
write.csv(sig_emb_time_terms, "emb_time_sig_terms.csv")

plot_genes_violin(cds_subset, color_by="embryo.time.bin") +
    theme(axis.text.x=element_text(angle=45, hjust=1))

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
emb_time_terms = emb_time_terms %>% mutate(q_value = p.adjust(p_value))
emb_time_terms %>% filter (q_value < 0.05) %>% select(gene_short_name, term, q_value, estimate)

### Comparing models

time_models = fit_models(cds_subset, model_formula_str = "~embryo.time", expression_family="negbinomial")
time_batch_models = fit_models(cds_subset, model_formula_str = "~embryo.time + batch", expression_family="negbinomial")

emb_model_lr_test = compare_models(time_batch_models, time_models)
write.csv(emb_model_lr_test, "emb_model_lr_test.csv")


cds_subset = cds[rowData(cds)$gene_short_name %in% ciliated_genes,
                 is.finite(colData(cds)$Pseudotime)]

gene_fits = fit_models(cds_subset, model_formula_str = "~splines::ns(Pseudotime, df=3)", verbose=TRUE)
fit_coefs = coefficient_table(gene_fits)
emb_time_terms = fit_coefs %>% filter(grepl("Pseudotime", term))
emb_time_terms = emb_time_terms %>% mutate(q_value = p.adjust(p_value))
emb_time_terms %>% filter (q_value < 0.05) %>% select(gene_short_name, term, q_value, estimate)

# Principal graph test:

ciliated_cds_pr_test_res = graph_test(cds, cores=4)
pr_deg_ids = row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))

gene_cluster_df = monocle3:::cluster_genes(cds[pr_deg_ids,], resolution=c(0,10^seq(-6,-1)))
cell_group_df = tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$cell.type)
pheatmap::pheatmap(aggregate_gene_expression(cds, gene_cluster_df, cell_group_df),
                    scale="column", clustering_method="ward.D2")

fData(cds)[gene_cluster_df %>% filter(cluster == 2) %>% pull(id),]

plot_cells(cds, genes=gene_cluster_df, color_by="cell.type", show_trajectory_graph=FALSE)

plot_cells(cds, genes=c("ctc-3", "nduo-4", "T05C7.4", "atp-6", "cutl-24"), label_branch_points=FALSE, label_roots=FALSE, label_leaves=FALSE)

plot_cells(cds, genes=c("F29C4.2", "grld-1", "csn-5", "vamp-7"), label_branch_points=FALSE, label_roots=FALSE, label_leaves=FALSE)


######## Branch analysis:

cds_subset = choose_cells(cds)

subset_pr_test_res = graph_test(cds_subset, cores=4)
pr_deg_ids = row.names(subset(subset_pr_test_res, q_value < 0.05))

gene_cluster_df = monocle3:::cluster_genes(cds_subset[pr_deg_ids,], resolution=c(0,10^seq(-6,-1)))
png("branch_modules.png", res=600, width=6, height=6, units="in")
plot_cells(cds_subset, genes=gene_cluster_df, show_trajectory_graph=FALSE, color_by="cell.type")
dev.off()
