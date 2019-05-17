rm(list = ls())  # Clear the environment
#options(warn=-1) # Turn off warning message globally
library(monocle3) # Load Monocle
library(tidymodels)
theme_set(theme_gray(base_size = 6))

expression_matrix = readRDS("L2_expression.rda")
cell_metadata = readRDS("L2_cell_metadata.rda")
gene_annotation = readRDS("L2_gene_metadata.rda")

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100, residual_model_formula_str = "~ plate")

## Step 2: Reduce the dimensionality of the data
cds <- reduce_dimension(cds)
plot_cells(cds, cell_size=0.35, color_by = "cao_cell_type")

## Step 3: (Optional) Cluster cells
cds <- cluster_cells(cds, resolution=c(10^seq(-6,-1)), verbose=TRUE)

plot_cells(cds)

## Step 4: Annotate cells by type:
#ceModel = readRDS(url("https://cole-trapnell-lab.github.io/garnett/classifiers/ceWhole"))
library(org.Ce.eg.db)
load(url("https://cole-trapnell-lab.github.io/garnett/classifiers/ceWhole"))

colData(cds)$garnett_cluster = get_clusters(cds)
cds <- classify_cells(cds, ceWhole,
                           db = org.Ce.eg.db,
                           cluster_extend = TRUE,
                           cds_gene_id_type = "ENSEMBL")

plot_cells(cds, color_by="cluster_ext_type")

pr_graph_test_res = graph_test(cds, neighbor_graph="knn")
pr_deg_ids = row.names(subset(pr_graph_test_res, q_value < 0.05))
gene_cluster_df = monocle3:::cluster_genes(cds[pr_deg_ids,], resolution=0.001)

text_df <- gene_cluster_df %>% dplyr::group_by(cluster) %>% dplyr::summarize(text_x = median(x = dim_1),
                                                                             text_y = median(x = dim_2))
ggplot2::qplot(dim_1, dim_2, color=cluster, data=gene_cluster_df) +
  ggplot2::geom_text(data=text_df, mapping = ggplot2::aes_string(x = "text_x", y = "text_y", label = "cluster"), color=I("black"),  size = 4)


png("module_graph.png", res=600, width=8, height=8, units="in")
monocle3:::plot_cells(cds, genes=gene_cluster_df, cell_size=0.5, show_backbone=FALSE, label_branch_points=FALSE, label_leaves=FALSE, label_roots=FALSE)
dev.off()


