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

plot_cells(cds, color_cells_by="cao_cell_type")  + ggsave("L2_umap_corrected_cao_type.png", width=5, height=4, dpi = 600)
plot_cells(cds, color_cells_by="cao_cell_type", label_groups_by_cluster=FALSE)  + ggsave("L2_umap_corrected_cao_type_no_cluster_label.png", width=5, height=4, dpi = 600)

## Step 4: Find genes expressed in each cluster or cell type
marker_test_res = top_markers(cds, group_cells_by="cluster", cores=8)

top_specific_markers = as.character(unique(marker_test_res %>% 
    filter(marker_test_q_value < 0.05) %>%
    group_by(cell_group) %>% 
    top_n(1, pseudo_R2) %>% pull(gene_id)))

plot_genes_by_group(cds, top_specific_markers, group_cells_by="cao_cell_type", ordering_type="maximal_on_diag")

## TODO: need to use choose_cells() here to pick some specific clusters:
pr_graph_test_res = graph_test(cds, neighbor_graph="knn", cores=8)
pr_deg_ids = row.names(subset(pr_graph_test_res, q_value < 0.05))

gene_cluster_df = cluster_genes(cds[pr_deg_ids,], resolution=c(0,10^seq(-6,-1)))

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
