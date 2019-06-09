rm(list = ls())  # Clear the environment
#options(warn=-1) # Turn off warning message globally
library(monocle3) # Load Monocle
library(tidymodels)
theme_set(theme_gray(base_size = 6))

# expression_matrix = readRDS("L2_expression.rda")
# cell_metadata = readRDS("L2_cell_metadata.rda")
# gene_annotation = readRDS("L2_gene_metadata.rda")

expression_matrix = readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_expression.rds"))
cell_metadata = readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_colData.rds"))
gene_annotation = readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_rowData.rds"))

cds = new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


## Step 1: Normalize and pre-process the data
cds = preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds) + ggsave("L2_pc_variance_explained.png", width=3, height=2, dpi = 600)

## Step 2: Reduce the dimensionality of the data
#### Without batch correction:
cds = reduce_dimension(cds)
plot_cells(cds) + ggsave("L2_umap_no_color.png", width=5, height=4, dpi = 600)
plot_cells(cds, color_cells_by="cao_cell_type") + ggsave("L2_umap_color_cells_by_cao_type.png", width=5, height=4, dpi = 600)
plot_cells(cds, color_cells_by="plate", label_cell_groups=FALSE) + ggsave("L2_umap_plate.png", width=5, height=4, dpi = 600)
plot_cells(cds, genes=c("cpna-2", "egl-21", "ram-2", "inos-1")) + ggsave("L2_umap_gene_markers.png", width=5, height=4, dpi = 600)


#### With batch correction:
cds = preprocess_cds(cds, num_dim = 100, residual_model_formula_str = "~ plate")
cds = reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1)
plot_cells(cds, color_cells_by="plate", label_cell_groups=FALSE) + ggsave("L2_umap_corrected_plate.png", width=5, height=4, dpi = 600)

# exmaple of using t-SNE:
cds = reduce_dimension(cds, reduction_method="tSNE")
plot_cells(cds, reduction_method="tSNE", color_cells_by="cao_cell_type") + ggsave("L2_tsne_corrected_cao_type.png", width=5, height=4, dpi = 600)


## Step 3: (Optional) Cluster cells
cds = cluster_cells(cds, resolution=c(10^seq(-6,-1)))
plot_cells(cds) + ggsave("L2_umap_color_cells_by_cluster.png", width=5, height=4, dpi = 600)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition") + ggsave("L2_umap_color_cells_by_partition.png", width=5, height=4, dpi = 600)

plot_cells(cds, color_cells_by="cao_cell_type")  + ggsave("L2_umap_corrected_cao_type.png", width=5, height=4, dpi = 600)
plot_cells(cds, color_cells_by="cao_cell_type", label_groups_by_cluster=FALSE) + ggsave("L2_umap_corrected_cao_type_no_cluster_label.png", width=5, height=4, dpi = 600)

## Step 4: Compare and contrast cells in different clusters
# Identify top marker genes for each cluster
marker_test_res = top_markers(cds, group_cells_by="partition", reference_cells=1000, cores=8)

top_specific_markers = marker_test_res %>%
    filter(fraction_expressing >= 0.10) %>%
    group_by(cell_group) %>%
    top_n(1, pseudo_R2)

top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=3) + ggsave("L2_plot_top_partition_marker.png", width=5, height=4, dpi = 600)

top_specific_markers = marker_test_res %>%
    filter(fraction_expressing >= 0.10) %>%
    group_by(cell_group) %>%
    top_n(3, pseudo_R2)

top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="cluster_row_col",
                    max.size=3) + ggsave("L2_plot_top3_partition_marker.png", width=5, height=6, dpi = 600)

## Cheat a bit here by automatically generating assignments that would stem from
# manual annotation. We'll do this so that when the clusters change we don't have
# manually go through them. Obviously this isn't possible when you don't already
# the answers!
type_partition_olap = as.matrix(prop.table(table(partitions(cds), colData(cds)$cao_cell_type), margin=1))
pheatmap::pheatmap(type_partition_olap, cluster_rows=FALSE, cluster_cols=FALSE)

paritition_identities = max.col(type_partition_olap,ties.method="first")
colnames(type_partition_olap)[paritition_identities]
for (i in (1:length(paritition_identities))){
    if (i < length(paritition_identities)){
        cat('"',i,'"="',colnames(type_partition_olap)[paritition_identities[i]], '",\n', sep="")
    } else {
       cat('"',i,'"="',colnames(type_partition_olap)[paritition_identities[i]], '"\n', sep="")
    }
}

## Step 5: Annotate cells by type:
colData(cds)$assigned_cell_type = as.character(partitions(cds))
colData(cds)$assigned_cell_type = dplyr::recode(colData(cds)$assigned_cell_type,
"1"="Body wall muscle",
"2"="Germline",
"3"="Unclassified neurons",
"4"="Seam cells",
"5"="Coelomocytes",
"6"="Pharyngeal epithelia",
"7"="Vulval precursors",
"8"="Non-seam hypodermis",
"9"="Intestinal/rectal muscle",
"10"="Touch receptor neurons",
"11"="Pharyngeal neurons",
"12"="Am/PH sheath cells",
"13"="NA",
"14"="flp-1(+) interneurons",
"15"="Canal associated neurons",
"16"="Pharyngeal gland",
"17"="Unclassified neurons",
"18"="Ciliated sensory neurons",
"19"="Ciliated sensory neurons",
"20"="Ciliated sensory neurons",
"21"="Ciliated sensory neurons",
"22"="Ciliated sensory neurons",
"23"="Ciliated sensory neurons",
"24"="Ciliated sensory neurons",
"25"="Oxygen sensory neurons",
"26"="Ciliated sensory neurons",
"27"="Unclassified neurons",
"28"="Pharyngeal gland",
"29"="Ciliated sensory neurons",
"30"="Ciliated sensory neurons",
"31"="Ciliated sensory neurons",
"32"="Ciliated sensory neurons",
"33"="Failed QC",
"34"="Pharyngeal muscle")
plot_cells(cds, group_cells_by="partition", color_cells_by="assigned_cell_type") + ggsave("L2_plot_cells_by_initial_annotation.png", width=5, height=4, dpi = 600)

# # Make a violin that
# plot_genes_violin(cds[rowData(cds)$gene_short_name %in% c("mec-7", "mec-17", "mig-6", "flp-1", "nlp-12"),
#                       colData(cds)$cao_cell_type %in% c("GABAergic neurons", "Ciliated sensory neurons", "Cholinergic neurons", "Pharyngeal neurons")],
#                       group_cells_by="cao_cell_type", ncol=1) +
#     theme(axis.text.x=element_text(angle=45, hjust=1)) + ggsave("L2_interneuron_violin.png", width=4, height=4, dpi = 600)


# Drilling into specific subsets of cells
# Draw bounding box around partition 7:
cds_subset = choose_cells(cds)

pr_graph_test_res = graph_test(cds_subset, neighbor_graph="knn", cores=8)
pr_deg_ids = row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))

gene_module_df = find_gene_modules(cds_subset[pr_deg_ids,], resolution=1e-3)
#c(0,10^seq(-6,-1))
plot_cells(cds_subset, genes=gene_module_df, show_trajectory_graph=FALSE, label_cell_groups=FALSE) +
    ggsave("L2_sex_partition_modules.png", width=5, height=4, dpi = 600)

plot_cells(cds_subset, color_cells_by="cluster") +
    ggsave("L2_sex_partition_color_by_cluster.png", width=5, height=4, dpi = 600)

## Cheat a bit here by automatically generating assignments that would stem from
# manual annotation. We'll do this so that when the clusters change we don't have
# manually go through them. Obviously this isn't possible when you don't already
# the answers!
type_partition_olap = as.matrix(prop.table(table(clusters(cds_subset), colData(cds_subset)$cao_cell_type), margin=1))
pheatmap::pheatmap(type_partition_olap, cluster_rows=FALSE, cluster_cols=FALSE)

paritition_identities = max.col(type_partition_olap,ties.method="first")
colnames(type_partition_olap)[paritition_identities]
for (i in (1:length(paritition_identities))){
    if (i < length(paritition_identities)){
        cat('"',i,'"="',colnames(type_partition_olap)[paritition_identities[i]], '",\n', sep="")
    } else {
       cat('"',i,'"="',colnames(type_partition_olap)[paritition_identities[i]], '"\n', sep="")
    }
}

colData(cds_subset)$assigned_cell_type = as.character(clusters(cds_subset)[colnames(cds_subset)])
colData(cds_subset)$assigned_cell_type = dplyr::recode(colData(cds_subset)$assigned_cell_type,
                                                "27"="Sex myoblasts",
                                                "28"="Vulval precursors",
                                                "31"="Somatic gonad precursors",
                                                "38"="Unclassified neurons",
                                                "42"="Sex myoblasts",
                                                "39"="Failed QC"
                                                )
plot_cells(cds_subset, group_cells_by="cluster", color_cells_by="assigned_cell_type")
colData(cds)[colnames(cds_subset),]$assigned_cell_type = colData(cds_subset)$assigned_cell_type
cds = cds[,colData(cds)$assigned_cell_type != "Failed QC" | is.na(colData(cds)$assigned_cell_type )]
plot_cells(cds, group_cells_by="partition", color_cells_by="assigned_cell_type", labels_per_group=5) + ggsave("L2_plot_cells_by_refined_annotation.png", width=5, height=4, dpi = 600)

# Plot the overlap between clusters and annotated cell types:
pheatmap::pheatmap(log(table(colData(cds)$assigned_cell_type, colData(cds)$cao_cell_type)+1),
                   clustering_method="ward.D2",
                   fontsize=6, width=5, height=5, filename="L2_cao_cell_type_by_manual_type.png")



#ceModel = readRDS(url("https://cole-trapnell-lab.github.io/garnett/classifiers/ceWhole"))

assigned_type_marker_test_res = top_markers(cds[,is.na(colData(cds)$assigned_cell_type) == FALSE & colData(cds)$assigned_cell_type != "Failed QC"],
                                            group_cells_by="assigned_cell_type",
                                            reference_cells=1000,
                                            cores=8)

garnett_markers = assigned_type_marker_test_res %>%
    filter(marker_test_q_value < 0.01 & specificity >= 0.5) %>%
    group_by(cell_group) %>%
    top_n(5, marker_score)
garnett_markers = garnett_markers %>% group_by(gene_short_name) %>% 
    filter(n() == 1) 

plot_genes_by_group(cds,
                    unique(as.character(garnett_markers$gene_id)),
                    group_cells_by="assigned_cell_type",
                    ordering_type="cluster_row_col",
                    max.size=3) + ggsave("L2_plot_top3_cao_type_markers.png", width=5, height=6, dpi = 600)


generate_garnett_marker_file(garnett_markers)


library(garnett)

worm_classifier <- train_cell_classifier(cds = cds,
                                         marker_file = "./marker_file.txt",
                                         db=org.Ce.eg.db::org.Ce.eg.db,
                                         cds_gene_id_type = "ENSEMBL",
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "SYMBOL",
                                         cores=8)

colData(cds)$garnett_cluster = clusters(cds)
cds = classify_cells(cds, worm_classifier,
                           db = org.Ce.eg.db::org.Ce.eg.db,
                           cluster_extend = TRUE,
                           cds_gene_id_type = "ENSEMBL")

plot_cells(cds, 
           group_cells_by="partition", 
           color_cells_by="cluster_ext_type") + ggsave("L2_umap_corrected_garnett_ext_type_our_annotations.png", width=5, height=4, dpi = 600)

# Plot the overlap between clusters and annotated cell types:
pheatmap::pheatmap(log(table(colData(cds)$cluster_ext_type, colData(cds)$cao_cell_type)+1),
                   clustering_method="ward.D2",
                   fontsize=6, width=5, height=5, filename="L2_cao_cell_type_by_garnett_ext_type_our_annotations.png")

load(url("https://cole-trapnell-lab.github.io/garnett/classifiers/ceWhole"))
colData(cds)$garnett_cluster = clusters(cds)
cds = classify_cells(cds, ceWhole,
                           db = org.Ce.eg.db::org.Ce.eg.db,
                           cluster_extend = TRUE,
                           cds_gene_id_type = "ENSEMBL")

plot_cells(cds, 
           group_cells_by="cluster", 
           color_cells_by="cluster_ext_type") + ggsave("L2_umap_corrected_garnett_ext_type_cao_annoutations.png", width=5, height=4, dpi = 600)

# Plot the overlap between clusters and annotated cell types:
pheatmap::pheatmap(prop.table(table(partitions(cds)[colnames(cds)], colData(cds)$cluster_ext_type), margin=1),
                   #clustering_method="ward.D2",
                   cluster_cols=FALSE,
                   cluster_rows=FALSE,
                   fontsize=6, width=5, height=8, filename="L2_cell_type_by_partition.png")



########
# Globally clustering genes into modules:

neurons_cds = cds[,colData(cds)$assigned_cell_type == "Neurons"]
plot_cells(neurons_cds, color_cells_by="partition") + ggsave("L2_umap_neurons.png", width=5, height=4, dpi = 600)

pr_graph_test_res = graph_test(neurons_cds, neighbor_graph="knn", cores=8)
pr_deg_ids = row.names(subset(pr_graph_test_res, q_value < 0.05))

gene_module_df = monocle3:::find_gene_modules(neurons_cds[pr_deg_ids,], resolution=1e-2)

cell_group_df = tibble::tibble(cell=row.names(colData(neurons_cds)), cell_group=partitions(cds)[colnames(neurons_cds)])
agg_mat = aggregate_gene_expression(neurons_cds, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) = stringr::str_c("Partition ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6, width=5, height=8,
                   file="L2_neuron_module_heatmap.png")

plot_cells(neurons_cds,
           genes=gene_module_df %>% filter(module %in% c(16,38,33,42)),
           cell_size=1,
           group_cells_by="partition",
           color_cells_by="partition",
           show_trajectory_graph=FALSE) + ggsave("L2_umap_neurons_selected_modules.png", width=5, height=4, dpi = 600)


# ########
# # Look at a specific set of clusters

# cds_subset = choose_cells(cds)

# pr_graph_test_res = graph_test(cds_subset, neighbor_graph="knn", cores=8)
# pr_deg_ids = row.names(subset(pr_graph_test_res, q_value < 0.05))

# gene_module_df = monocle3:::find_gene_modules(cds[pr_deg_ids,], resolution=c(0,10^seq(-6,-1)), verbose=TRUE)
# cell_group_df = tibble::tibble(cell=row.names(colData(cds)), cell_group=clusters(cds))
# agg_mat = aggregate_gene_expression(cds, gene_module_df, cell_group_df)
# #agg_mat
# pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
#                     scale="column", clustering_method="ward.D2", fontsize=6)

# png("branch_modules.png", res=600, width=6, height=6, units="in")
# plot_cells(cds_subset, genes=gene_module_df, show_trajectory_graph=FALSE, label_cell_groups=FALSE)
# dev.off()


# mod_df = data.frame(resolution= c(0,10^seq(-6,-1)),
#                     modularity = c(0.0463080808623619, 0.0463080808623619, 0.519517905711304, 0.722441186959141, 0.893506647690171, 0.926503643050994, 0.861890178533872),
#                     clusters = c(5, 5, 6, 8, 18, 60, 179))

# qplot(clusters, modularity, data=mod_df)

# qplot(diff(mod_df$modularity), diff(mod_df$clusters))

# qplot(10^seq(-6,-1), diff(mod_df$modularity)/diff(mod_df$clusters), log="x")
