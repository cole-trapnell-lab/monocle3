# This script creates the plots used by the clustering.html page
# and stores the new files in the directory given by the bash
# variable 'manual_images_dir', as defined in globals.R.
#
# Call this script from R using the command
#   source('clustering.R', local=TRUE, echo=TRUE)
#
# There is information about using ggsave in the blog
# page at
#
#   URL: https://minimaxir.com/2017/08/ggplot2-web/
#
# A plot template looks like
#
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

# You can load the data into Monocle 3 like this: 
expression_matrix <- readRDS(url("https://staff.washington.edu/hpliner/data/cao_l2_expression.rds"))
cell_metadata <- readRDS(url("https://staff.washington.edu/hpliner/data/cao_l2_colData.rds"))
gene_annotation <- readRDS(url("https://staff.washington.edu/hpliner/data/cao_l2_rowData.rds"))

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


# When using PCA, you should specify the number of principal components you want Monocle to compute.
cds <- preprocess_cds(cds, num_dim = 100)

# It's a good idea to check that you're using enough PCs to capture most of the variation in gene expression across all the cells in the data set. You can look at the fraction of variation explained by each PC using plot_pc_variance_explained(): 
plot_file_name <- 'L2_pc_variance_explained.png'
plot_cmd <- plot_pc_variance_explained(cds)
ggplot_cells_png(plot_cmd, plot_extra='', manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# To reduce the dimensionality of the data down into the X, Y plane so we can plot it easily, call reduce_dimension(): 
cds <- reduce_dimension(cds)

# To plot the data, use Monocle's main plotting function, plot_cells():
plot_file_name <- 'L2_umap_no_color.png'
plot_cmd <- plot_cells(cds)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# We can color the cells in the UMAP plot by the authors' original annotations using the color_cells_by argument to plot_cells(). 
plot_file_name <- 'L2_umap_color_by_cao_type.png'
plot_cmd <- plot_cells(cds, color_cells_by="cao_cell_type")
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# You can also color your cells according to how much of a gene or set of genes they express:
plot_file_name <- 'L2_umap_gene_markers.png'
plot_cmd <- plot_cells(cds, genes=c("cpna-2", "egl-21", "ram-2", "inos-1"))
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# If you want, you can also use t-SNE to visualize your data. First, call reduce_dimension with reduction_method="tSNE". 
cds <- reduce_dimension(cds, reduction_method="tSNE")

# Then, when you call plot_cells(), pass reduction_method="tSNE" to it as well:
plot_file_name <- 'L2_tsne_corrected_cao_type.png'
plot_cmd <- plot_cells(cds, reduction_method="tSNE", color_cells_by="cao_cell_type")
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

#  Then you can simply color the cells by batch. Cao & Packer et al included a "plate" annotation in their data, which specifies which sci-RNA-seq plate each cell originated from. Coloring the UMAP by plate reveals:
plot_file_name <- 'L2_umap_plate.png'
plot_cmd <- plot_cells(cds, color_cells_by="plate", label_cell_groups=FALSE)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# Nevertheless, we can try and remove what batch effect is by running the align_cds() function:
cds <- align_cds(cds, num_dim = 100, alignment_group = "plate")
cds <- reduce_dimension(cds)

plot_file_name <- 'L2_umap_corrected_plate.png'
plot_cmd <- plot_cells(cds, color_cells_by="plate", label_cell_groups=FALSE)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# You can cluster your cells using the cluster_cells() function, like this: 
cds <- cluster_cells(cds, resolution=1e-5)
plot_file_name <- 'L2_umap_color_cells_by_cluster.png'
plot_cmd <- plot_cells(cds)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# You can visualize these partitions like this: 
plot_file_name <- 'L2_umap_color_cells_by_partition.png'
plot_cmd <- plot_cells(cds, color_cells_by="partition", group_cells_by="partition")
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# For example, the call below colors the cells according to their cell type annotation, and each cluster is labeled according the most common annotation within it
plot_file_name <- 'L2_umap_corrected_cao_type.png'
plot_cmd <- plot_cells(cds, color_cells_by="cao_cell_type")
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# Finally, you can disable this labeling policy, making plot_cells() behave like it did before we called cluster_cells(), like this
plot_file_name <- 'L2_umap_corrected_cao_type_no_cluster_label.png'
plot_cmd <- plot_cells(cds, color_cells_by="cao_cell_type", label_groups_by_cluster=FALSE)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# To do that, start by calling the top_markers() function:
marker_test_res <- top_markers(cds, group_cells_by="partition", reference_cells=1000, cores=4)

# For example, pseudo_R2 is one such measure. We can rank markers according to pseudo_R2 like this
top_specific_markers <- marker_test_res %>%
                            filter(fraction_expressing >= 0.10) %>%
                            group_by(cell_group) %>%
                            top_n(1, pseudo_R2)
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

# Now, we can plot the expression and fraction of cells that express each marker in each group with the plot_genes_by_group function:
plot_file_name <- 'L2_plot_top_partition_marker.png'
plot_cmd <- plot_genes_by_group(cds, top_specific_marker_ids, group_cells_by="partition", ordering_type="maximal_on_diag", max.size=3)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# It's often informative to look at more than one marker, which you can do just by changing the first argument to top_n(): 
top_specific_markers <- marker_test_res %>% filter(fraction_expressing >= 0.10) %>% group_by(cell_group) %>% top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_file_name <- 'L2_plot_top3_partition_marker.png'
plot_cmd <- plot_genes_by_group(cds, top_specific_marker_ids, group_cells_by="partition", ordering_type="cluster_row_col", max.size=3)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)


# consensus <- function(x) {
#   uniqx <- unique(na.omit(x))
#   uniqx[which.max(tabulate(match(x, uniqx)))]
# }
# 
# lpartition <- unique(partitions(cds))
# l_cell_type <- list()
# for(ipartition in lpartition) {
#   l_cell_type[as.character(ipartition)] <- consensus(colData(cds)[partitions(cds)==ipartition,][['cao_cell_type']])
# }
# colData(cds)[['assigned_cell_type']] <- as.character(l_cell_type[as.character(partitions(cds))])

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

# Let's see how the new annotations look:
plot_file_name <- 'L2_plot_cells_by_initial_annotation.png'
plot_cmd <- plot_cells(cds, group_cells_by="partition", color_cells_by="assigned_cell_type")
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# So we can isolate it with the choose_cells() function for further analysis: 
cds_subset <- choose_cells(cds)

# We can use graph_test() to identify genes that are differentially expressed in different subsets of cells from this partition: 
pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="knn", cores=4)
pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))

# We can take all the genes that vary across this set of cells and group those that have similar patterns of expression into modules: 
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=1e-3)

# Plotting these modules' aggregate expression values reveals which cells express which modules. 
plot_file_name <- 'L2_sex_partition_modules.png'
plot_cmd <- plot_cells(cds_subset, genes=gene_module_df, show_trajectory_graph=FALSE, label_cell_groups=FALSE)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# Let's recluster the cells at finer resolution and then see how they overlap with the clusters in the partition: 
cds_subset <- cluster_cells(cds_subset, resolution=1e-2)
plot_file_name <- 'L2_sex_partition_color_by_cluster.png'
plot_cmd <- plot_cells(cds_subset, color_cells_by="cluster")
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# Based on how the patterns line up, we'll make the following assignments:
# No plot png file.

# colData(cds_subset)$assigned_cell_type <- as.character(clusters(cds_subset)[colnames(cds_subset)])
# colData(cds_subset)$assigned_cell_type <- dplyr::recode(colData(cds_subset)$assigned_cell_type, "1"="Somatic gonad precursors", "2"="Somatic gonad precursors", "3"="Vulval precursors", "4"="Sex myoblasts", "5"="Sex myoblasts", "6"="Vulval precursors", "7"="Failed QC", "8"="Vulval precursors", "10"="Unclassified neurons", "11"="Distal tip cells")
# plot_cells(cds_subset, group_cells_by="cluster", color_cells_by="assigned_cell_type")
#plot_file_name <- xxxx
#message('plot ', plot_file_name)
#png(file.path(manual_images_dir, plot_file_name))
#plot_cmd <- plot_cells(cds_subset, group_cells_by="cluster", color_cells_by="assigned_cell_type")
#plot_cmd
#dev.off()

# consensus <- function(x) {
#   uniqx <- unique(na.omit(x))
#   uniqx[which.max(tabulate(match(x, uniqx)))]
# }

# l_cell_type <- list()
# colData(cds_subset)[['cluster']] <- as.character(clusters(cds_subset))
# lcluster <- unique(colData(cds_subset)[['cluster']])
# for(icluster in lcluster) {
#   l_cell_type[as.character(icluster)] <- consensus(colData(cds_subset)[colData(cds_subset)[['cluster']]==icluster,][['cao_cell_type']])
# }
# colData(cds_subset)[['assigned_cell_type']] <- as.character(clusters(cds_subset)[colnames(cds_subset)])
# colData(cds_subset)[['assigned_cell_type']] <- as.character(l_cell_type[as.character(colData(cds_subset)[['cluster']])])

colData(cds_subset)$assigned_cell_type <- as.character(clusters(cds_subset)[colnames(cds_subset)])
colData(cds_subset)$assigned_cell_type <- dplyr::recode(colData(cds_subset)$assigned_cell_type,
                                                        "1"="Sex myoblasts",
                                                        "2"="Somatic gonad precursors",
                                                        "3"="Vulval precursors",
                                                        "4"="Sex myoblasts",
                                                        "5"="Vulval precursors",
                                                        "6"="Somatic gonad precursors",
                                                        "7"="Sex myoblasts",
                                                        "8"="Sex myoblasts",
                                                        "9"="Ciliated neurons",
                                                        "10"="Vulval precursors",
                                                        "11"="Somatic gonad precursor",
                                                        "12"="Distal tip cells",
                                                        "13"="Somatic gonad precursor",
                                                        "14"="Sex myoblasts",
                                                        "15"="Vulval precursors")

# Now we can transfer the annotations from the cds_subset object back to the full dataset. We'll also filter out low-quality cells at this stage 
colData(cds)[colnames(cds_subset),]$assigned_cell_type <- colData(cds_subset)$assigned_cell_type
cds <- cds[,colData(cds)$assigned_cell_type != "Failed QC" & !is.na(colData(cds)$assigned_cell_type )]
plot_file_name <- 'L2_plot_cells_by_refined_annotation.png'
plot_cmd <- plot_cells(cds, group_cells_by="partition", color_cells_by="assigned_cell_type", labels_per_group=5)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# To generate a Garnett file, first find the top markers that each annotated cell type expresses:
assigned_type_marker_test_res <- top_markers(cds, group_cells_by="assigned_cell_type", reference_cells=1000, cores=4)

# Next, filter these markers according to how stringent you want to be:
# Require that markers have at least JS specificty score > 0.5 and
# be significant in the logistic test for identifying their cell type:
garnett_markers <- assigned_type_marker_test_res %>% filter(marker_test_q_value < 0.01 & specificity >= 0.5) %>% group_by(cell_group) %>% top_n(5, marker_score)
# Exclude genes that are good markers for more than one cell type:
garnett_markers <- garnett_markers %>% group_by(gene_short_name) %>% filter(n() == 1)
 
# Then call generate_garnett_marker_file:
generate_garnett_marker_file(garnett_markers, file="./marker_file.txt")

# will produce a text file like this:
# 
# > Cell type Ciliated sensory neurons
# expressed: che-3, scd-2, C33A12.4, R102.2, F27C1.11
# 
# > Cell type Non-seam hypodermis
# expressed: col-14, col-180, F11E6.3, grsp-1, C06A8.3
# 
# > Cell type Seam cells
# expressed: col-65, col-77, col-107, ram-2, Y47D7A.13 

# When you're ready run Garnett, load the package:
## Install the monocle3 branch of garnett
BiocManager::install(c("org.Mm.eg.db", "org.Hs.eg.db"))
devtools::install_github("cole-trapnell-lab/garnett", ref="monocle3")

library(garnett)
# install gene database for worm
BiocManager::install("org.Ce.eg.db")

# Now train a Garnett classifier based on your marker file like this:
colData(cds)$garnett_cluster <- clusters(cds)
worm_classifier <- train_cell_classifier(cds = cds, marker_file = "./marker_file.txt", db=org.Ce.eg.db::org.Ce.eg.db, cds_gene_id_type = "ENSEMBL", num_unknown = 50, marker_file_gene_id_type = "SYMBOL", cores=4)

# Now that we've trained a classifier worm_classifier, we can use it to annotate the L2 cells according to type:
cds <- classify_cells(cds, worm_classifier, db = org.Ce.eg.db::org.Ce.eg.db, cluster_extend = TRUE, cds_gene_id_type = "ENSEMBL")

# Here's how Garnett annotated the cells:
plot_file_name <- 'L2_umap_corrected_garnett_ext_type.png'
plot_cmd <- plot_cells(cds, group_cells_by="partition", color_cells_by="cluster_ext_type")
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# You can classify cells with it by first downloading and then passing it to the classify_cells() function
ceWhole <- readRDS(url("https://cole-trapnell-lab.github.io/garnett/classifiers/ceWhole_20191017.RDS"))
cds <- classify_cells(cds, ceWhole, db = org.Ce.eg.db, cluster_extend = TRUE, cds_gene_id_type = "ENSEMBL")

