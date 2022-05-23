# This script creates the plots used by the trajectories.html page
# and stores the new files in the directory given by the bash
# variable 'manual_images_dir', as defined in globals.R.
#
# Call this script from R using the command
#   source('trajectories.R', local=TRUE, echo=TRUE)
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

# Overlaying the manual annotations on the UMAP reveals that these branches are principally occupied by one cell type.
plot_file_name <- 'embryo_umap_packer_cell_type.png'
plot_cmd <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "cell.type")
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# Let's look at some genes with interesting patterns of expression in ciliated neurons:
ciliated_genes <- c("che-1", "hlh-17", "nhr-6", "dmd-6", "ceh-36", "ham-1")
plot_file_name <- 'embryo_ciliated_markers.png'
plot_cmd <- plot_cells(cds, genes=ciliated_genes, label_cell_groups=FALSE, show_trajectory_graph=FALSE)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# When you are learning trajectories, each partition will eventually become a separate trajectory. We run cluster_cells()as before.
cds <- cluster_cells(cds)
plot_file_name <- 'embryo_umap_partition.png'
plot_cmd <- plot_cells(cds, color_cells_by = "partition")
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# Next, we will fit a principal graph within each partition using the learn_graph() function:
cds <- learn_graph(cds)
plot_file_name <- 'embryo_pr_graph_packer_cell_type.png'
plot_cmd <- plot_cells(cds, color_cells_by = "cell.type", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# In time series experiments, this can usually be accomplished by finding spots in the UMAP space that are occupied by cells from early time points:
plot_file_name <- 'embryo_pr_graph_by_time.png'
plot_cmd <- plot_cells(cds, color_cells_by = "embryo.time.bin", label_cell_groups=FALSE, label_leaves=TRUE, label_branch_points=TRUE, graph_label_size=1.5)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# If you don't provide them as an argument, it will launch a graphical user interface for selecting one or more root nodes.
cds <- order_cells(cds)

# Plotting the cells and coloring them by pseudotime shows how they were ordered:
plot_file_name <- 'embryo_pr_graph_by_pseudotime.png'
plot_cmd <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# Then it picks the node that is most heavily occupied by early cells and returns that as the root.
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
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

# Passing the programatically selected root node to order_cells() via the root_pr_nodeargument yields:
plot_file_name <- 'embryo_pr_graph_by_pseudotime_programmatically_ordered.png'
plot_cmd <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

# It is often useful to subset cells based on their branch in the trajectory. The function choose_graph_segments allows you to do so interactively.
cds_sub <- choose_graph_segments(cds)

# Working with 3D trajectories
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))
plot_file_name <- 'emb_3d_by_partition.html'
plot_cmd <- plot_cells_3d(cds_3d, color_cells_by="partition")
htmlwidgets::saveWidget(plot_cmd, file.path(manual_images_dir, plot_file_name), selfcontained = F, libdir = "lib")

