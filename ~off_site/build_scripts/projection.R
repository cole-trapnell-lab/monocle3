# This script creates the plots used by the projection.html page
# and stores the new files in the directory given by the bash
# variable 'manual_images_dir', as defined in globals.R.
#
# Call this script from R using the command
#   source('projection.R', local=TRUE, echo=TRUE)
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

library(Matrix)


# We begin by loading the reference and query data sets into Monocle3.
matrix_ref <- readMM(gzcon(url("https://staff.washington.edu/bge/monocle3/projection/data/cao.mouse_embryo.sample.mtx.gz")))
cell_ann_ref <- read.csv(gzcon(url("https://staff.washington.edu/bge/monocle3/projection/data/cao.mouse_embryo.sample.coldata.txt.gz"), text=TRUE), sep='\t')
gene_ann_ref <- read.csv(gzcon(url("https://staff.washington.edu/bge/monocle3/projection/data/cao.mouse_embryo.sample.rowdata.txt.gz"), text=TRUE), sep='\t')

cds_ref <- new_cell_data_set(matrix_ref,
                             cell_metadata = cell_ann_ref,
                             gene_metadata = gene_ann_ref)

matrix_qry <- readMM(gzcon(url("https://staff.washington.edu/bge/monocle3/projection/data/srivatsan.mouse_embryo_scispace.sample.mtx.gz")))
cell_ann_qry <- read.csv(gzcon(url("https://staff.washington.edu/bge/monocle3/projection/data/srivatsan.mouse_embryo_scispace.sample.coldata.txt.gz"), text=TRUE), sep='\t')
gene_ann_qry <- read.csv(gzcon(url("https://staff.washington.edu/bge/monocle3/projection/data/srivatsan.mouse_embryo_scispace.sample.rowdata.txt.gz"), text=TRUE), sep='\t')

cds_qry <- new_cell_data_set(matrix_qry,
                             cell_metadata = cell_ann_qry,
                             gene_metadata = gene_ann_qry)


# It is essential that the query cds has the same genes in the same order as the reference cds, so we identify the shared genes and sort them.
genes_ref <- row.names(cds_ref)
genes_qry <- row.names(cds_qry)
genes_shared <- intersect(genes_ref, genes_qry)
cds_ref <- cds_ref[genes_shared,]
cds_qry <- cds_qry[genes_shared,]


# After applying the gene and UMI filters, we re-calculate the size factors of both data sets.
cds_ref <- estimate_size_factors(cds_ref)
cds_qry <- estimate_size_factors(cds_qry)


# After processing, we save the transform models and nearest neighbor index so that we can use them to transform the query data set into the reference space.
cds_ref <- preprocess_cds(cds_ref, num_dim=100)
cds_ref <- reduce_dimension(cds_ref, build_nn_index=TRUE)
save_transform_models(cds_ref, 'cds_ref_test_models')


# The load_transform_models() function loads the reference transform models into the query cds where they are used by preprocess_transform() and reduce_dimension_transform().
cds_qry <- load_transform_models(cds_qry, 'cds_ref_test_models')
cds_qry <- preprocess_transform(cds_qry)
cds_qry <- reduce_dimension_transform(cds_qry)

# First we plot the reference and query cells in UMAP space.
plot_file_name <- 'L2_projection_reference.png'
plot_cmd <- plot_cells(cds_ref)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)

plot_file_name <- 'L2_projection_query.png'
plot_cmd <- plot_cells(cds_qry)
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)


# Now we label the cells in the reference and query cdses, combine the cdses, and plot the combined cells.
colData(cds_ref)[['data_set']] <- 'reference'
colData(cds_qry)[['data_set']] <- 'query'
cds_combined <- combine_cds(list(cds_ref, cds_qry),  keep_all_genes=TRUE, cell_names_unique=TRUE, keep_reduced_dims=TRUE)
plot_file_name <- 'L2_projection_combined.png'
plot_cmd <- plot_cells(cds_combined, color_cells_by='data_set')
ggplot_cells_png(plot_cmd, manual_images_dir=manual_images_dir, plot_file_name=plot_file_name)


