# context("test-subset_cells")
#
# cds <- load_a549()
#
# test_that("test subset_along_path error messages work", {
#   expect_error(cds <- subset_along_path(cds),
#                "No dimensionality reduction for UMAP calculated. Please run reduce_dimension with reduction_method = UMAP and partition_cells before running learn_graph.")
#   cds <- preprocess_cds(cds)
#   expect_error(cds <- subset_along_path(cds),
#                "No dimensionality reduction for UMAP calculated. Please run reduce_dimension with reduction_method = UMAP and partition_cells before running learn_graph.")
#   cds <- reduce_dimension(cds)
#   expect_error(cds <- subset_along_path(cds),
#                "No cell partition for UMAP calculated. Please run partition_cells with reduction_method = UMAP before running learn_graph.")
#   cds <- partition_cells(cds)
#   expect_error(cds <- subset_along_path(cds),
#                "No principal_graph for UMAP calculated. Please run learn_graph with reduction_method = UMAP before running subset_along_path")
#
#   #expect_error(cds <- learn_graph(cds, learn_graph_control = list(FALSE)), "")
#   #expect_error(cds <- subset_along_path(cds, learn_graph_control = list(prune = FALSE)), "Unknown variable in learn_graph_control")
# })
#
# cds <- preprocess_cds(cds)
# cds <- reduce_dimension(cds)
# cds <- reduce_dimension(cds)
# cds <- partition_cells(cds)
# cds <- learn_graph(cds)
#
# # This is a helper function to find the pr graph node that has the highest concentration of vehicle cells
# find_vehicle_pr_node = function(cds){
#   cell_ids <- which(colData(cds)[, "vehicle"])
#
#   closest_vertex <-
#     cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
#   closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
#   root_pr_nodes <-
#     igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
#                                                               (which.max(table(closest_vertex[cell_ids,]))))]
#
#   root_pr_nodes
# }
#
# cds = order_cells(cds, root_pr_nodes = find_vehicle_pr_node(cds))
#
# plot_cell_trajectory(cds)
