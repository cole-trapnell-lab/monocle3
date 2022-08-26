context("test-cluster_genes")

test_that('test cluster_genes', {
  cds <- monocle3:::load_worm_embryo()
  set.seed(42)
  #
  afb_cells <- rownames(colData(cds)[!is.na(colData(cds)$cell.type)&colData(cds)$cell.type=='AFD',])
  afb_cds <- cds[,afb_cells]

  afb_cds <- preprocess_cds(afb_cds, num_dim = 100)
  afb_cds <- align_cds(afb_cds, alignment_group='batch', residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")
  afb_cds <- reduce_dimension(afb_cds, umap.fast_sgd=FALSE, cores=1)
  afb_cds <- cluster_cells(afb_cds)
  afb_cds <- learn_graph(afb_cds, learn_graph_control=list(ncenter=1000), close_loop=TRUE)
  #
  pr_graph_test_res <- graph_test(afb_cds)
  pr_deg_ids <- subset(pr_graph_test_res, q_value < 0.05)$id
  gene_module_df <- find_gene_modules(afb_cds[pr_deg_ids,], resolution=1e-2)
  cell_group_df <- tibble::tibble(cell=row.names(colData(afb_cds)),
                                  cell_group=partitions(afb_cds)[colnames(afb_cds)])
  agg_mat <- aggregate_gene_expression(afb_cds, gene_module_df, cell_group_df)

  # gene_cluster_df = monocle3:::cluster_genes(test_cds[pr_deg_ids,], resolution=0.001)
  #
  #text_df <- gene_cluster_df %>% dplyr::group_by(cluster) %>% dplyr::summarize(text_x = stats::median(x = dim_1),
  #                                                                             text_y = stats::median(x = dim_2))
  # ggplot2::qplot(dim_1, dim_2, color=cluster, data=gene_cluster_df) +
  #   ggplot2::geom_text(data=text_df, mapping = ggplot2::aes_string(x = "text_x", y = "text_y", label = "cluster"), color=I("black"),  size = 4)
  # png("module_graph.png", res=600, width=8, height=8, units="in")
  # monocle3:::plot_cells(test_cds, genes=gene_cluster_df, cell_size=0.5, show_backbone=FALSE, label_branch_points=FALSE, label_leaves=FALSE, label_roots=FALSE)
  # dev.off()
} )
