# == Cell Ranger 3.0 data

context("test-io")
skip_not_travis <- function ()
{
  if (identical(Sys.getenv("TRAVIS"), "true")) {
    return(invisible(TRUE))
  }
  skip("Not on Travis")
}


transform_models <- 'transform_models.01'
monocle_objects <- 'monocle_objects.01'

#
# Expected values.
#
# pdir <- "../testdata/cr3.0/outs/filtered_feature_bc_matrix"
# sum( assay( cds ) )
cr3_assay_sum <- 80147

# as.vector(rownames(assay(cds)))
cr3_matrix_rowname <- c( "ENSG00000243485", "ENSG00000237613", "ENSG00000268674", "CD3_GCCTGACTAGATCCA", "CD19_CGTGCAACACTCGTA" )
# fData(cds)$V2
cr3_feature_V2 <- c( "RP11-34P13.3", "FAM138A", "FAM231B", "CD3", "CD19" )
# fData(cds)$V3
cr3_feature_V3 <- c( "Gene Expression", "Gene Expression", "Gene Expression", "Antibody Capture", "Antibody Capture" )

# as.vector(pData(assay(cds)))
cr3_matrix_colname <- c( "AAAGTAGCACAGTCGC-1", "AAATGCCCACCCAGTG-1", "AACCGCGCAGGCGATA-1", "AACTCAGAGAACTCGG-1", "AATCCAGCAGTAACGG-1" )


#
# Test.
#
test_that("load_mm_data: load cell ranger matrix 3.0", {
  pdir <- "../testdata/cr3.0/outs/filtered_feature_bc_matrix/"
  pmat <- paste0( pdir, "matrix.mtx.gz" )
  prow <- paste0( pdir, "features.tsv.gz" )
  pcol <- paste0( pdir, "barcodes.tsv.gz" )
  cds <- load_mm_data( pmat, prow, pcol, header=FALSE )
  expect_equal( sum( assay( cds ) ), cr3_assay_sum )
  expect_true( all( as.vector( rownames( assay( cds ) ) ) == cr3_matrix_rowname ) )
  expect_true( all( as.vector( colnames( assay( cds ) ) ) == cr3_matrix_colname ) )
  expect_true( all( fData( cds )$V2 == cr3_feature_V2 ) )
  expect_true( all( fData( cds )$V3 == cr3_feature_V3 ) )
} )


# == Cell Ranger 2.0 data

#
# Expected values.
#
# pdir <- "../testdata/cr2.0/outs/filtered_gene_bc_matrices/hg19/"
# sum( assay( cds ) )
cr2_assay_sum <- 22265

# as.vector(rownames(assay(cds)))[1:5]
cr2_matrix_rowname <- c( "ENSGXXXX00", "ENSGXXXX01", "ENSGXXXX02", "ENSGXXXX03", "ENSGXXXX04" )
# fData(cds)$V2[1:5]
cr2_feature_V2 <- c( "MS4A1", "CD79B", "CD79A", "HLA-DRA", "TCL1A" )

# as.vector(colnames(assay(cds)))[1:5]
cr2_matrix_colname <- c( "GAACCTGATGAACC-1", "TGACTGGATTCTCA-1", "AGTCAGACTGCACA-1", "AATGTTGACAGTCA-1", "AGAGATGATCTCGC-1" )

#
# Test.
#
test_that("load_mm_data: load cell ranger matrix 2.0", {
  pdir <- "../testdata/cr2.0/outs/filtered_gene_bc_matrices/hg19/"
  pmat <- paste0( pdir, "matrix.mtx" )
  prow <- paste0( pdir, "genes.tsv" )
  pcol <- paste0( pdir, "barcodes.tsv" )
  cds <- load_mm_data( pmat, prow, pcol, header=FALSE )
  expect_equal( sum( assay( cds ) ), cr2_assay_sum )
  expect_true( all( as.vector( rownames( assay( cds ) ) )[1:5] == cr2_matrix_rowname ) )
  expect_true( all( as.vector( colnames( assay( cds ) ) )[1:5] == cr2_matrix_colname ) )
  expect_true( all( fData( cds )$V2[1:5] == cr2_feature_V2 ) )
} )


# == Matrix Market data derived from Cell Ranger

# Expected values.
#
# pdir <- "../testdata/MatrixMarket"
# as.vector( rownames( assay( cds ) ) )
mm_matrix_rowname <- c( "ENSG00000243485", "ENSG00000237613", "ENSG00000268674", "CD3_GCCTGACTAGATCCA", "CD19_CGTGCAACACTCGTA" )
# as.vector( fData( cds )$gene_short_name )
mm_feature_gene_short_name <- c( "RP11-34P13.3", "FAM138A", "FAM231B", "CD3", "CD19" )
mm_feature_gene_short_name_wquote <- c( "RP1\'1-34P13.3", "FAM\'138A", "FAM\'231B", "CD3", "CD19" )
# as.vector( fData( cds )$source )
mm_feature_source <- c( "Gene_Expression", "Gene_Expression", "Gene_Expression", "Antibody_Capture", "Antibody_Capture" )

# as.vector( colnames( assay( cds ) ) )
mm_matrix_colname <-c( "AAAGTAGCACAGTCGC-1", "AAATGCCCACCCAGTG-1", "AACCGCGCAGGCGATA-1", "AACTCAGAGAACTCGG-1", "AATCCAGCAGTAACGG-1" )
# as.vector( pData( cds )$cell_number )
mm_cell_cell_number <- c( "cell1", "cell2",  "cell3",  "cell4",  "cell5" )
# as.vector( pData( cds )$umi_token )
mm_cell_umi_token <- c( "umi=10",  "umi=11",  "umi=12",  "umi=13",  "umi=14" )

# sum( assay( cds ) )
mm_assay_sum <- 80147


#
# Tests.
#
test_that("load_mm_data: load MatrixMarket with annotations file ncol=1:nheader=0:header=FALSE", {
  pdir <- "../testdata/MatrixMarket/"
  pmat <- paste0( pdir, "matrix.mtx.gz" )
  prow <- paste0( pdir, "features_c1h0.txt" )
  pcol <- paste0( pdir, "barcodes_c1h0.txt" )
  cds <- load_mm_data( pmat, prow, pcol, header=FALSE, sep="" )
  expect_equal( sum( assay( cds ) ), mm_assay_sum )
  expect_true( all( as.vector( rownames( assay( cds ) ) ) == mm_matrix_rowname ) )
  expect_true( all( as.vector( colnames( assay( cds ) ) ) == mm_matrix_colname ) )
} )


test_that("load_mm_data: load MatrixMarket with annotations file ncol=1:nheader=1:header=TRUE", {
  pdir <- "../testdata/MatrixMarket/"
  pmat <- paste0( pdir, "matrix.mtx.gz" )
  prow <- paste0( pdir, "features_c1h1.txt" )
  pcol <- paste0( pdir, "barcodes_c1h1.txt" )
  cds <- load_mm_data( pmat, prow, pcol, header=TRUE, sep="" )
  expect_equal( sum( assay( cds ) ), mm_assay_sum )
  expect_true( all( as.vector( rownames( assay( cds ) ) ) == mm_matrix_rowname ) )
  expect_true( all( as.vector( colnames( assay( cds ) ) ) == mm_matrix_colname ) )
  expect_true( all( as.vector( pData( cds )$cell_number ) == mm_cell_cell_number ) )
  expect_true( all( as.vector( pData( cds )$umi_token ) == mm_cell_umi_token ) )
} )


test_that("load_mm_data: load MatrixMarket with annotations file ncol=3:nheader=0:header=FALSE", {
  pdir <- "../testdata/MatrixMarket/"
  pmat <- paste0( pdir, "matrix.mtx.gz" )
  prow <- paste0( pdir, "features_c3h0.txt" )
  pcol <- paste0( pdir, "barcodes_c3h0.txt" )
  cds <- load_mm_data( pmat, prow, pcol, header=FALSE, sep="" )
  expect_equal( sum( assay( cds ) ), mm_assay_sum )
  expect_true( all( as.vector( rownames( assay( cds ) ) ) == mm_matrix_rowname ) )
  expect_true( all( as.vector( colnames( assay( cds ) ) ) == mm_matrix_colname ) )
  expect_true( all( fData( cds )$V2 == mm_feature_gene_short_name ) )
  expect_true( all( fData( cds )$V3 == mm_feature_source ) )
  expect_true( all( as.vector( pData( cds )$cell_number ) == mm_cell_cell_number ) )
  expect_true( all( as.vector( pData( cds )$umi_token ) == mm_cell_umi_token ) )
} )


test_that("load_mm_data: load MatrixMarket with annotations file ncol=3:nheader=2:header=TRUE", {
  pdir <- "../testdata/MatrixMarket/"
  pmat <- paste0( pdir, "matrix.mtx.gz" )
  prow <- paste0( pdir, "features_c3h2.txt" )
  pcol <- paste0( pdir, "barcodes_c3h2.txt" )
  cds <- load_mm_data( pmat, prow, pcol, header=TRUE, sep="" )
  expect_equal( sum( assay( cds ) ), mm_assay_sum )
  expect_true( all( as.vector( rownames( assay( cds ) ) ) == mm_matrix_rowname ) )
  expect_true( all( as.vector( colnames( assay( cds ) ) ) == mm_matrix_colname ) )
  expect_true( all( fData( cds )$gene_short_name == mm_feature_gene_short_name ) )
  expect_true( all( fData( cds )$source == mm_feature_source ) )
  expect_true( all( as.vector( pData( cds )$cell_number ) == mm_cell_cell_number ) )
  expect_true( all( as.vector( pData( cds )$umi_token ) == mm_cell_umi_token ) )
} )


test_that("load_mm_data: load MatrixMarket with annotations file ncol=3:nheader=3:header=TRUE", {
  pdir <- "../testdata/MatrixMarket/"
  pmat <- paste0( pdir, "matrix.mtx.gz" )
  prow <- paste0( pdir, "features_c3h3.txt" )
  pcol <- paste0( pdir, "barcodes_c3h3.txt" )
  cds <- load_mm_data( pmat, prow, pcol, header=TRUE, sep="" )
  expect_equal( sum( assay( cds ) ), mm_assay_sum )
  expect_true( all( as.vector( rownames( assay( cds ) ) ) == mm_matrix_rowname ) )
  expect_true( all( as.vector( colnames( assay( cds ) ) ) == mm_matrix_colname ) )
  expect_true( all( fData( cds )$gene_short_name == mm_feature_gene_short_name ) )
  expect_true( all( fData( cds )$source == mm_feature_source ) )
  expect_true( all( as.vector( pData( cds )$cell_number ) == mm_cell_cell_number ) )
  expect_true( all( as.vector( pData( cds )$umi_token ) == mm_cell_umi_token ) )
} )


test_that("load_mm_data: load MatrixMarket with annotations file ncol=3:nheader=3:header=TRUE and replace metadata dimension labels", {
  pdir <- "../testdata/MatrixMarket/"
  pmat <- paste0( pdir, "matrix.mtx.gz" )
  prow <- paste0( pdir, "features_c3h3.txt" )
  pcol <- paste0( pdir, "barcodes_c3h3.txt" )
  cds <- load_mm_data( pmat, prow, pcol,
                       feature_metadata_column_names=c('features1','features2'),
                       cell_metadata_column_names=c('cells1','cells2'),
                       header=TRUE, sep="" )
  expect_equal( sum( assay( cds ) ), mm_assay_sum )
  expect_true( all( as.vector( rownames( assay( cds ) ) ) == mm_matrix_rowname ) )
  expect_true( all( as.vector( colnames( assay( cds ) ) ) == mm_matrix_colname ) )
  expect_true( all( fData( cds )$features1 == mm_feature_gene_short_name ) )
  expect_true( all( fData( cds )$features2 == mm_feature_source ) )
  expect_true( all( as.vector( pData( cds )$cells1 ) == mm_cell_cell_number ) )
  expect_true( all( as.vector( pData( cds )$cells2 ) == mm_cell_umi_token ) )
} )


test_that("load_mm_data: load MatrixMarket with annotations file ncol=2:nheader=1:header=TRUE", {
  pdir <- "../testdata/MatrixMarket/"
  pmat <- paste0( pdir, "matrix.mtx.gz" )
  prow <- paste0( pdir, "features_c2h1.txt" )
  pcol <- paste0( pdir, "barcodes_c3h3.txt" )
  cds <- load_mm_data( pmat, prow, pcol,
                       header=TRUE, sep="" )
  expect_equal( sum( assay( cds ) ), mm_assay_sum )
  expect_true( all( as.vector( rownames( assay( cds ) ) ) == mm_matrix_rowname ) )
  expect_true( all( as.vector( colnames( assay( cds ) ) ) == mm_matrix_colname ) )
  expect_true( all( fData( cds )$gene_short_name == mm_feature_gene_short_name ) )
  expect_true( all( as.vector( pData( cds )$cells1 ) == mm_cell_cell_number ) )
  expect_true( all( as.vector( pData( cds )$cells2 ) == mm_cell_umi_token ) )
} )


test_that("load_mm_data: load MatrixMarket with annotations file ncol=2:nheader=1:header=FALSE and there is a single quote in the gene names", {
  pdir <- "../testdata/MatrixMarket/"
  pmat <- paste0( pdir, "matrix.mtx.gz" )
  prow <- paste0( pdir, "features_c3h2q.txt" )
  pcol <- paste0( pdir, "barcodes_c3h3.txt" )
  cds <- load_mm_data( pmat, prow, pcol,
                       quote="\"",
                       header=TRUE, sep="" )
  expect_equal( sum( assay( cds ) ), mm_assay_sum )
  expect_true( all( as.vector( rownames( assay( cds ) ) ) == mm_matrix_rowname ) )
  expect_true( all( as.vector( colnames( assay( cds ) ) ) == mm_matrix_colname ) )
  expect_true( all( fData( cds )$gene_short_name == mm_feature_gene_short_name_wquote ) )
  expect_true( all( as.vector( pData( cds )$cells1 ) == mm_cell_cell_number ) )
  expect_true( all( as.vector( pData( cds )$cells2 ) == mm_cell_umi_token ) )
} )


test_that("save_transform_models and load_transform_models", {
  skip_on_travis()
  # set up a cds with nearest neighbor indices and transform models
  cds <- load_a549()
  cds <- preprocess_cds(cds, build_nn_index=TRUE, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  cds <- align_cds(cds, preprocess_method='PCA', residual_model_formula_str='~n.umi', build_nn_index=TRUE, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  cds <- reduce_dimension(cds, preprocess_method='Aligned', build_nn_index=TRUE, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  save_transform_models(cds, directory_path=transform_models)

  # load cds, load transform models, and process cds with transform models
  rm(cds)
  cds <- load_a549()
  cds <- load_transform_models(cds, directory_path=transform_models)
  cds <- preprocess_transform(cds)
  cds <- align_beta_transform(cds, preprocess_method='PCA')
  cds <- reduce_dimension_transform(cds, preprocess_method='Aligned', reduction_method='UMAP')

  nn_index <- get_cds_nn_index(cds, reduction_method='PCA', nn_control=list(method='annoy'))
  expect_equal(as.character(class(nn_index[['ann']])), 'Rcpp_AnnoyEuclidean')

  # check PCA reduced dims matrix and nearest neighbors
  expect_equivalent(ncol(reducedDims(cds)[['PCA']]), 50)
  expect_equivalent(nrow(reducedDims(cds)[['PCA']]), 500)
  expect_equivalent(reducedDims(cds)[['PCA']][[1,1]], 2.420739, tol=1e-5)
  nn_res <- search_nn_index(query_matrix=reducedDims(cds)[['PCA']], nn_index=get_cds_nn_index(cds, reduction_method='PCA', nn_control=list(method='annoy', metric='euclidean', n_trees=50)), k=5, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  expect_equivalent(nn_res[['nn.idx']][[1]], 1)
  expect_equivalent(nn_res[['nn.dists']][[1]], 0)

  # check Aligned reduced dims matrix and nearest neighbors
  expect_equivalent(ncol(reducedDims(cds)[['Aligned']]), 50)
  expect_equivalent(nrow(reducedDims(cds)[['Aligned']]), 500)
  expect_equivalent(reducedDims(cds)[['Aligned']][[1,1]], 3.870306, tol=1e-5)
  nn_res <- search_nn_index(query_matrix=reducedDims(cds)[['Aligned']], nn_index=get_cds_nn_index(cds, reduction_method='Aligned', nn_control=list(method='annoy', metric='euclidean', n_trees=50)), k=5, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  expect_equivalent(nn_res[['nn.idx']][[1]], 1)
  expect_equivalent(nn_res[['nn.dists']][[1]], 0)

  # check UMAP reduced dims matrix and nearest neighbors
  expect_equivalent(ncol(reducedDims(cds)[['UMAP']]), 2)
  expect_equivalent(nrow(reducedDims(cds)[['UMAP']]), 500)
  expect_equivalent(reducedDims(cds)[['UMAP']][[1,1]], 1.96, tol=1e-2)
  nn_res <- search_nn_index(query_matrix=reducedDims(cds)[['UMAP']], nn_index=get_cds_nn_index(cds, reduction_method='UMAP', nn_control=list(method='annoy', metric='euclidean', n_trees=50)), k=5, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  expect_equivalent(nn_res[['nn.idx']][[1]], 1)
  expect_equivalent(nn_res[['nn.dists']][[1]], 0)
  system(paste0('rm -r ', transform_models))

  rm(cds)
  rm(nn_index)
  cds <- load_a549()
  cds <- preprocess_cds(cds, build_nn_index=TRUE, nn_control=list(method='hnsw', metric='euclidean', n_trees=50))
  save_transform_models(cds, directory_path=transform_models)
  rm(cds)
  cds <- load_a549()
  cds <- preprocess_cds(cds) 
  cds <- load_transform_models(cds, directory_path=transform_models)
  nn_index <- get_cds_nn_index(cds, reduction_method='PCA', nn_control=list(method='hnsw'))
  expect_equal(as.character(class(nn_index)), 'Rcpp_HnswL2')
  system(paste0('rm -r ', transform_models))
} )


test_that("save_monocle_objects and load_monocle_objects", {
  skip_on_travis()
  # set up a cds with nearest neighbor indices and monocle_objects
  cds <- load_a549()
  cds <- preprocess_cds(cds, build_nn_index=TRUE, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  cds <- align_cds(cds, preprocess_method='PCA', residual_model_formula_str='~n.umi', build_nn_index=TRUE, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  cds <- reduce_dimension(cds, preprocess_method='Aligned', build_nn_index=TRUE, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  save_monocle_objects(cds, directory_path=monocle_objects)

  # load monocle objects models
  rm(cds)
  cds <- load_monocle_objects(directory_path=monocle_objects)

  nn_index <- get_cds_nn_index(cds, reduction_method='PCA', nn_control=list(method='annoy'))
  expect_equal(as.character(class(nn_index[['ann']])), 'Rcpp_AnnoyEuclidean')

  # check count matrix
  expect_equivalent(ncol(counts(cds)), 500)
  expect_equivalent(nrow(counts(cds)), 198)
  expect_equivalent(Matrix::colSums(counts(cds))[[1]], 19)

  # check PCA reduced dims matrix and nearest neighbors
  expect_equivalent(ncol(reducedDims(cds)[['PCA']]), 50)
  expect_equivalent(nrow(reducedDims(cds)[['PCA']]), 500)
  expect_equivalent(reducedDims(cds)[['PCA']][[1,1]], 2.420739, tol=1e-5)
  nn_res <- search_nn_index(query_matrix=reducedDims(cds)[['PCA']], nn_index=get_cds_nn_index(cds, reduction_method='PCA', nn_control=list(method='annoy', metric='euclidean', n_trees=50)), k=5, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  expect_equivalent(nn_res[['nn.idx']][[1]], 1)
  expect_equivalent(nn_res[['nn.dists']][[1]], 0)

  # check Aligned reduced dims matrix and nearest neighbors
  expect_equivalent(ncol(reducedDims(cds)[['Aligned']]), 50)
  expect_equivalent(nrow(reducedDims(cds)[['Aligned']]), 500)
  expect_equivalent(reducedDims(cds)[['Aligned']][[1,1]], 3.870306, tol=1e-5)
  nn_res <- search_nn_index(query_matrix=reducedDims(cds)[['Aligned']], nn_index=get_cds_nn_index(cds, reduction_method='Aligned', nn_control=list(method='annoy', metric='euclidean', n_trees=50)), k=5, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  expect_equivalent(nn_res[['nn.idx']][[1]], 1)
  expect_equivalent(nn_res[['nn.dists']][[1]], 0)

  # check UMAP reduced dims matrix and nearest neighbors
  expect_equivalent(ncol(reducedDims(cds)[['UMAP']]), 2)
  expect_equivalent(nrow(reducedDims(cds)[['UMAP']]), 500)
  expect_equivalent(reducedDims(cds)[['UMAP']][[1,1]], 1.96, tol=1e-2)
  nn_res <- search_nn_index(query_matrix=reducedDims(cds)[['UMAP']], nn_index=get_cds_nn_index(cds, reduction_method='UMAP', nn_control=list(method='annoy', metric='euclidean', n_trees=50)), k=5, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  expect_equivalent(nn_res[['nn.idx']][[1]], 1)
  expect_equivalent(nn_res[['nn.dists']][[1]], 0)
  system(paste0('rm -r ', monocle_objects))

  rm(cds)
  rm(nn_index)
  cds <- load_a549()
  cds <- preprocess_cds(cds, build_nn_index=TRUE, nn_control=list(method='hnsw', metric='euclidean', n_trees=50))
  save_monocle_objects(cds, directory_path=monocle_objects)
  rm(cds)
  cds <- load_a549()
  cds <- load_monocle_objects(directory_path=monocle_objects)
  nn_index <- get_cds_nn_index(cds, reduction_method='PCA', nn_control=list(method='hnsw'))
  expect_equal(as.character(class(nn_index)), 'Rcpp_HnswL2')
  system(paste0('rm -r ', monocle_objects))
} )


test_that("save_transform_models and load_transform_models", {
  skip_not_travis()
  # set up a cds with nearest neighbor indices and transform models
  cds <- load_a549()
  cds <- preprocess_cds(cds, build_nn_index=TRUE, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  cds <- align_cds(cds, preprocess_method='PCA', residual_model_formula_str='~n.umi', build_nn_index=TRUE, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  cds <- reduce_dimension(cds, preprocess_method='Aligned', build_nn_index=TRUE, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  save_transform_models(cds, directory_path=transform_models)

  # load cds, load transform models, and process cds with transform models
  rm(cds)
  cds <- load_a549()
  cds <- load_transform_models(cds, directory_path=transform_models)
  cds <- preprocess_transform(cds)
  cds <- align_beta_transform(cds, preprocess_method='PCA')
  cds <- reduce_dimension_transform(cds, preprocess_method='Aligned', reduction_method='UMAP')

  nn_index <- get_cds_nn_index(cds, reduction_method='PCA', nn_control=list(method='annoy'))
  expect_equal(as.character(class(nn_index[['ann']])), 'Rcpp_AnnoyEuclidean')

  # check PCA reduced dims matrix and nearest neighbors
  expect_equivalent(ncol(reducedDims(cds)[['PCA']]), 50)
  expect_equivalent(nrow(reducedDims(cds)[['PCA']]), 500)
  expect_equivalent(reducedDims(cds)[['PCA']][[1,1]], 2.420739, tol=1e-5)
  nn_res <- search_nn_index(query_matrix=reducedDims(cds)[['PCA']], nn_index=get_cds_nn_index(cds, reduction_method='PCA', nn_control=list(method='annoy', metric='euclidean', n_trees=50)), k=5, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  expect_equivalent(nn_res[['nn.idx']][[1]], 1)
  expect_equivalent(nn_res[['nn.dists']][[1]], 0)

  # check Aligned reduced dims matrix and nearest neighbors
  expect_equivalent(ncol(reducedDims(cds)[['Aligned']]), 50)
  expect_equivalent(nrow(reducedDims(cds)[['Aligned']]), 500)
  expect_equivalent(reducedDims(cds)[['Aligned']][[1,1]], 3.870306, tol=1e-5)
  nn_res <- search_nn_index(query_matrix=reducedDims(cds)[['Aligned']], nn_index=get_cds_nn_index(cds, reduction_method='Aligned', nn_control=list(method='annoy', metric='euclidean', n_trees=50)), k=5, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  expect_equivalent(nn_res[['nn.idx']][[1]], 1)
  expect_equivalent(nn_res[['nn.dists']][[1]], 0)

  # check UMAP reduced dims matrix and nearest neighbors
  expect_equivalent(ncol(reducedDims(cds)[['UMAP']]), 2)
  expect_equivalent(nrow(reducedDims(cds)[['UMAP']]), 500)
  expect_equivalent(reducedDims(cds)[['UMAP']][[1,1]], 1.96, tol=1e-2)
  nn_res <- search_nn_index(query_matrix=reducedDims(cds)[['UMAP']], nn_index=get_cds_nn_index(cds, reduction_method='UMAP', nn_control=list(method='annoy', metric='euclidean', n_trees=50)), k=5, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  expect_equivalent(nn_res[['nn.idx']][[1]], 1)
  expect_equivalent(nn_res[['nn.dists']][[1]], 0)
  system(paste0('rm -r ', transform_models))

  rm(cds)
  rm(nn_index)
  cds <- load_a549()
  cds <- preprocess_cds(cds, build_nn_index=TRUE, nn_control=list(method='hnsw', metric='euclidean', n_trees=50))
  save_transform_models(cds, directory_path=transform_models)
  rm(cds)
  cds <- load_a549()
  cds <- preprocess_cds(cds)
  cds <- load_transform_models(cds, directory_path=transform_models)
  nn_index <- get_cds_nn_index(cds, reduction_method='PCA', nn_control=list(method='hnsw'))
  expect_equal(as.character(class(nn_index)), 'Rcpp_HnswL2')
  system(paste0('rm -r ', transform_models))
} )


test_that("save_monocle_objects and load_monocle_objects", {
  skip_not_travis()
  # set up a cds with nearest neighbor indices and monocle_objects
  cds <- load_a549()
  cds <- preprocess_cds(cds, build_nn_index=TRUE, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  cds <- align_cds(cds, preprocess_method='PCA', residual_model_formula_str='~n.umi', build_nn_index=TRUE, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  cds <- reduce_dimension(cds, preprocess_method='Aligned', build_nn_index=TRUE, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  save_monocle_objects(cds, directory_path=monocle_objects)

  # load monocle objects models
  rm(cds)
  cds <- load_monocle_objects(directory_path=monocle_objects)

  nn_index <- get_cds_nn_index(cds, reduction_method='PCA', nn_control=list(method='annoy'))
  expect_equal(as.character(class(nn_index[['ann']])), 'Rcpp_AnnoyEuclidean')

  # check count matrix
  expect_equivalent(ncol(counts(cds)), 500)
  expect_equivalent(nrow(counts(cds)), 198)
  expect_equivalent(Matrix::colSums(counts(cds))[[1]], 19)

  # check PCA reduced dims matrix and nearest neighbors
  expect_equivalent(ncol(reducedDims(cds)[['PCA']]), 50)
  expect_equivalent(nrow(reducedDims(cds)[['PCA']]), 500)
  expect_equivalent(reducedDims(cds)[['PCA']][[1,1]], 2.420739, tol=1e-5)
  nn_res <- search_nn_index(query_matrix=reducedDims(cds)[['PCA']], nn_index=get_cds_nn_index(cds, reduction_method='PCA', nn_control=list(method='annoy', metric='euclidean', n_trees=50)), k=5, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  expect_equivalent(nn_res[['nn.idx']][[1]], 1)
  expect_equivalent(nn_res[['nn.dists']][[1]], 0)

  # check Aligned reduced dims matrix and nearest neighbors
  expect_equivalent(ncol(reducedDims(cds)[['Aligned']]), 50)
  expect_equivalent(nrow(reducedDims(cds)[['Aligned']]), 500)
  expect_equivalent(reducedDims(cds)[['Aligned']][[1,1]], 3.870306, tol=1e-5)
  nn_res <- search_nn_index(query_matrix=reducedDims(cds)[['Aligned']], nn_index=get_cds_nn_index(cds, reduction_method='Aligned', nn_control=list(method='annoy', metric='euclidean', n_trees=50)), k=5, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  expect_equivalent(nn_res[['nn.idx']][[1]], 1)
  expect_equivalent(nn_res[['nn.dists']][[1]], 0)

  # check UMAP reduced dims matrix and nearest neighbors
  expect_equivalent(ncol(reducedDims(cds)[['UMAP']]), 2)
  expect_equivalent(nrow(reducedDims(cds)[['UMAP']]), 500)
  expect_equivalent(reducedDims(cds)[['UMAP']][[1,1]], 1.96, tol=1e-2)
  nn_res <- search_nn_index(query_matrix=reducedDims(cds)[['UMAP']], nn_index=get_cds_nn_index(cds, reduction_method='UMAP', nn_control=list(method='annoy', metric='euclidean', n_trees=50)), k=5, nn_control=list(method='annoy', metric='euclidean', n_trees=50))
  expect_equivalent(nn_res[['nn.idx']][[1]], 1)
  expect_equivalent(nn_res[['nn.dists']][[1]], 0)
  system(paste0('rm -r ', monocle_objects))

  rm(cds)
  rm(nn_index)
  cds <- load_a549()
  cds <- preprocess_cds(cds, build_nn_index=TRUE, nn_control=list(method='hnsw', metric='euclidean', n_trees=50))
  save_monocle_objects(cds, directory_path=monocle_objects)
  rm(cds)
  cds <- load_a549()
  cds <- load_monocle_objects(directory_path=monocle_objects)
  nn_index <- get_cds_nn_index(cds, reduction_method='PCA', nn_control=list(method='hnsw'))
  expect_equal(as.character(class(nn_index)), 'Rcpp_HnswL2')
  system(paste0('rm -r ', monocle_objects))
} )


