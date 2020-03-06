# == Cell Ranger data

#
# Expected values.
#
# pdir <- "../testdata/cr3.0/outs/filtered_feature_bc_matrix"
# sum( assay( cds ) )
cr_assay_sum <- 80147

# as.vector(rownames(assay(cds)))
cr_matrix_rowname <- c( "ENSG00000243485", "ENSG00000237613", "ENSG00000268674", "CD3_GCCTGACTAGATCCA", "CD19_CGTGCAACACTCGTA" )
# fData(cds)$V2
cr_feature_V2 <- c( "RP11-34P13.3", "FAM138A", "FAM231B", "CD3", "CD19" )
# fData(cds)$V3
cr_feature_V3 <- c( "Gene Expression", "Gene Expression", "Gene Expression", "Antibody Capture", "Antibody Capture" )

# as.vector(pData(assay(cds)))
cr_matrix_colname <- c( "AAAGTAGCACAGTCGC-1", "AAATGCCCACCCAGTG-1", "AACCGCGCAGGCGATA-1", "AACTCAGAGAACTCGG-1", "AATCCAGCAGTAACGG-1" )


#
# Test.
#
test_that( "load cell ranger matrix 3.0", {
  pdir <- "../testdata/cr3.0/outs/filtered_feature_bc_matrix/"
  pmat <- paste0( pdir, "matrix.mtx.gz" )
  prow <- paste0( pdir, "features.tsv.gz" )
  pcol <- paste0( pdir, "barcodes.tsv.gz" )
  cds <- load_mm_data( pmat, prow, pcol, header=FALSE )
  expect_equal( sum( assay( cds ) ), cr_assay_sum )
  expect_true( all( as.vector( rownames( assay( cds ) ) ) == cr_matrix_rowname ) )
  expect_true( all( as.vector( colnames( assay( cds ) ) ) == cr_matrix_colname ) )
  expect_true( all( fData( cds )$V2 == cr_feature_V2 ) )
  expect_true( all( fData( cds )$V3 == cr_feature_V3 ) )
} )



# == Matrix Market data derived from Cell Ranger

# Expected values.
#
# pdir <- "../testdata/MatrixMarket"
# as.vector( rownames( assay( cds ) ) )
mm_matrix_rowname <- c( "ENSG00000243485", "ENSG00000237613", "ENSG00000268674", "CD3_GCCTGACTAGATCCA", "CD19_CGTGCAACACTCGTA" )
# as.vector( fData( cds )$gene_short_name )
mm_feature_gene_short_name <- c( "RP11-34P13.3", "FAM138A", "FAM231B", "CD3", "CD19" )
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
test_that( "load MatrixMarket with annotations file ncol=1:nheader=0:header=FALSE", {
  pdir <- "../testdata/MatrixMarket/"
  pmat <- paste0( pdir, "matrix.mtx.gz" )
  prow <- paste0( pdir, "features_c1h0.txt" )
  pcol <- paste0( pdir, "barcodes_c1h0.txt" )
  cds <- load_mm_data( pmat, prow, pcol, header=FALSE, sep="" )
  expect_equal( sum( assay( cds ) ), mm_assay_sum )
  expect_true( all( as.vector( rownames( assay( cds ) ) ) == mm_matrix_rowname ) )
  expect_true( all( as.vector( colnames( assay( cds ) ) ) == mm_matrix_colname ) )
} )


test_that( "load MatrixMarket with annotations file ncol=1:nheader=1:header=TRUE", {
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


test_that( "load MatrixMarket with annotations file ncol=3:nheader=0:header=FALSE", {
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


test_that( "load MatrixMarket with annotations file ncol=3:nheader=2:header=TRUE", {
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


test_that( "load MatrixMarket with annotations file ncol=3:nheader=3:header=TRUE", {
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


test_that( "load MatrixMarket with annotations file ncol=3:nheader=3:header=TRUE and replace metadata dimension labels", {
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

