context("test-identity")

skip_not_travis <- function ()
{
  if (identical(Sys.getenv("TRAVIS"), "true")) {
    return(invisible(TRUE))
  }
  skip("Not on Travis")
}

cds <- load_a549()
mtxfile <- 'expr.mtx'
pdtfile <- 'pdat.txt'
fdtfile <- 'fdat.txt'

transform_models <- 'transform_models.02'

Matrix::writeMM(counts(cds), file=mtxfile)
write.table(pData(cds), file=pdtfile)
write.table(fData(cds), file=fdtfile)

test_that("identity strings", {
  cds <- load_mm_data(mtxfile,fdtfile, pdtfile,header=TRUE)
  cds <- preprocess_cds(cds, num_dim=25)
  cds <- align_cds(cds, preprocess_method='PCA', residual_model_formula_str='~n.umi')
  cds <- reduce_dimension(cds, preprocess_method='Aligned')

  expect_equal(cds@int_metadata[['counts_metadata']][['identity']][['matrix_type']], mtxfile)

  expect_equal(cds@int_metadata[['reduce_dim_metadata']][['PCA']][['identity']][['matrix_type']], 'matrix:PCA')
  expect_equal(cds@int_metadata[['reduce_dim_metadata']][['PCA']][['identity']][['matrix_id']], cds@int_metadata[['reduce_dim_metadata']][['PCA']][['identity']][['model_id']])
  expect_equal(cds@int_metadata[['reduce_dim_metadata']][['PCA']][['identity']][['prev_matrix_type']], 'expr.mtx')
  expect_equal(cds@int_metadata[['reduce_dim_metadata']][['PCA']][['identity']][['prev_matrix_id']], cds@int_metadata[['counts_metadata']][['identity']][['matrix_id']])
  expect_equal(cds@int_metadata[['reduce_dim_metadata']][['PCA']][['identity']][['model_type']], 'matrix:PCA')
  expect_equal(cds@int_metadata[['reduce_dim_metadata']][['PCA']][['identity']][['model_id']], cds@reduce_dim_aux[['PCA']][['model']][['identity']][['model_id']])

  expect_equal(cds@int_metadata[['reduce_dim_metadata']][['Aligned']][['identity']][['matrix_type']], 'matrix:Aligned')
  expect_equal(cds@int_metadata[['reduce_dim_metadata']][['Aligned']][['identity']][['matrix_id']], cds@int_metadata[['reduce_dim_metadata']][['Aligned']][['identity']][['model_id']])
  expect_equal(cds@int_metadata[['reduce_dim_metadata']][['Aligned']][['identity']][['prev_matrix_type']], 'matrix:PCA')
  expect_equal(cds@int_metadata[['reduce_dim_metadata']][['Aligned']][['identity']][['prev_matrix_id']], cds@int_metadata[['reduce_dim_metadata']][['PCA']][['identity']][['matrix_id']])
  expect_equal(cds@int_metadata[['reduce_dim_metadata']][['Aligned']][['identity']][['model_type']], 'matrix:Aligned')
  expect_equal(cds@int_metadata[['reduce_dim_metadata']][['Aligned']][['identity']][['model_id']], cds@reduce_dim_aux[['Aligned']][['model']][['identity']][['model_id']])

  expect_equal(cds@int_metadata[['reduce_dim_metadata']][['UMAP']][['identity']][['matrix_type']], 'matrix:UMAP')
  expect_equal(cds@int_metadata[['reduce_dim_metadata']][['UMAP']][['identity']][['matrix_id']], cds@int_metadata[['reduce_dim_metadata']][['UMAP']][['identity']][['model_id']])
  expect_equal(cds@int_metadata[['reduce_dim_metadata']][['UMAP']][['identity']][['prev_matrix_type']], 'matrix:Aligned')
  expect_equal(cds@int_metadata[['reduce_dim_metadata']][['UMAP']][['identity']][['prev_matrix_id']], cds@int_metadata[['reduce_dim_metadata']][['Aligned']][['identity']][['matrix_id']])
  expect_equal(cds@int_metadata[['reduce_dim_metadata']][['UMAP']][['identity']][['model_type']], 'matrix:UMAP')
  expect_equal(cds@int_metadata[['reduce_dim_metadata']][['UMAP']][['identity']][['model_id']], cds@reduce_dim_aux[['UMAP']][['model']][['identity']][['model_id']])

  expect_equal(cds@reduce_dim_aux[['PCA']][['model']][['identity']][['model_type']], 'matrix:PCA')
  
  expect_equal(cds@reduce_dim_aux[['Aligned']][['model']][['identity']][['model_type']], 'matrix:Aligned')
  expect_equal(cds@reduce_dim_aux[['Aligned']][['model']][['identity']][['prev_model_type']], 'matrix:PCA')
  expect_equal(cds@reduce_dim_aux[['Aligned']][['model']][['identity']][['prev_model_id']], cds@reduce_dim_aux[['PCA']][['model']][['identity']][['model_id']])

  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['identity']][['model_type']], 'matrix:UMAP')
  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['identity']][['prev_model_type']], 'matrix:Aligned')
  expect_equal(cds@reduce_dim_aux[['UMAP']][['model']][['identity']][['prev_model_id']], cds@reduce_dim_aux[['Aligned']][['model']][['identity']][['model_id']])

  save_transform_models(cds, directory_path=transform_models)

  cds2 <- load_mm_data(mtxfile,fdtfile, pdtfile,header=TRUE)
  cds2 <- load_transform_models(cds2, directory_path=transform_models)
  cds2 <- preprocess_transform(cds2)
  cds2 <- align_beta_transform(cds2)
  cds2 <- reduce_dimension_transform(cds2)

  expect_equal(cds2@int_metadata[['reduce_dim_metadata']][['PCA']][['identity']][['matrix_type']], 'matrix:PCA')
  expect_equal(cds2@int_metadata[['reduce_dim_metadata']][['PCA']][['identity']][['prev_matrix_type']], 'expr.mtx')
  expect_equal(cds2@int_metadata[['reduce_dim_metadata']][['PCA']][['identity']][['prev_matrix_id']], cds2@int_metadata[['counts_metadata']][['identity']][['matrix_id']])
  expect_equal(cds2@int_metadata[['reduce_dim_metadata']][['PCA']][['identity']][['model_type']], 'matrix:PCA')
  expect_equal(cds2@int_metadata[['reduce_dim_metadata']][['PCA']][['identity']][['model_id']], cds@int_metadata[['reduce_dim_metadata']][['PCA']][['identity']][['model_id']])

  expect_equal(cds2@int_metadata[['reduce_dim_metadata']][['Aligned']][['identity']][['matrix_type']], 'matrix:Aligned')
  expect_equal(cds2@int_metadata[['reduce_dim_metadata']][['Aligned']][['identity']][['prev_matrix_type']], 'matrix:PCA')
  expect_equal(cds2@int_metadata[['reduce_dim_metadata']][['Aligned']][['identity']][['prev_matrix_id']], cds2@int_metadata[['reduce_dim_metadata']][['PCA']][['identity']][['matrix_id']])
  expect_equal(cds2@int_metadata[['reduce_dim_metadata']][['Aligned']][['identity']][['model_type']], 'matrix:Aligned')
  expect_equal(cds2@int_metadata[['reduce_dim_metadata']][['Aligned']][['identity']][['model_id']], cds@int_metadata[['reduce_dim_metadata']][['Aligned']][['identity']][['model_id']])

  expect_equal(cds2@int_metadata[['reduce_dim_metadata']][['UMAP']][['identity']][['matrix_type']], 'matrix:UMAP')
  expect_equal(cds2@int_metadata[['reduce_dim_metadata']][['UMAP']][['identity']][['prev_matrix_type']], 'matrix:Aligned')
  expect_equal(cds2@int_metadata[['reduce_dim_metadata']][['UMAP']][['identity']][['prev_matrix_id']], cds2@int_metadata[['reduce_dim_metadata']][['Aligned']][['identity']][['matrix_id']])
  expect_equal(cds2@int_metadata[['reduce_dim_metadata']][['UMAP']][['identity']][['model_type']], 'matrix:UMAP')
  expect_equal(cds2@int_metadata[['reduce_dim_metadata']][['UMAP']][['identity']][['model_id']], cds@int_metadata[['reduce_dim_metadata']][['UMAP']][['identity']][['model_id']])


  expect_equal(cds2@int_metadata[['reduce_dim_metadata']][['PCA']][['identity']][['model_id']], cds@int_metadata[['reduce_dim_metadata']][['PCA']][['identity']][['model_id']])
  expect_equal(cds2@int_metadata[['reduce_dim_metadata']][['Aligned']][['identity']][['model_id']], cds@int_metadata[['reduce_dim_metadata']][['Aligned']][['identity']][['model_id']])
  expect_equal(cds2@int_metadata[['reduce_dim_metadata']][['UMAP']][['identity']][['model_id']], cds@int_metadata[['reduce_dim_metadata']][['UMAP']][['identity']][['model_id']])

  expect_equal(cds2@reduce_dim_aux[['PCA']][['model']][['identity']][['model_type']], 'matrix:PCA')  
  expect_equal(cds2@reduce_dim_aux[['PCA']][['model']][['identity']][['model_id']], cds@reduce_dim_aux[['PCA']][['model']][['identity']][['model_id']])
  expect_equal(cds2@reduce_dim_aux[['PCA']][['model']][['identity']][['model_path']], normalizePath(transform_models, mustWork=FALSE))  

  expect_equal(cds2@reduce_dim_aux[['Aligned']][['model']][['identity']][['model_type']], 'matrix:Aligned')
  expect_equal(cds2@reduce_dim_aux[['Aligned']][['model']][['identity']][['model_id']], cds@reduce_dim_aux[['Aligned']][['model']][['identity']][['model_id']])
  expect_equal(cds2@reduce_dim_aux[['Aligned']][['model']][['identity']][['prev_model_type']], cds@reduce_dim_aux[['Aligned']][['model']][['identity']][['prev_model_type']])
  expect_equal(cds2@reduce_dim_aux[['Aligned']][['model']][['identity']][['prev_model_id']], cds@reduce_dim_aux[['Aligned']][['model']][['identity']][['prev_model_id']])
  expect_equal(cds2@reduce_dim_aux[['Aligned']][['model']][['identity']][['model_path']], normalizePath(transform_models, mustWork=FALSE))

  expect_equal(cds2@reduce_dim_aux[['UMAP']][['model']][['identity']][['model_type']], 'matrix:UMAP')
  expect_equal(cds2@reduce_dim_aux[['UMAP']][['model']][['identity']][['model_id']], cds@reduce_dim_aux[['UMAP']][['model']][['identity']][['model_id']])
  expect_equal(cds2@reduce_dim_aux[['UMAP']][['model']][['identity']][['prev_model_type']], cds@reduce_dim_aux[['UMAP']][['model']][['identity']][['prev_model_type']])
  expect_equal(cds2@reduce_dim_aux[['UMAP']][['model']][['identity']][['prev_model_id']], cds@reduce_dim_aux[['UMAP']][['model']][['identity']][['prev_model_id']])
  expect_equal(cds2@reduce_dim_aux[['UMAP']][['model']][['identity']][['model_path']], normalizePath(transform_models, mustWork=FALSE))

  system(paste0('rm -r ', mtxfile, ' ', pdtfile, ' ', fdtfile, ' ', transform_models))
})
