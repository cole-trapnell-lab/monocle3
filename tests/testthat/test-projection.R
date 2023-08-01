context("test-projection")

skip_not_travis <- function ()
{
  if (identical(Sys.getenv("TRAVIS"), "true")) {
    return(invisible(TRUE))
  }
  skip("Not on Travis")
}


transform_models <- 'transform_models.01'

test_that("preprocess_transform using dgCMatrix", {
  nn_control <- list(method='annoy', metric='euclidean')
  cds1 <- load_a549()
  cds1 <- preprocess_cds(cds1, num_dim=50, nn_control=nn_control, build_nn_index=TRUE)
  cds1 <- preprocess_cds(cds1, method='LSI', num_dim=50, nn_control=nn_control, build_nn_index=TRUE)
  cds1 <- reduce_dimension(cds1, preprocess_method='PCA', nn_control=nn_control, build_nn_index=TRUE)
  save_transform_models(cds1, directory_path=transform_models)

  cds2 <- load_a549()
  cds2 <- load_transform_models(cds2, directory_path=transform_models)
  cds2 <- preprocess_transform(cds2)
  cds2 <- preprocess_transform(cds2, reduction_method='LSI')
  cds2 <- reduce_dimension_transform(cds2)

  expect_equal(reducedDims(cds1)[['PCA']][[1,1]], reducedDims(cds2)[['PCA']][[1,1]], tol=1e-6)
  expect_equal(reducedDims(cds1)[['PCA']][[2,1]], reducedDims(cds2)[['PCA']][[2,1]], tol=1e-6)
  expect_equal(reducedDims(cds1)[['PCA']][[1,2]], reducedDims(cds2)[['PCA']][[1,2]], tol=1e-6)

  expect_equal(reducedDims(cds1)[['LSI']][[1,1]], reducedDims(cds2)[['LSI']][[1,1]], tol=1e-6)
  expect_equal(reducedDims(cds1)[['LSI']][[2,1]], reducedDims(cds2)[['LSI']][[2,1]], tol=1e-6)
  expect_equal(reducedDims(cds1)[['LSI']][[1,2]], reducedDims(cds2)[['LSI']][[1,2]], tol=1e-6)

  expect_equal(reducedDims(cds1)[['UMAP']][[1,1]], reducedDims(cds2)[['UMAP']][[1,1]], tol=1e-6)
  expect_equal(reducedDims(cds1)[['UMAP']][[2,1]], reducedDims(cds2)[['UMAP']][[2,1]], tol=1e-6)
  expect_equal(reducedDims(cds1)[['UMAP']][[1,2]], reducedDims(cds2)[['UMAP']][[1,2]], tol=1e-6)

  system(paste0('rm -rf ', transform_models))
})


test_that("preprocess_transform using BPCells matrix", {
  nn_control <- list(method='annoy', metric='euclidean')
  cds1 <- load_a549(matrix_control=list(matrix_class='BPCells'))
  cds1 <- preprocess_cds(cds1, num_dim=50, nn_control=nn_control, build_nn_index=TRUE)
  cds1 <- preprocess_cds(cds1, method='LSI', num_dim=50, nn_control=nn_control, build_nn_index=TRUE)
  cds1 <- reduce_dimension(cds1, preprocess_method='PCA', nn_control=nn_control, build_nn_index=TRUE)
  save_transform_models(cds1, directory_path=transform_models)

  cds2 <- load_a549(matrix_control=list(matrix_class='BPCells'))
  cds2 <- load_transform_models(cds2, directory_path=transform_models)
  cds2 <- preprocess_transform(cds2)
  cds2 <- preprocess_transform(cds2, reduction_method='LSI')
  cds2 <- reduce_dimension_transform(cds2)

  expect_equal(reducedDims(cds1)[['PCA']][[1,1]], reducedDims(cds2)[['PCA']][[1,1]], tol=1e-6)
  expect_equal(reducedDims(cds1)[['PCA']][[2,1]], reducedDims(cds2)[['PCA']][[2,1]], tol=1e-6)
  expect_equal(reducedDims(cds1)[['PCA']][[1,2]], reducedDims(cds2)[['PCA']][[1,2]], tol=1e-6)

  expect_equal(reducedDims(cds1)[['LSI']][[1,1]], reducedDims(cds2)[['LSI']][[1,1]], tol=1e-6)
  expect_equal(reducedDims(cds1)[['LSI']][[2,1]], reducedDims(cds2)[['LSI']][[2,1]], tol=1e-6)
  expect_equal(reducedDims(cds1)[['LSI']][[1,2]], reducedDims(cds2)[['LSI']][[1,2]], tol=1e-6)

  expect_equal(reducedDims(cds1)[['UMAP']][[1,1]], reducedDims(cds2)[['UMAP']][[1,1]], tol=1e-6)
  expect_equal(reducedDims(cds1)[['UMAP']][[2,1]], reducedDims(cds2)[['UMAP']][[2,1]], tol=1e-6)
  expect_equal(reducedDims(cds1)[['UMAP']][[1,2]], reducedDims(cds2)[['UMAP']][[1,2]], tol=1e-6)

  system(paste0('rm -rf ', transform_models))
})
