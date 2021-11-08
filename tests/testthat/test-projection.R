context("test-projection")

skip_not_travis <- function ()
{
  if (identical(Sys.getenv("TRAVIS"), "true")) {
    return(invisible(TRUE))
  }
  skip("Not on Travis")
}


transform_models <- 'transform_models.01'

test_that("preprocess_transform", {
  nn_control <- list(method='annoy', metric='euclidean')
  cds <- load_a549()
  cds <- preprocess_cds(cds, num_dim=50, nn_control=nn_control, build_nn_index=TRUE)
  cds <- preprocess_cds(cds, method='LSI', num_dim=50, nn_control=nn_control, build_nn_index=TRUE)
  cds <- reduce_dimension(cds, preprocess_method='PCA', nn_control=nn_control, build_nn_index=TRUE)
  save_transform_models(cds, directory_path=transform_models)

  cds2 <- load_a549()
  cds2 <- load_transform_models(cds2, directory_path=transform_models)
  cds2 <- preprocess_transform(cds2)
  cds2 <- preprocess_transform(cds2, reduction_method='LSI')
  cds2 <- reduce_dimension_transform(cds2)

  expect_equal(reducedDims(cds)[['PCA']][[1,1]], reducedDims(cds2)[['PCA']][[1,1]], tol=1e-6)
  expect_equal(reducedDims(cds)[['PCA']][[2,1]], reducedDims(cds2)[['PCA']][[2,1]], tol=1e-6)
  expect_equal(reducedDims(cds)[['PCA']][[1,2]], reducedDims(cds2)[['PCA']][[1,2]], tol=1e-6)

  expect_equal(reducedDims(cds)[['LSI']][[1,1]], reducedDims(cds2)[['LSI']][[1,1]], tol=1e-6)
  expect_equal(reducedDims(cds)[['LSI']][[2,1]], reducedDims(cds2)[['LSI']][[2,1]], tol=1e-6)
  expect_equal(reducedDims(cds)[['LSI']][[1,2]], reducedDims(cds2)[['LSI']][[1,2]], tol=1e-6)

  expect_equal(reducedDims(cds)[['UMAP']][[1,1]], reducedDims(cds2)[['UMAP']][[1,1]], tol=1e-6)
  expect_equal(reducedDims(cds)[['UMAP']][[2,1]], reducedDims(cds2)[['UMAP']][[2,1]], tol=1e-6)
  expect_equal(reducedDims(cds)[['UMAP']][[1,2]], reducedDims(cds2)[['UMAP']][[1,2]], tol=1e-6)

  system(paste0('rm -rf ', transform_models))
})

