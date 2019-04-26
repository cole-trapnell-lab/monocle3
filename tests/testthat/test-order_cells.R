context("test-order_cells")

cds <- load_a549()
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = 20)
cds <- reduce_dimension(cds)
cds <- partition_cells(cds)
cds <- learn_graph(cds)


test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
