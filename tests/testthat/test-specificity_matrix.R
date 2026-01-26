context("test-specificity_matrix")

makeprobsvec <- getFromNamespace("makeprobsvec", "monocle3")
JSdistVec <- getFromNamespace("JSdistVec", "monocle3")
specificity_matrix <- getFromNamespace("specificity_matrix", "monocle3")

ref_specificity_matrix <- function(agg_expr_matrix) {
  if (!is.matrix(agg_expr_matrix)) {
    agg_expr_matrix <- as.matrix(agg_expr_matrix)
  }
  out <- lapply(row.names(agg_expr_matrix), function(x) {
    agg_exprs <- as.numeric(agg_expr_matrix[x, ])
    agg_exprs <- makeprobsvec(agg_exprs)
    perfect_spec_matrix <- diag(ncol(agg_expr_matrix))
    sapply(1:ncol(agg_expr_matrix), function(col_idx) {
      1 - JSdistVec(agg_exprs, perfect_spec_matrix[, col_idx])
    })
  })
  out <- do.call(rbind, out)
  colnames(out) <- colnames(agg_expr_matrix)
  row.names(out) <- row.names(agg_expr_matrix)
  out
}

test_that("specificity_matrix matches JS reference on small input", {
  mat <- matrix(c(1, 2, 3,
                  4, 5, 6),
                nrow = 2, byrow = TRUE)
  row.names(mat) <- c("g1", "g2")
  colnames(mat) <- c("c1", "c2", "c3")

  ref <- ref_specificity_matrix(mat)
  res <- specificity_matrix(mat, cores = 1)
  expect_equal(res, ref, tolerance = 1e-12)
})

test_that("specificity_matrix handles NA and negative values consistently", {
  mat <- matrix(c(1, NA, 3,
                  -1, 2, 3),
                nrow = 2, byrow = TRUE)
  row.names(mat) <- c("g1", "g2")
  colnames(mat) <- c("c1", "c2", "c3")

  ref <- ref_specificity_matrix(mat)
  res <- specificity_matrix(mat, cores = 1)
  expect_equal(res, ref, tolerance = 1e-12)
})
