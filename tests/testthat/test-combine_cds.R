testthat::skip_if_offline()
cds <- load_a549()
cds2 <- monocle3:::load_worm_embryo()

test_that("combine_cds works", {
  # all genes shared
  comb <- combine_cds(list(cds, cds),
                      keep_all_genes = TRUE,
                      cell_names_unique = FALSE)
  testthat::expect_equal(nrow(exprs(comb)),
                         nrow(exprs(cds)))
  testthat::expect_equal(ncol(exprs(comb)),
                         ncol(exprs(cds)) + ncol(exprs(cds)))

  comb <- combine_cds(list(cds, cds),
                      keep_all_genes = FALSE,
                      cell_names_unique = FALSE)
  testthat::expect_equal(nrow(exprs(comb)),
                         nrow(exprs(cds)))
  testthat::expect_equal(ncol(exprs(comb)),
                         ncol(exprs(cds)) + ncol(exprs(cds)))

  testthat::expect_error(comb <- combine_cds(list(cds, cds),
                                             keep_all_genes = TRUE,
                                             cell_names_unique = TRUE),
                         "Cell names are not unique across CDSs - cell_names_unique must be TRUE.")
  testthat::expect_error(comb <- combine_cds(list(cds, cds),
                                             keep_all_genes = FALSE,
                                             cell_names_unique = TRUE),
                         "Cell names are not unique across CDSs - cell_names_unique must be TRUE.")

  # some genes shared
  cds3 <- cds[1:100,]
  cds3 <- cds3[,Matrix::rowSums(exprs(cds3)) != 0]
  expect_warning(comb <- combine_cds(list(cds, cds3),
                                     keep_all_genes = TRUE,
                                     cell_names_unique = FALSE))
  testthat::expect_equal(nrow(exprs(comb)),
                         nrow(exprs(cds)))
  testthat::expect_equal(ncol(exprs(comb)),
                         ncol(exprs(cds)) + ncol(exprs(cds3)))

  expect_warning(comb <- combine_cds(list(cds, cds3),
                                     keep_all_genes = FALSE,
                                     cell_names_unique = FALSE))
  testthat::expect_equal(nrow(exprs(comb)),
                         nrow(exprs(cds3)))
  testthat::expect_equal(ncol(exprs(comb)),
                         ncol(exprs(cds)) + ncol(exprs(cds3)))

  testthat::expect_error(comb <- combine_cds(list(cds, cds3),
                                             keep_all_genes = TRUE,
                                             cell_names_unique = TRUE),
                         "Cell names are not unique across CDSs - cell_names_unique must be TRUE.")
  testthat::expect_error(comb <- combine_cds(list(cds, cds3),
                                             keep_all_genes = FALSE,
                                             cell_names_unique = TRUE),
                         "Cell names are not unique across CDSs - cell_names_unique must be TRUE.")


  # no genes shared
  expect_warning(comb <- combine_cds(list(cds, cds2),
                                     keep_all_genes = TRUE,
                                     cell_names_unique = FALSE),
                 "No genes are shared amongst all the CDS objects.")
  testthat::expect_equal(nrow(exprs(comb)),
                         nrow(exprs(cds)) + nrow(exprs(cds2)))
  testthat::expect_equal(ncol(exprs(comb)),
                         ncol(exprs(cds)) + ncol(exprs(cds2)))

  expect_error(comb <- combine_cds(list(cds, cds2),
                                   keep_all_genes = FALSE,
                                   cell_names_unique = FALSE),
               paste("No genes are shared amongst all the CDS objects. To generate a",
                     "combined CDS with all genes, use keep_all_genes = TRUE"))

  expect_warning(comb <- combine_cds(list(cds, cds2),
                                     keep_all_genes = TRUE,
                                     cell_names_unique = TRUE),
                 "No genes are shared amongst all the CDS objects.")

  testthat::expect_equal(nrow(exprs(comb)),
                         nrow(exprs(cds)) + nrow(exprs(cds2)))
  testthat::expect_equal(ncol(exprs(comb)),
                         ncol(exprs(cds)) + ncol(exprs(cds2)))

  expect_error(comb <- combine_cds(list(cds, cds2),
                                   keep_all_genes = FALSE,
                                   cell_names_unique = TRUE),
               paste("No genes are shared amongst all the CDS objects. To generate a",
                     "combined CDS with all genes, use keep_all_genes = TRUE"))

  # triples
  expect_warning(expect_error(comb <- combine_cds(list(cds, cds2, cds3),
                                                  keep_all_genes = TRUE,
                                                  cell_names_unique = TRUE),
                              "Cell names are not unique across CDSs - cell_names_unique must be TRUE."),
                 "No genes are shared amongst all the CDS objects.")

  expect_error(comb <- combine_cds(list(cds, cds2, cds3),
                                   keep_all_genes = FALSE,
                                   cell_names_unique = TRUE),
               paste("No genes are shared amongst all the CDS objects. To generate a",
                     "combined CDS with all genes, use keep_all_genes = TRUE"))

  expect_warning(comb <- combine_cds(list(cds, cds2, cds3),
                                     keep_all_genes = TRUE,
                                     cell_names_unique = FALSE),
                 "No genes are shared amongst all the CDS objects.")
  testthat::expect_equal(nrow(exprs(comb)),
                         length(unique(c(row.names(exprs(cds)),
                                         row.names(exprs(cds2)),
                                         row.names(exprs(cds3))))))
  testthat::expect_equal(ncol(exprs(comb)),
                         ncol(exprs(cds)) + ncol(exprs(cds2)) + ncol(exprs(cds3)))

  expect_error(comb <- combine_cds(list(cds, cds2, cds3),
                                   keep_all_genes = FALSE,
                                   cell_names_unique = FALSE),
               paste("No genes are shared amongst all the CDS objects. To generate a",
                     "combined CDS with all genes, use keep_all_genes = TRUE"))



})
