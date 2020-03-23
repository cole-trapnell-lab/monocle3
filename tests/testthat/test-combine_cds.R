testthat::skip_if_offline()
cds <- load_a549()
cds2 <- monocle3:::load_worm_embryo()

test_that("combine_cds works", {
  check_comb <- function(cds1, cds2, comb, keep_all_genes) {

    if (keep_all_genes) {
    testthat::expect_equal(nrow(exprs(comb)),
                           length(unique(c(row.names(exprs(cds1)),
                                           row.names(exprs(cds2))))))

      gene_order <- unique(c(row.names(exprs(cds1)), row.names(exprs(cds2))))
      temp1 <- Matrix::rowSums(exprs(cds1))[gene_order]
      names(temp1) <- gene_order
      temp1[is.na(temp1)] <- 0
      temp2 <- Matrix::rowSums(exprs(cds2))[gene_order]
      temp2[is.na(temp2)] <- 0
      names(temp2) <- gene_order
      testthat::expect_equal(temp1 + temp2,
                             Matrix::rowSums(exprs(comb))[gene_order])
    } else {

        testthat::expect_equal(nrow(exprs(comb)),
                               length(intersect(row.names(exprs(cds1)),
                                               row.names(exprs(cds2)))))

      gene_order <- intersect(row.names(exprs(cds1)), row.names(exprs(cds2)))
      testthat::expect_equal(Matrix::rowSums(exprs(cds1))[gene_order] +
                               Matrix::rowSums(exprs(cds2))[gene_order],
                             Matrix::rowSums(exprs(comb))[gene_order])

    }
    testthat::expect_equal(ncol(exprs(comb)),
                           ncol(exprs(cds1)) + ncol(exprs(cds2)))
    if ("ENSG00000260917.1" %in% row.names(fData(cds1)) & "ENSG00000260917.1" %in% row.names(fData(cds2))) {
      testthat::expect_equal(sum(exprs(cds1)["ENSG00000260917.1",]) +
                               sum(exprs(cds2)["ENSG00000260917.1",]),
                             sum(exprs(comb)["ENSG00000260917.1",]))
    }



  }

  # all genes shared
  testthat::expect_warning(combine_cds(list(cds, cds),
                                       keep_all_genes = TRUE,
                                       cell_names_unique = FALSE),
                           paste0("By default, the combine_cds function adds ",
                                  "a column called'sample' which indicates ",
                                  "which initial cds a cell came from. One ",
                                  "or more of your input cds objects ",
                                  "contains a 'sample' column, which will be ",
                                  "overwritten. We recommend you rename this ",
                                  "column or provide an alternative column ",
                                  "name using the 'sample_col_name' ",
                                  "parameter."))
  names(pData(cds))[names(pData(cds)) == "sample"] <- "sample_type"
  names(pData(cds2))[names(pData(cds2)) == "sample"] <- "sample_type"

  comb_h <- combine_cds(list(cds, cds),
                      keep_all_genes = TRUE,
                      cell_names_unique = FALSE,
                      sample_col_name = "hannah")

  comb <- combine_cds(list(cds, cds),
                      keep_all_genes = TRUE,
                      cell_names_unique = FALSE)

  testthat::expect_true("hannah" %in% names(colData(comb_h)))
  testthat::expect_true(all(colData(comb)$sample == colData(comb_h)$hannah))

  check_comb(cds, cds, comb, keep_all_genes = TRUE)

  comb <- combine_cds(list(cds, cds),
                      keep_all_genes = FALSE,
                      cell_names_unique = FALSE)

  check_comb(cds, cds, comb, keep_all_genes = FALSE)

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
  cds3 <- detect_genes(cds3)
  testthat::expect_warning(comb2 <- combine_cds(list(cds, cds3),
                                                keep_all_genes = TRUE,
                                                cell_names_unique = FALSE))
  testthat::expect_equal(sum(fData(comb2)$num_cells_expressed == "conf"), 100)
  check_comb(cds, cds3, comb2, keep_all_genes = TRUE)

  testthat::expect_warning(comb3 <- combine_cds(list(cds, cds3),
                                                keep_all_genes = FALSE,
                                                cell_names_unique = FALSE))
  check_comb(cds, cds3, comb3, keep_all_genes = FALSE)

  testthat::expect_error(comb <- combine_cds(list(cds, cds3),
                                             keep_all_genes = TRUE,
                                             cell_names_unique = TRUE),
                         "Cell names are not unique across CDSs - cell_names_unique must be TRUE.")
  testthat::expect_error(comb <- combine_cds(list(cds, cds3),
                                             keep_all_genes = FALSE,
                                             cell_names_unique = TRUE),
                         "Cell names are not unique across CDSs - cell_names_unique must be TRUE.")

  # some genes shared mixed order
  cds3 <- cds[c(50:100, 1:49),]
  cds3 <- cds3[,Matrix::rowSums(exprs(cds3)) != 0]

  temp <- exprs(cds3)
  row.names(temp)[1] <- "HAP00000240216.7"
  temp2 <- fData(cds3)
  row.names(temp2)[1] <- "HAP00000240216.7"

  cds3 <- suppressWarnings(new_cell_data_set(temp, pData(cds3), temp2))
  cds3 <- cds3[,Matrix::rowSums(exprs(cds3)) != 0]

  testthat::expect_warning(comb4 <- combine_cds(list(cds, cds3),
                                                keep_all_genes = TRUE,
                                                cell_names_unique = FALSE))

  check_comb(cds, cds3, comb4, keep_all_genes = TRUE)

  testthat::expect_warning(comb5 <- combine_cds(list(cds, cds3),
                                                keep_all_genes = FALSE,
                                                cell_names_unique = FALSE))

  check_comb(cds, cds3, comb5, keep_all_genes = FALSE)

  testthat::expect_error(comb <- combine_cds(list(cds, cds3),
                                             keep_all_genes = TRUE,
                                             cell_names_unique = TRUE),
                         "Cell names are not unique across CDSs - cell_names_unique must be TRUE.")
  testthat::expect_error(comb <- combine_cds(list(cds, cds3),
                                             keep_all_genes = FALSE,
                                             cell_names_unique = TRUE),
                         "Cell names are not unique across CDSs - cell_names_unique must be TRUE.")

  # no genes shared
  testthat::expect_warning(comb6 <- combine_cds(list(cds, cds2),
                                     keep_all_genes = TRUE,
                                     cell_names_unique = FALSE),
                 "No genes are shared amongst all the CDS objects.")

  check_comb(cds, cds2, comb6, keep_all_genes = TRUE)
  testthat::expect_error(comb <- combine_cds(list(cds, cds2),
                                   keep_all_genes = FALSE,
                                   cell_names_unique = FALSE),
               paste("No genes are shared amongst all the CDS objects. To generate a",
                     "combined CDS with all genes, use keep_all_genes = TRUE"))

  testthat::expect_warning(comb7 <- combine_cds(list(cds, cds2),
                                     keep_all_genes = TRUE,
                                     cell_names_unique = TRUE),
                 "No genes are shared amongst all the CDS objects.")

  check_comb(cds, cds2, comb7, keep_all_genes = TRUE)

  testthat::expect_error(comb <- combine_cds(list(cds, cds2),
                                   keep_all_genes = FALSE,
                                   cell_names_unique = TRUE),
               paste("No genes are shared amongst all the CDS objects. To generate a",
                     "combined CDS with all genes, use keep_all_genes = TRUE"))

  # triples
  testthat::expect_warning(expect_error(comb <- combine_cds(list(cds, cds2, cds3),
                                                  keep_all_genes = TRUE,
                                                  cell_names_unique = TRUE),
                              "Cell names are not unique across CDSs - cell_names_unique must be TRUE."),
                 "No genes are shared amongst all the CDS objects.")

  testthat::expect_error(comb <- combine_cds(list(cds, cds2, cds3),
                                   keep_all_genes = FALSE,
                                   cell_names_unique = TRUE),
               paste("No genes are shared amongst all the CDS objects. To generate a",
                     "combined CDS with all genes, use keep_all_genes = TRUE"))

  testthat::expect_warning(comb6 <- combine_cds(list(cds, cds2, cds3),
                                     keep_all_genes = TRUE,
                                     cell_names_unique = FALSE),
                 "No genes are shared amongst all the CDS objects.")
  testthat::expect_equal(nrow(exprs(comb6)),
                         length(unique(c(row.names(exprs(cds)),
                                         row.names(exprs(cds2)),
                                         row.names(exprs(cds3))))))
  testthat::expect_equal(ncol(exprs(comb6)),
                         ncol(exprs(cds)) + ncol(exprs(cds2)) +
                           ncol(exprs(cds3)))
  testthat::expect_equal(sum(exprs(cds)["ENSG00000260917.1",]) +
                           sum(exprs(cds3)["ENSG00000260917.1",]),
                         sum(exprs(comb6)["ENSG00000260917.1",]))

  testthat::expect_error(comb <- combine_cds(list(cds, cds2, cds3),
                                   keep_all_genes = FALSE,
                                   cell_names_unique = FALSE),
               paste("No genes are shared amongst all the CDS objects. To generate a",
                     "combined CDS with all genes, use keep_all_genes = TRUE"))



})
