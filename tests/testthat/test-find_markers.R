context("test-find_markers")

library(dplyr)
set.seed(1000)
cds <- load_a549()
assigned_type_marker_test_res = top_markers(cds,
                                            group_cells_by="top_oligo",
                                            reference_cells=1000,
                                            cores=8)

garnett_markers <- assigned_type_marker_test_res %>%
  filter(marker_test_p_value < 0.01 & specificity >= 0.5) %>%
  group_by(cell_group) %>%
  top_n(5, marker_score)


test_that("generate_garnett_marker_file_works", {
  temp <- tempdir()
  generate_garnett_marker_file(garnett_markers, file=paste0(temp,"/marker_file.txt"))
  marker <- readLines(paste0(temp,"/marker_file.txt"))
  expect_equal(length(marker), 42)
  expect_equal(marker[1], "> Cell type Dex_0_AD04")

  garnett_markers$cell_group[1] <- "test_(_:_)_#_,_>"
  expect_warning(generate_garnett_marker_file(garnett_markers, file=paste0(temp,"/marker_file.txt")))
  marker <- readLines(paste0(temp,"/marker_file.txt"))
  expect_equal(length(marker), 45)
  expect_equal(marker[1], "> Cell type test_._._._._._.")
  unlink(temp)
})
