context("test-io")
skip_not_travis <- function ()
{
  if (identical(Sys.getenv("TRAVIS"), "true")) {
    return(invisible(TRUE))
  }
  skip("Not on Travis")
}

# functions
#   select_matrix_parameter_value <- function(parameter, matrix_control, matrix_control_default, default_value)
#   check_matrix_control <- function(matrix_control=list(), control_type=c('unrestricted', 'pca'), check_conditional=FALSE)
#   set_matrix_control <- function(matrix_control=list(), matrix_control_default=list(), control_type=c('unrestricted', 'pca'))
#   show_matrix_control <- function(matrix_control, label=NULL)
#   push_matrix_path <- function(mat)
#   bpcells_find_base_matrix <- function(mat)
#   get_matrix_class <- function(mat)
#   get_matrix_info <- function(mat)
#   show_matrix_info <- function(matrix_info, indent='')
#   set_matrix_class <- function(mat, matrix_control=list())
#   rm_bpcells_dir <- function(mat)
#   set_cds_row_order_matrix <- function(cds)
#   convert_counts_matrix <- function(cds, matrix_control=list(matrix_class='BPCells'))
#

# functions that need additional review
#   (done) check_matrix_control <- function(matrix_control=list(), control_type=c('unrestricted', 'pca'), check_conditional=FALSE)
#   (done) set_matrix_control <- function(matrix_control=list(), matrix_control_default=list(), control_type=c('unrestricted', 'pca'))
#   (done) set_matrix_class <- function(mat, matrix_control=list())
#   (done) convert_counts_matrix <- function(cds, matrix_control=list(matrix_class='BPCells'))
#   (done) get_matrix_info <- function(mat)
#   (done) show_matrix_info <- function(matrix_info, indent='')
#   (done) rm_bpcells_dir <- function(mat)
#   (done) set_cds_row_order_matrix <- function(cds)

# (exported) functions that have a matrix_control parameter. Check that they have a call to set_matrix_control/check_matrix_control
#   io.R
#     (OK) load_a549()
#     (OK) load_worm_embryo()
#     (OK) load_worm_l2()
#     (OK) load_mm_data()
#     (OK) load_mtx_data()
#     (OK) load_monocle_objects()
#
#   load_cellranger_data.R
#     (OK) load_cellranger_data()
#
#   matrix.R
#     (OK) check_matrix_control()
#     (OK) set_matrix_control()
#     (OK) set_matrix_class()
#     (OK) convert_counts_matrix()
#     (OK) load_bpcells_matrix_dir()
#
#   pca.R
#     (OK) set_pca_matrix_control()
#     (OK) bpcells_prcomp_irlba()
#
#   projection.R
#     (OK) preprocess_transform()
# 
#   utils.R
#     (OK) combine_cds()
#

# (all) functions that have a matrix_control parameter: is matrix_control used only for set_matrix_class?
#   io.R
#     (yes) load_a549()
#     (yes) load_worm_embryo()
#     (yes) load_worm_l2()
#     (yes) load_mm_data()
#     (yes) load_mtx_data()
#     (yes) load_bpcells_matrix_dir()
#     (yes) load_monocle_objects() # used only in call to load_bpcells_matrix_dir
#
#   load_cellranger_data.R
#     (yes) load_cellranger_data()
#
#   matrix.R
#     (NA) check_matrix_control() # matrix_control is checked
#     (NA) set_matrix_control() # matrix_control defaults set
#     (NA) set_matrix_class()
#     (NA) compare_matrix_control() # compares matrix_info to matrix_control
#     (yes) convert_counts_matrix()
#
#   pca.R
#     (NA) set_pca_matrix_control() # set matrix_control defaults: used only it projection.R (check)
#     (yes) bpcells_prcomp_irlba()
#
#   projection.R
#     (yes) preprocess_transform()
#
#   utils.R
#     (yes) combine_cds()
#

# # Test priorities.
# 
# (done) check_matrix_control
# (done) set_matrix_control
# (done) set_matrix_class
# 
# (done) convert_counts_matrix
# 
# (done) save_monocle_objects
# (done) load_monocle_objects
# 

#
# These tests are not exhaustive but may check more frequently used and
# important combinations.
#

test_that("check_matrix_control", {

  # matrix_class = bad_class
  matrix_control <- list(matrix_class='bad_class')
  testthat::expect_error(check_matrix_control(matrix_control=matrix_control))

  # matrix_class = dgCMatrix
  matrix_control <- list(matrix_class='dgCMatrix')
  testthat::expect_true(check_matrix_control(matrix_control=matrix_control))
  testthat::expect_true(check_matrix_control(matrix_control=matrix_control, control_type='unrestricted', check_conditional=FALSE))
  testthat::expect_true(check_matrix_control(matrix_control=matrix_control, control_type='unrestricted', check_conditional=TRUE))
  testthat::expect_true(check_matrix_control(matrix_control=matrix_control, control_type='pca', check_conditional=TRUE))


  # matrix_class = BPCells
  matrix_control <- list(matrix_class='BPCells')
  testthat::expect_true(check_matrix_control(matrix_control=matrix_control))

  testthat::expect_true(check_matrix_control(matrix_control=matrix_control, control_type='unrestricted', check_conditional=FALSE))
  testthat::expect_true(check_matrix_control(matrix_control=matrix_control, control_type='pca', check_conditional=FALSE))

  testthat::expect_error(check_matrix_control(matrix_control=matrix_control, control_type='unrestricted', check_conditional=TRUE))
  testthat::expect_error(check_matrix_control(matrix_control=matrix_control, control_type='pca', check_conditional=TRUE))

  matrix_control <- list(matrix_class='BPCells', matrix_mode='mem')
  testthat::expect_error(check_matrix_control(matrix_control=matrix_control, control_type='unrestricted', check_conditional=FALSE))
  testthat::expect_error(check_matrix_control(matrix_control=matrix_control, control_type='pca', check_conditional=FALSE))

  matrix_control <- list(matrix_class='BPCells', matrix_mode='dir')
  testthat::expect_true(check_matrix_control(matrix_control=matrix_control, control_type='unrestricted', check_conditional=FALSE))
  testthat::expect_true(check_matrix_control(matrix_control=matrix_control, control_type='pca', check_conditional=FALSE))

  matrix_control <- list(matrix_class='BPCells', matrix_mode='dir')
  testthat::expect_error(check_matrix_control(matrix_control=matrix_control, control_type='unrestricted', check_conditional=TRUE))
  testthat::expect_error(check_matrix_control(matrix_control=matrix_control, control_type='pca', check_conditional=TRUE))

  matrix_control <- list(matrix_class='BPCells', matrix_mode='dir', matrix_type='double', matrix_compress=FALSE, matrix_path='.', matrix_buffer_size=8192L, matrix_bpcells_copy=TRUE)
  testthat::expect_true(check_matrix_control(matrix_control=matrix_control, control_type='unrestricted', check_conditional=TRUE))
  testthat::expect_true(check_matrix_control(matrix_control=matrix_control, control_type='pca', check_conditional=TRUE))

  matrix_control <- list(matrix_class='BPCells', matrix_mode='dir', matrix_type='float', matrix_compress=FALSE, matrix_path='.', matrix_buffer_size=8192L, matrix_bpcells_copy=TRUE)
  testthat::expect_true(check_matrix_control(matrix_control=matrix_control, control_type='unrestricted', check_conditional=TRUE))
  testthat::expect_true(check_matrix_control(matrix_control=matrix_control, control_type='pca', check_conditional=TRUE))

  matrix_control <- list(matrix_class='BPCells', matrix_mode='dir', matrix_type='unint32_t', matrix_compress=FALSE, matrix_path='.', matrix_buffer_size=8192L, matrix_bpcells_copy=TRUE)
  testthat::expect_error(check_matrix_control(matrix_control=matrix_control, control_type='unrestricted', check_conditional=TRUE))
  testthat::expect_error(check_matrix_control(matrix_control=matrix_control, control_type='pca', check_conditional=TRUE))

} )


test_that("set_matrix_control", {

  # Check defaults are set.
  matrix_control <- list()
  matrix_control_default <- get_global_variable('matrix_control_csparsematrix_unrestricted')
  matrix_control_res <- set_matrix_control(matrix_control=matrix_control, matrix_control_default=matrix_control_default, control_type='unrestricted')
  testthat::expect_true(all(matrix_control_res %in% matrix_control_default))

  # Check dgCMatrix defaults are set.
  matrix_control <- list(matrix_class='dgCMatrix')
  matrix_control_default <- get_global_variable('matrix_control_csparsematrix_unrestricted')
  matrix_control_res <- set_matrix_control(matrix_control=matrix_control, matrix_control_default=matrix_control_default, control_type='unrestricted')
  testthat::expect_true(all(matrix_control_res %in% matrix_control_default))

  # Check BPCells defaults are set.
  matrix_control <- list(matrix_class='BPCells')
  matrix_control_default <- get_global_variable('matrix_control_bpcells_unrestricted')
  matrix_control_res <- set_matrix_control(matrix_control=matrix_control, matrix_control_default=matrix_control_default, control_type='unrestricted')
  testthat::expect_true(all(matrix_control_res %in% matrix_control_default))

  # Check BPCells non-defaults are set.
  matrix_control <- get_global_variable('matrix_control_bpcells_unrestricted')
  matrix_control[['matrix_type']] <- 'float'
  matrix_control_default <- get_global_variable('matrix_control_bpcells_unrestricted')
  matrix_control_res <- set_matrix_control(matrix_control=matrix_control, matrix_control_default=matrix_control_default, control_type='unrestricted')
  testthat::expect_true(all(matrix_control_res %in% matrix_control))

  matrix_control <- get_global_variable('matrix_control_bpcells_unrestricted')
  matrix_control[['matrix_compress']] <- TRUE
  matrix_control_default <- get_global_variable('matrix_control_bpcells_unrestricted')
  matrix_control_res <- set_matrix_control(matrix_control=matrix_control, matrix_control_default=matrix_control_default, control_type='unrestricted')
  testthat::expect_true(all(matrix_control_res %in% matrix_control))

  matrix_control <- get_global_variable('matrix_control_bpcells_unrestricted')
  matrix_control[['matrix_path']] <- 'new_dir'
  matrix_control_default <- get_global_variable('matrix_control_bpcells_unrestricted')
  matrix_control_res <- set_matrix_control(matrix_control=matrix_control, matrix_control_default=matrix_control_default, control_type='unrestricted')
  testthat::expect_true(all(matrix_control_res %in% matrix_control))

  matrix_control <- get_global_variable('matrix_control_bpcells_unrestricted')
  matrix_control[['matrix_bpcells_copy']] <- FALSE
  matrix_control_default <- get_global_variable('matrix_control_bpcells_unrestricted')
  matrix_control_res <- set_matrix_control(matrix_control=matrix_control, matrix_control_default=matrix_control_default, control_type='unrestricted')
  testthat::expect_true(all(matrix_control_res %in% matrix_control))
} )


test_that("set_matrix_class", {

  # Check dgCMatrix -> dgCMatrix matrix
  matrix_control <- get_global_variable('matrix_control_csparsematrix_unrestricted')
  mat1 <- counts(load_a549())
  mat2 <- set_matrix_class(mat1, matrix_control=matrix_control) 
  testthat::expect_true(is(mat2, 'dgCMatrix'))

  # Check BPCells -> BPCells matrix with copy
  matrix_control <- get_global_variable('matrix_control_bpcells_unrestricted')
  mat1 <- counts(load_a549(matrix_control=list(matrix_class='BPCells')))
  mat2 <- set_matrix_class(mat1, matrix_control=matrix_control)
  testthat::expect_true(is(mat2, 'IterableMatrix'))
  testthat::expect_false(mat2@dir == mat1@dir)
  unlink(c(mat1@dir, mat2@dir), recursive=TRUE)
  rm(mat1, mat2)

  # Check BPCells -> BPCells matrix without copy
  matrix_control <- get_global_variable('matrix_control_bpcells_unrestricted')
  matrix_control[['matrix_bpcells_copy']] <- FALSE
  mat1 <- counts(load_a549(matrix_control=list(matrix_class='BPCells')))
  mat2 <- set_matrix_class(mat1, matrix_control=matrix_control)
  testthat::expect_true(is(mat2, 'IterableMatrix'))
  testthat::expect_true(mat2@dir == mat1@dir)
  unlink(c(mat1@dir, mat2@dir), recursive=TRUE)
  rm(mat1, mat2)

  # Check dgCMatrix -> BPCells matrix
  matrix_control <- get_global_variable('matrix_control_bpcells_unrestricted')
  mat1 <- counts(load_a549(matrix_control=list(matrix_class='dgCMatrix')))
  mat2 <- set_matrix_class(mat1, matrix_control=matrix_control)
  testthat::expect_true(is(mat2, 'IterableMatrix'))
  unlink(mat2@dir, recursive=TRUE)
  rm(mat1, mat2)

  # Check BPCells -> dgCMatrix matrix
  matrix_control <- get_global_variable('matrix_control_csparsematrix_unrestricted')
  mat1 <- counts(load_a549(matrix_control=list(matrix_class='BPCells')))
  mat2 <- set_matrix_class(mat1, matrix_control=matrix_control)
  testthat::expect_true(is(mat2, 'dgCMatrix'))
  unlink(mat1@dir, recursive=TRUE)
  rm(mat1, mat2)

} )


test_that("convert_counts_matrix", {
  cds1 <- load_a549()
  testthat::expect_true(is(counts(cds1), 'dgCMatrix'))
  matrix_control <- list(matrix_class='BPCells')
  cds2 <- convert_counts_matrix(cds1, matrix_control=matrix_control)
  testthat::expect_true(is(counts(cds2), 'IterableMatrix'))
} )


test_that("save_monocle_objects and load_monocle_objects", {
  cds1 <- load_a549(matrix_control=list(matrix_class='BPCells'))
  save_monocle_objects(cds1, directory_path='monocle_objects_test.tmp')
  cds2 <- load_monocle_objects(directory_path='monocle_objects_test.tmp')
  testthat::expect_true(is(counts(cds2), 'IterableMatrix'))
  testthat::expect_true(compare_matrix_control(get_matrix_info(mat=counts(cds1)), get_matrix_info(mat=counts(cds2)), compare_matrix_path_flag=FALSE))
  unlink('monocle_objects_test.tmp', recursive=TRUE)
} )

