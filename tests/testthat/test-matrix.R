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
#     (OK) set_matrix_control_pca()
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
#     (NA) set_matrix_control_pca() # set matrix_control defaults: used only it projection.R (check)
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

  # Check that missing matrix_control[['matrix_class']] throws an error.
  testthat::expect_error(set_matrix_control(matrix_control=list(matrix_path='uhoh'), matrix_control_default=list()))
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
  # Convert dgCMatrix to dgCMatrix matrix.
  cds1 <- load_a549()
  testthat::expect_true(is(counts(cds1), 'dgCMatrix'))
  cds2 <- convert_counts_matrix(cds1, matrix_control=list(matrix_class='dgCMatrix'))
  testthat::expect_true(is(counts(cds2), 'dgCMatrix'))

  # Convert dgCMatrix to BPCells matrix.
  cds1 <- load_a549()
  testthat::expect_true(is(counts(cds1), 'dgCMatrix'))
  cds2 <- convert_counts_matrix(cds1, matrix_control=list(matrix_class='BPCells'))
  testthat::expect_true(is(counts(cds2), 'IterableMatrix'))

  # Convert BPCells to dgCMatrix matrix.
  cds1 <- load_a549(matrix_control=list(matrix_class='BPCells'))
  testthat::expect_true(is(counts(cds1), 'IterableMatrix'))
  cds2 <- convert_counts_matrix(cds1, matrix_control=list(matrix_class='dgCMatrix'))
  testthat::expect_true(is(counts(cds2), 'dgCMatrix'))

  # Convert BPCells to dgCMatrix matrix.
  cds1 <- load_a549(matrix_control=list(matrix_class='BPCells'))
  testthat::expect_true(is(counts(cds1), 'IterableMatrix'))
  cds2 <- convert_counts_matrix(cds1, matrix_control=list(matrix_class='BPCells'))
  testthat::expect_true(is(counts(cds2), 'IterableMatrix'))

  # Check that missing matrix_control[['matrix_class']] throws an error.
  testthat::expect_error(convert_counts_matrix(cds1, matrix_control=list(matrix_path='uhoh')))
} )


test_that("save_monocle_objects and load_monocle_objects", {
  cds1 <- load_a549(matrix_control=list(matrix_class='BPCells'))
  save_monocle_objects(cds1, directory_path='monocle_objects_test.tmp')
  cds2 <- load_monocle_objects(directory_path='monocle_objects_test.tmp')
  testthat::expect_true(is(counts(cds2), 'IterableMatrix'))
  testthat::expect_true(compare_matrix_control(get_matrix_info(mat=counts(cds1)), get_matrix_info(mat=counts(cds2)), compare_matrix_path_flag=FALSE))
  unlink('monocle_objects_test.tmp', recursive=TRUE)
} )


test_that("set_matrix_control_default", {
  # Set default matrix_class to dgCMatrix matrix.
  monocle3:::set_global_variable('matrix_class_default', 'dgCMatrix')

  # Set default matrix_control for dgCMatrix matrix.
  set_global_variable('matrix_control_csparsematrix_unrestricted', list(matrix_class='dgCMatrix'))
  set_global_variable('matrix_control_csparsematrix_pca', list(matrix_class='dgCMatrix'))

  # Load cds with default matrix.
  cds <- load_a549()
  testthat::expect_true(is(counts(cds), 'dgCMatrix'))

  # Load cds with requested dgCMatrix matrix.
  cds <- load_a549(matrix_control=list(matrix_class='dgCMatrix'))
  testthat::expect_true(is(counts(cds), 'dgCMatrix'))

  # Load cds with requested BPCells matrix.
  cds <- load_a549(matrix_control=list(matrix_class='BPCells'))
  testthat::expect_true(is(counts(cds), 'IterableMatrix'))


  # Set default matrix_class to BPCells matrix.
  monocle3:::set_global_variable('matrix_class_default', 'BPCells')

  # Set default matrix_control for BPCells matrix.
  monocle3:::set_global_variable('matrix_control_bpcells_unrestricted', list(matrix_class='BPCells', matrix_mode='dir', matrix_type='double', matrix_compress=FALSE, matrix_path='~/git/monocle3/bpcells_dir_tmp', matrix_buffer_size=8192L, matrix_bpcells_copy=TRUE))
  monocle3:::set_global_variable('matrix_control_bpcells_pca', list(matrix_class='BPCells', matrix_mode='dir', matrix_type='double', matrix_compress=FALSE, matrix_path='~/git/monocle3/bpcells_dir_tmp', matrix_buffer_size=8192L, matrix_bpcells_copy=TRUE))

  # Load cds with default matrix.
  cds <- load_a549()
  testthat::expect_true(is(counts(cds), 'IterableMatrix'))

  # Load cds with requested dgCMatrix matrix.
  cds <- load_a549(matrix_control=list(matrix_class='dgCMatrix'))
  testthat::expect_true(is(counts(cds), 'dgCMatrix'))

  # Load cds with requested BPCells matrix.
  cds <- load_a549(matrix_control=list(matrix_class='BPCells'))
  testthat::expect_true(is(counts(cds), 'IterableMatrix'))

  # Load cds with requested matrix_class and matrix_path.
  cds <- load_a549(matrix_control=list(matrix_class='BPCells', matrix_path='./bpcells_matrix_dir_test'))
  testthat::expect_true(length(grep('bpcells_matrix_dir_test', x=counts(cds)@dir, ignore.case=FALSE, perl=TRUE, value=TRUE)) > 0)
  unlink('bpcells_matrix_dir_test', recursive=TRUE)

  # Load cds with requested matrix_class and matrix_type='float'.
  cds <- load_a549(matrix_control=list(matrix_class='BPCells', matrix_type='float'))
  testthat::expect_true(counts(cds)@type == 'float')

  # Load cds with requested matrix_class and matrix_compress=TRUE.
  cds <- load_a549(matrix_control=list(matrix_class='BPCells', matrix_compress=TRUE))
  testthat::expect_true(counts(cds)@compressed == TRUE)

  # Check that missing matrix_control[['matrix_class']] throws an error.
  testthat::expect_error(load_a549(matrix_control=list(matrix_path='uhoh')))

  # Restore matrix_class_default to dgCMatrix for downstream
  # testing.
  monocle3:::set_global_variable('matrix_class_default', 'dgCMatrix')
} )


test_that("set_matrix_control_pca", {
  # Load and preprocess cds with dgCMatrix matrix.
  cds <- load_a549(matrix_control=list(matrix_class='dgCMatrix'))
  testthat::expect_message(preprocess_cds(cds, verbose=TRUE), regexp='pca: sparse_prcomp_irlba: matrix class: dgCMatrix')

  # Load and preprocess cds with BPCells matrix.
  cds <- load_a549(matrix_control=list(matrix_class='BPCells'))
  testthat::expect_message(preprocess_cds(cds, verbose=TRUE), regexp='pca: bpcells_prcomp_irlba: matrix class: TransformScaleShift')

  # Check that set_matrix_control_pca adjusts matrix_compress=TRUE but not matrix_type='float'.
  cds <- load_a549(matrix_control=list(matrix_class='BPCells', matrix_type='float', matrix_compress=TRUE))
  testthat::expect_message(preprocess_cds(cds, verbose=TRUE), 'type:        float')
  testthat::expect_message(preprocess_cds(cds, verbose=TRUE), 'compress:    FALSE')

  # Check preprocess_transform use of set_matrix_control_pca with dgCMatrix cds_qry.
  cds_ref <- load_a549()
  cds_qry <- load_a549()

  cds_ref <- preprocess_cds(cds_ref)
  save_transform_models(cds=cds_ref, directory_path='monocle_transform_models')
  cds_qry <- load_transform_models(cds_qry, directory_path='monocle_transform_models')
  testthat::expect_message(preprocess_transform(cds_qry, verbose=TRUE), 'projection: sparse_apply_transform: matrix class: dgCMatrix')

  # Check preprocess_transform use of set_matrix_control_pca with BPCells cds_qry.
  cds_qry <- load_a549(matrix_control=list(matrix_class='BPCells'))
  cds_qry <- load_transform_models(cds_qry, directory_path='monocle_transform_models')
  testthat::expect_message(preprocess_transform(cds_qry, verbose=TRUE), 'projection: bpcells_apply_transform: matrix class: TransformScaleShift')
} )


test_that("set_matrix_control_combine_cds", {
  cds_bpc <- monocle3:::load_worm_embryo(matrix_control=list(matrix_class='dgCMatrix'))

  # Combine two dgCMatrix matrix cdses.
  cds_bpc1 <- cds_bpc[,1:1000]
  cds_bpc2 <- cds_bpc[,1001:6188]  
  cds_combined <- combine_cds(list(cds_bpc1, cds_bpc2), cell_names_unique=TRUE)
  testthat::expect_true(is(counts(cds_combined), 'dgCMatrix'))

  # Combine one dgCMatrix and one BPCells  matrix cdses.
  cds_bpc1 <- convert_counts_matrix(cds_bpc1, matrix_control=list(matrix_class='BPCells'))
  cds_combined <- combine_cds(list(cds_bpc1, cds_bpc2), cell_names_unique=TRUE)
  testthat::expect_true(is(counts(cds_combined), 'IterableMatrix'))

  # Combine two BPCells matrix cdses.
  cds_bpc2 <- convert_counts_matrix(cds_bpc2, matrix_control=list(matrix_class='BPCells'))
  cds_combined <- combine_cds(list(cds_bpc1, cds_bpc2), cell_names_unique=TRUE)
  testthat::expect_true(is(counts(cds_combined), 'IterableMatrix'))

  # Check that missing matrix_control[['matrix_class']] throws an error.
  testthat::expect_error(combine_cds(list(cds_bpc1, cds_bpc2), cell_names_unique=TRUE, matrix_control=list(matrix_path='uhoh')))

} )

