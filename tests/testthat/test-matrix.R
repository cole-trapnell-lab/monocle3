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
#   check_matrix_control <- function(matrix_control=list(), control_type=c('unrestricted', 'counts', 'mm', 'pca'), check_conditional=FALSE)
#   set_matrix_control <- function(matrix_control=list(), matrix_control_default=list(), control_type=c('unrestricted', 'counts', 'mm', 'pca'))
#   show_matrix_control <- function(matrix_control, label=NULL)
#   push_matrix_path <- function(mat)
#   bpcells_find_base_matrix <- function(mat)jA
#   get_matrix_class <- function(mat)
#   get_matrix_info <- function(mat)
#   show_matrix_info <- function(matrix_info, indent='')
#   set_matrix_class <- function(mat, matrix_control=list())
#   rm_bpcells_dir <- function(mat)
#   set_cds_row_order_matrix <- function(cds)
#   convert_counts_matrix <- function(cds, matrix_control=list(matrix_class='BPCells'))
#

# functions that need additional review
#   check_matrix_control <- function(matrix_control=list(), control_type=c('unrestricted', 'counts', 'mm', 'pca'), check_conditional=FALSE)
#   set_matrix_control <- function(matrix_control=list(), matrix_control_default=list(), control_type=c('unrestricted', 'counts', 'mm', 'pca'))
#   set_matrix_class <- function(mat, matrix_control=list())
#   convert_counts_matrix <- function(cds, matrix_control=list(matrix_class='BPCells'))
#   (done) get_matrix_info <- function(mat)
#   (done) show_matrix_info <- function(matrix_info, indent='')
#   (done) rm_bpcells_dir <- function(mat)
#   (done) set_cds_row_order_matrix <- function(cds)

# (exported) files that have functions that have a matrix_control parameter. Check that they have a call to set_matrix_control/check_matrix_control
# grep 'matrix_control=lis' *.R
#   io.R
#   load_cellranger_data.R
#   matrix.R
#   pca.R
#   projection.R
#
# grep 'matrix_control=' *.R
#   io.R
#   load_cellranger_data.R
#   matrix.R
#   pca.R
#   projection.R
#   utils.R

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


