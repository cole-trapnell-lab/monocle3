#
# Set up a global-variable-like environment.
#

set_global_variable <- function(variable_name, value) {
  assign(variable_name, value, envir=._._global_variable_env_._.)
}


get_global_variable <- function(variable_name) {
  value <- tryCatch({
                      get('guard_element', envir=._._global_variable_env_._.) 
                    }, error=function(msg) {
                      message('Global variable storage is compromised.')
                      return(NA)
                    })

  value <- tryCatch({
                      get(variable_name, envir=._._global_variable_env_._.)
                    }, error=function(msg) {
                      message(paste0("\'variable_name\'", ' is not a global variable.'))
                      return(NA)
                    })
  return(value)
}


# The ..global_variable_env.. stores globally accessible
# objects accessible using the set_global_variable and
# get_global_variable functions.
._._global_variable_env_._. <- new.env()

# Define some global variables.
.onLoad <- function(libname, pkgname) {
  set_global_variable('guard_element', 'sanity_check')
  set_global_variable('id_count', 1)
  set_global_variable('reduce_dim_pca_model_version', 1)
  set_global_variable('reduce_dim_lsi_model_version', 1)
  set_global_variable('reduce_dim_aligned_model_version', 1)
  set_global_variable('reduce_dim_tsne_model_version', 1)
  set_global_variable('reduce_dim_umap_model_version', 1)
  set_global_variable('monocle_objects_version', 1)
  set_global_variable('transform_models_version', 1)

  # Watching preprocess_cds() it appears that R uses OMP_NUM_THREADS
  # threads if OMP_NUM_THREADS > 1 and OPENBLAS_NUM_THREADS is NA.
  # After setting, RhpcBLASctl::blas_set_num_threads(1L) and then
  # RhpcBLASctl::blas_set_num_threads(NA), preprocess_cds() uses
  # one thread. So I assume that the number of desired threads for
  # both OMP_NUM_THREADS and OPENBLAS_NUM_THREADS is the greater
  # of the two.
  omp_num_threads <- as.numeric(Sys.getenv('OMP_NUM_THREADS'))
  blas_num_threads <- as.numeric(Sys.getenv('OPENBLAS_NUM_THREADS'))

  if(is.na(omp_num_threads))
    omp_num_threads <- 1L
  if(is.na(blas_num_threads))
    blas_num_threads <- 1L

  if(omp_num_threads < blas_num_threads)
    omp_num_threads <- blas_num_threads
  else
    blas_num_threads <- omp_num_threads
  
  set_global_variable('omp_num_threads', omp_num_threads)
  set_global_variable('blas_num_threads', blas_num_threads)
}

