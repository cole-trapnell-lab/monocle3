#
# Set global options.
#
options("sp_evolution_status"=2)

#
# Set up a global-variable-like environment.
#

set_global_variable <- function(variable_name, value) {
  assign(variable_name, value, envir=._._global_variable_env_._.)
}


# Return value of variable_name. If variable_name is NULL, return a list
# of all global variables.
get_global_variable <- function(variable_name=NULL) {
  value <- tryCatch({
                      v <- get('guard_element', envir=._._global_variable_env_._.) 
                      if(v != 'sanity_check') stop()
                      v
                    }, error=function(msg) {
                      message('Global variable storage is compromised.')
                      return(NA)
                    })

  if(is.na(value)) {
    return(NA)
  }

  if(!is.null(variable_name)) {
    value <- tryCatch({
                        get(variable_name, envir=._._global_variable_env_._.)
                      }, error=function(msg) {
                        message("\'", variable_name, "\'", ' is not a global variable.')
                        return(NA)
                      })
  }
  else {
    value=list()
    variable_names <- ls(envir=._._global_variable_env_._.)
    for(variable_name in variable_names) {
      value[[variable_name]] <- get(variable_name, envir=._._global_variable_env_._.)
    }
  }

  return(value)
}


# The ._._global_variable_env_._. environment stores global
# objects. Use the set_global_variable and get_global_variable
# functions to access them.
._._global_variable_env_._. <- new.env(parent=emptyenv())

# Define some global variables.
.onLoad <- function(libname, pkgname) {
  # A value used to ensure that this is the Monocle3
  # global variables.
  set_global_variable('guard_element', 'sanity_check')

  # Counter used to ensure that matrix time stamps generated
  # by the get_unique_id() function are distinct. This value
  # is incremented each time that get_unique_id() is called
  # during an R session.
  set_global_variable('id_count', 1)

  # Object version numbers. These are used by the Monocle3
  # to recognize objects so we suggest that you not change
  # them.
  set_global_variable('reduce_dim_pca_model_version', 1)
  set_global_variable('reduce_dim_lsi_model_version', 1)
  set_global_variable('reduce_dim_aligned_model_version', 1)
  set_global_variable('reduce_dim_tsne_model_version', 1)
  set_global_variable('reduce_dim_umap_model_version', 1)
  set_global_variable('monocle_objects_version', 1)
  set_global_variable('transform_models_version', 1)
  set_global_variable('monocle3_annoy_index_version', 2)
  set_global_variable('monocle3_hnsw_index_version', 1)
  set_global_variable('monocle3_timer_t0', 0)
  set_global_variable('monocle3_timer_msg', "")
  set_global_variable('monocle_gc_matrix_path', list())

  # Default nn_control list for functions that do not need
  # an index, which is all but the label transfer functions.
  set_global_variable('nn_control_annoy_euclidean', list(method='annoy', metric='euclidean', n_trees=50, M=48, ef_construction=200, ef=150, grain_size=1, cores=1))

  # Default nn_control list for functions that need an index,
  # which are the label transfer functions.
  set_global_variable('nn_control_annoy_cosine', list(method='annoy', metric='cosine', n_trees=50, M=48, ef_construction=200, ef=150, grain_size=1, cores=1))

  # Default matrix_control list for any.
  set_global_variable('matrix_control_csparsematrix_unrestricted', list(matrix_class='dgCMatrix'))
  set_global_variable('matrix_control_bpcells_unrestricted', list(matrix_class='BPCells', matrix_mode='dir', matrix_type='double', matrix_compress=FALSE, matrix_path='.', matrix_buffer_size=8192L, matrix_bpcells_copy=TRUE))

  # Default matrix_control list for pca.
   set_global_variable('matrix_control_csparsematrix_pca', list(matrix_class='dgCMatrix'))
   set_global_variable('matrix_control_bpcells_pca', list(matrix_class='BPCells', matrix_mode='dir', matrix_type='double', matrix_compress=FALSE, matrix_path='.', matrix_buffer_size=8192L, matrix_bpcells_copy=TRUE))


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

  # for travis
  Sys.setenv('TESTTHAT_MAX_FAILS' = Inf)
}

#
# Try to clean up any temporary matrix files and directories on exiting.
#
._._gc_matrix_object_remove_._. <- function(env) {
  matrix_path_list <- get_global_variable('monocle_gc_matrix_path')
  for(matrix_path in matrix_path_list) {
    if(file.exists(matrix_path) || dir.exists(matrix_path)) {
      unlink(matrix_path, recursive=TRUE)
    }
  }
}

reg.finalizer(._._global_variable_env_._., ._._gc_matrix_object_remove_._., onexit=TRUE)

