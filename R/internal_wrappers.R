bpcells_available <- function() {
  requireNamespace("BPCells", quietly = TRUE)
}

require_bpcells <- function(fn_name = NULL) {
  if(!bpcells_available()) {
    if(is.null(fn_name) || !nzchar(fn_name)) {
      msg <- "BPCells must be installed to use this functionality."
    }
    else {
      msg <- paste0("BPCells must be installed to use ", fn_name, ".")
    }
    stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}

is_iterable_matrix <- function(x) {
  bpcells_available() && methods::is(x, "IterableMatrix")
}

bpcells_linear_operator <- function(x) {
  require_bpcells("bpcells_linear_operator")
  getFromNamespace("linear_operator", "BPCells")(x)
}

bpcells_matrix_inputs <- function(x) {
  require_bpcells("bpcells_matrix_inputs")
  getFromNamespace("matrix_inputs", "BPCells")(x)
}

delayedarray_set_verbose_block_processing <- function(value) {
  getFromNamespace("set_verbose_block_processing", "DelayedArray")(value)
}
