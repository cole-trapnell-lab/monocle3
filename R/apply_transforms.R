sparse_apply_transform <- function(FM, rotation_matrix, vcenter=vcenter, vscale=vscale, block_size=NULL, cores=1, verbose=FALSE) {
  if(verbose) {
    message('sparse_apply_transform: start')
  }

  if(!is.null(block_size)) {
    block_size0 <- DelayedArray::getAutoBlockSize()
    DelayedArray::setAutoBlockSize(block_size)
  }

  # Thank you Maddy.
  intersect_genes <- intersect(rownames(rotation_matrix), rownames(FM))
  intersect_indices <- match(intersect_genes, rownames(rotation_matrix))

  if(any(is.na(intersect_indices))) {
    stop('gene sets differ: genes in the subject matrix are not in the reference matrix')
  }

  if((length(intersect_indices)/length(rownames(FM))) < 0.5) {
    warning('fewer than half the genes in the subject matrix are also in the reference matrix: are the matrices prepared using the same gene set?')
  } 

  # [intersect_genes,] orders FM rows by intersect_genes
  FM <- FM[intersect_genes,]
  
  xt <- Matrix::t(FM)
  xtda <- DelayedArray::DelayedArray(xt)

  vcenter <- vcenter[intersect_indices]
  vscale <- vscale[intersect_indices]

  xtdasc <- t(xtda) - vcenter
  xtdasc <- t(xtdasc / vscale)

  irlba_res <- list()
#  irlba_res$x <- xtdasc %*% rotation_matrix[intersect_indices,]
  irlba_res$x <- matrix_multiply_multicore(mat_a=xtdasc,
                                         mat_b=rotation_matrix[intersect_indices,],
                                         cores)
  irlba_res$x <- as.matrix(irlba_res$x)
  class(irlba_res) <- c('irlba_prcomp', 'prcomp')

  if(verbose) {
    message('sparse_apply_transform: finish')
  }

  return(irlba_res)
}


bpcells_apply_transform <- function(FM, rotation_matrix, vcenter=vcenter, vscale=vscale, pca_control=list(), verbose=FALSE) {
  if(verbose) {
    message('bpcells_apply_transform: start')
  }

  # Thank you Maddy.
  intersect_genes <- intersect(rownames(rotation_matrix), rownames(FM))
  intersect_indices <- match(intersect_genes, rownames(rotation_matrix))

  if(any(is.na(intersect_indices))) {
    stop('gene sets differ: genes in the subject matrix are not in the reference matrix')
  }

  if((length(intersect_indices)/length(rownames(FM))) < 0.5) {
    warning('fewer than half the genes in the subject matrix are also in the reference matrix: are the matrices prepared using the same gene set?')
  }

  # bge
  # [intersect_genes,] orders FM rows by intersect_genes 
  # This almost certainly requires rewriting the matrix.
  FM <- FM[intersect_genes,]
 
  xt <- BPCells::t(FM)

  vcenter <- vcenter[intersect_indices]
  vscale <- vscale[intersect_indices]

  xtsc <- BPCells::t(xt) - vcenter
  xtsc <- BPCells::t(xtsc / vscale)

  # make intermediate matrix.

  irlba_res <- list()
  irlba_res$x <- xtsc %*% rotation_matrix[intersect_indices,]

  irlba_res$x <- as.matrix(irlba_res$x)
  class(irlba_res) <- c('irlba_prcomp', 'prcomp')

  if(verbose) {
    message('bpcells_apply_transform: finish')
  }

  return(irlba_res)
}


