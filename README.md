MONOCLE 3
=======================

Monocle 3 is an analysis toolkit for single-cell RNA-Seq experiments.  To use this package, you will need the R statistical computing environment (version 3.0 or later) and several packages available through Bioconductor and CRAN.

Details on how to install and use Monocle 3 are available on our website:

http://cole-trapnell-lab.github.io/monocle3/

## Monocle3 with BPCells count matrix support

This development branch version of Monocle3 adds the ability to store the count matrix on-disk using the BPCells package. By default, Monocle3 stores the counts matrix in memory as a sparse matrix, as in previous versions. In order to store the matrix on-disk, you must set the matrix_control list value 'matrix_class="BPCells"' in the affected commands. For example, to load a MatrixMarket file as an on-disk matrix, use the command

```
cds <- load_mm_data(mat_path=<path_to_mtx_file>,
                    feature_anno_path=<path_to_feature_anno_file>,
                    cell_anno_path=<path_to_cell_anno_file>,
                    matrix_control=list(matrix_class='BPCells'))
```

### Install Monocle3 with BPCells

You must install BPCells from Github before you can install the Monocle3 version, and BPCells requires an HDF5 object library for installation. After installing the HDF5 library, you can install BPCells using the command

```
remotes::install_github("bnprks/BPCells")
```

The [BPCells Github site](https://github.com/bnprks/BPCells)  has additional information.

Some Linux distributions provide the HDF5 library as an option. The
BPCells site has information about installing the HDF5 library on various operating systems.

I used Homebrew to install an HDF5 library on MacOS. I seemed to need to install the pkg-config package as well, and add a pkg-config configuration file for HDF5. Homebrew installed pkg-config in '/opt/homebrew' so I added the hdf5.pc file in

> /opt/homebrew/lib/pkgconfig/hdf5.pc

with the contents

```
prefix=/opt/homebrew/Cellar/hdf5/1.12.2_2
exec_prefix=${prefix}
includedir=${prefix}/include
libdir=/opt/homebrew/Cellar/hdf5/1.12.2_2/lib
  
Name: hdf5
Description: HDF5
URL: xx
Version: 1.12.2_2
Cflags: -I${includedir}
Libs: -L${libdir} -lhdf5
```

You may need to update the version strings in your hdf5.pc file.

### Notes

- The function `convert_counts_matrix()` converts the counts matrix in an existing CDS. For example, `cds <- convert_counts_matrix(cds, matrix_control=list(matrix_class="BPCells"))` converts an in-memory dgCMatrix that is stored in the cds to an on-disk BPCells matrix.
- The functions `load_mm_data()`, `load_mtx_data()`, `load_a549()`, `load_worm_embryo()`, and `load_worm_l2()` store the count matrix as an on-disk BPCells matrix when the matrix_control parameter is set as `matrix_control=list(matrix_class="BPCells")`.
- The function `new_cell_data_set()` accepts a BPCells on-disk count matrix.
- The functions `save_monocle_objects()` and `load_monocle_objects()` store and load BPCells on-disk matrices, when the cds count matrix is an on-disk BPCells matrix.
- The function `combine_cds()` combines cdses with mixes of dgCMatrix and BPCells on-disk count matrices into a BPCells on-disk count matrix. When called with the `matrix_control=list(matrix_class="BPCells")` parameter, `combine_cds()` combines cdses with all dgCMatrix count matrices into a BPCells on-disk count matrix.
- By default, the BPCells on-disk matrix is stored in a directory that is created where R is started. You can change the directory location using the `matrix_path` value in the `matrix_control` parameter. For more information, see the help document for the function `set_matrix_control()`.
- Note that when the count matrix is stored as a BPCells on-disk matrix, the `new_cell_data_set()` function stores a second BPCells on-disk copy of the matrix in the cds assays slot with the name `counts_row_order`. The `counts_row_order` matrix is used by Monocle3 when count matrices are accessed intensively by row. The reason is that R stores matrices as a one-dimensional vector in column-major order, and BPCells does as well by default. As a result, column access is fast and row access is slow. We use BPCell's ability to store and access matrices in row-major order, which allows fast row access. However, this means that the two copies of the count matrices must be kept up to date. If you replace the cds's count matrix using the `counts(cds)` setter method, you must update the `counts_row_order` matrix, which you can do using the function `set_cds_row_order_matrix()`.
- I tested this version using BPCells counts matrices on the examples in the Monocle3 documentation although I did not try all of the plotting functions.
- The BPCells matrix in-memory storage is not supported by Monocle3; that is, `matrix_control=list(matrix_class="BPCells", matrix_mode="mem")`. Do not use it. The default is `matrix_mode="dir"`, which is on-disk storage.

