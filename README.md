MONOCLE 3
=======================

Monocle 3 is an analysis toolkit for single-cell RNA-Seq experiments.  To use this package, you will need the R statistical computing environment (version 3.0 or later) and several packages available through Bioconductor and CRAN.

Details on how to install and use Monocle 3 are available on our website:

http://cole-trapnell-lab.github.io/monocle3/

## Monocle3 with BPCells counts matrix support

This development branch version of Monocle3 adds the ability to store the counts matrix on-disk using the BPCells package. By default, Monocle3 stores the counts matrix in-memory as a sparse matrix, as in previous versions. In order to store the matrix on-disk, you must set the matrix_control list value `matrix_class="BPCells"` in the affected commands. For example, to load a MatrixMarket file as an on-disk matrix, use the command

```
cds <- load_mm_data(mat_path=<path_to_mtx_file>,
                    feature_anno_path=<path_to_feature_anno_file>,
                    cell_anno_path=<path_to_cell_anno_file>,
                    matrix_control=list(matrix_class='BPCells'))
```

### Install Monocle3 with BPCells

You must install BPCells from Github before you can install this Monocle3 version, and BPCells requires an HDF5 object library for installation. After installing the HDF5 library, you install BPCells using the command

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

- Monocle3 can use the BPCells package to store the feature-cell counts matrix on-disk rather than in-memory, which enables analysis of considerably larger data sets than before. By default, Monocle3 stores the counts matrix in-memory as a sparse matrix, as it has in the past. To store the counts matrix on-disk, use the parameter `matrix_control=list(matrix_class="BPCells")` when you make the CDS or convert the counts matrix using one of the functions
  - `load_mm_data()`
  - `load_mtx_data()`
  - `load_cellranger_data()`
  - `load_a549()`
  - `load_worm_embryo()`
  - `load_worm_l2()`
  - `convert_counts_matrix()`

  For example, to convert a dgCMatrix counts matrix to a BPCells on-disk matrix in an existing CDS, use the command `cds <- convert_counts_matrix(cds, matrix_control=list(matrix_class="BPCells"))`.
- The method `new_cell_data_set()` accepts a BPCells on-disk counts matrix.
- The functions `save_monocle_objects()` and `load_monocle_objects()` store and load BPCells on-disk matrices when the CDS counts matrix is an on-disk BPCells matrix.
- The Monocle3 `saveRDS()` function warns the user to use `save_monocle_objects()` when saving a CDS with a BPCells on-disk counts matrix. If you insist on using the `saveRDS()` function, the BPCells on-disk matrix directory will not be stored and you will be unable to load it with the `readRDS()` function.
- The function `combine_cds()` combines CDSes with mixes of dgCMatrix and BPCells on-disk counts matrices into a BPCells on-disk counts matrix. When called with the `matrix_control=list(matrix_class="BPCells")` parameter, `combine_cds()` combines CDSes with all dgCMatrix counts matrices into a BPCells on-disk counts matrix.
- Note that when the counts matrix is stored as a BPCells on-disk matrix, the `new_cell_data_set()` method stores a second BPCells on-disk copy of the matrix in the CDS assays slot with the name `counts_row_order`. The `counts_row_order` matrix is used by Monocle3 when the counts matrix is accessed intensively by row. The reason is that, by default, BPCells stores the matrix as a one-dimensional vector in column-major order, as does R. As a result, column access is fast and row access is slow. We use BPCell's ability to also store and access matrices in row-major order, which gives fast row access. However, this means that the two copies of the counts matrix must have the same count values. If you replace or change the CDS's counts matrix, you must also update the `counts_row_order` matrix, which you can do using the function `set_cds_row_order_matrix()`.
  - The CDS assays slot is a named list where the standard, column-major order, matrix is called `counts` and the BPCells row-major order matrix is called `counts_row_order`.
  - The `counts` matrix getter and setter methods are `counts(cds)` and `counts(cds)<-`.
  - The Monocle3 setter warns about re-setting the BPCells `counts_row_order` matrix, unless called with the parameter `bpcells_warn=FALSE`.
  - The `counts_row_order` setter method is called `counts_row_order`.
  - There is no corresponding `counts_row_order` setter method
- By default, the BPCells on-disk matrix is stored in a directory that is created where R is started. You can change the directory location using the `matrix_path` value in the `matrix_control` parameter.
- For more information about the `matrix_control` values, see the help document for the function `set_matrix_control()`.
- I tested this version using BPCells counts matrices on the examples in the Monocle3 documentation although I did not try all of the plotting functions.

