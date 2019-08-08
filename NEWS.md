# monocle3 0.1.3

### Changes
* Added bootstrap bars to plot_percent_cells_positive().
* Added a UMI cutoff for load_cellranger_data().
* Added clear_cds_slots() to clear slots from a cds.
* Added cell_stroke parameter to plot_cells.

### Bug fixes
* Fixed #154, which actually resulted from a problem in how ncenter was calculated for each partition.
* Fixed an issue that prevents plotting cells by pseudotime when root_cells are passed to order_cells().
* Fixed #167 - NAs from plot_genes_by_group.

# monocle3 0.1.2

### Changes
* Added a `NEWS.md` file to track changes to the package.
* Added a combine_cds function that will combine cell_data_sets .
* Added a message when duplicate markers are present in generate_garnett_marker_file and option to exclude them.

### Bug fixes
* Fixed a bug that made expression matrix the wrong class with load_cellranger.
