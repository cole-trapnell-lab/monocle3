
# monocle3 0.2.0

### Major changes
* Added mutual-nearest-neighbor batch correction (MNNCorrect).
* Switched to leiden-based clustering, dropped reticulate/python dependency.
* Added a mechanism to get the citations used during an analysis get_citations().

### Other Changes
* Added non-standard color options for plot_cell_3d.
* Added norm_method = 'none' option for importing pre-normalized data.

### Bug fixes
* Fixed a bug that effected cell size in plot_cells_3d.
* Added a check for illegal characters in generate_garnett_marker_file.
* Fixed a bug in the alpha parameter in plot_cells.
* Fixed multiple minor reported bugs.

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
