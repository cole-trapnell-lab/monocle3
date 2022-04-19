== WWW site files

  The document root is the top directory in the Monocle3 gh-pages directory and
  the initial page is the conventional index.html.
  Most html files are in the gh-pages/docs directory.
  Most image files are in gh-pages/images/manual_images.
  Most table csv files are in _data/de_tests.

== image and table build scripts

  R scripts that generate most of the image PNG files and table CSV files are in the directory
  gh-pages/\~off_site/build_scripts. Run them as follows

    # Edit gh-pages/\~off_site/build_scripts/globals.R to set destination directories for the image and table files

    cd gh-pages/\~off_site/build_scripts
    R
    > devtools::load_all('<monocle3 directory>')
    > source('<build script>', local=TRUE, echo=TRUE)  # <build script> is clustering.R, differential.R, projection.R, or trajectories.R.

    # When you are satisfied with the files, move them to the images and _data directories.

  The cluster.html page has a section called 'Annotate your cells according to type', which
  includes two tables that assign cell type names to clusters. The cell type names are based
  on Jonathan Packers re-analysis of Jun Cao's cell type assignments. The L2 worm cell metadata
  downloaded from http://staff.washington.edu/hpliner/data has only Jun Cao's cell type
  assignments. I found a file of Jonathan's cell type assignments for L2 worm and I transferred
  those labels to Jun's data set and created a file that has Jonathan's cell type assignment
  for each of Jun's cells. This work is in gh_pages/\~off_site/celegans_l2. The file NOTES
  describes how I transferred the labels, and how to use the file with Jonathan's cell types
  to label clusters.


== retired site files

  Old site files are moved to gh-pages/\~off_site/retired.gh_pages/<n>
  where <n> is a date.

== additional documentation files

  o  Jekyll setup and running
  o  building GIF animation files
  o  add html pages
  o  additional site information
  o  header information in R build scripts
