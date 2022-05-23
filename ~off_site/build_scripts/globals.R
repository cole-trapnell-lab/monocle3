#
# Global definitions used by *.R files that test and generate png
# files for Monocle3 documentation. This file is sourced by the
# test/generate files.
#

# Monocle3 is loaded from here.
monocle3_git_dir <- '/home/brent/git/monocle3'

# Images and tables are stored within this root directory.
monocle3_git_gh_dir <- '/home/brent/git_gh_pages/monocle3/\~off_site'

# Images are store here.
images_new_dir <- file.path(monocle3_git_gh_dir, 'images.new')
manual_images_dir <- file.path(images_new_dir, 'manual_images')

# CSV file tables are stored here.
data_new_dir <- file.path(monocle3_git_gh_dir, '_data,new')
data_de_dir <- file.path(data_new_dir, 'de_tests')


check_images_new_dir <- function(images_new_dir) {
  if(!file.exists(images_new_dir)) {
    dir.create(images_new_dir)
  }
  if(!file.exists(manual_images_dir)) {
    dir.create(manual_images_dir)
  }
}

check_data_new_dir <- function(data_new_dir) {
  if(!file.exists(data_new_dir)) {
    dir.create(data_new_dir)
  }
  if(!file.exists(data_de_dir)) {
    dir.create(data_de_dir)
  }
}


ggplot_cells_png <- function( plot_cmd, plot_extra=theme(text=element_text(size=6)), manual_images_dir, plot_file_name) {
  message('plot ', plot_file_name)
  if(plot_extra != '') {
    plot_cmd + plot_extra
  }
  else
  {
    plot_cmd
  }
  ggsave(file.path(manual_images_dir, plot_file_name), units='in', width=5, height=4, dpi=600, device='png')
}


check_images_new_dir(images_new_dir)
check_data_new_dir(data_new_dir)

