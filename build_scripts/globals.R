#
# Global definitions used by *.R files that test and generate png
# files for Monocle3 documentation. This file is sourced by the
# test/generate files.
#

monocle3_git_dir <- '/home/brent/git/monocle3'
monocle3_git_gh_dir <- '/home/brent/git_gh_pages/monocle3'
images_new_dir <- file.path(monocle3_git_gh_dir, 'images.new')
manual_images_dir <- file.path(images_new_dir, 'manual_images')
data_new_dir <- file.path(monocle3_git_gh_dir, '_data_new')
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

check_images_new_dir(images_new_dir)
check_data_new_dir(data_new_dir)
