#!/bin/bash

# Copy image files that are expected to not change.
# Note that the images files will be copied to the
# directories ../images.new and ../images.new/manual_images.

# Copy images for initial page.
if [ ! -d ../images.new ]
then
  echo 'Images directory does not exist...creating it now.'
  mkdir ../images.new
fi

lfil='worm_trajectory.png worm_umap.png worm_violin.png worm_classification.png top_markers.png pseudotime_dynamics.png branch_module.png choose_graph_segments.gif favicon.png monocle3_new_workflow.png'

for fil in $lfil
do
  echo "copy $fil"
  cp ../images_warehouse/$fil ../images.new
done


if [ ! -d ../images.new/manual_images ]
then
  echo 'Images.new/manual_images directory does not exist...creating it now.'
  mkdir ../images.new/manual_images
fi

lfil='manual_images/choose_cells_recording.gif manual_images/choose_root_recording.gif manual_images/choose_branch_recording.gif'

for fil in $lfil
do
  echo "copy $fil"
  cp ../images_warehouse/$fil ../images.new/manual_images
done

