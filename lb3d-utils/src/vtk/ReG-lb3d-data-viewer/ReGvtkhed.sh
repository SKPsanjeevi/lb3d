#!/bin/bash
# this script calls the viewer script in order to generate a hedgehog
# image of the surfactant dirctions. The script takes a desnity and a
# directions file as arguments.
 
#Make sure the LBEINST/VIEWERPATH and LD_LIBRARY_PATH variables are set
#properly. If you set LBEINST in your shell profile, ReGvtkhed.sh should
#find what it needs. 
#Original version by A. Porter, changes by Jens.


if [ "$LBEINST" == "" ] ; then
  LBEINST=/home/jens/codes/lbe/branches/version5/lbe
fi
export PATH=/usr/local/VTK/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/VTK/lib:/usr/local/lib/vtk/$LD_LIBRARY_PATH
VIEWERPATH=$LBEINST/utils/vtk/ReG-lb3d-data-viewer
vtk $VIEWERPATH/viewer_component.tcl $VIEWERPATH/viewer_component_ui.tcl 1 hh $1 $2

