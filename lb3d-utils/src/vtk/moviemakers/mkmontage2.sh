#!/bin/bash

# This is an example script which takes already rendered images from 4
# different directories, adds a string for the timestep and combines all 4
# images into a single one. Output is written with and without timestamp.

# Jens, 01.03.07

# The original version was used for VTK output of 850x330.

DIR1="Ftot1.6_128cub_sh0.1fg0"
DIR2="Ftot1.6_128cub_sh0.1fg0.1"
DIR3="Ftot1.6_128cub_sh0.1fg0.2"
DIR4="Ftot1.6_128cub_sh0.1fg0.3"
OUTDIR="movie_128cub_sh0.1fgvar-2x2-tilted"

mkdir -p $OUTDIR
mkdir -p $OUTDIR-notimestamp

LSDIR1=`ls $DIR1/colo*jpg`

for i in $LSDIR1 ; do
TIMESTEP=`echo $i | sed 's/.*t\([0-9][0-9][0-9][0-9][0-9][0-9]\).*/\1/'`
echo $TIMESTEP
echo $i
j=`ls  $DIR2/colour_Ftot*$TIMESTEP.vtk.jpg`
k=`ls  $DIR3/colour_Ftot*$TIMESTEP.vtk.jpg`
l=`ls  $DIR4/colour_Ftot*$TIMESTEP.vtk.jpg`
montage -geometry 400x -tile 2x2 $i $j $k $l $OUTDIR/$TIMESTEP.jpg
cp $OUTDIR/$TIMESTEP.jpg $OUTDIR-notimestamp/$TIMESTEP.jpg 
mogrify  -font helvetica -fill blue -pointsize 20 -draw "text 05,20 't=$TIMESTEP'" $OUTDIR/$TIMESTEP.jpg
done

# mogrify -geometry 352x -border 0x100 -bordercolor white *jpg
