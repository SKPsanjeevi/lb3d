#!/bin/bash

# This is an example script which takes already rendered images from 3
# different directories, adds a string for the timestep and combines all 3
# images into a single one.


DIR1="Ftot1.6_256_sh0fg0"
DIR2="Ftot1.6_256_sh0fg0.15"
DIR3="Ftot1.6_256_sh0fg0.3"
OUTDIR="movie_256sh0fgvar-new"

mkdir -p $OUTDIR

LSDIR1=`ls $DIR1/colo*jpg`

for i in $LSDIR1 ; do
TIMESTEP=`echo $i | sed s/_t0/_t/ | sed 's/.*t\([0-9][0-9][0-9][0-9][0-9]\).*/\1/'`
echo $TIMESTEP
echo $i
j=`ls  $DIR2/colour_Ftot*$TIMESTEP.vtk.jpg`
k=`ls  $DIR3/colour_Ftot*$TIMESTEP.vtk.jpg`
montage -geometry 300x -tile 3x1 $i $j $k $OUTDIR/$TIMESTEP.jpg
mogrify  -font helvetica -fill blue -pointsize 30 -draw "text 05,30 't=$TIMESTEP'" $OUTDIR/$TIMESTEP.jpg
done

