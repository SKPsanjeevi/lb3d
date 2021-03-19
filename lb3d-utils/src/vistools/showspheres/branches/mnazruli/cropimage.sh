#!/bin/bash
mkdir -p cropped/small

FILES=`ls Streamline.*.png`

for i in $FILES ; do

TIMESTEP=`echo $i | sed 's/.*.\([0-9][0-9][0-9][0-9][0-9][0-9]\).*/\1/'`
echo $TIMESTEP
convert -crop 1867x1400+510+600 -font helvetica -fill blue -pointsize 100\
  -draw "text 1800,750 't=$TIMESTEP'" $i cropped/$TIMESTEP.jpg
convert -geometry 768x576 cropped/$TIMESTEP.jpg cropped/small/$TIMESTEP.small.jpg
done
