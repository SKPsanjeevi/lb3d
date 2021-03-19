#!/bin/bash
# Small example script that shows how to get 1D plots out of a lot of data
# sets in
# different directories. It goes through all the directories given at the
# command line and converts some of the files.
# 
# THIS SCRIPT IS NOT MEANT TO BE USED AS IT IS. IT HAS TO BE ADAPTED.
# Jens, 30.10.02

CUT1D=../../utils/ftools/cut1d_hdf
FILES="od*h5 vel*h5" # Which files do you want to 'cut' ?
FTYPE=2          # 1=double, 2=float
TYPE=1           # 1=scalar, 2=2scalar
VDIR="X Y Z"     # Variable directions
AB="8,8"       # Fixed points in other directions

for j in $* ; do
 echo $j
 cd $j
  for k in $VDIR ; do
   for i in `ls $FILES` ; do
      echo $i >/tmp/cut1d.$$
      echo "$FTYPE" >>/tmp/cut1d.$$
      echo "$TYPE" >>/tmp/cut1d.$$
      echo "$k" >>/tmp/cut1d.$$
      echo "$AB" >>/tmp/cut1d.$$
      cat /tmp/cut1d.$$| $CUT1D  >>/dev/null
      rm -f /tmp/cut1d.$$
   done
  done
 cd ..
done
