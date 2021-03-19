#!/bin/bash
# Small example script that shows how to convert a lot of datasets in
# different directories. It goes through all the directories given at the
# command line and converts some of the files. You have to set the lattice
# size and data type.
# THIS SCRIPT IS NOT MEANT TO BE USED AS IT IS. IT HAS TO BE ADAPTED.
# Jens, 09.10.02

NX=64
NY=64
NZ=256
DBL2SGL=/volta/scratch/jens/lbeV4/utils/ftools/rotxdr
for j in $* ; do
 echo $j
 cd $j
  for i in `ls colour*.xdr` ; do
   echo $i >/tmp/dblsgl.$$
   echo "$NX,$NY,$NZ" >>/tmp/dblsgl.$$
   cat /tmp/dblsgl.$$| $DBL2SGL
   rm -f /tmp/dblsgl.$$
   rm $i
  done
 cd ..
done
