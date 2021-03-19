#!/bin/bash
# Small example script that shows how to convert a lot of datasets in
# different directories. It goes through all the directories given at the
# command line and converts some of the files. You have to set the lattice
# size and data type.
# THIS SCRIPT IS NOT MEANT TO BE USED AS IT IS. IT HAS TO BE ADAPTED.
# Jens, 09.10.02

DBL2SGL=/volta/scratch/jens/lbe/utils/dbl2sgl/dbl2sgl
for j in $* ; do
 echo $j
 cd $j
 rm -f *.fld
  for i in `ls dir*0.xdr flo*0.xdr` ; do
   echo $i >/tmp/dblsgl.$$
   echo "64,64,64" >>/tmp/dblsgl.$$
   echo "4" >>/tmp/dblsgl.$$
   cat /tmp/dblsgl.$$| $DBL2SGL
   rm -f /tmp/dblsgl.$$
   rm $i
  done
  for i in `ls sur*0.xdr colour*0.xdr vel*0.xdr` ; do
   echo $i >/tmp/dblsgl.$$
   echo "64,64,64" >>/tmp/dblsgl.$$
   echo "1" >>/tmp/dblsgl.$$
   cat /tmp/dblsgl.$$| $DBL2SGL
   rm -f /tmp/dblsgl.$$
   rm $i
  done
  for i in `ls owd*0.xdr` ; do
   echo $i >/tmp/dblsgl.$$
   echo "64,64,64" >>/tmp/dblsgl.$$
   echo "2" >>/tmp/dblsgl.$$
   cat /tmp/dblsgl.$$| $DBL2SGL
   rm -f /tmp/dblsgl.$$
   rm $i
  done
 cd ..
done
