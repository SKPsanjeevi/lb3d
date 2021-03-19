#!/bin/bash
# Small example script that shows how to rename a lot of datasets in
# different directories. It goes through all the directories given at the
# command line and renames some of the files.
# THIS SCRIPT IS NOT MEANT TO BE USED AS IT IS. IT HAS TO BE ADAPTED.
# Jens, 09.10.02

for j in $* ; do
 echo $j
 cd $j
  for i in `ls *2.xdr` ; do
  NEW=`echo $i|sed s/.2.xdr/.xdr/g`
  mv $i $NEW
  done
  for i in `ls *2.fld` ; do
  NEW=`echo $i|sed s/.2.fld/.fld/g`
  mv $i $NEW
  done
 cd ..
done
