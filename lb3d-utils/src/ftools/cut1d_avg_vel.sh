#!/usr/freeware/bin/bash
# /bin/bash   # use this for linux, with !
# Small example script that shows how to get 1D plots out of a lot of data
# sets in
# different directories. It goes through all the directories given at the
# command line and converts some of the files.
# 
# THIS SCRIPT IS NOT MEANT TO BE USED AS IT IS. IT HAS TO BE ADAPTED.
# Jens, 30.10.02
#
#
DIRECTORY=/sanhp/uqngs/lbe/lbe/studies/shear_steady/output
CUT1D=$DIRECTORY/cut1d_avg
FILES="vel*.xdr" # Which files do you want to 'cut' ?
FTYPE=2          # 1=double, 2=float
DIMS="128,128,128"  # Lattice size
TYPE=1           # 1=scalar, 2=2scalar
VDIR="X"         # Variable directions: DELTA_? only work for VDIR='X' !!!
AVG3D=1          # Want to do 3D-avg (1/0)?
DELTA_Z=16       # How much region z\in(0,delta)U(NZ-delta,NZ) want to exclude?
                 # Max three digits for now.
DELTA_X=8

#
# Next line caters for old _t?????? format
# and new _t??????-?????????? fingerprint format
# for filenames
#
STATF="vel_*_t*_YZavg.X.Stats"
STATF2="vel_"$*"_DX"$DELTA_X"_DZ"$DELTA_Z"_YZavg.X.Stats"

for j in $* ; do
 echo 'Directory to read from is:'
 echo $j
 cd $j
  pwd
  for k in $VDIR ; do
   echo 'Files to read from are:'
   ls $FILES
   for i in `ls $FILES` ; do
      echo $i >/tmp/cut1d.$$
      echo "$FTYPE" >>/tmp/cut1d.$$
      echo "$DIMS" >>/tmp/cut1d.$$
      echo "$TYPE" >>/tmp/cut1d.$$
      echo "$k" >>/tmp/cut1d.$$
      echo "$AVG3D" >> /tmp/cut1d.$$
      echo "$DELTA_Z" >>/tmp/cut1d.$$
      echo "$DELTA_X" >>/tmp/cut1d.$$
      cat /tmp/cut1d.$$
      cat /tmp/cut1d.$$| $CUT1D > log
 #     echo '	** TMP FILES TO DELETE'
 #     find /tmp -name cut1d.$$ | xargs
 #     rm -f /tmp/cut1d.$$
   done
 # Delete previously written STATF file:
 #  find $DIRECTORY/$* -name '$STATF2' | xargs
   rm -f $DIRECTORY/$*/$STATF2
   echo 'Now appending .Stat files into:' 
   echo $DIRECTORY/$*/$STATF2
   for i in `ls $STATF` ; do
 #$STATF` ; do
 #     This filename doesn't specify sampling (cut) axis yet
      echo 'Writing '$i' to ' $DIRECTORY/$*/$STATF2
      cat $i >> $DIRECTORY/$*/$STATF2
   done
   echo '	** TEMPORAL STATS FILES TO DELETE'
 #   find $DIRECTORY/$* -name '$STATF'
   rm -f $DIRECTORY/$*/$STATF
  done
 cd ..
done

