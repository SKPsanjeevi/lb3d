#!/bin/bash

TF="permeability.txt"

echo "# time permeability massflux" > $TF

for f in profile-z*asc
do
  T=`echo $f | sed 's/.*t\([0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]\)-.*/\1/'`
  MF=`grep -ve "^#" $f | awk '{sum+=$7;} END {printf("%E", sum/NR);}'`
  PERM=`grep -ve "^#" $f | awk '{sum+=$11;} END {printf("%E", sum/NR);}'`
  echo "$T $PERM $MF" >> $TF
done

