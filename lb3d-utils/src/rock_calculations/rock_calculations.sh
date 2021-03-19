#!/bin/bash

# This script calls ./rock_calculations for every timestep in a particular run
# File name formats expected are od_<name>_tnnnnnnnn-<code>.h5, vel_<name>_tnnnnnnnn-<code>.h5 and rock_<name>_t00000000-<code>.h5
# Paramfile is expected to be an ASCII file containing wall thickness, iolet length, tau, delta_x, delta_t, all on new lines
# Without optional flags output is written to perm_<name>_tnnnnnnnn-<code>.txt files
# If --dumppermeability is specified, output files only contain the calculated permeability. Additionally, the file perm_<name>-<code>.txt contains data in a two-column format <t> <perm>

echo "Usage $0 <name> <code> <paramfile> [--dumppermeability] "

rm perm_$1-$2.txt -f

for file in od_$1_*$2*
do
  t=`echo $file | sed "s/od_$1_t//" | sed "s/-$2.h5//" `
  echo "Executing for t = $t"
  cat $3 | ./rock_calculations od_$1_t$t-$2.h5 vel_$1_t$t-$2.h5 rock_$1_t00000000-$2.h5 > perm_$1_t$t-$2.txt
  if [ $# -gt 3 ]
  then
    if [ $4 == "--dumppermeability" ]
    then
      cat $3 | ./rock_calculations od_$1_t$t-$2.h5 vel_$1_t$t-$2.h5 rock_$1_t00000000-$2.h5 $4 > perm_$1_t$t-$2.dump.txt
  
      p=`cat perm_$1_t$t-$2.dump.txt`
      echo "$t $p" >> perm_$1-$2.txt
      rm perm_$1_t$t-$2.dump.txt
    fi
  fi
done
