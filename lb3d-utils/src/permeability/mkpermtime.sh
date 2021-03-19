#!/bin/bash
# This script takes the flow output files of mflowscript.sh and greps for
# permeability and timestep.
#
# Jens, 20.03.08

rm -f permtime.dat
for i in flow*.dat ; do
echo -n "`echo $i|sed 's/^.*_t// ; s/-[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].*$//'` "
grep sample_perm_tot $i|sed s/"!sample_perm_tot= "//
done

