#!/bin/bash

if [ $# -ne 2 ]
then
    echo "usage: trajectories.sh <particle-uid> <path to asc md-cfg file>"
    echo
    echo "Writes all md-cfg data related to particle with uid <particle-uid>"
    echo "from files that differ from <path to asc md-cfg file> only in the"
    echo "time step to standard output. One line is written per time step"
    echo "and the time step is prepended to each line."
    exit -1
fi

particle=`printf '%.10i' $1`

dir=`dirname $2`
basename=`basename $2 .asc`
uid=${basename#*_t????????-}
name=${basename%_t????????-??????????}

ls $dir \
    | grep "^${name}_t[0-9]\{8\}-${uid}\.asc$" \
    | while read file
  do
  timestep=${file%-??????????.asc}
  echo -n "${timestep##*_t} "
  grep $particle ${dir}/$file
done
