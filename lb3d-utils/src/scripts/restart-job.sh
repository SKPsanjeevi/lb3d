#!/bin/bash

if [ $# -ne 2 ]
then
    echo "usage: $0 <new final timestep> <directory>"
    echo
    echo "Creates/modifies new lb3d (version5-md) input files 'restart'"
    echo "and 'restart.md' in the directory <directory> that---when run---"
    echo "evolve the simulation from the latest checkpoint to time step"
    echo "<new final timestep>. The source input files are restart[.md] or"
    echo "---if not available---input-file[.md]."
    exit
fi
end=$1
dir=$2

cd $dir

lastif=`ls input-file restart |sort |tail -n 1`
if [ "$lastif" != "restart" ]
then
    cp input-file restart
    cp input-file.md restart.md
fi

chkuid=`ls checkparam*xdr |tail -n 1`
chkuid=${chkuid%.*}
chkuid=${chkuid##*_}

sed "s/n_iteration[ 	]*=[ 	]*[0-9][0-9]*/n_iteration = $end/" restart \
    |sed \
    "s/restore_string[ 	]*=[ 	]*\"[-t0-9]*\"/restore_string = \"$chkuid\"/" \
    |sed "s/init_cond[ 	]*=[ 	]*[-0-9]*/init_cond = 7/" \
    >restart.tmp

mv restart.tmp restart
