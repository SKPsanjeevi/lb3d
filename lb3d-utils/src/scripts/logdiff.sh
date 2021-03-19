#!/bin/bash

echo "Usage $0 <file1> <file2> "
echo "This will attempt to read a pair of LB3D stdout files, remove timestamps, and run a diff on it."

if [ $# -ne 2 ] ; then
        echo "ERROR: logdiff requires two hdf5 files, aborting..."
        exit 0
fi

diff <(cat $1 | grep -v -i "allocating" | sed 's/[0-9][0-9]:[0-9][0-9]:[0-9][0-9] - \(.*\)/\1/') <(cat $2 | grep -v -i "allocating" | sed 's/[0-9][0-9]:[0-9][0-9]:[0-9][0-9] - \(.*\)/\1/')

