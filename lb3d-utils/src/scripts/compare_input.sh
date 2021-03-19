#!/bin/bash

echo "Usage $0 <hdf5file1> <hdf5file2> "
echo "This will attempt to read a pair of HDF5 input files in LB3D format, extract the metadata, and run a diff on it."

if [ $# -ne 2 ] ; then
	echo "ERROR: compare_input requires two hdf5 files, aborting..."
	exit 0
fi

# Use h5dump to get plaintext attributes (lb3d metadata), use grep to select only lines with (nnn) at the start.
# Then use sed to remove the numbers, remove the whitespace, and the first and last line (Which contain the headers).
# Finally, remove leading and trailing whitespace. Everything is then shoved into diff.

diff <(h5dump -A $1 | grep '([0-9]*)' | sed 's/: "\(.*\)",/\1/' | sed 's/([0-9]*)\(.*\)/\1/' | sed 's/^[ \t]*//' | sed '1d' | sed '$d' | sed 's/^[ \t]*//;s/[ \t]*$//') <(h5dump -A $2 | grep '([0-9]*)' | sed 's/: "\(.*\)",/\1/' | sed 's/([0-9]*)\(.*\)/\1/' | sed 's/^[ \t]*//' | sed '1d' | sed '$d' | sed 's/^[ \t]*//;s/[ \t]*$//')

