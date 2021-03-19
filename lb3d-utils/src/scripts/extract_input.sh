#!/bin/bash

echo "Usage $0 <inputfile> <outputfile> "
echo "This will attempt to read an HDF5 input file in LB3D format, extract the metadata, and reproduce the input-file, input-file-diff and input-file.md to outputfile, outputfile-diff and outputfile.md, respectively."

if [ $# -ne 2 ] ; then
        echo "ERROR: extract_input requires input and output file name, aborting..."
        exit 0
fi

# Use h5dump to get plaintext attributes (lb3d metadata), use grep to select only lines with (nnn) at the start, use awk to select the right sections
# Then use sed to remove the numbers, remove the whitespace, and the first and last line (Which contain the headers).
# Finally, remove leading and trailing whitespace.
h5dump -A $1 | grep '([0-9]*)' | awk '/\( Start input file \)/,/\( End input file \)/' | sed 's/: "\(.*\)",/\1/' | sed 's/([0-9]*)\(.*\)/\1/' | sed 's/^[ \t]*//' | sed '1d' | sed '$d' |  sed 's/^[ \t]*//;s/[ \t]*$//' > $2
h5dump -A $1 | grep '([0-9]*)' | awk '/\( Start input file MD \)/,/\( End input file MD \)/' | sed 's/: "\(.*\)",/\1/' | sed 's/([0-9]*)\(.*\)/\1/' | sed 's/^[ \t]*//' | sed '1d' | sed '$d' | sed 's/^[ \t]*//;s/[ \t]*$//' > $2.md
h5dump -A $1 | grep '([0-9]*)' | awk '/\( Start differential input file \)/,/\( End differential input file \)/' | sed 's/: "\(.*\)",/\1/' | sed 's/([0-9]*)\(.*\)/\1/' | sed 's/^[ \t]*//' | sed '1d' | sed '$d' | sed 's/^[ \t]*//;s/[ \t]*$//' > $2-diff

echo "Done."

