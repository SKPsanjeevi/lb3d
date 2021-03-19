#!/bin/bash

# h5repack path
H5REPACK=h5repack

# Level of nice
NICELVL="-19"

# Requested repack options
RQOPTS="-f SHUF -f GZIP=9"

# Requested filter output
RQFILTERS="FILTERS { PREPROCESSING SHUFFLE COMPRESSION DEFLATE { LEVEL 9 } }"

# Override defaults
while getopts "g:n:" flag
do
case "$flag" in
        "g")
                RQOPTS="-f SHUF -f GZIP=$OPTARG"
                RQFILTERS="FILTERS { PREPROCESSING SHUFFLE COMPRESSION DEFLATE { LEVEL $OPTARG } }"
                ;;
        "n")
                NICELVL="-$OPTARG"
                ;;
esac
done

echo "=== STARTING BATCH HDF5 REPACKING ==="
echo
echo "Executing nice $NICELVL $H5REPACK $RQOPTS on HDF5 files."

# LOOP STARTS HERE

ERR=0
DERR=0

for f in `find . -name "*.h5"` ; do
  TS=`date +%H:%M:%S`
  echo "$TS - $f"

  HDF5FILTERS=`h5dump -pH $f | sed 's/[ \t]*//' | awk 'BEGIN {f=0;} {if ($1 == "FILTERS") {f=1;} if (f==1) {printf("%s ",$0)}; if ($1 == "}") {f=0;} }'`
  # Get rid of linebreak?
  HDF5FILTERS=`echo $HDF5FILTERS`

  if [[ $HDF5FILTERS == $RQFILTERS ]]
  then
    echo "  Already compressed."
  else
    echo "  Repacking..."
    nice $NICELVL $H5REPACK $RQOPTS -i $f -o $f.shuf.gz9
    if [ $? -eq 0 ] ; then
      HDF5INFO=`h5dump -pH $f.shuf.gz9 | sed 's/[ \t]*//' | awk 'BEGIN {f=0;} {if ($1 == "FILTERS" || $1 == "STORAGE_LAYOUT") {f=1;} if (f==1) {printf("%s ",$0)}; if ($1 == "}") {f=0;} }'`
      echo "  No error: $HDF5INFO"
      mv $f.shuf.gz9 $f
    else
      let ERR=ERR+1
      let DERR=DERR+1
      echo "  Failure #$ERR"
    fi
  fi
done # file loop

echo
echo "Completed with $ERR errors."
echo
echo "=== DONE! === "
echo

