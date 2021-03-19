#!/bin/bash

SFILE=$1
SCALE=$2

echo "Converting $SFILE to block tool format"
BFILE="$1.block"
SIZETEXT=`./spherepacking_to_blocktool.py $SFILE $BFILE -s $SCALE | grep "System size:"`
echo "Converting $BFILE to xdr"
./block $BFILE | grep -ve "^sphere"
rm "$BFILE.txt"

NX=`echo $SIZETEXT | awk '{print $3}'`
NY=`echo $SIZETEXT | awk '{print $4}'`
NZ=`echo $SIZETEXT | awk '{print $5}'`

XFILE="$BFILE.xdr"
XOUTFILE="$XFILE-out.h5"

../h5converter/h5converter -i $XFILE xdr -o $XOUTFILE hdf -s $NX $NY $NZ -d OutArray 

HFILE=`echo $XOUTFILE | sed 's/block.xdr-out.//'`

rm -f $BFILE
rm -f $XFILE

mv $XOUTFILE $HFILE

