#!/bin/bash

if [ $# -lt 3 ]
then
    echo "usage: $0 <tnx> <tny> <tnz> [xdr file[s ...]]"
    echo
    echo "Converts zero or more scalar or vector xdr files related to lb3d to"
    echo "vtk format by prepending an appropriate vtk header. The dimensions"
    echo "of the input files must be given as tnx, tny, and tnz (resembling"
    echo "nx, ny, and nz in the respective lb3d input file."
    exit -1
fi

tnx=${1##0}
tny=${2##0}
tnz=${3##0}
shift 3

npoints=$[tnx*tny*tnz]

for path in $*
do
    file=${path##*/}
    case ${file%%_*} in
	$file)
            # this avoids that all obstacle files have to start with "rock_"
	    echo "assuming $file to be a rock file..."
	    data="SCALARS rock_state float 1
LOOKUP_TABLE default"
	    ;;
	od | pxx | pxy | pxz | pyy | pyz | pzz | v1 | v2 | v3 | rock | vel | sum-lbe-force)
	    data="SCALARS OutArray float 1
LOOKUP_TABLE default"
	    ;;
        arr)
	    data="VECTORS OutArray float"
            ;;
	checkparams | checkpoint | md-checkpoint)
	    echo "skipping checkpoint file: $file"
	    continue
	    ;;
	*)
	    echo "unknown kind of file: $file has kind ${file%%_*}"
	    exit
	    ;;
    esac
    (
	cat <<EOF
# vtk DataFile Version 2.0
data
BINARY
DATASET STRUCTURED_POINTS
DIMENSIONS $tnx $tny $tnz
ORIGIN 1 1 1
SPACING 1 1 1
POINT_DATA $npoints
$data
EOF
	cat $path
	)>${path%.xdr}.vtk
done
