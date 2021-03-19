#!/bin/bash
#
#
#  Takes particle coordinates as given by MD-SRD and creates a discretized
#  rock file for lb3d
#
#  Jens, 20.03.08
#


# Set variables here:

# TOOLPATH contains the directory where voxels and createtunnel reside
TOOLPATH=/data/data2/hecht/lbe/branches/version5/lbe/utils/permeability

#Particle file as obtained from MD-SRD
PARTICLEFILE=$1
#Name of final output (without suffix)
OUTFILE=$2

NUMSPHERES=`wc -l $PARTICLEFILE | awk '{print $1}'` 

RADIUS=$3

DIMX=`awk 'BEGIN{maxx=0.;minx=1.;}
{if($1>maxx){maxx=$1;}if($1<minx){minx=$1;}}
END{print (maxx - minx)}' $PARTICLEFILE`

DIMY=`awk 'BEGIN{maxy=0.;miny=1.;}
{if($2>maxy){maxy=$2;}if($2<miny){miny=$2;}}
END{print (maxy - miny)}' $PARTICLEFILE`

DIMZ=`awk -v r=$RADIUS 'BEGIN{maxz=0.;minz=1.;}
{if($3>maxz){maxz=$3;}if($3<minz){minz=$3;}}
END{print (maxz - minz + 2.0*r)}' $PARTICLEFILE`

OFFSETX=`awk 'BEGIN{maxx=0.;minx=1.;}
{if($1>maxx){maxx=$1;}if($1<minx){minx=$1;}}
END{print (minx)}' $PARTICLEFILE`

OFFSETY=`awk 'BEGIN{maxy=0.;miny=1.;}
{if($2>maxy){maxy=$2;}if($2<miny){miny=$2;}}
END{print (miny)}' $PARTICLEFILE`

OFFSETZ=`awk -v r=$RADIUS 'BEGIN{maxz=0.;minz=1.;}
{if($3>maxz){maxz=$3;}if($3<minz){minz=$3;}}
END{print (minz - r)}' $PARTICLEFILE`

echo "DIMX=$DIMX"
echo "DIMY=$DIMY"
echo "DIMZ=$DIMZ"
echo "OFFSETX=$OFFSETX"
echo "OFFSETY=$OFFSETY"
echo "OFFSETZ=$OFFSETZ"
NUMVOXELS=100

ROCKDIMX=128
ROCKDIMY=128
ROCKDIMZ=128


##############################################################

TMPFILE=particles2rock.$$
TMPFILE2=particles2rock2.$$
TMPFILE3=particles2rock3.$$
#Convert particle positions

DIM=`echo | awk -v x=$DIMX -v y=$DIMY '{if(x>y){print x}else{print y}}'`

awk "{print (\$1 - $OFFSETX) \" \"  (\$2 - $OFFSETY) \" \"  (\$3 - $OFFSETZ) \" $RADIUS\"}" $PARTICLEFILE >$TMPFILE
#cp $TMPFILE $TMPFILE.$$
#Run Bibhu's voxel tool
cat > $TMPFILE2 << EOF
$TMPFILE
$TMPFILE3
$NUMSPHERES
$NUMVOXELS
$DIM
EOF

cat $TMPFILE2 | $TOOLPATH/voxels >/dev/null
#cp $TMPFILE3 $TMPFILE3.$$
#Threshold:
awk '{if($1>108){print 1}else{print 0}}' $TMPFILE3 > $TMPFILE

#Create lb3d rock file
$TOOLPATH/createtunnel $TMPFILE $NUMVOXELS $NUMVOXELS $NUMVOXELS $ROCKDIMX $ROCKDIMY $ROCKDIMZ u2m >/dev/null

FILENAME=$TMPFILE"_"$ROCKDIMX"_"$ROCKDIMY"_"$ROCKDIMZ"_r7_u2m"
OUTFILE=$OUTFILE"_"$ROCKDIMX"_"$ROCKDIMY"_"$ROCKDIMZ
echo
echo Writing $OUTFILE.xdr and $OUTFILE.vtk.
mv $FILENAME.xdr $OUTFILE.xdr
mv $FILENAME.vtk $OUTFILE.vtk

rm -f $TMPFILE $TMPFILE2 $TMPFILE3
