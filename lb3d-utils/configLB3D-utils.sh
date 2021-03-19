#!/bin/bash

# This script sets up a Makefile suitable for the platform the code
# is to be compiled on.
#
# If you try to compile on some unsupported platform foo, then create
# a file called "defines.foo" with the appropriate flags for the compilers,
# then type "./configLB3D.sh foo".
#

######################################################################
# Variable definition
######################################################################

MYDIR=$(pwd)
ME=$0

LIBDIRS="lib/xdrf"
CODEDIRS="src/borders src/cpost src/demos src/dropletTools src/ellipsoidPickeringTools src/ftools src/h5converter src/local_tau src/misc src/moreDropletTools src/permeability src/post src/rock_calculations src/rock_embedding src/rock_gen_spheres src/scripts src/struct_xdr src/tests src/vtk src/xdrtools src/XMTutils"

######################################################################
# Function definition
######################################################################

function usage()
{
echo Syntax is $0 OPTION PLATFORM
echo
echo OPTION is one of:
echo CONFIG
echo MAKE
echo CLEAN
echo
echo PLATFORM is one of:
for i in defines.* ; do
        new=`echo $i|sed s/defines.//`
        echo -n $new " "
done
echo 
exit -1

}

######################################################################
# Check commandline arguments
######################################################################

GLOBAL_ERR=0

if  expr $# '<' 1 > /dev/null  ; then
        usage
fi

# Check if option is valid

OPTION=$1
shift

case ${OPTION} in
    CONFIG)
        ;;
    MAKE)
        ;;
    CLEAN)
        ;;
    *)
        echo OPTION contains an unknown value.
        exit -1
esac

######################################################################
# Execute CONFIG
######################################################################

if [ "${OPTION}" = "CONFIG" ] ; then

# Check if platform is valid

PLATFORM=$1
shift
DEFFILE=defines.$PLATFORM

        if [ ! -f $DEFFILE ]; then
	        echo I do not know about platform \"$PLATFORM\".
	        echo
                usage
        fi

echo
echo "Creating Makefiles ..."
echo

for i in $LIBDIRS; do
    echo $i
    cat $MYDIR/$DEFFILE > $MYDIR/$i/Makefile
    echo -e "MAKEFFLAGS=$MAKEFFLAGS" >> $MYDIR/$i/Makefile
    echo -e "PLATFORM=$PLATFORM" >> $MYDIR/$i/Makefile
    echo "XDRLIB=-L${MYDIR}/lib/xdrf" >> $MYDIR/$i/Makefile
    cat $MYDIR/$i/Makefile.STUB >> $MYDIR/$i/Makefile
done

for i in $CODEDIRS ; do
    echo $i
    EXECUTABLE=lbe
    TEMPLATE=$MYDIR/$i/Makefile.template
    OUTFILE=$MYDIR/$i/Makefile
    TMPFILE=$MYDIR/$i/Makefile.tmp

    # change command line syntax for preprocessor switches if LB3D is
    # to be compiled by an IBM compiler
    if echo $MAKEFFLAGS |grep '[^ ]' > /dev/null
    then
        PLATFORMNF=`echo $PLATFORM | sed 's/\.NOFLAGS//'` # Compatibility fix for testbuilds.sh 
	case $PLATFORMNF in
	    JUQUEEN | HUYGENS | JUGENE | SP2MPI)
		CODEMAKEFFLAGS=-WF,`echo $MAKEFFLAGS |sed 's/ /,/g'`
		;;
	    *)
		CODEMAKEFFLAGS=$MAKEFFLAGS
		;;
	esac
    fi

    echo "MAKEFFLAGS=$CODEMAKEFFLAGS" > $TMPFILE
    echo -e "PLATFORM=$PLATFORM" >> $TMPFILE
    echo "XDRLIB=-L${MYDIR}/lib/xdrf" >> $TMPFILE
    cat $MYDIR/$DEFFILE >> $TMPFILE

    cd $MYDIR/$i

    if [ -f "Makefile.STUB" ]; then
	cat Makefile.STUB >> $TMPFILE
    fi

    mv $TMPFILE $OUTFILE
    cd $MYDIR
done

fi

######################################################################
# Execute MAKE
######################################################################

case ${OPTION} in
    MAKE)

echo
echo "Starting compilation ..."
echo

MAKEDIRS="$LIBDIRS $CODEDIRS" 

for i in $MAKEDIRS; do
    echo --------- $i ---------
    cd $MYDIR/$i
    make
    ERR=$?
    if [ "$ERR" -ne "0" ]; then
        echo "Exit code $ERR returned during compilation of '$i', but continuing to next util ..."
        GLOBAL_ERR=$ERR
    else
        echo "Succesfully compiled $i !"
    fi
    cd $MYDIR
done

;;

esac

######################################################################
# Execute CLEAN
######################################################################

case ${OPTION} in 
    CLEAN)

echo
echo "Cleaning up ..."
echo

CLEANDIRS="$LIBDIRS $CODEDIRS"

for i in $CLEANDIRS; do
    echo --------- $i ---------
    cd $MYDIR/$i
    make clean
    cd $MYDIR
done

;;

esac

exit $GLOBAL_ERR

