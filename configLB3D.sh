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
CODEDIRS="src"
DOCDIRS="doc doc/doxygen"

######################################################################
# Function definition
######################################################################

function getFlags()
{
for CODEDIR in ${CODEDIRS}; do
        cd ${CODEDIR}
        FLAGLIST=$(ls *.F90 *.h *.c | xargs grep '^#e\?l\?ifn\?\(def\)' | sed 's/.*#e\?l\?ifn\?\(def\) \(.*\)/\2/' | sed 's/[ \t]*$//' | sort | uniq)
        for FLAG in ${FLAGLIST}; do
                echo -D${FLAG}
        done
        cd ..
done
}

function usage()
{
echo Syntax is $0 OPTION PLATFORM COMPILERFLAGS
echo
echo OPTION is one of:
echo CONFIG
echo MAKE
echo MAKE-NODOC
echo MAKE-ONLYDOC
echo MKLBE
echo CLEAN
echo CLEANLBE
echo
echo PLATFORM is one of:
for i in defines.* ; do
        new=`echo $i|sed s/defines.//`
        echo -n $new " "
done
echo 
echo
echo  COMPILERFLAGS is one or more of:
getFlags
exit -1

}

######################################################################
# Check commandline arguments 
######################################################################

if  expr $# '<' 1 > /dev/null  ; then
        usage
fi

# Check if option is valid

OPTION=$1
shift

case ${OPTION} in
    CONFIG)
        ;;
    MAKE | MAKE-NODOC | MAKE-ONLYDOC | MKLBE)
        ;;
    CLEAN | CLEANLBE)
        ;;
    *)
        echo OPTION contains an unknown value.
        usage
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

# Check if compilerflags are valid

VALIDFLAGS=$(getFlags)
while [[ $1 ]]; do
ISVALID=0
for VALIDFLAG in ${VALIDFLAGS};do
        if [ $1 == ${VALIDFLAG} ]; then
                ISVALID=1
                MAKEFFLAGS="$MAKEFFLAGS $1"
		if [ $1 == "-DP3M" ]; then
		    FFTWFLAGS="-ldfftw"
		fi
        fi
done
        if [ ${ISVALID} == 0 ]; then
                echo I do not know about compilerflag \"$1\"
                echo 
                usage
        fi
        shift
done

echo
echo "Creating Makefiles ..."
echo

for i in $LIBDIRS $DOCDIRS; do
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
		CODEMAKECFLAGS=$MAKEFFLAGS
                ;;
            *)
	        CODEMAKEFFLAGS=$MAKEFFLAGS
		CODEMAKECFLAGS=$MAKEFFLAGS
                ;;
        esac
    fi

    echo "MAKEFFLAGS=$CODEMAKEFFLAGS" > $TEMPLATE
    echo "MAKECFLAGS=$CODEMAKECFLAGS" >> $TEMPLATE
    echo -e "PLATFORM=$PLATFORM" >> $TEMPLATE
    echo "XDRLIB=-L${MYDIR}/lib/xdrf" >> $TEMPLATE
    echo "FFTWFLAGS=${FFTWFLAGS}" >> $TEMPLATE
    cat $MYDIR/$DEFFILE >> $TEMPLATE

    cd $MYDIR/$i
    $MYDIR/mkmf.pl -m $OUTFILE -t $TEMPLATE -p $EXECUTABLE
    # Adding the revision script
    cat $OUTFILE > $TMPFILE
    $MYDIR/get-version.sh -m >> $TMPFILE
    # Change the Makefile so that 'make' will also recreate the revision.h file.
    cat $TMPFILE | sed "s/all: $EXECUTABLE/all: version $EXECUTABLE/" > $OUTFILE
    rm -rf $TMPFILE
    cd $MYDIR
done

fi

######################################################################
# Execute different variations of MAKE
######################################################################

case ${OPTION} in
    MAKE | MAKE-NODOC | MAKE-ONLYDOC | MKLBE)

echo
echo "Starting compilation ..."
echo

case ${OPTION} in
    "MAKE")
        MAKEDIRS="$LIBDIRS $CODEDIRS $DOCDIRS"
        ;;
    "MAKE-NODOC")
        MAKEDIRS="$LIBDIRS $CODEDIRS"
        ;;
    "MAKE-ONLYDOC")
        MAKEDIRS="$DOCDIRS"
        ;;
    "MKLBE")
        MAKEDIRS="$CODEDIRS"
        ;;
esac

for i in $MAKEDIRS; do
    echo --------- $i ---------
    cd $MYDIR/$i
    make  -j 12
    ERR=$?
    if [ "$ERR" -ne "0" ]; then
        echo "Exit code $ERR returned during compilation of '$i', aborting $0 ..."
        exit $ERR
    else
        echo "Succesfully compiled $i !"
    fi
    cd $MYDIR
done

;;

esac

######################################################################
# Execute different variations of CLEAN
######################################################################

case ${OPTION} in 
    CLEAN | CLEANLBE)

echo
echo "Cleaning up ..."
echo

case ${OPTION} in
    CLEAN)
        CLEANDIRS="$LIBDIRS $CODEDIRS $DOCDIRS"
	;;
    CLEANLBE)
        CLEANDIRS="$CODEDIRS"
	;;
esac

for i in $CLEANDIRS; do
    echo --------- $i ---------
    cd $MYDIR/$i
    make clean
    cd $MYDIR
done

# Harsher cleanup for CLEAN
case ${OPTION} in
    CLEAN)
	# Clean Makefiles etc.
        for i in $CLEANDIRS; do
            cd $MYDIR/$i
            rm -fv Makefile Makefile.template lbe_version.h
            cd $MYDIR
        done
        # Clean up after testbuilds.sh
        rm -fv defines.*.NOFLAGS
        rm -fv tests/*.txt
        ;;

esac

;;

esac

