#!/bin/bash

if [ $# -lt 1 ]
then
    echo "usage: clean-old-checkpoints.sh [dir[s ...]]"
    echo
    echo "Deletes all but the latest lb3d checkpoint files in each of the"
    echo "directories passed as argument. Does nothing if the latest"
    echo "checkpoint is incomplete. Relies on new rXDR format (including"
    echo "checktopo_*.xdr$ file) for completeness check. Use on your own"
    echo "risk!"
    exit -1
fi

for d in $*
do
    oldd=`pwd`
    echo -n "$d/..."
    cd $d

    if [ `ls |grep ^md- |wc -l` -gt 0 ]
	then
	echo -n ' [MD]'
	with_MD='true'
    else
	with_MD='false'
    fi

    if [ `ls |grep _IBM_ |wc -l` -gt 0 ]
	then
	echo -n ' [IBM]'
	with_IBM='true'
    else
	with_IBM='false'
    fi

    echo -ne "\n"

    # find last checkpoint
    chkuid=`ls checkparam*xdr |tail -n 1`
    chkuid=${chkuid%.*}
    chkuid=${chkuid##*_}

    # check for completeness of this checkpoint
    lastcomplete='true'

    topofile=checktopo_*_$chkuid.xdr
    if [ -f $topofile ]
    then
	cdx=`dump-xdrints.py 1 $topofile`
	cdy=`dump-xdrints.py 2 $topofile`
	cdz=`dump-xdrints.py 3 $topofile`
	nprocs=$[cdx*cdy*cdz]
    else
	nprocs=-1
	echo "$d/: missing file $topofile"
	lastcomplete='false'
    fi

    ncheckpointfiles=`ls \
	    |grep \
	    ^checkpoint_.\*_${chkuid}_p[0-9][0-9][0-9][0-9][0-9][0-9]\\.xdr \
            |grep -v ^checkpoint_IBM_ \
	    |wc -l`
    if [ "$nprocs" -ne "$ncheckpointfiles" ]
	then
	echo "$d/: missing LB checkpoint file(s) for $chkuid \
(found only $ncheckpointfiles chunks instead of $nprocs)"
	lastcomplete='false'
    fi

    if [ "$with_IBM" = 'true' ]
	then
	nibmcheckpointfiles=`ls \
	    |grep \
	    ^checkpoint_IBM_.\*_${chkuid}_p[0-9][0-9][0-9][0-9][0-9][0-9]\\.xdr \
	    |wc -l`
	if [ "$nprocs" -ne "$nibmcheckpointfiles" ]
	    then
	    echo "$d/: missing IBM checkpoint file(s) for $chkuid \
(found only $nibmcheckpointfiles chunks instead of $nprocs)"
	    lastcomplete='false'
	fi
    fi

    if [ "$with_MD" = 'true' ]
	then
	mdfile=md-checkpoint_*_$chkuid.xdr
	if [ ! -f $mdfile ]
	    then
	    echo "$d/: missing file $mdfile"
	    lastcomplete='false'
	fi
    fi

    # delete all checkpoints but the last if the last is complete
    if [ "$lastcomplete" = 'true' ]
    then
	ls |grep '^[^_]*check.*\.\(h5\|xdr\)$' |grep -v $chkuid |while read f
	do
	    rm $f
	done
	echo "$d/: latest checkpoint is $chkuid"
    else
	echo "$d/: latest checkpoint $chkuid is incomplete---directory skipped"
    fi

    cd $oldd
done
