#!/bin/bash

if [ $# -lt 1 ]
then
    echo "usage: clean-duplicates.sh [dir[s ...]]"
    echo
    echo "Searches the current directory for files with different uids"
    echo "belonging to the same time step. If such files are found, all but"
    echo "the one with the uid belonging to the latest run are removed. This"
    echo "way, a continuous and unique time series of data is achieved even if"
    echo "for the simulation setup, restoring from a checkpoint should give"
    echo "rise to some deviation of the results."
    echo
    echo "Don't use this without prior testing for your directory structure,"
    echo "file name convention, etc. !"
    exit -1
fi

O=CD.txt

for d in $*
do
    oldd=`pwd`
    echo -n "$d: "
    cd $d

    NPRINT=100

    files=`ls -rt \
	|grep -v '^check' \
	|grep '.*_t[0-9]\{8\}-[0-9]\{10\}\.\(asc\|xdr\|h5\)$'`

    tsuids=`echo "$files" |grep -o '[0-9]\{8\}-[0-9]\{10\}' |sort |uniq`

    nuids=`echo "$tsuids" |grep -o '[0-9]\{10\}' |sort |uniq |wc -l`
    nuid=0

    ndelete=0

    while [ -n "$tsuids" ]
    do
	luid=`echo "$tsuids" |tail -n 1`
	luid=${luid##*-}

	ofiles=`echo "$files" |grep -v $luid`

	nuid=$[nuid+1]
	echo -e "\nchecking for duplicates of $luid... [uid $nuid/$nuids]" >>$O

	lfiles=`echo "$files" |grep "_t[0-9]\{8\}-$luid\.\(asc\|xdr\|h5\)\$"`

	nlfiles=`echo "$lfiles" |wc -l`
	nlfile=0

	while read lfile
	do
	    nlfile=$[nlfile+1]
	    if [ "$[nlfile%NPRINT]" -eq 0 ]
	    then
		echo -e "\n$lfile [file $nlfile/$nlfiles] [uid $nuid/$nuids]" >>$O
	    else
		echo -n "." >>$O
	    fi

	    uptots=${lfile%-[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].*}
	    ext=${lfile##*.}

	    oofiles=`echo "$ofiles" |grep "^$uptots-[0-9]\{10\}.$ext\$"`

	    if [ -n "$oofiles" ]
	    then
		if [ -f $lfile ]
		then
		    while read rfile
		    do
			echo -e "\nrm -f $rfile" >>$O
			rm -f $rfile
			ndelete=$[ndelete+1]
		    done < <( echo "$oofiles" )
		fi
	    fi
	done < <( echo "$lfiles" )

	end=`echo "$tsuids" |grep -nm 1 $luid`
	end=${end%:*}
	tsuids=`echo "$tsuids" |head -n $[end-1]`
    done
    echo -e "\n" >>$O

    echo "deleted $ndelete files" >>$O
    echo "deleted $ndelete files"

    cd $oldd
done
