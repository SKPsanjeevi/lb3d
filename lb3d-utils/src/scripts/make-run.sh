################################################################################
#!/bin/bash

# /---------------------------------------------------------------------------\
# |                                                                           |
# |     HELP ON CUSTOMIZATION:                                                |
# |                                                                           |
# |     To use this script, personalize the variables below.                  |
# |     - EMAIL will be used to set notify_user in LL directives              |
# |     - BGSIZE will be used to set bgsize in LL directives                  |
# |     - LBEFILE is the name of the executable, and will be used to generate |
# |       jobnames                                                            |
# |     - INPUTBASE will be passed to LB3D as -f <INPUTBASE>                  |
# |     - INPUTDIFF will be passed to LB3D as -d <INPUTDIFF>                  |
# |     - CODEDIR specifies where the original LB3D binary can be found       |
# |     - ROCKDIR specifies where the rock files can be found                 |
# |                                                                           |
# |     WARNING: One should always make sure that BGSIZEDEF, INPUTBASEDEF,    |
# |       LBEFILEDEF match BGSIZE, INPUTBASE, LBEFILE, respectively.          |
# |     WARNING: Use double quotes in the definition of LBEFILE and           |
# |       INPUTBASE .                                                         |
# |                                                                           |
# |     HELP ON USE:                                                          |
# |                                                                           |
# |     Run ./make-run.sh (possibly with optional arguments):                 |
# |       -b <BGSIZE>    : overrides default BGSIZE                           |
# |       -d <INPUTDIFF> : overrides default INPUTDIFF                        |
# |       -f <INPUTBASE> : overrides default INPUTBASE                        |
# |       -l <LBEFILE>   : overrides default LBEFILE                          |
# |     This will generate a file 'run-<JOBNAME>.ll', which can then be       |
# |     submitted to the LoadLeveler.                                         |
# |                                                                           |
# \---------------------------------------------------------------------------/

# LoadLeveler variables
EMAIL='my@email.address'
BGSIZE=512

# Directories / files
LBEFILE="lbe"
INPUTBASE="input-file"
INPUTDIFF=
CODEDIR='/PATH/TO/EXECUTABLE'
ROCKDIR='/PATH/TO/ROCKS'

# Defaults
BGSIZEDEF=512
INPUTBASEDEF="input-file"
LBEFILEDEF="lbe"

# Should not need to change
WAITSTEP=60
MINDELTA=3 # Minimum allowed time between two runs

# Automatic:
INPUTDIR=$( pwd )
SIMDIR=$INPUTDIR
SELF=$0
SUFFIX=${SELF:(-3)}
NP=$[8*BGSIZE]

# /---------------------------------------------------------------------------\
# |                                                                           |
# |     Helper functions                                                      |
# |                                                                           |
# \---------------------------------------------------------------------------/

function checkdir
{
    DIRPROMPT=$1
    DIRNAME=$2
    if [[ -d $DIRNAME ]]; then
        echo "  Checked $DIRPROMPT = '$DIRNAME'. "
    else
        if [[ $DIRNAME == $SIMDIR ]]; then
            mkdir -p $SIMDIR
            echo "  NOTIFY: $DIRPROMPT = '$DIRNAME' was created. "
        else
            mailmessage "ERROR" "Directory '$DIRPROMPT' is missing. "
            exit -1
        fi
    fi
}

function mailmessage
{
    echo "  $1 : $2"
    echo "  $LOADL_STEP_ID : $2 " | mail -s "$1 : $LOADL_JOB_NAME " '$EMAIL'
    echo "    *** Mail sent to '$EMAIL'. "
}

function checkfile
{
    FILEPROMPT=$1
    FILENAME=$2
    if [[ -e $FILENAME ]]; then
        echo "  Checked $FILEPROMPT = '$FILENAME'."
    else
        mailmessage "ERROR" "File '$FILEPROMPT' is missing. "
        exit -1
    fi
}

function runcount
{
    if [[ ! -e RUNS ]]; then
        RUNCOUNT=0
    else
        RUNCOUNT=$( cat RUNS )
    fi
    let RUNCOUNT=RUNCOUNT+1
    echo $RUNCOUNT > RUNS
    RUNCOUNTLAST=${RUNCOUNT:(-1)}
    case $RUNCOUNTLAST in
        1) RUNCOUNTORDINAL=`echo -n $RUNCOUNT; echo -n 'st'`;;
        2) RUNCOUNTORDINAL=`echo -n $RUNCOUNT; echo -n 'nd'`;;
        3) RUNCOUNTORDINAL=`echo -n $RUNCOUNT; echo -n 'rd'`;;
        *) RUNCOUNTORDINAL=`echo -n $RUNCOUNT; echo -n 'th'`;;
    esac
    echo "  Running simulation for the $RUNCOUNTORDINAL time..."
}

function submit
{
    runcount
    if [[ $( ls $OUTPUTFOLDERCP | grep checkparams ) ]]; then
        RESTORE=$( ls -tr $OUTPUTFOLDERCP/checkparams*|tail -1| sed 's/.*\(t[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]-[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]\).*/\1/' )
        echo "  Restarting simulation from checkpoint $RESTORE... "
        if [ -n "$INPUTDIFF" ]; then
            runjob --exe ./$LBEFILE --args -f $INPUTBASE --args -d $INPUTDIFF --args -r $RESTORE --ranks-per-node 8 & JOBID="$!"
        else
            runjob --exe ./$LBEFILE --args -f $INPUTBASE --args -r $RESTORE --ranks-per-node 8 & JOBID="$!"
        fi
    else
        echo "  Starting new simulation... "
        if [ -n "$INPUTDIFF" ]; then
            runjob --exe ./$LBEFILE --args -f $INPUTBASE --args -d $INPUTDIFF --ranks-per-node 8 & JOBID="$!"
        else
            runjob --exe ./$LBEFILE --args -f $INPUTBASE --ranks-per-node 8 & JOBID="$!"
        fi
    fi
    echo $TIMESTAMP > LASTRUN
}

function checkfinish
{
    FINISH=1
    ps -p $JOBID 1> /dev/null 2> /dev/null
    if [[ $? -eq 0 ]]; then
        FINISH=0
        #echo -n '.'
        let WAITSTEPS=WAITSTEPS+1
        #if [[ $(( WAITSTEPS%60 )) == 0 ]] ; then
        #    echo -ne "\n "
        #fi
    else
        echo -ne "\n  Process $JOBID is not running... "
        FINISH=1
    fi
}

function waitforfinish
{
    echo -ne "  Waiting for process $JOBID to finish...\n "
    echo
    while true
    do
        sleep $WAITSTEP
        checkfinish
        if [[ $FINISH -gt 0 ]]; then
            echo "  Exiting main."
            return 0
        fi
    done
}

function resubmit
{
    ME=$( basename $0 )
    echo "  Resubmitting $INPUTDIR/$ME..."
    llsubmit $INPUTDIR/$ME
}

function postprocess
{
    echo "  Postprocess (dummy function). "
}

function checkexit
{
    ERRFILE=`echo $LOADL_STEP_ERR | sed 's/\.2\.out/\.1\.out/'` # both stdout and stderr of jobstep 2 are written to .out file
    echo "  Attempting to find exit status of simulation step by analyzing '$ERRFILE'. "
    if [ -e "$ERRFILE" ]; then
        ECODE=`grep -e "ibm.runjob.LogSignalInfo: received signal [0-9]*" $ERRFILE | sed 's/.*ibm.runjob.LogSignalInfo: received signal[ \t]*\([0-9]*\).*/\1/' `
        case "$ECODE" in
            0)
                echo "    ??? (usually means simulation complete). "
                ;;
            24)
                echo "    Received signal 24 (usually means wallclock limit has been reached). "
                ;;
            134)
                mailmessage "ERROR" "Simulation step ended with exit status 134. Aborting..."
                exit -1
                ;;
            143)
                mailmessage "ERROR" "Simulation step ended with exit status 143. Aborting... "
                exit -1
                ;;
            *)
                mailmessage "WARNING" "Simulation step ended with unknown exit status '$ECODE'. "
                ;;
        esac
    else
        mailmessage "WARNING" "Could not find error file '$ERRFILE'. "
    fi
}

function cleancp
{
    for i in `ls -tr $OUTPUTFOLDERCP/checkparams* | head -n $KEEPCP`
    do
        CPID=`echo $i | sed 's/.*\(t[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]-[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]\).*/\1/'`
        echo "  Cleaning checkpoint files for checkpoint '$CPID'... "
        for f in $OUTPUTFOLDERCP/*check*$CPID*
        do
            rm $f
        done
    done

    PREVIOUSWD=`pwd`
    cd $OUTPUTFOLDERCP
    for f in checkpoint*
    do
        CPARAM=`echo $f | sed -e 's/checkpoint/checkparams/' | sed -e 's/_p[0-9][0-9][0-9][0-9][0-9][0-9]\.xdr/\.xdr/'`
        if [ ! -f $CPARAM ]; then
            mailmessage "WARNING" "Found checkpoint file $OUTPUTFOLDERCP/$f without matching checkparams file, some manual cleanup might be required. No further warnings on this subject will be generated. "
            #rm $f
            break
        fi
    done
    cd $PREVIOUSWD
}

#function checkstatus
#{
#    if [[ $( ls $OUTPUTFOLDER | grep .h5 ) ]]; then
#        LASTTIMESTEP=$( ls -rt $OUTPUTFOLDER/*h5 | tail -1 | sed 's/.*t\([0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]\).*/\1/' | awk '{printf "%d\n", $1}')
#        if [ -e "$INPUTDIR/$INPUTBASE" ]; then
#            ASPIREDTIMESTEP=$( cat $INPUTDIR/$INPUTBASE | grep -e '^n_iteration' | cut -d '=' -f 2 | sed -e s/" "//g )
#            if [ -e "$INPUTDIR/$INPUTDIFF" ]; then
#                ASPIREDTIMESTEPDIFF=$( cat $INPUTDIR/$INPUTDIFF | grep -e '^n_iteration' | cut -d '=' -f 2 | sed -e s/" "//g )
#                if [ -n "$ASPIREDTIMESTEPDIFF" ]; then
#                    ASPIREDTIMESTEP=$ASPIREDTIMESTEPDIFF
#                fi
#            fi
#            echo "  Timestep $LASTTIMESTEP was reached, target timestep is $ASPIREDTIMESTEP..."
#            if [[ $LASTTIMESTEP == $ASPIREDTIMESTEP ]]; then
#                echo "  Simulation $SIMNAME is done."
#                DONEFILE="DONE.$JOBNAME"
#                echo $TIMESTAMP > $DONEFILE
#            fi
#        else
#            mailmessage "WARNING" "Checkstatus found no input file - can't grep n_iteration. "
#        fi
#    else
#        mailmessage "WARNING" "Checkstatus found no output files. "
#    fi
#

function checkstatus
{
    if [ -f "$INPUTDIR/$INPUTBASE" ]; then
        ASPIREDTIMESTEP=$( cat $INPUTDIR/$INPUTBASE | grep -e '^n_iteration' | cut -d '=' -f 2 | sed -e s/" "//g )
        if [ -f "$INPUTDIR/$INPUTDIFF" ]; then
            ASPIREDTIMESTEPDIFF=$( cat $INPUTDIR/$INPUTDIFF | grep -e '^n_iteration' | cut -d '=' -f 2 | sed -e s/" "//g )
            if [ -n "$ASPIREDTIMESTEPDIFF" ]; then
                ASPIREDTIMESTEP=$ASPIREDTIMESTEPDIFF
            fi
        fi
        echo "  Target timestep is $ASPIREDTIMESTEP..."
    else
        mailmessage "WARNING" "Checkstatus found no input file - can't grep n_iteration, continuing run. Manual abort might be required. "
    fi

    if [ -n "$ASPIREDTIMESTEP" ]; then
        TSSTR=`printf "t%08d" $ASPIREDTIMESTEP`
	TSCNT=`find $OUTPUTFOLDER -name "*$TSSTR*" -print -quit 2>/dev/null |wc -l`
        if [ "$TSCNT" -gt 0 ]; then
            echo "  Found output files for target timestep $ASPIREDTIMESTEP. Simulation is complete. "
            mailmessage "FINISHED" "Simulation has completed successfully. "
            DONEFILE="DONE.$JOBNAME"
            echo $TIMESTAMP > $DONEFILE
        else
            echo "  Output files for target timestep $ASPIREDTIMESTEP not found. Continuing... "
        fi
    else
        mailmessage "WARNING" "ASPIREDTIMESTEP not detected, continuing run. Manual abort might be required. "
    fi
}

function checktimedelta
{
    LASTRUN=$( cat LASTRUN )
    DELTA=$(( $TIMESTAMP-$LASTRUN ))
    echo "  Checking time difference: $DELTA = $TIMESTAMP - $LASTRUN ."
    if [[ $DELTA -lt $MINDELTA ]]; then
        mailmessage "ERROR" "Last run started only $DELTA seconds ago... exiting queue. "
        exit -1
    fi
}

function printheader
{
    printhr
    echo    "  STARTING JOBSTEP..."
    echo
    echo -e "  Job: $JOBNAME\t Step: $LOADL_STEP_NAME \t Start: $STARTTIME"
}

function printhr
{
    echo
    echo ' -------------------------------------------------------------------------------'
    echo
}


# /---------------------------------------------------------------------------\
# |                                                                           |
# |     Some pre-processing                                                   |
# |                                                                           |
# \---------------------------------------------------------------------------/

printhr

echo "  Parsing command line arguments... "
while getopts b:d:f:l: op
do case "$op" in
        b)
            echo "    Setting BGSIZE to $OPTARG from command line."
            BGSIZE=$OPTARG
            ;;
        d)
            echo "    Setting INPUTDIFF to $OPTARG from command line."
            INPUTDIFF=$OPTARG
            ;;
        f)
            echo "    Setting INPUTBASE to $OPTARG from command line."
            INPUTBASE=$OPTARG
            ;;
        l)
            echo "    Setting LBEFILE to $OPTARG from command line."
            LBEFILE=$OPTARG
            ;;
    esac
done

echo
echo "  Pre-check..."

if [[ -e $INPUTBASE ]]; then
    echo "  Checked INPUTBASE = '$INPUTBASE'."
else
    echo "  ERROR: INPUTBASE file '$INPUTBASE' is missing, cannot read parameters. "
    exit -1
fi

ROCKFILE=$( cat $INPUTDIR/$INPUTBASE | grep '^obs_file' | sed "s/\"/'/g" | cut -d "'" -f 2 )
OUTPUTFOLDER=$( cat $INPUTDIR/$INPUTBASE | grep '^folder' | sed "s/\"/'/g" | cut -d "'" -f 2 )
OUTPUTFOLDERCP=$( cat $INPUTDIR/$INPUTBASE | grep '^cpfolder' | sed "s/\"/'/g" | cut -d "'" -f 2 )
KEEPCP=$( cat $INPUTDIR/$INPUTBASE | grep '^num_chkp_files' | cut -d "=" -f 2 )
CDX=$( cat $INPUTDIR/$INPUTBASE | grep '^cdx' | cut -d "=" -f 2 )
CDY=$( cat $INPUTDIR/$INPUTBASE | grep '^cdy' | cut -d "=" -f 2 )
CDZ=$( cat $INPUTDIR/$INPUTBASE | grep '^cdz' | cut -d "=" -f 2 )
SIMNAME=$( cat $INPUTDIR/$INPUTBASE | grep '^gr_out_file' | cut -d "=" -f 2 | sed "s/'//" | sed "s/ //" )
if [ -n "$INPUTDIFF" ]; then

    if [[ -e $INPUTDIFF ]]; then
        echo "  Checked INPUTDIFF = '$INPUTDIFF'."
    else
        echo "  ERROR: INPUTDIFF file '$INPUTDIFF' is missing, cannot read parameters. "
        exit -1
    fi

    SIMNAME=$( echo $INPUTDIFF | sed -e s/input-// | sed -e s/-diff// )

    ROCKFILEDIFF=$( cat $INPUTDIR/$INPUTDIFF | grep '^obs_file' | sed "s/\"/'/g" | cut -d "'" -f 2 )
    if [ -n "$ROCKFILEDIFF" ]; then
        ROCKFILE=$ROCKFILEDIFF
    fi
    OUTPUTFOLDERDIFF=$( cat $INPUTDIR/$INPUTDIFF | grep '^folder' | sed "s/\"/'/g" | cut -d "'" -f 2 )
    if [ -n "$OUTPUTFOLDERDIFF" ]; then
        OUTPUTFOLDER=$OUTPUTFOLDERDIFF
    fi
    OUTPUTFOLDERCPDIFF=$( cat $INPUTDIR/$INPUTDIFF | grep '^cpfolder' | sed "s/\"/'/g" | cut -d "'" -f 2 )
    if [ -n "$OUTPUTFOLDERCPDIFF" ]; then
        OUTPUTFOLDERCP=$OUTPUTFOLDERCPDIFF
    fi
    KEEPCPDIFF=$( cat $INPUTDIR/$INPUTDIFF | grep '^num_chkp_files' | cut -d "=" -f 2 )
    if [ -n "$KEEPCPDIFF" ]; then
        KEEPCP=$KEEPCPDIFF
    fi
    CDXDIFF=$( cat $INPUTDIR/$INPUTDIFF | grep '^cdx' | cut -d "=" -f 2 )
    if [ -n "$CDXDIFF" ]; then
        CDX=$CDXDIFF
    fi
    CDYDIFF=$( cat $INPUTDIR/$INPUTDIFF | grep '^cdy' | cut -d "=" -f 2 )
    if [ -n "$CDYDIFF" ]; then
        CDY=$CDYDIFF
    fi
    CDZDIFF=$( cat $INPUTDIR/$INPUTDIFF | grep '^cdz' | cut -d "=" -f 2 )
    if [ -n "$CDZDIFF" ]; then
        CDZ=$CDZDIFF
    fi
fi

KEEPCP=$[(-1)*KEEPCP]

OUTPUTFOLDER="$INPUTDIR/$OUTPUTFOLDER"
OUTPUTFOLDERCP="$OUTPUTFOLDER/$OUTPUTFOLDERCP"

JOBNAME="$SIMNAME.$LBEFILE"

#echo
#echo "  DBG: ROCKFILE       = $ROCKFILE"
#echo "  DBG: OUTPUTFOLDER   = $OUTPUTFOLDER"
#echo "  DBG: OUTPUTFOLDERCP = $OUTPUTFOLDERCP"
#echo "  DBG: SIMNAME        = $SIMNAME"
#echo "  DBG: EMAIL          = $EMAIL"
#echo "  DBG: BGSIZE         = $BGSIZE"
#echo "  DBG: KEEPCP         = $KEEPCP"

TODAY=$( date +'%R %Y/%m/%d' )
STARTTIME=$( date +%R )
TIMESTAMP=$( date +%s )

RUNCOUNT=0
WAITSTEPS=0

# /---------------------------------------------------------------------------\
# |                                                                           |
# |     LoadLeveler script                                                    |
# |                                                                           |
# \---------------------------------------------------------------------------/
#
# @shell = /bin/bash
# @environment = COPY_ALL
# @job_name = $JOBNAME
# @notify_user = $EMAIL
#
#
# @wall_clock_limit = 00:05:00
# @job_type = serial
# @notification = error
# @output = $(job_name).$(jobid).$(stepid).out
# @error  = $(job_name).$(jobid).$(stepid).out
# @step_name = pre
# @queue
#
#
# @wall_clock_limit = 05:50:00
# @job_type = bluegene
# @bg_size = $BGSIZE
# @notification = error
# @output = $(job_name).$(jobid).$(stepid).out
# @error  = $(job_name).$(jobid).$(stepid).err
# @step_name = main
# @ bg_connectivity = torus
# make sure this step is only executed when 'pre' is finished
# @dependency = (pre == 42)
# @queue
#
#
# @wall_clock_limit = 00:45:00
# @job_type = serial
# @notification = error
# @output = $(job_name).$(jobid).$(stepid).out
# @error  = $(job_name).$(jobid).$(stepid).out
# @step_name = post
# make sure this step is only executed when 'main' is finished
# @dependency = (pre == 42 && main != 100)
# @queue
#

# /---------------------------------------------------------------------------\
# |                                                                           |
# |     Script generation                                                     |
# |                                                                           |
# \---------------------------------------------------------------------------/

if [ $SUFFIX == ".sh" ]; then

    SIMNAME=$( cat $INPUTDIR/$INPUTBASE | grep '^gr_out_file' | cut -d "=" -f 2 | sed "s/'//g" | sed "s/ //")    
    if [ -n "$INPUTDIFF" ]; then
         SIMNAME=$( echo $INPUTDIFF | sed -e s/input-// | sed -e s/-diff// )
    fi
    JOBNAME="$SIMNAME.$LBEFILE"
    LLFILE="run-$JOBNAME.ll"

    echo
    echo "  Hardcoding parameters... "
    echo "    BGSIZE      = $BGSIZE "
    echo "    INPUTDIFF   = $INPUTDIFF "
    echo "    INPUTBASE   = $INPUTBASE "
    echo "    LBEFILE     = $LBEFILE "

    cat $SELF | sed -e 's/^BGSIZE='$BGSIZEDEF'/BGSIZE='$BGSIZE'/' | sed -e 's/^INPUTDIFF=.*/INPUTDIFF=\"'$INPUTDIFF'\"/' | \
                sed -e 's/^INPUTBASE=\"'$INPUTBASEDEF'\"/INPUTBASE=\"'$INPUTBASE'\"/' | sed -e 's/^LBEFILE=\"'$LBEFILEDEF'\"/LBEFILE=\"'$LBEFILE'\"/' > $LLFILE.tmp 

    echo
    echo "  Hardcoding LL variables... "
    echo "    notify_user = $EMAIL "
    echo "    jobname     = $JOBNAME "
    echo "    bg_size     = $BGSIZE "

    cat $LLFILE.tmp | sed -e 's/$EMAIL/'$EMAIL'/' | sed -e 's/$JOBNAME/'$JOBNAME'/' | sed -e 's/$BGSIZE/'$BGSIZE'/' > $LLFILE
    rm $LLFILE.tmp
    chmod 755 $LLFILE

    echo
    echo "  Wrote LoadLeveler script to $LLFILE ."
    echo "    Use 'llsubmit $LLFILE' to start the process ."
    echo

    printhr

else

# /---------------------------------------------------------------------------\
# |                                                                           |
# |     Pre-processing                                                        |
# |                                                                           |
# \---------------------------------------------------------------------------/

case $LOADL_STEP_NAME in
    pre)

        printheader

        printhr

        checkdir  'INPUTDIR       ' $INPUTDIR
        checkdir  'SIMDIR         ' $SIMDIR
        checkdir  'OUTPUTFOLDER   ' $OUTPUTFOLDER
        checkdir  'OUTPUTFOLDERCP ' $OUTPUTFOLDERCP

        printhr

        checkfile 'INPUTBASE ' $INPUTDIR/$INPUTBASE
        checkfile 'INPUTDIFF ' $INPUTDIR/$INPUTDIFF

        if [[ -n "$CDX" && -n "$CDY" && -n "$CDZ" ]]; then
            CNP=$[CDX*CDY*CDZ]
            if [ "$CNP" -ne "$NP" ]; then
                mailmessage "ERROR" "Incorrect decomposition: $CDX * $CDY * $CDZ != $NP"
                exit -1
            fi
        else
            mailmessage "WARNING" "Not all three dimensions specified for cartesian decomposition. This might cause problems when restoring."
        fi

        cd $SIMDIR
        printhr

        echo "  Finished pre-processing step. "
        echo
        ;;

# /---------------------------------------------------------------------------\
# |                                                                           |
# |     Simulation                                                            |
# |                                                                           |
# \---------------------------------------------------------------------------/

    main)

        printheader

        printhr

        cd $SIMDIR
        submit
        waitforfinish

        printhr

        echo "  Finished submit step."
        echo
        ;;

# /---------------------------------------------------------------------------\
# |                                                                           |
# |     Post-processing                                                       |
# |                                                                           |
# \---------------------------------------------------------------------------/

    post)

        printheader

        printhr

        cd $SIMDIR
        checkexit
        checkstatus
        cleancp

        DONEFILE="DONE.$JOBNAME"

        if [[ ! -e $DONEFILE ]]; then
            checktimedelta
            resubmit
        else
            postprocess
        fi

        printhr

        echo "  Finished post-processing step."
        echo
        ;;

esac

exit 42

fi
