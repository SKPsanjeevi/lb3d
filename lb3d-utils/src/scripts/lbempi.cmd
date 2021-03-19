###########################################################################
#
# This is the condor command script to run lbe with Condor and MPICH2.
# It calls "lbecondormp2script" as an executable.
# "lbecondormp2script" in turn sets up mpd and launches lbe on the
# requested number of CPUs. If checkpoints are available, we will
# automatically restart. If not, we start at t=0. 
#
# - Make sure "lbecondormp2script" is in your "initialdir".
# - Set inititaldir to your lbe/code directory or to where your jobs 
#   should start.
# - "machine_count" sets the number of CPUs to use. So far, jobs are 
#   restricted to a single machine, i.e. you should not use more CPUs 
#   than cores available in a single SMP box (2-4).
# - Set inputfile to your lbe input-file.
#
# Note: Condor MPI jobs have to be submitted on the dedicated scheduler.
# Currently, this is "nippon". Jobs will not start if submitted 
# somewhere else.
#
# Jens, 24.10.07
#
###########################################################################

universe = parallel
inputfile = input-file
initialdir = /data/data1/jens/lbestuff/lbe/branches/version5/lbe/code
machine_count = 2
#output = lbecondor.$$.out
output = $(inputfile).$(cluster).out
error = $(inputfile).$(cluster).err
log = $(inputfile).$(cluster).log

# You should not need to change anything below.

executable = lbecondormp2script
arguments = lbe
should_transfer_files = yes
when_to_transfer_output = on_exit
+WantParallelSchedulingGroups = True
+WantIOProxy=True
environment = INITDIR=$(initialdir);INPUTFILE=$(inputfile)
queue
