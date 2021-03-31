#!/bin/bash
#
# Abvove line specifies which shell should process the 
# shell commands
# SLURM directives start with #SBATCH

# In this job an MPI program is started.
# This job will run on two nodes, 16 processes on each node, 
# making a total of 32 processes.

#SBATCH -N 16 
#SBATCH --ntasks-per-node=24
#
# Job will take at most 1 hour, 10 minutes and 20 seconds wallclock time

#SBATCH -t 18:59:59

cd $(pwd)

# and call the MPI program
# Note that the srun command is aware of the number of processes to start
# These are defined above with the -N and --ntasks-per-node flags
#

rm rockconf ../rocks/*
# sh angle.sh 45 
# ./createblock.x
mkdir -p Production
rm Production/*
srun $HOME/Desktop/lb3d/src/lbe -f input-default
tempvar=${PWD##*/}
cp in* Production/
h5tovtk Production/vel_*.h5
N=00000; for X in Production/vel_*vtk*; do mv $X Production/vel_$N.vtk; N=$(($N+100)); done
