#!/bin/bash
#
# This jobfile is to be used with executables compiled
# using the defines.CARTESIUS file (using Intel MPI).
#
# Abvove line specifies which shell should process the 
# shell commands
# SLURM directives start with #SBATCH
# See 'man sbatch' for a complete list.
#
# In this job an MPI program is started.
# This job will run on 1024 CPUs, on thin nodes (p = normal). 
#
#SBATCH -p normal
#SBATCH -n 1024 
#
# Job will take at most 1 hour, 10 minutes and 20 seconds wallclock time
#
#SBATCH -t 1:10:20
#
# Specify a job name
#
#SBATCH -J CARTESIUS_job 
#
# Get mail if anything happens
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.lastname@your.site.edu
#
# Here the shell script starts. 
# Note that the srun command is aware of the number of processes to start
# This is defined above with the -n flag.
#
# Use 'sbatch <jobfile>' to submit.

module load fftw2/intel/2.1.5

srun ./lbe -f input-file 

