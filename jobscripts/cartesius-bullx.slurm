#!/bin/bash
#
# This jobfile is to be used with executables compiled
# using the defines.CARTESIUS-BULLX file (using bullx MPI).
# This configuration is recommended for large runs (n > 1024).
#
# Abvove line specifies which shell should process the 
# shell commands
# SLURM directives start with #SBATCH
# See 'man sbatch' for a complete list.
#
# In this job an MPI program is started.
# This job will run on 4096 CPUs, on thin nodes (p = normal). 
#
#SBATCH -p normal
#SBATCH -n 4096 
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

module unload mpi/impi
module load mpi/xmpi

srun --mpi=openmpi ./lbe -f input-file 

