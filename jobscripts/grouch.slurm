#!/bin/bash
#SBATCH -o job.%j.%N.out
#SBATCH -J TestJob
#SBATCH --get-user-env
#SBATCH -n 4
#SBATCH --time=08:00:00

module load gcc
module load openmpi
module load hdf5

mpirun lbe -f input-param
