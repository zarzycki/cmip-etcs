#!/bin/bash -l

#PBS -N forcing-gen 
#PBS -A UPSU0032
#PBS -l select=1:ncpus=4:mpiprocs=4:mem=380GB
#PBS -l walltime=24:00:00
#PBS -q casper

module load ncarenv
module load python/3.7.12
ncar_pylib

time python ptype-driver2.py


