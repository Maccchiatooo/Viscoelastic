#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=200:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=2   # number of nodes
#SBATCH --ntasks=32   # 16 processor core(s) per node 
#SBATCH --job-name="complete"
#SBATCH --mem=64GB

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load intel/2021.1
#compile the code
mpiifort cart_mpi.f90 lbm.f90 main.f90 -o czhao29
# RUN the work
mpirun -n 32 ./czhao29
