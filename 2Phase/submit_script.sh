#!/bin/bash
#SBATCH --account=nn8079k

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=1:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=8   # number of nodes
#SBATCH --ntasks=128   # 16 processor core(s) per node 
#SBATCH --job-name="complete"

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load intel/2021b
#compile the code
mpiifort cart_mpi.f90 lbm.f90 main.f90 -o czhao29
# RUN the work
mpirun -n 128 ./czhao29
