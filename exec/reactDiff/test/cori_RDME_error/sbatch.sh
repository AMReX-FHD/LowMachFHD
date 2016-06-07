#!/bin/bash -l

#SBATCH -p debug 
#SBATCH -N 1
#SBATCH -t 00:3:00 
#SBATCH -J TEST   

srun -n 32 ../main.Linux.Intel.mpi.exe inputs 
