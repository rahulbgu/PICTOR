#!/bin/bash -l
#SBATCH -N 2
#SBATCH -p debug
#SBATCH -C knl,quad,cache
#SBATCH -t 30:00
#SBATCH -L SCRATCH
 
export OMP_NUM_THREADS=1 # only needed for hybrid MPI/OpenMP codes built with "-qopenmp" flag

srun -n 136 ./PICTOR