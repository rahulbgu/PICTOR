#!/bin/bash -l
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --license=SCRATCH  #note: specify license need for the file systems your job needs, such as SCRATCH,project
 
srun -n 16 ./PIC