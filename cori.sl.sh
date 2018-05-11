#!/bin/bash -l

#SBATCH --nodes=2
#SBATCH --time=00:30:00
#SBATCH --partition=regular
#SBATCH --license=SCRATCH   #note: specify license need for the file systems your job needs, such as SCRATCH,project
#SBATCH --constraint=haswell
#SBATCH --mail-user=rahuliitk@gmail.com
 
srun -n 32 -c 4 ./PIC
