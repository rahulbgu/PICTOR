#!/bin/bash
#SBATCH -N 19  #Node count
#SBATCH -n 512 #core count
#SBATCH --mem-per-cpu=4000
#SBATCH -t 24:00:00
 
module load openmpi/gcc/1.10.2/64
#module load openmpi
#cd /tigress/rk10/Projects/AlfvenTurb/1D/run1
 
srun ./PIC
