#!/bin/bash  
#SBATCH --job-name="PIC"  
#SBATCH --output="PIC.%j.%N.out"  
#SBATCH --partition=compute  
#SBATCH --nodes=72  
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=5000  
#SBATCH --export=ALL  
#SBATCH -t 48:00:00  

#This job runs with 2 nodes, 24 cores per node for a total of 48 cores.  
#ibrun in verbose mode will give binding detail  

ibrun -v ../PIC