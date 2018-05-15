#PBS -S /bin/bash
#PBS -l size=2048,walltime=24:00:00

cd $PBS_O_WORKDIR
aprun -n 2048 ./pictor.x
