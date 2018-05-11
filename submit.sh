#! /bin/bash
#$ -V
#$ -cwd
#$ -N test.sh 
#$ -pe rmpi 40 
#$ -j y
#$ -e test.e
#$ -o test.o
export MPIRUN=/storage/openmpi-1.5_openib/bin/mpirun
export PATH=/storage/openmpi-1.5_openib/bin:$PATH
export LD_LIBRARY_PATH=/storage/openmpi-1.5_openib/lib:$LD_LIBRARY_PATH

machinefile=$TMPDIR/nodes
awk '{ for (i=0;i<$2;++i) {print $1} }' $PE_HOSTFILE >> $machinefile

$MPIRUN --hostfile $machinefile -np $NSLOTS PIC
