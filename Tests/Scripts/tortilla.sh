#!/bin/sh
#BSUB -o siesta.out # output file
#BSUB -e siesta.err # error file
#BSUB -q extra   # queue
#BSUB -J prueba  # name of the job
#BSUB -n 4       # number of processors

echo Running on `hostname`
##echo This jobs runs on the following processors:
##echo `cat $PBS_NODEFILE`
##NPROCS=`wc -l < $PBS_NODEFILE`
##echo This job has allocated $NPROCS nodes
#

##cd $PBS_O_WORKDIR

prog=../../../siesta

make check tests="h2o_libomm" SIESTA="mpirun $prog" > Out.dat


## Maybe:    
## make SIESTA="mpirun --gm-kill 600 $prog"



