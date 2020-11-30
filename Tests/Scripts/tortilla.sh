#!/bin/sh
#BSUB -o siesta.out # output file
#BSUB -e siesta.err # error file
#BSUB -q extra   # queue
#BSUB -J prueba  # name of the job
#BSUB -n 4       # number of processors

prog=../../../siesta

make check tests="h2o_libomm" SIESTA="mpirun $prog" > Out.dat




