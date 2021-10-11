#!/bin/sh
##BSUB -o siesta.out # output file
##BSUB -e siesta.err # error file
##BSUB -q extra   # queue
##BSUB -J prueba  # name of the job
##BSUB -n 4       # number of processors
#SBATCH --ntasks=4 ##Number of processes, in this case 20. Maximum 48.
#SBATCH --time=01:00:00 ## Running time in HH:MM:SS in case you want to limit it. If this flag is removed, the wall time will be set to UNLIMITED.
#SBATCH --cpus-per-task=1 ##Number of cores per process. Better left untouched so that ntask sets the number of cores.
#SBATCH --mem-per-cpu=3Gb ##Memory per core. Given the machine, optimal would be somewhere in between 3.5 and 3.9 Gb, but it only accepts integers. If you want something else, use Mb instead. Maximum allocatable memory per job: 177Gb
#SBATCH --job-name=siesta  ## Job name

prog=../../../siesta

make check tests="h2o_blomm_dense" SIESTA="mpirun $prog" > Out.dat




