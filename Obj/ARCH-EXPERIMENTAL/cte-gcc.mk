#
# Defs for CTE-Power at BSC. gcc 0.6.4 compilers
#
# Need to:
#
#              ml load scalapack lapack netcdf fftw
#
NETCDF_LIBS=-L$(NETCDF_HOME)/lib -lnetcdff -lnetcdf -L/apps/HDF5/1.8.20/GCC/OPENMPI/lib -lhdf5_hl -lhdf5 -lcurl -lz
NETCDF_INCFLAGS=-I$(NETCDF_HOME)/include
#
SCALAPACK_LIBS=-lscalapack
LAPACK_LIBS=-llapack -lblas
FFTW_ROOT=/apps/FFTW/3.3.8/GCC/OPENMPI/
# Needed for PEXSI (ELSI) support
LIBS_CPLUS=-lstdc++ -lmpi_cxx
#--------------------------------------------------------
#
# Define compiler names and flags
#
FC_PARALLEL=mpif90
FC_SERIAL=gfortran
FPP = $(FC_SERIAL) -E -P -x c
FFLAGS = -O2 
FFLAGS_DEBUG= -g -O0
RANLIB=echo
