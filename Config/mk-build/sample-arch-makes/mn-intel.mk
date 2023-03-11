# Example fortran.mk include file
# MareNostrum IV at BSC, Intel suite
#
SIESTA_ARCH=prace-intel
#
# You need to have the appropriate modules
# (intel/2017.4 environment)
# ml netcdf
# ml flook
# ml elpa
# ml fftw
#
#--------------------------------------------------------
# Use these symbols to request particular features
# To turn on, set '=1'.
#--------------
# These are set in build.mk
# WITH_PSML=1
# WITH_GRIDXC=1
#-------------
#
WITH_CMAKE_LIB_IS_LIB64=1
WITH_ELPA=1
WITH_FLOOK=1
WITH_MPI=1
WITH_NETCDF=1
WITH_NCDF=1
WITH_LEGACY_GRIDXC_INSTALL=0
WITH_GRID_SP=0
#
#===========================================================
# Make sure you have the appropriate library symbols
# (Either explicitly here, or through shell variables, perhaps
#  set by a module system)
# Define also compiler names and flags
#--------------------------------------------------------
#XMLF90_ROOT=
#PSML_ROOT=
#GRIDXC_ROOT=
#ELSI_ROOT=
#ELPA_ROOT=
#ELPA_INCLUDE_DIRECTORY=
#FLOOK_ROOT=
#-------------------------------------------------------
# For netCDF, they have all the dependencies (hdf5, etc) encoded in
# the fortran library, so it is only necessary to specify NETCDF_ROOT (with 'ml netcdf')
# BUT the spec has to be -L$(NETCDF_ROOT)/lib -lnetcdff
# An explicit (static) library load will not work
#
SCALAPACK_LIBS=-L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
LAPACK_LIBS=-L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
FFTW_ROOT=/apps/FFTW/3.3.6/INTEL/IMPI
# Needed for PEXSI (ELSI) support
LIBS_CPLUS=-lstdc++ 
#--------------------------------------------------------
#
# Define compiler names and flags
#
FC_PARALLEL=mpiifort
FC_SERIAL=ifort
#
FPP = $(FC_SERIAL) -E -P -x c
#
# (add -qopenmp to all three for OpenMP support)
#
# Warning: Some of these 'safety' options might impact performance
# 
FFLAGS = -g -traceback -O2 -prec-div -prec-sqrt -fp-model source
FFLAGS_DEBUG= -g -traceback -O1 -prec-div -prec-sqrt -fp-model source
LDFLAGS= 
#
# Delicate files with Intel compiler ----------------------------
#
atom.o: atom.F
	$(FC) -c $(FFLAGS_DEBUG) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) $<
state_analysis.o: state_analysis.F
	$(FC) -c $(FFLAGS_DEBUG) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) $<
create_Sparsity_SC.o: create_Sparsity_SC.F90
	$(FC) -c $(FFLAGS_DEBUG) $(INCFLAGS) $(FPPFLAGS)  $<
#
RANLIB=echo
# ----------------------------------------------------------
