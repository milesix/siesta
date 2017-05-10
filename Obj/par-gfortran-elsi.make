# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#
SIESTA_ARCH=gfortran-openmpi-elsi
# The only thing you should change is the location of the libraries
# on your computer
#
FC=mpif90
FC_SERIAL=gfortran
#
FC_ASIS=$(FC)
#
FFLAGS= -g -O2 # -Wall
FFLAGS_CHECKS= -O0 -g -fcheck=all
#FFLAGS= -O0 -g -fcheck=all
FFLAGS_DEBUG= -g -O0
RANLIB=echo
COMP_LIBS=
#
NETCDF_ROOT=/usr/local
NETCDF_INCFLAGS=-I$(NETCDF_ROOT)/include
FPPFLAGS_CDF=-DCDF
NETCDF_LIBS= -L$(NETCDF_ROOT)/lib -lnetcdff
#
# Note: ELSI library compiled without PEXSI for this version
#
ELSI_HOME=$(HOME)/code/ELSI/elsi-interface
ELSI_INCFLAGS = -I${ELSI_HOME}/include
ELSI_LIB = -L${ELSI_HOME}/lib -lelsi -lOMM -lMatrixSwitch \
                -lpspblas -lelpa -ltomato
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=.      # Note . for no-op
FPPFLAGS_MPI=-DMPI

LIBS= $(ELSI_LIB) -L/opt/scalapack/openmpi-1.6.1-gfortran/lib \
        -lscalapack -ltmg -lreflapack -lrefblas \
      $(NETCDF_LIBS)

SYS=nag
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI) -DF2003 
#
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#
