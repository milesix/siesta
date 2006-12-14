# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996-2006.
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
SIESTA_ARCH=macosx-xlf
#
# IBM XL Fortran Compiler on MacOS X  (beta version)
#
# Issues addressed: 
#
#      - Case-insensitive filesystem
# Other issues:
#
#      ONE HAS TO BE EXTRA CAREFUL WHEN DELETING FILES, DUE TO THE
#      CASE-INSENSITIVENESS OF THE FILE SYSTEM ON THE MAC.
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90
#
FC=xlf
FC_ASIS=$(FC)
FREE_F90=-qsuffix=f=f90 -qfree=f90
FREE_F90_CPP=-qsuffix=cpp=F90 -qfree=f90
#
RANLIB=ranlib
#
FFLAGS=-O3 -qarch=auto -qtune=auto -qcache=auto -qnolm
FFLAGS_DEBUG= -g -O0
LDFLAGS=
COMP_LIBS=
#
NETCDF_LIBS=          #-L/opt/lib -lnetcdf
NETCDF_INTERFACE=     # libnetcdf_f90.a
DEFS_CDF=             # -WF,-DCDF
FFLAGS_NETCDF= -qsuffix=f=f90:cpp=F90 
#
MPI_INTERFACE=
MPI_INCLUDE=
DEFS_MPI=
#
DEFS_EXTRA= -WF,-DFC_HAVE_ABORT
#
LINEAR_ALGEBRA_LIBS=-Wl,-framework -Wl,vecLib
#
LIBS=  $(MPI_LIBS)  $(LINEAR_ALGEBRA_LIBS) $(NETCDF_LIBS) 
COMP_LIBS=dc_lapack.a
SYS=xlf
DEFS= $(DEFS_CDF) $(DEFS_MPI) $(DEFS_EXTRA)
#
# Beta compiler crashes on these, or takes a long time...
#
parsing.o:
	$(FC) -c $(FFLAGS_DEBUG)  parsing.f
siesta.o:
	$(FC) -c  $(FFLAGS_DEBUG)  $(DEFS)  siesta.F
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FREE_F90_CPP) $(FFLAGS) $(INCFLAGS)  $(DEFS) $<
.f90.o:
	$(FC) -c $(FREE_F90) $(FFLAGS) $(INCFLAGS)   $<
#









