# Siesta/Util/SiestaSubroutine/SimpleTest/Makefile
# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
# Makefile for Client-Server example
#
# JM Soler, A Garcia, E Anglada
#
#  FIXME: Pending appropriate treament of libfdict and libncdf
#
TOPDIR=.
MAIN_OBJDIR=.
#
.SUFFIXES:
.SUFFIXES: .c .f .F .o .a .f90 .F90
#
#
#  simple_pipes_serial     : Simple serial test with pipes version
#  simple_pipes_parallel   : Simple parallel test with pipes version
#  simple_sockets_serial   : Simple serial test with sockets version
#  simple_sockets_parallel : Simple parallel test with sockets version
#  simple_mpi_serial       : Simple serial test with MPI version
#  simple_mpi_parallel     : Simple parallel test with MPI version
#
#  Do not move default list from here, or else a default target from arch.make can
#  interfere
all: default
#
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
include $(MAIN_OBJDIR)/check_for_build_mk.mk
#
VPATH:=$(TOPDIR)/Util/SiestaSubroutine/SimpleTest/Src:$(TOPDIR)/Src
#
PROGS:= simple_pipes_serial simple_sockets_serial simple_mpi_serial
ifeq ($(WITH_MPI),1)
  PROGS+= simple_pipes_parallel simple_sockets_parallel simple_mpi_parallel
  PROGS+= phonons
endif
default: $(PROGS)
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
# Routines used by pipes and sockets versions (VPATH knows where they are)
#
FSIESTA_PIPES= posix_calls.o fsiesta_pipes.o pxf.o 
FSIESTA_SOCKETS= posix_calls.o fsockets.o sockets.o fsiesta_sockets.o 
#
# Siesta libraries used by MPI version
#
SIESTA_LIB=$(MAIN_OBJDIR)/Src/libSiestaForces.a

XMLPARSER=$(XMLF90_LIBS)
XC=$(GRIDXC_LIBS)
PSML=$(PSML_LIBS)
#
INCFLAGS+=$(NETCDF_INCFLAGS) $(FDF_INCFLAGS) $(MS_INCFLAGS) \
          $(NCPS_INCFLAGS) $(PSOP_INCFLAGS) $(PSML_INCFLAGS) \
          $(XMLF90_INCFLAGS) $(GRIDXC_INCFLAGS)
INCFLAGS+= $(MPI_WRAPPERS_INCFLAGS)
INCFLAGS+= $(DFTD3_INCFLAGS)
INCFLAGS+= $(SIESTA_LIB_INCFLAGS)
INCFLAGS+= $(LIBSYS_INCFLAGS)
#
ALL_LIBS= $(SIESTA_LIB) $(FDF_LIBS) $(XC) $(LIBSYS) \
	$(NCPS) $(PSML) $(XMLF90_LIBS) $(PSOP) $(MS) \
	$(MPI_INTERFACE) $(COMP_LIBS) $(DFTD3_LIBS) $(LIBS)
#
# Pipes version
#
simple_pipes_serial: $(FSIESTA_PIPES) simple_serial.o
	$(FC) $(LDFLAGS) -o simple_pipes_serial $(FSIESTA_PIPES) simple_serial.o 
#
simple_pipes_parallel: $(FSIESTA_PIPES) simple_parallel.o
	$(FC) $(LDFLAGS) -o simple_pipes_parallel \
	      $(FSIESTA_PIPES) simple_parallel.o
#
# Sockets version
#
simple_sockets_serial: $(FSIESTA_SOCKETS) simple_serial.o
	$(FC) $(LDFLAGS) -o simple_sockets_serial \
	      $(FSIESTA_SOCKETS) simple_serial.o 
#
simple_sockets_parallel: $(FSIESTA_SOCKETS) simple_parallel.o
	$(FC) $(LDFLAGS) -o simple_sockets_parallel \
	      $(FSIESTA_SOCKETS) simple_parallel.o
#
# MPI version
#
phonons: $(SIESTA_LIB) phonons.o
	$(FC) $(LDFLAGS) -o phonons phonons.o  $(ALL_LIBS)
#
simple_mpi_serial: $(SIESTA_LIB) simple_serial.o 
	$(FC) $(LDFLAGS) -o simple_mpi_serial simple_serial.o  $(ALL_LIBS)
#
simple_mpi_parallel: $(SIESTA_LIB) simple_mpi_parallel.o 
	$(FC) $(LDFLAGS) -o simple_mpi_parallel \
	      simple_mpi_parallel.o  $(ALL_LIBS)
#
intermediate: 
	@echo "==> Cleaning intermediate files ONLY"
	@echo "==> Use make clean to remove executable files"
	rm -f *.a *.o *.mod 
#
pristine: clean
clean: 
	@echo "==> Cleaning all intermediate and executable files"
	rm -f simple_pipes_serial simple_pipes_parallel 
	rm -f simple_sockets_serial simple_sockets_parallel 
	rm -f simple_mpi_serial simple_mpi_parallel
	rm -f phonons
	rm -f *.a *.o *.mod

install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin

