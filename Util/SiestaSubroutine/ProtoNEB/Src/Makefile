# Siesta/Util/SiestaSubroutine/SimpleTest/Makefile
# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
# Makefile for Siesta as Subroutine example
#
# JM Soler, J Nara, A Garcia, E Anglada
#
TOPDIR=.
MAIN_OBJDIR=.
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90
#
all: default
#
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
include $(MAIN_OBJDIR)/check_for_build_mk.mk
#
VPATH:=$(TOPDIR)/Util/SiestaSubroutine/ProtoNEB/Src  #:$(TOPDIR)/Src
#
PROGS:= protoNEB
default: $(PROGS)
#
XC=$(GRIDXC_LIBS)
PSML=$(PSML_LIBS)
#
# For compilation of local sources
#
INCFLAGS+= $(FDF_INCFLAGS) 
INCFLAGS+= $(MPI_WRAPPERS_INCFLAGS)
INCFLAGS+= $(DFTD3_INCFLAGS)
INCFLAGS+= $(SIESTA_LIB_INCFLAGS)
#
# The "Siesta library" is used by this MPI-based driver
SIESTA_LIB=$(MAIN_OBJDIR)/Src/libSiestaForces.a
#
# Link-time dependencies (could be handled by SIESTA_LIB logic...)
#
ALL_LIBS= $(SIESTA_LIB) $(FDF_LIBS)  $(XC) \
	$(NCPS) $(PSML) $(XMLF90_LIBS) $(PSOP) $(MS)\
	$(MPI_INTERFACE) $(COMP_LIBS) $(DFTD3_LIBS) $(LIBS)
#
protoNEB: $(SIESTA_LIB) protoNEB.o handlers.o
	$(FC) $(LDFLAGS) -o protoNEB \
	      protoNEB.o handlers.o $(ALL_LIBS)
#
intermediate: 
	@echo "==> Cleaning intermediate files ONLY"
	@echo "==> Use make clean to remove executable files"
	rm -f *.a *.o *.mod
#
pristine: clean
clean:
	@echo "==> Cleaning all intermediate and executable files"
	rm -f protoNEB
	rm -f *.a *.o *.mod 
#
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin

