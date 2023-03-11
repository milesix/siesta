# 
# Copyright (C) 1996-2021       The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#
# Makefile for PLSTM 
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a .f90 .F90

TOPDIR=.
MAIN_OBJDIR=.

PROGS:= plstm plsts
all: $(PROGS)

VPATH=$(TOPDIR)/Util/STM/simple-stm:$(TOPDIR)/Src:$(TOPDIR)/Util/Grid

#
# The m_gridfunc module (a dependency) has netCDF support, which is
# currently not used by plstm.  Compiling with netCDF in some systems
# forces loading the parallel libraries, which sometimes prevents
# interactive use of the program.  For now we deactivate netCDF in
# m_gridfunc in this Makefile.  When netCDF support is added to plstm,
# the following statement can be removed
#
override WITH_NETCDF=
override WITH_MPI=
override WITH_FLOOK=

ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
include $(MAIN_OBJDIR)/check_for_build_mk.mk


FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

plstm.o: m_gridfunc.o m_getopts.o
plsts.o: m_gridfunc.o
plstm: m_getopts.o plstm.o m_gridfunc.o
	$(FC) $(LDFLAGS) -o $@  plstm.o m_gridfunc.o m_getopts.o $(LIBS)

plsts: plsts.o m_gridfunc.o 
	$(FC) $(LDFLAGS) -o $@  plsts.o m_gridfunc.o $(LIBS)


clean: 
	rm -f *.o *.mod plstm plsts

install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin
