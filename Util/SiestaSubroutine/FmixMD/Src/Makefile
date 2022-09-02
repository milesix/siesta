# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
# Makefile for Client-Server example
#
# JM Soler, A Garcia
#
TOPDIR=.
MAIN_OBJDIR=.
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90
#
dummy: default
#
override WITH_NETCDF=
override WITH_MPI=
#
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
include $(MAIN_OBJDIR)/check_for_build_mk.mk
#
VPATH:=$(TOPDIR)/Util/SiestaSubroutine/FmixMD/Src:$(TOPDIR)/Src
#
#       driver: Mixed molecular dynamics
#       simple: Simple test
#       para  : Test driver for parallel calculation
#
PROGS:= driver simple para   # Consider keeping only the first one
default: $(PROGS)
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
driver.o: fmixmd.o sample.o
fmixmd.o: sample.o
fsiesta_pipes.o: posix_calls.o
#
OBJS= fsiesta_pipes.o sample.o fmixmd.o driver.o pxf.o posix_calls.o
driver : $(OBJS)
	$(FC) -o driver  $(LDFLAGS) $(OBJS) 
#
simple: fsiesta_pipes.o simple.o  pxf.o posix_calls.o
	$(FC) $(LDFLAGS) -o simple simple.o fsiesta_pipes.o pxf.o posix_calls.o
#
para: fsiesta_pipes.o para.o pxf.o posix_calls.o
	$(FC) $(LDFLAGS) -o para para.o fsiesta_pipes.o pxf.o posix_calls.o
#
clean: 
	@echo "==> Cleaning object and executable files"
	rm -f $(PROGS) *.o *.mod

install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin



