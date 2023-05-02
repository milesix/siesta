# 
# Copyright (C) 1996-2021       The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
# Makefile for WFS utilities
#
# The VPATH directive below allows to re-use m_getopts.f90 from the top Src.
# Other points to note, until we switch to a better building system:
#
#  The arch.make file is supposed to be in $(OBJDIR). This is normally
#  the top Obj, but if you are using architecture-dependent build directories
#  you might want to change this. (If you do not understand this, you do not
#  need to change anything. Power users can do "make OBJDIR=Whatever".)
#
#  If your main Siesta build used an mpi compiler, you might need to
#  define an FC_SERIAL symbol in your top arch.make, to avoid linking
#  in the mpi libraries even if we explicitly undefine MPI below.
#  
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
dummy: default
#
override WITH_MPI=

TOPDIR=.
MAIN_OBJDIR=.

VPATH=$(TOPDIR)/Util/WFS:$(TOPDIR)/Src

include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

PROGS:= readwf readwfx wfs2wfsx wfsx2wfs
ifeq ($(WITH_NETCDF),1)
 PROGS+=  wfsnc2wfsx
endif
default: $(PROGS)

INCFLAGS += $(UNITS_INCFLAGS)
#
# This is needed on some systems to avoid loading the parallel libraries, which
# sometimes force running on queuing systems
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#

readwf: readwf.o
	$(FC) $(LDFLAGS) -o $@ $<  
readwfx: readwfx.o m_getopts.o
	$(FC) $(LDFLAGS) -o $@ $<  m_getopts.o
wfs2wfsx: wfs2wfsx.o
	$(FC) $(LDFLAGS) -o $@ $<  
wfsx2wfs: wfsx2wfs.o
	$(FC) $(LDFLAGS) -o $@ $<  

wfsnc2wfsx.o: precision.o 
wfsnc2wfsx: wfsnc2wfsx.o precision.o $(UNITS)
	$(FC) $(LDFLAGS) -o $@ $< $(NETCDF_LIBS) $(UNITS)

clean: 
	rm -f *.o  *.*d $(PROGS)
#
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin

readwfx.o: m_getopts.o











