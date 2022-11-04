# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#
# Makefile for Contour utilities
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a .f90 .F90
#
TOPDIR=.
MAIN_OBJDIR=.

VPATH=$(TOPDIR)/Util/Contour:$(TOPDIR)/Src

PROGS= grid1d grid2d
default: $(PROGS)

include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

FC:=$(FC_SERIAL)         # Make it non-recursive

grid1d: grid1d.o
	$(FC) $(LDFLAGS) -o $@  grid1d.o
grid2d: grid2d.o
	$(FC) $(LDFLAGS) -o $@  grid2d.o

clean: 
	rm -f *.o $(PROGS)
#
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin
