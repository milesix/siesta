# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#
# Makefile for Optical utilities
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a .f90 .F90

TOPDIR=.
MAIN_OBJDIR=.

VPATH=$(TOPDIR)/Util/Optical

PROGS= optical_input optical
default: $(PROGS)

include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

FC:=$(FC_SERIAL)         # Make it non-recursive

optical_input: optical_input.o
	$(FC) $(LDFLAGS) -o $@ optical_input.o

optical: optical.o
	$(FC) $(LDFLAGS) -o $@ optical.o

clean: 
	rm -f *.o $(PROGS)
#
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin

