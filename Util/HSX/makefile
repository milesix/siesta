# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#
# Makefile for HSX utilities
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
TOPDIR=.
MAIN_OBJDIR=.

PROGS:= hsx2hs hs2hsx
all: $(PROGS)

VPATH=$(TOPDIR)/Util/HSX

include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
FFLAGS=$(FFLAGS_DEBUG)
#
#------------------
hsx2hs: hsx_m.o hsx2hs.o
	$(FC) $(LDFLAGS) -o $@ hsx_m.o  hsx2hs.o
#
hs2hsx: hsx_m.o hs2hsx.o
	$(FC) $(LDFLAGS) -o $@ hsx_m.o  hs2hsx.o
#
#------------------
clean: 
	rm -f *.o hsx2hs hs2hsx *.o *.*d
#
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin
	cp -p hsx_m.mod $(SIESTA_INSTALL_DIRECTORY)/include
	cp -p hsx_m.o $(SIESTA_INSTALL_DIRECTORY)/lib




