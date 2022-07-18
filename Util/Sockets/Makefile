# ---
# Copyright (C) 1996-2022      The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for socket communication example
#

.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90 .c
#
TOPDIR=.
MAIN_OBJDIR=.

PROGS:= f2fmaster f2fslave
default: $(PROGS)

override WITH_MPI=

VPATH=$(TOPDIR)/Util/Sockets:$(TOPDIR)/Src

ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
include $(MAIN_OBJDIR)/check_for_build_mk.mk
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
SYSOBJ=$(SYS).o
#
#------------------------------------------------------------------------
#
OTHER_OBJS=fsockets.o sockets.o

f2fmaster: f2fmaster.o $(OTHER_OBJS)
	$(FC) -o $@ \
	       $(LDFLAGS) $(OTHER_OBJS) f2fmaster.o
#
f2fslave: f2fslave.o $(OTHER_OBJS)
	$(FC) -o $@ \
	       $(LDFLAGS) $(OTHER_OBJS) f2fslave.o
#
f2fmaster.o f2fslave.o: fsockets.o sockets.o
fsockets.o: sockets.o
#
clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f $(PROGS)
	rm -f *.mod *.o

install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin

