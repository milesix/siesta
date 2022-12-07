# ---
# Copyright (C) 1996-2021       The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# The VPATH directive below allows to re-use m_getopts.f90 from the top Src.

.SUFFIXES:
.SUFFIXES: .f .f90 .F .F90 .o

TOPDIR=.
MAIN_OBJDIR=.

PROGS:= eigfat2plot gnubands
all: $(PROGS)

override WITH_MPI=

VPATH=$(TOPDIR)/Util/Bands:$(TOPDIR)/Src

ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
include $(MAIN_OBJDIR)/check_for_build_mk.mk

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

OBJS=m_getopts.o

gnubands: $(OBJS) gnubands.o
	$(FC) -o $@ $(LDFLAGS) $(OBJS) gnubands.o

eigfat2plot: $(OBJS) eigfat2plot.o
	$(FC) -o $@ $(LDFLAGS) $(OBJS) eigfat2plot.o

clean:
	rm -f *.o *.mod eigfat2plot gnubands

install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin

dep:
	-sfmakedepend --depend=obj --modext=o \
		*.f90 *.F90 *.f *.F

# DO NOT DELETE THIS LINE - used by make depend
