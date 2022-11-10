# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
TOPDIR=.
MAIN_OBJDIR=.

VPATH=$(TOPDIR)/Util/Projections

PROGS= orbmol_proj
default: $(PROGS)

include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

FC:=$(FC_SERIAL)         # Make it non-recursive
#------------------
orbmol_proj: orbmol_proj.o
	$(FC) $(LDFLAGS) -o $@  orbmol_proj.o
#------------------
clean: 
	rm -f *.o orbmol_proj *.o *.*d
#
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin











