# 
# Copyright (C) 1996-2021       The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
# Makefile for PDOS XML processing
#

.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
TOPDIR=.
MAIN_OBJDIR=.

VPATH=$(TOPDIR)/Util/pdosxml

PROGS= pdosxml
default: $(PROGS)

override WITH_MPI=

include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

FC:=$(FC_SERIAL)         # Make it non-recursive

INCFLAGS+= $(XMLF90_INCFLAGS)
#
OBJS=  m_orbital_chooser.o m_pdos.o pdosxml.o
#     
pdosxml:$(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(XMLF90_LIBS)
#
clean: 
	rm -f *.o pdosxml *.o *.*d
#
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin










