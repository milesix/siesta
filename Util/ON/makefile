# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#
# Makefile for LWF to netCDF converter
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
TOPDIR=.
MAIN_OBJDIR=.

VPATH=$(TOPDIR)/Util/ON

include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

PROGS=
ifeq ($(WITH_NETCDF),1)
 PROGS+= lwf2cdf
 COPY="cp -p"
else
 COPY=echo
endif
default: $(PROGS)

FC:=$(FC_SERIAL)         # Make it non-recursive

INCFLAGS:= $(NETCDF_INCFLAGS) $(INCFLAGS)
#
lwf2cdf: lwf2cdf.o
	$(FC) $(LDFLAGS) -o $@ lwf2cdf.o  $(NETCDF_LIBS)
#
clean: 
	rm -f *.o *.*d
	rm -f $(PROGS)
#
install: $(PROGS)
	$(COPY) $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin





