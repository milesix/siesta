# 
# Copyright (C) 1996-2021       The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#
# Makefile for DM to DM-netCDF converter
# and other utilities
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
#
dummy: default
#
TOPDIR=.
MAIN_OBJDIR=.

VPATH=$(TOPDIR)/Util/DensityMatrix:$(TOPDIR)/Src

include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

PROGS:=dmfilter dmbs2dm dmUnblock dm_noncol_sign_flip4
ifeq ($(WITH_NETCDF),1)
 PROGS+= dm2cdf cdf2dm 
endif
default: $(PROGS)

OPT_OBJS=m_getopts.o
#
dmbs2dm.o: $(OPT_OBJS)
dmbs2dm: dmbs2dm.o $(OPT_OBJS)
	$(FC) $(LDFLAGS) -o $@ dmbs2dm.o $(OPT_OBJS)
#
dmUnblock.o: $(OPT_OBJS)
dmUnblock: dmUnblock.o $(OPT_OBJS)
	$(FC) $(LDFLAGS) -o $@ dmUnblock.o $(OPT_OBJS)
#
dm2cdf: dm2cdf.o
	$(FC) $(LDFLAGS) -o $@ dm2cdf.o  $(NETCDF_LIBS)
#
cdf2dm: cdf2dm.o
	$(FC) $(LDFLAGS) -o $@ cdf2dm.o  $(NETCDF_LIBS)
#
dmfilter: dmfilter.o
	$(FC) $(LDFLAGS) -o $@ dmfilter.o
#
dm_noncol_sign_flip4: dm_noncol_sign_flip4.o
	$(FC) $(LDFLAGS) -o $@ dm_noncol_sign_flip4.o
#
clean: 
	rm -f *.o *.*d
	rm -f $(PROGS)
#
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin
#





