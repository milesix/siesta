# 
# Copyright (C) 1996-2014	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#
# Makefile for Grid utilities
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
VPATH:=$(shell pwd)/../../Src
OBJDIR=Obj
#
all: grid2cdf cdf2xsf cdf2grid grid2val grid2cube
#
include ../../$(OBJDIR)/arch.make
#
# This is needed on some systems to avoid loading the parallel libraries, which
# sometimes force running on queuing systems
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
##FFLAGS=$(FFLAGS_DEBUG)
#
# We assume that the top arch.make defines the relevant NetCDF symbols,
# including the INCFLAGS line. If not, uncomment the following lines, 
# using appropriate information
#
# NETCDF_ROOT=/some/path
# NETCDF_INCFLAGS=-I$(NETCDF_ROOT)/include
# NETCDF_LIBS=-L$(NETCDF_ROOT)/lib -lnetcdf
#
INCFLAGS:= $(NETCDF_INCFLAGS) $(INCFLAGS)
#------------------
grid2cdf: grid2cdf.o
	$(FC) $(LDFLAGS) -o $@  grid2cdf.o  $(NETCDF_LIBS)
#
cdf2grid: cdf2grid.o
	$(FC) $(LDFLAGS) -o $@  cdf2grid.o  $(NETCDF_LIBS)
#
cdf2xsf: cdf2xsf.o
	$(FC) $(LDFLAGS) -o $@  cdf2xsf.o  $(NETCDF_LIBS)
#
cdfget: cdfget.o
	$(FC) $(LDFLAGS) -o $@  cdfget.o  $(NETCDF_LIBS)
cdf_diff: m_grid.o cdf_diff.o 
	$(FC) $(LDFLAGS) -o $@  m_grid.o cdf_diff.o  $(NETCDF_LIBS)
cdf_diff.o: m_grid.o
#------------------
cdf_fft: cdf_fft.o fft3d.o m_fft_gpfa.o local_die.o
	$(FC) $(LDFLAGS) -o $@  cdf_fft.o \
               local_die.o fft3d.o m_fft_gpfa.o $(NETCDF_LIBS)
#
fft3d.o: m_fft_gpfa.o
#------------------
grid2val.o: f2kcli.o
grid2val: grid2val.o
	$(FC) $(LDFLAGS) -o $@ f2kcli.o grid2val.o 
#
#------------------
grid2cube: grid2cube.o
	$(FC) $(LDFLAGS) -o $@ grid2cube.o 
#
#------------------
clean: 
	rm -f *.o grid2cdf grid2val grid2cube cdf2grid cdf2xsf *.o *.*d
#










