# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#-------------------------------------------------------------------
# arch.make file for gfortran compiler.
# To use this arch.make file you should rename it to
#   arch.make
# or make a sym-link.
# For an explanation of the flags see DOCUMENTED-TEMPLATE.make

.SUFFIXES:
.SUFFIXES: .f .F .o .c .a .f90 .F90

SIESTA_ARCH = unknown

CC = gcc
FPP = $(FC) -E -P -x c
FC = gfortran
FC_SERIAL = gfortran

FFLAGS = -O2 -fPIC -ftree-vectorize

AR = ar
RANLIB = ranlib

SYS = nag

LDFLAGS =

#
# Required libraries.
# Please edit to match your installation directories
#
PSML_ROOT=$(HOME)/lib/gfortran-5.2.0/libpsml-0.9.4
XMLF90_ROOT=$(HOME)/lib/gfortran-5.2.0/xmlf90-1.5.0
GRIDXC_ROOT=$(HOME)/lib/gfortran-5.2.0/gridxc-0.6.1

#
# These are non-optimized libraries. You should
# use optimized versions for production runs.
#
COMP_LIBS = libsiestaLAPACK.a libsiestaBLAS.a

FPPFLAGS := $(DEFS_PREFIX)-DFC_HAVE_ABORT $(FPPFLAGS)
FPPFLAGS := $(DEFS_PREFIX)-D__NO_PROC_POINTERS__ $(FPPFLAGS)

LIBS = $(COMP_LIBS)

# Dependency rules ---------

FFLAGS_DEBUG = -g -O1   # your appropriate flags here...

# The atom.f code is very vulnerable. Particularly the Intel compiler
# will make an erroneous compilation of atom.f with high optimization
# levels.
atom.o: atom.F
	$(FC) -c $(FFLAGS_DEBUG) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) $< 

.c.o:
	$(CC) -c $(CFLAGS) $(INCFLAGS) $(CPPFLAGS) $< 
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F)  $< 
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_free_F90) $< 
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_fixed_f)  $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_free_f90)  $<
