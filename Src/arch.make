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
FPPFLAGS= -fbounds-check -g -fcheck=all -DFC_HAVE_FLUSH -DFC_HAVE_ABORT -fbacktrace


AR = ar
RANLIB = ranlib

SYS = nag

SP_KIND = 4
DP_KIND = 8
KINDS = $(SP_KIND) $(DP_KIND)

LDFLAGS =

COMP_LIBS = libsiestaLAPACK.a libsiestaBLAS.a

#FPPFLAGS = $(DEFS_PREFIX)-DFC_HAVE_ABORT

LIBS = $(COMP_LIBS)

###ADD LIBXC STUFF
#libxc_root=points to your libxc installation path
LIBXC_ROOT=/home/likewise-open/ICN/sillera/Desktop/LIBXC/LinresLibxc/libxc-2.2.3/INSTAL
LIBXC_INCFLAGS= -I $(LIBXC_ROOT)/include
LIBXC_LIBS=-L$(LIBXC_ROOT)/lib -l xcf90 -l xc


#LIBXC_INCFLAGS= 
USE_LIBXC=1


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

