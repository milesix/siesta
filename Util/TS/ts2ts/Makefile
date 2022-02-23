# ---
# Copyright (C) 1996-2021	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for stand-alone get_chem_labels

OBJDIR=Obj

.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

EXE = ts2ts
default: $(EXE)

VPATH:=$(shell pwd)/../../../Src
MAIN_OBJDIR:=$(shell cd -P -- ../../../\$(OBJDIR) && pwd -P)
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)

dep:
	sfmakedepend --depend=obj  --modext=o \
          $(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 $(VPATH)/*.T90
	@sed -i -e 's/\.T90\.o:/.o:/g' Makefile

override WITH_MPI=

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

DEFS:=$(DEFS_PREFIX) -UMPI $(DEFS) $(DEFS_PREFIX) -UMPI
FPPFLAGS:=$(DEFS_PREFIX) -UMPI $(FPPFLAGS) $(DEFS_PREFIX) -UMPI

INCFLAGS:=$(NETCDF_INCFLAGS) $(FDF_INCFLAGS) $(INCFLAGS)

# Uncomment the following line for debugging support
#FFLAGS=$(FFLAGS_DEBUG)

SYSOBJ=$(SYS).o

# Note that machine-specific files are now in top Src directory.
DEP_OBJS = precision.o parallel.o \
          alloc.o m_io.o pxf.o units.o sys.o $(SYSOBJ)

LOCAL_OBJS += ts2ts.o handlers.o

OBJS = $(DEP_OBJS) $(LOCAL_OBJS)

$(EXE): $(FDF_LIBS) $(OBJS)
	$(FC) -o $(EXE) \
	       $(LDFLAGS) $(OBJS) $(FDF_LIBS)

clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f $(EXE) *.o  *.a *.mod

# Dependencies
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.f) $(DEP_OBJS:.o=.f90)) \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.F) $(DEP_OBJS:.o=.F90)) \
		$(LOCAL_OBJS:.o=.f90) $(LOCAL_OBJS:.o=.F90) \
		$(LOCAL_OBJS:.o=.f) $(LOCAL_OBJS:.o=.F) || true

# DO NOT DELETE THIS LINE - used by make depend
m_io.o: sys.o
units.o: precision.o
ts2ts.o: precision.o units.o
