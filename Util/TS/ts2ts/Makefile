# ---
# Copyright (C) 1996-2021	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for stand-alone get_chem_labels

.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

EXE = ts2ts
default: $(EXE)

override WITH_MPI=

TOPDIR=.
MAIN_OBJDIR=.

VPATH=$(TOPDIR)/Util/TS/ts2ts:$(TOPDIR)/Src
TOPSRC=$(TOPDIR)/Src

ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
include $(MAIN_OBJDIR)/check_for_build_mk.mk

dep:
	sfmakedepend --depend=obj  --modext=o \
          $(TOPSRC)/*.f $(TOPSRC)/*.f90 $(TOPSRC)/*.F $(TOPSRC)/*.F90 $(TOPSRC)/*.T90
	@sed -i -e 's/\.T90\.o:/.o:/g' Makefile


FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

INCFLAGS:=$(NETCDF_INCFLAGS) $(FDF_INCFLAGS) $(INCFLAGS)
INCFLAGS+=$(LIBSYS_INCFLAGS) 

# Uncomment the following line for debugging support
#FFLAGS=$(FFLAGS_DEBUG)

SYSOBJ=$(SYS).o

# Note that machine-specific files are now in top Src directory.
DEP_OBJS = precision.o parallel.o \
          alloc.o m_io.o units.o $(SYSOBJ)

LOCAL_OBJS += ts2ts.o

OBJS = $(DEP_OBJS) $(LOCAL_OBJS)

$(EXE): $(OBJS) $(LIBSYS)
	$(FC) -o $(EXE) \
	       $(LDFLAGS) $(OBJS) $(FDF_LIBS) $(LIBSYS)

clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f $(EXE) *.o  *.a *.mod

PROGS:= ts2ts
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin

# Dependencies
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o \
		$(addprefix $(TOPSRC)/,$(DEP_OBJS:.o=.f) $(DEP_OBJS:.o=.f90)) \
		$(addprefix $(TOPSRC)/,$(DEP_OBJS:.o=.F) $(DEP_OBJS:.o=.F90)) \
		$(LOCAL_OBJS:.o=.f90) $(LOCAL_OBJS:.o=.F90) \
		$(LOCAL_OBJS:.o=.f) $(LOCAL_OBJS:.o=.F) || true

# DO NOT DELETE THIS LINE - used by make depend
units.o: precision.o
ts2ts.o: precision.o units.o
