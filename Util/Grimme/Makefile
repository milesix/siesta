# ---
# Copyright (C) 1996-2021	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---

.SUFFIXES: .f .F .o .a  .f90 .F90

EXE = fdf2grimme
default: $(EXE)

override WITH_MPI=
override WITH_NETCDF=

TOPDIR=.
MAIN_OBJDIR=.

VPATH=$(TOPDIR)/Util/Grimme:$(TOPDIR)/Src

include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

INCFLAGS+= $(FDF_INCFLAGS) 
INCFLAGS+= $(LIBSYS_INCFLAGS) 

# Note that machine-specific files are now in top Src directory.
DEP_OBJS = precision.o parallel.o \
          chemical.o io.o m_io.o  \
          periodic_table.o \
          atom_options.o units.o

LOCAL_OBJS= grimme.o
OBJS=$(LOCAL_OBJS) $(DEP_OBJS)

$(EXE): $(OBJS) $(LIBSYS)
	$(FC) -o $(EXE) \
	       $(LDFLAGS) $(OBJS) $(FDF_LIBS) $(LIBSYS)

dep:
	sfmakedepend --depend=obj --modext=o \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.f) $(DEP_OBJS:.o=.f90)) \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.F) $(DEP_OBJS:.o=.F90)) \
		$(LOCAL_OBJS:.o=.f90) $(LOCAL_OBJS:.o=.F90) \
		$(LOCAL_OBJS:.o=.f) $(LOCAL_OBJS:.o=.F) \
                || true
clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f $(EXE) *.o  *.a *.mod

PROGS:= fdf2grimme
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin
#

# DO NOT DELETE THIS LINE - used by make depend
chemical.o: parallel.o precision.o
io.o: m_io.o
grimme.o: chemical.o periodic_table.o precision.o units.o
