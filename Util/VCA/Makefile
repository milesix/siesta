# ---
# Copyright (C) 1996-2021	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for stand-alone mixps and fractional
#
# The idea is to use the code in the top Src directory as much as possible.
# This is achieved by the VPATH directive below.
# Other points to note, until we switch to a better building system:
#
#  The arch.make file is supposed to be in $(OBJDIR). This is normally
#  the top Obj, but if you are using architecture-dependent build directories
#  you might want to change this. (If you do not understand this, you do not
#  need to change anything. Power users can do "make OBJDIR=Whatever".)
#
#  If your main Siesta build used an mpi compiler, you might need to
#  define an FC_SERIAL symbol in your top arch.make, to avoid linking
#  in the mpi libraries even if we explicitly undefine MPI below.
#  
#  The dependency list at the end is overly large, but harmless
#
OBJDIR=Obj
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90

override WITH_MPI=
override WITH_NETCDF=

default: mixps fractional

TOPDIR=.
MAIN_OBJDIR=.

VPATH=$(TOPDIR)/Util/VCA:$(TOPDIR)/Src

include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk


FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive


INCFLAGS:=$(NETCDF_INCFLAGS)\
          $(NCPS_INCFLAGS) \
          $(PSML_INCFLAGS) \
          $(FDF_INCFLAGS) \
          $(INCFLAGS)

SYSOBJ=$(SYS).o

#------------------------------------------------------------------------

MIXPS_DEPS= interpolation.o atom_options.o \
            m_io.o io.o alloc.o periodic_table.o  precision.o
MIXPS_LOCAL_OBJS= mixps.o handlers.o local_sys.o
MIXPS_OBJS= $(MIXPS_DEPS) $(MIXPS_LOCAL_OBJS)
#

mixps:  $(FDF_LIBS) $(NCPS) $(MIXPS_OBJS)
	$(FC) -o mixps \
	       $(LDFLAGS) $(MIXPS_OBJS) $(NCPS) $(FDF_LIBS) $(PSML_LIBS) $(XMLF90_LIBS)

#----------------------------------------------------------------------------
FRACTIONAL_DEPS= m_io.o alloc.o interpolation.o \
                 atom_options.o io.o periodic_table.o precision.o
FRACTIONAL_LOCAL_OBJS= fractional.o handlers.o local_sys.o
FRACTIONAL_OBJS= $(FRACTIONAL_DEPS) $(FRACTIONAL_LOCAL_OBJS)

fractional: $(FDF_LIBS)  $(NCPS) $(FRACTIONAL_OBJS)
	$(FC) -o fractional \
	       $(LDFLAGS) $(FRACTIONAL_OBJS) $(NCPS) $(PSML) $(FDF_LIBS)\
               $(PSML_LIBS) $(XMLF90_LIBS)
clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f mixps fractional *.o  *.a
	rm -f *.mod

PROGS:= mixps fractional
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin

DEP_OBJS= $(MIXPS_DEPS) $(FRACTIONAL_DEPS)
LOCAL_OBJS=$(FRACTIONAL_LOCAL_OBJS) $(MIXPS_LOCAL_OBJS)
dep:
	sfmakedepend --depend=obj --modext=o \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.f) $(DEP_OBJS:.o=.f90)) \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.F) $(DEP_OBJS:.o=.F90)) \
		$(LOCAL_OBJS:.o=.f90) $(LOCAL_OBJS:.o=.F90) \
		$(LOCAL_OBJS:.o=.f) $(LOCAL_OBJS:.o=.F) \
                || true

# DO NOT DELETE THIS LINE - used by make depend
atom_options.o: local_sys.o
io.o: m_io.o
m_io.o: local_sys.o
fractional.o: periodic_table.o precision.o
mixps.o: interpolation.o periodic_table.o precision.o
sys.o: local_sys.o
