# ---
# Copyright (C) 1996-2021	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for stand-alone get_chem_labels

.SUFFIXES: .f .F .o .a  .f90 .F90

TOPDIR=.
MAIN_OBJDIR=.

default: get_chem_labels
#
override WITH_MPI=
override WITH_NETCDF=
override WITH_FLOOK=

VPATH=$(TOPDIR)/Util/Helpers:$(TOPDIR)/Src

ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
include $(MAIN_OBJDIR)/check_for_build_mk.mk


FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

INCFLAGS+=$(FDF_INCFLAGS)
#
SYSOBJ=$(SYS).o
#
# Note that machine-specific files are now in top Src directory.
#
SIESTA_OBJS=precision.o parallel.o \
          chemical.o io.o m_io.o alloc.o \
          memory_log.o memory.o periodic_table.o\
          pxf.o atom_options.o  \
          $(SYSOBJ)

LOCAL_OBJS=get_chem_labels.o handlers.o local_sys.o
OBJS_GET_CHEM_LABELS=$(SIESTA_OBJS) $(LOCAL_OBJS)
#
#
get_chem_labels: $(FDF_LIBS)  $(OBJS_GET_CHEM_LABELS)
	$(FC) -o get_chem_labels \
	       $(LDFLAGS) $(OBJS_GET_CHEM_LABELS) $(FDF_LIBS) $(LIBS)
#
clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f get_chem_labels  *.o  *.a
	rm -f *.mod .tmp_fdflog

PROGS:= get_chem_labels
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin

# Dependencies
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o \
		$(addprefix $(VPATH)/,$(SIESTA_OBJS:.o=.f) $(SIESTA_OBJS:.o=.f90)) \
		$(addprefix $(VPATH)/,$(SIESTA_OBJS:.o=.F) $(SIESTA_OBJS:.o=.F90)) \
		$(LOCAL_OBJS:.o=.f90) $(LOCAL_OBJS:.o=.F90) \
		$(LOCAL_OBJS:.o=.f) $(LOCAL_OBJS:.o=.F) || true

# DO NOT DELETE THIS LINE - used by make depend
atom_options.o: local_sys.o
chemical.o: local_sys.o parallel.o precision.o
io.o: m_io.o
m_io.o: local_sys.o
memory.o: memory_log.o
memory_log.o: m_io.o parallel.o precision.o
get_chem_labels.o: chemical.o
sys.o: local_sys.o
