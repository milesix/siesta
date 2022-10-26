# ---
# Copyright (C) 1996-2021       The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
# Makefile for macroave
#
.SUFFIXES: .f .F .o .a  .f90 .F90
#
default: macroave
#

TOPDIR=.
MAIN_OBJDIR=.

VPATH=$(TOPDIR)/Util/Macroave/Src:$(TOPDIR)/Src

include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
TOP_SRCS:= interpolation.f90 m_fft_gpfa.F
LOCAL_SRCS:= defs_basis.f90 defs_common.f90  io.f \
  iorho.f  macroave.f thetaft.f surpla.f \
  volcel.f hdr_io.f90 local_die.f dpnint1.f90
#
ALL_SRCS:=$(TOP_SRCS) $(LOCAL_SRCS)
# Generate local .o file list
ALL_OBJS:=$(call change_extensions, $(ALL_SRCS))
#
#
macroave: $(ALL_OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o macroave $(ALL_OBJS) 
#
permute:   permute.o 
	$(FC) $(FFLAGS) $(LDFLAGS) -o permute permute.o
#
# Dependencies
# Generate sources list with proper placement for sfmakedepend
TOP_SRCS:= $(addprefix $(TOPDIR)/Src/,$(TOP_SRCS))
LOCAL_SRCS:= $(addprefix $(TOPDIR)/Util/Macroave/Src/,$(LOCAL_SRCS))
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o $(TOP_SRCS) $(LOCAL_SRCS) || true

clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f macroave permute *.o  *.a *.pcl *.pc *.mod
	rm -f *.mod
#
PROGS:= macroave permute
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin


# DO NOT DELETE THIS LINE - used by make depend
defs_common.o: defs_basis.o
hdr_io.o: defs_basis.o defs_common.o
macroave.o: defs_basis.o defs_common.o interpolation.o m_fft_gpfa.o
thetaft.o: defs_basis.o
gpfa_core_dp.o: m_fft_gpfa.o
gpfa_core_sp.o: m_fft_gpfa.o
gpfa_fft.o: m_fft_gpfa.o
