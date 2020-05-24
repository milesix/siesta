# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for matel-test program and utils
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

OBJDIR=Obj

.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90

VPATH:=$(shell pwd)/../../Src

default: matel-test

override WITH_MPI=

MAIN_OBJDIR:=$(shell cd -P -- ../../\$(OBJDIR) && pwd -P)
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive


DEFS:= $(DEFS) $(DEFS_PREFIX) -UMPI
FPPFLAGS:=$(FPPFLAGS) $(DEFS_PREFIX) -UMPI
INCFLAGS:=$(NETCDF_INCFLAGS) $(INCFLAGS)

SYSOBJ=$(SYS).o

#
# Note that machine-specific files are now in top Src directory.
#
DEPS_MATEL-TEST=io.o m_io.o precision.o  parallel.o \
          atm_types.o chkdim.o propor.o \
          matel_mod.o matel_params.o matel_table.o matel_ylm.o \
          radial.o interpolation.o xml.o alloc.o \
          radfft.o  bessph.o m_fft_gpfa.o \
          pxf.o dot.o register_rfs.o matel_registry.o  \
          m_trialorbitalclass.o units.o $(SYSOBJ) spher_harm.o

#
MATEL-TEST_LOCAL_OBJS=local_sys.o handlers.o m_ion_io.o matel-test.o 
OBJS_MATEL-TEST=$(MATEL-TEST_LOCAL_OBJS) $(DEPS_MATEL-TEST) 

#
matel-test: $(OBJS_MATEL-TEST) 
	$(FC) -o matel-test \
	       $(LDFLAGS) $(OBJS_MATEL-TEST) $(NETCDF_LIBS) $(COMP_LIBS) $(LIBS)

m_ion_io.o: atm_types.o
matel-test.o: matel_mod.o m_ion_io.o

clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f matel-test *.o  *.a
	rm -f *.mod
	rm -f _tmp_deps deps.list  protomake*

ALL_OBJS=$(OBJS_MATEL-TEST)

DEP_OBJS= $(DEPS_MATEL-TEST) 
LOCAL_OBJS=$(MATEL_TEST_LOCAL_OBJS)

dep:
	sfmakedepend --depend=obj --modext=o \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.f) $(DEP_OBJS:.o=.f90)) \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.F) $(DEP_OBJS:.o=.F90)) \
		$(LOCAL_OBJS:.o=.f90) $(LOCAL_OBJS:.o=.F90) \
		$(LOCAL_OBJS:.o=.f) $(LOCAL_OBJS:.o=.F) \
                || true

# DO NOT DELETE THIS LINE - used by make depend
atm_types.o: precision.o radial.o
atm_types.o: precision.o radial.o
bessph.o: precision.o
bessph.o: precision.o
io.o: m_io.o
io.o: m_io.o
m_trialorbitalclass.o: precision.o units.o
m_trialorbitalclass.o: precision.o units.o
matel_mod.o: alloc.o atm_types.o matel_registry.o matel_table.o matel_ylm.o
matel_mod.o: precision.o radfft.o
matel_mod.o: alloc.o atm_types.o matel_registry.o matel_table.o matel_ylm.o
matel_mod.o: precision.o radfft.o
matel_registry.o: m_trialorbitalclass.o precision.o radial.o spher_harm.o
matel_registry.o: m_trialorbitalclass.o precision.o radial.o spher_harm.o
matel_table.o: alloc.o interpolation.o matel_params.o matel_registry.o
matel_table.o: matel_ylm.o parallel.o precision.o radfft.o spher_harm.o
matel_table.o: alloc.o interpolation.o matel_params.o matel_registry.o
matel_table.o: matel_ylm.o parallel.o precision.o radfft.o spher_harm.o
matel_ylm.o: alloc.o matel_params.o matel_registry.o parallel.o precision.o
matel_ylm.o: radfft.o spher_harm.o
matel_ylm.o: alloc.o matel_params.o matel_registry.o parallel.o precision.o
matel_ylm.o: radfft.o spher_harm.o
propor.o: precision.o
propor.o: precision.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
radial.o: alloc.o interpolation.o precision.o xml.o
radial.o: alloc.o interpolation.o precision.o xml.o
register_rfs.o: atm_types.o matel_registry.o parallel.o precision.o radial.o
register_rfs.o: atm_types.o matel_registry.o parallel.o precision.o radial.o
spher_harm.o: precision.o
spher_harm.o: precision.o
units.o: precision.o
units.o: precision.o
xml.o: precision.o
xml.o: precision.o
gpfa_core_dp.o: m_fft_gpfa.o
gpfa_core_sp.o: m_fft_gpfa.o
gpfa_fft.o: m_fft_gpfa.o
m_bessph.o: bessph.o
m_matel_registry.o: matel_registry.o
m_radfft.o: radfft.o
trialorbitalclass.o: m_trialorbitalclass.o
