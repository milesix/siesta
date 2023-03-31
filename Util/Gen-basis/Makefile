# ---
# Copyright (C) 1996-2021	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for stand-alone Gen-basis and ioncat
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90

TOPDIR=.
MAIN_OBJDIR=.

default: gen-basis ioncat

override WITH_MPI=
override WITH_NCDF=
override WITH_FLOOK=

VPATH=$(TOPDIR)/Util/Gen-basis:$(TOPDIR)/Src

ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
include $(MAIN_OBJDIR)/check_for_build_mk.mk

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

INCFLAGS:=$(NETCDF_INCFLAGS) $(FDF_INCFLAGS) \
          $(NCPS_INCFLAGS) $(PSOP_INCFLAGS) $(PSML_INCFLAGS) \
          $(INCFLAGS)
SYSOBJ=$(SYS).o

#
# Note that machine-specific files are now in top Src directory.
#
DEPS_GEN-BASIS=m_io.o io.o alloc.o basis_types.o precision.o  parallel.o \
          parsing.o basis_io.o chemical.o atm_types.o \
          atmparams.o radial.o interpolation.o \
	  bessph.o radfft.o sorting.o m_filter.o \
          basis_specs.o dftu_specs.o atom.o periodic_table.o\
          pxf.o dot.o  atom_options.o arw.o \
          xml.o m_walltime.o read_xc_info.o  $(SYSOBJ) \
          m_fft_gpfa.o m_cite.o version.o posix_calls.o reclat.o units.o \
          hamann.o files.o sys.o m_spin.o m_vee_integrals.o psml_die.o

GEN-BASIS_LOCAL_OBJS=gen-basis.o handlers.o 
OBJS_GEN-BASIS=$(DEPS_GEN-BASIS) $(GEN-BASIS_LOCAL_OBJS)

DEPS_IONCAT=m_getopts.o basis_types.o precision.o  parallel.o \
          parsing.o m_io.o io.o alloc.o atom_options.o basis_io.o atm_types.o \
          atmparams.o radial.o  chkdim.o m_filter.o \
	  bessph.o radfft.o m_fft_gpfa.o chemical.o reclat.o \
          basis_specs.o periodic_table.o pxf.o \
          atom.o atmfuncs.o spher_harm.o interpolation.o arw.o \
          xml.o m_walltime.o dftu_specs.o m_spin.o m_vee_integrals.o \
          hamann.o m_cite.o version.o posix_calls.o units.o files.o sys.o \
          psml_die.o $(SYSOBJ)

IONCAT_LOCAL_OBJS=ioncat.o handlers.o 
OBJS_IONCAT=$(DEPS_IONCAT) $(IONCAT_LOCAL_OBJS)

#
XC = $(GRIDXC_LIBS)
INCFLAGS:= $(INCFLAGS) $(GRIDXC_INCFLAGS)

gen-basis: $(NCPS) $(PSOP) $(FDF_LIBS) $(OBJS_GEN-BASIS) $(COMP_LIBS)
	$(FC) -o gen-basis \
	       $(LDFLAGS) $(OBJS_GEN-BASIS) $(FDF_LIBS) $(XC) \
              $(NCPS) $(PSOP) $(PSML_LIBS) $(XMLF90_LIBS) $(NETCDF_LIBS) $(COMP_LIBS) $(LIBS)
ioncat:  $(NCPS) $(PSOP) $(FDF_LIBS) $(OBJS_IONCAT) $(COMP_LIBS)
	$(FC) -o ioncat \
	       $(LDFLAGS) $(OBJS_IONCAT) $(FDF_LIBS) \
              $(NCPS) $(PSOP) $(XC) $(PSML_LIBS) $(XMLF90_LIBS) $(NETCDF_LIBS) $(COMP_LIBS) $(LIBS)
clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f gen-basis ioncat *.o  *.a
	rm -f *.mod
	rm -f _tmp_deps deps.list  protomake*

PROGS:= gen-basis ioncat
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin
#
ALL_OBJS=$(OBJS_IONCAT) $(OBJS_GEN-BASIS)

DEP_OBJS= $(DEPS_GEN-BASIS) $(DEPS_IONCAT)
LOCAL_OBJS=$(IONCAT_LOCAL_OBJS) $(GEN-BASIS_LOCAL_OBJS)

dep:
	sfmakedepend --depend=obj --modext=o \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.f) $(DEP_OBJS:.o=.f90)) \
		$(addprefix $(VPATH)/,$(DEP_OBJS:.o=.F) $(DEP_OBJS:.o=.F90)) \
		$(LOCAL_OBJS:.o=.f90) $(LOCAL_OBJS:.o=.F90) \
		$(LOCAL_OBJS:.o=.f) $(LOCAL_OBJS:.o=.F) \
                || true

# DO NOT DELETE THIS LINE - used by make depend
arw.o: alloc.o parallel.o precision.o sys.o
atm_types.o: precision.o radial.o
atmfuncs.o: atm_types.o precision.o radial.o spher_harm.o sys.o
atom.o: atm_types.o atmparams.o atom_options.o basis_specs.o basis_types.o
atom.o: hamann.o interpolation.o m_filter.o periodic_table.o precision.o
atom.o: radial.o sys.o
atom_options.o: sys.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: precision.o radial.o sys.o xml.o
basis_specs.o: atom_options.o basis_types.o chemical.o periodic_table.o
basis_specs.o: precision.o sys.o
basis_types.o: alloc.o atmparams.o precision.o sys.o
bessph.o: precision.o sys.o
chemical.o: parallel.o precision.o sys.o
chkdim.o: sys.o
dftu_specs.o: alloc.o atm_types.o atmparams.o atom.o atom_options.o
dftu_specs.o: basis_specs.o basis_types.o interpolation.o m_cite.o m_spin.o
dftu_specs.o: m_vee_integrals.o parallel.o precision.o radial.o sys.o units.o
io.o: m_io.o
m_filter.o: bessph.o precision.o radfft.o sys.o
m_io.o: sys.o
m_spin.o: precision.o
m_vee_integrals.o: precision.o sys.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
radial.o: alloc.o interpolation.o precision.o xml.o
read_xc_info.o: parallel.o precision.o sys.o
spher_harm.o: alloc.o precision.o sys.o
units.o: precision.o
xml.o: precision.o
gen-basis.o: atm_types.o atom.o atom_options.o basis_io.o basis_specs.o
gen-basis.o: basis_types.o chemical.o dftu_specs.o sys.o
ioncat.o: atm_types.o atmfuncs.o basis_io.o m_getopts.o precision.o
aux_proj.o: atom.o
gpfa_core_dp.o: m_fft_gpfa.o
gpfa_core_sp.o: m_fft_gpfa.o
gpfa_fft.o: m_fft_gpfa.o
m_bessph.o: bessph.o
m_hamann.o: hamann.o
m_radfft.o: radfft.o
t_spin.o: m_spin.o
version_info.o: version.o
