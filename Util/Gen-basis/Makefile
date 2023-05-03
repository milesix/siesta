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
INCFLAGS += $(UNITS_INCFLAGS)

#
# Note that machine-specific files are now in top Src directory.
#
GEN_BASIS_TOP_SRCS=    alloc.F90  arw.f  atm_types.f  atmparams.f  atom.F \
    atom_options.F90    basis_io.F    basis_specs.f    basis_types.f    bessph.f \
    cellsubs.f    chemical.f    dftu_specs.f    dot.f    files.f    get_kpoints_scale.f90 \
    hamann.f90    interpolation.f90    io.f    minvec.f    m_cite.F90    m_fft_gpfa.F \
    m_filter.f90    m_io.f   m_walltime.f90    memory.F    memory_log.F90 \
    nag.f    parallel.F    parsing.f    periodic_table.f    precision.F    pxf.F90 \
    radfft.f    radial.f    reclat.f    read_xc_info.F    siesta_geom.F90    sorting.f \
    sys.F    xml.f    m_spin.F90    m_vee_integrals.F90

GEN_BASIS_LOCAL_SRCS=gen-basis.F handlers.f

DEPS_GEN-BASIS= $(call change_extensions, $(GEN_BASIS_TOP_SRCS))
GEN-BASIS_LOCAL_OBJS = $(call change_extensions, $(GEN_BASIS_LOCAL_SRCS))

OBJS_GEN-BASIS=$(DEPS_GEN-BASIS) $(GEN-BASIS_LOCAL_OBJS)

IONCAT_TOP_SRCS=  m_getopts.f90  alloc.F90  arw.f  atm_types.f  atmparams.f  atom.F atmfuncs.f \
    atom_options.F90    basis_io.F    basis_specs.f    basis_types.f    bessph.f \
    cellsubs.f   chkdim.f chemical.f    dftu_specs.f    dot.f    files.f    get_kpoints_scale.f90 \
    hamann.f90    interpolation.f90    io.f    minvec.f    m_cite.F90    m_fft_gpfa.F \
    m_filter.f90    m_io.f  m_walltime.f90    memory.F    memory_log.F90 \
    nag.f    parallel.F    parsing.f    periodic_table.f    precision.F    pxf.F90 \
    radfft.f    radial.f    reclat.f    read_xc_info.F    siesta_geom.F90    sorting.f \
    sys.F    xml.f    m_spin.F90    posix_calls.f90 m_vee_integrals.F90 version.F90 spher_harm.F

IONCAT_LOCAL_SRCS=ioncat.f handlers.f

DEPS_IONCAT= $(call change_extensions, $(IONCAT_TOP_SRCS))
IONCAT_LOCAL_OBJS= $(call change_extensions, $(IONCAT_LOCAL_SRCS))

OBJS_IONCAT=$(DEPS_IONCAT) $(IONCAT_LOCAL_OBJS)

#
XC = $(GRIDXC_LIBS)
INCFLAGS += $(GRIDXC_INCFLAGS)

gen-basis: $(NCPS) $(PSOP) $(UNITS) $(OBJS_GEN-BASIS) $(COMP_LIBS)
	$(FC) -o gen-basis \
	       $(LDFLAGS) $(OBJS_GEN-BASIS) $(FDF_LIBS) $(XC) \
              $(NCPS) $(PSOP) $(UNITS) $(PSML_LIBS) $(XMLF90_LIBS) $(NETCDF_LIBS) $(COMP_LIBS) $(LIBS)
ioncat:  $(NCPS) $(PSOP) $(UNITS) $(OBJS_IONCAT) $(COMP_LIBS)
	$(FC) -o ioncat \
	       $(LDFLAGS) $(OBJS_IONCAT) $(FDF_LIBS) \
              $(NCPS) $(PSOP) $(UNITS) $(XC) $(PSML_LIBS) $(XMLF90_LIBS) $(NETCDF_LIBS) $(COMP_LIBS) $(LIBS)
clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f gen-basis ioncat *.o  *.a
	rm -f *.mod

PROGS:= gen-basis ioncat
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin
#
# Dependencies
GEN_BASIS_ABS_TOP_SRCS:= $(addprefix $(TOPDIR)/Src/,$(GEN_BASIS_TOP_SRCS))
GEN_BASIS_ABS_LOCAL_SRCS:= $(addprefix $(TOPDIR)/Util/Genbasis/,$(GEN_BASIS_LOCAL_SRCS))
IONCAT_ABS_TOP_SRCS:= $(addprefix $(TOPDIR)/Src/,$(IONCAT_TOP_SRCS))
IONCAT_ABS_LOCAL_SRCS:= $(addprefix $(TOPDIR)/Util/Genbasis/,$(IONCAT_LOCAL_SRCS))
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o $(GEN_BASIS_ABS_TOP_SRCS) $(GEN_BASIS_ABS_LOCAL_SRCS) \
	    $(IONCAT_ABS_TOP_SRCS) $(IONCAT_ABS_LOCAL_SRCS) || true

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
cellsubs.o: precision.o
chemical.o: parallel.o precision.o sys.o
chkdim.o: sys.o
dftu_specs.o: alloc.o atm_types.o atmparams.o atom.o atom_options.o
dftu_specs.o: basis_specs.o basis_types.o interpolation.o m_cite.o m_spin.o
dftu_specs.o: m_vee_integrals.o parallel.o precision.o radial.o sys.o
get_kpoints_scale.o: parallel.o siesta_geom.o
io.o: m_io.o
m_filter.o: bessph.o precision.o radfft.o sys.o
m_io.o: sys.o
m_spin.o: precision.o
m_vee_integrals.o: precision.o sys.o
memory.o: memory_log.o
memory_log.o: m_io.o parallel.o precision.o
minvec.o: cellsubs.o precision.o sorting.o sys.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
radial.o: alloc.o interpolation.o precision.o xml.o
read_xc_info.o: parallel.o precision.o sys.o
siesta_geom.o: precision.o
spher_harm.o: alloc.o precision.o sys.o
xml.o: precision.o
aux_proj.o: atom.o
gpfa_core_dp.o: m_fft_gpfa.o
gpfa_core_sp.o: m_fft_gpfa.o
gpfa_fft.o: m_fft_gpfa.o
m_bessph.o: bessph.o
m_get_kpoints_scale.o: get_kpoints_scale.o
m_hamann.o: hamann.o
m_minvec.o: minvec.o
m_radfft.o: radfft.o
t_spin.o: m_spin.o
version_info.o: version.o
