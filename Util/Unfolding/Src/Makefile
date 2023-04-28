# ---
# Copyright (C) 1996-2021	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---

# Makefile for unfold

.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

TOPDIR=.
MAIN_OBJDIR=.

.PHONY: default what dep version  clean

override WITH_NCDF=
override WITH_FLOOK=

default: what unfold 

VPATH=$(TOPDIR)/Util/Unfolding/Src:$(TOPDIR)/Src:$(TOPDIR)/Util/HSX

ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
include $(MAIN_OBJDIR)/check_for_build_mk.mk

FC_DEFAULT := $(FC)
FC_SERIAL  ?= $(FC_DEFAULT)

# Filter out non-needed libraries
COMP_LIBS := $(filter-out libncdf.a, $(COMP_LIBS))
COMP_LIBS := $(filter-out libfdict.a, $(COMP_LIBS))

# For removing any -ipo compilation in sub directories...
# Later compilers are also having this enabled, and it lets you
# specify what to not be used in libraries.
# This is so because inter-procedural compilation should NOT be performed in
# libraries (it removes names when a routine can be fully moved into others)
IPO_FLAG?= -ipo

SIESTA_SYS=$(SIESTA_ARCH)

what:
	@echo
	@echo "Compilation architecture to be used: $(SIESTA_SYS)"
	@echo "If that is not what you want, give the correct"
	@echo "value to the variable SIESTA_SYS in your shell environment."
	@echo $(COMMENTS)
	@echo

INCFLAGS += $(NETCDF_INCFLAGS) $(FDF_INCFLAGS) \
            $(NCPS_INCFLAGS) $(PSML_INCFLAGS) $(PSOP_INCFLAGS)
INCFLAGS += $(MPI_WRAPPERS_INCFLAGS)
INCFLAGS += $(UNITS_INCFLAGS)

XC = $(GRIDXC_LIBS)
INCFLAGS += $(GRIDXC_INCFLAGS)

SIESTA_SRCS=   precision.F  broadcast_fdf_struct.F90  sys.F  m_io.f  alloc.F90 \
  parallel.F  memory_log.F90  memory.F  pxf.F90  moreParallelSubs.F90  m_mpi_utils.F \
  interpolation.f90  xml.f  radial.f  atm_types.f  spher_harm.f  atmfuncs.f \
  reclat.f  idiag.f  cellsubs.f  sorting.f  minvec.f  io.f  files.f  timestamp.f90 \
  m_walltime.f90  m_wallclock.f90  diag_option.F90  diag.F90  bloch_unfold.F90  debugmpi.F \
  find_kgrid.F  posix_calls.f90  atmparams.f  basis_types.f  atom_options.F90  chemical.f \
  basis_io.F  siesta_geom.F90  get_kpoints_scale.f90  bessph.f  m_fft_gpfa.F \
  radfft.f  periodic_table.f  basis_specs.f  hamann.f90  m_filter.f90  atom.F \
  m_cite.F90  dftu_specs.f  broadcast_basis.F  m_cell.f  zmatrix.F  coor.F \
  chkdim.f  volcel.f  digcel.f  arw.f  redcel.F  printmatrix.F  schecomm.F \
  class_Geometry.F90  class_Data2D.F90  class_OrbitalDistribution.F90  class_Sparsity.F90 \
  class_SpData2D.F90  class_Pair_Geometry_SpData2D.F90  class_Fstack_Pair_Geometry_SpData2D.F90 \
  class_Data1D.F90  class_SpData1D.F90  sparse_matrices.F90  domain_decom.F  spatial.F \
  parallelsubs.F  pspltm1.F  mmio.F  qsort.F  timer_tree.f90  m_timer.F90  timer.F90 \
  m_uuid.f90  object_debug.F90  nag.f  m_spin.F90  m_vee_integrals.F90 

SIESTA_OBJS= $(call change_extensions, $(SIESTA_SRCS))

DEP_OBJS = $(SIESTA_OBJS) $(TS_OBJS) $(SYSOBJ) $(EXTRA)

UNFOLD_SRCS = array.F90 unfold.F90 handlers.f
HSX_SRCS= $(TOPDIR)/Util/HSX/hsx_m.f90

UNFOLD_OBJS  = array.o hsx_m.o
UNFOLD_OBJS += unfold.o 

# Handlers and resolvers of last resort
UNFOLD_OBJS += handlers.o
#

#######################

$(DEP_OBJS): $(COMP_LIBS) 

unfold: $(MPI_INTERFACE) $(COMP_LIBS) $(NCPS) $(PSOP) $(UNITS) $(DEP_OBJS) $(UNFOLD_OBJS) 
	$(FC) $(FFLAGS) $(LDFLAGS) -o unfold \
              $(DEP_OBJS) $(UNFOLD_OBJS) $(MPI_INTERFACE) \
              $(XC) $(NCPS) $(PSOP) $(UNITS) $(PSML_LIBS) $(XMLF90_LIBS) $(FDF_LIBS) $(COMP_LIBS) $(LIBS)

.PHONY: clean
clean: 
	@echo "==> Cleaning object, library, and executable files"
	-rm -f unfold *.o *.a *.mod

PROGS:= unfold
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin

# Dependencies
ABS_TOP_SRCS:= $(addprefix $(TOPDIR)/Src/,$(SIESTA_SRCS))
ABS_LOCAL_SRCS:= $(addprefix $(TOPDIR)/Util/Unfolding/Src/,$(UNFOLD_SRCS))
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o $(ABS_TOP_SRCS) $(ABS_LOCAL_SRCS) $(HSX_SRCS) || true

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
broadcast_basis.o: atm_types.o chemical.o dftu_specs.o m_mpi_utils.o m_spin.o
broadcast_basis.o: parallel.o radial.o
cellsubs.o: precision.o
chemical.o: parallel.o precision.o sys.o
chkdim.o: sys.o
class_Fstack_Pair_Geometry_SpData2D.o: class_Pair_Geometry_SpData2D.o
class_Geometry.o: alloc.o
class_Pair_Geometry_SpData2D.o: class_Geometry.o class_SpData2D.o
class_SpData1D.o: class_Data1D.o class_Data1D.o class_Data1D.o class_Data1D.o
class_SpData1D.o: class_Data1D.o
class_SpData2D.o: class_Data2D.o class_Data2D.o class_Data2D.o
class_Sparsity.o: alloc.o
coor.o: alloc.o atmfuncs.o parallel.o precision.o siesta_geom.o sys.o zmatrix.o
debugmpi.o: parallel.o
dftu_specs.o: alloc.o atm_types.o atmparams.o atom.o atom_options.o
dftu_specs.o: basis_specs.o basis_types.o interpolation.o m_cite.o m_spin.o
dftu_specs.o: m_vee_integrals.o parallel.o precision.o radial.o sys.o
diag.o: alloc.o diag_option.o parallel.o precision.o sys.o
diag_option.o: parallel.o precision.o
domain_decom.o: alloc.o parallel.o precision.o printmatrix.o schecomm.o
domain_decom.o: sparse_matrices.o sys.o
find_kgrid.o: alloc.o minvec.o parallel.o precision.o
get_kpoints_scale.o: parallel.o siesta_geom.o
idiag.o: parallel.o sys.o
io.o: m_io.o
m_cell.o: precision.o siesta_geom.o
m_filter.o: bessph.o precision.o radfft.o sys.o
m_io.o: sys.o
m_mpi_utils.o: precision.o sys.o
m_spin.o: precision.o
m_timer.o: m_io.o m_walltime.o moreParallelSubs.o parallel.o precision.o sys.o
m_vee_integrals.o: precision.o sys.o
m_wallclock.o: m_walltime.o
memory.o: memory_log.o
memory_log.o: m_io.o parallel.o precision.o
minvec.o: cellsubs.o precision.o sorting.o sys.o
moreParallelSubs.o: alloc.o m_io.o parallel.o precision.o sys.o
parallelsubs.o: domain_decom.o parallel.o precision.o spatial.o sys.o
printmatrix.o: alloc.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
radial.o: alloc.o interpolation.o precision.o xml.o
redcel.o: parallel.o precision.o sys.o
schecomm.o: alloc.o
siesta_geom.o: precision.o
sparse_matrices.o: alloc.o class_Fstack_Pair_Geometry_SpData2D.o
sparse_matrices.o: class_OrbitalDistribution.o class_SpData1D.o
sparse_matrices.o: class_SpData2D.o class_SpData2D.o class_Sparsity.o
sparse_matrices.o: precision.o
spatial.o: precision.o
spher_harm.o: alloc.o precision.o sys.o
timer.o: m_timer.o parallel.o sys.o timer_tree.o
timer_tree.o: m_walltime.o
xml.o: precision.o
zmatrix.o: alloc.o m_cell.o parallel.o precision.o siesta_geom.o sys.o
array.o: precision.o sys.o
unfold.o: alloc.o array.o atmfuncs.o basis_io.o basis_types.o cellsubs.o
unfold.o: diag_option.o get_kpoints_scale.o hsx_m.o m_io.o m_mpi_utils.o
unfold.o: m_timer.o memory_log.o parallel.o parallelsubs.o precision.o radfft.o
unfold.o: siesta_geom.o spher_harm.o sys.o
aux_proj.o: atom.o
bloch_unfold_m.o: bloch_unfold.o
class_ddata1d.o: class_Data1D.o
class_ddata2d.o: class_Data2D.o
class_dspdata1d.o: class_SpData1D.o
class_dspdata2d.o: class_SpData2D.o
class_fstack_pair_geometry_dspdata2d.o: class_Fstack_Pair_Geometry_SpData2D.o
class_gdata1d.o: class_Data1D.o
class_gspdata1d.o: class_SpData1D.o
class_idata1d.o: class_Data1D.o
class_idata2d.o: class_Data2D.o
class_ispdata1d.o: class_SpData1D.o
class_ispdata2d.o: class_SpData2D.o
class_ldata1d.o: class_Data1D.o
class_lspdata1d.o: class_SpData1D.o
class_pair_geometry_dspdata2d.o: class_Pair_Geometry_SpData2D.o
class_sdata1d.o: class_Data1D.o
class_sspdata1d.o: class_SpData1D.o
class_zdata1d.o: class_Data1D.o
class_zdata2d.o: class_Data2D.o
class_zspdata1d.o: class_SpData1D.o
class_zspdata2d.o: class_SpData2D.o
gpfa_core_dp.o: m_fft_gpfa.o
gpfa_core_sp.o: m_fft_gpfa.o
gpfa_fft.o: m_fft_gpfa.o
m_array.o: array.o
m_bessph.o: bessph.o
m_diag.o: diag.o
m_diag_option.o: diag_option.o
m_find_kgrid.o: find_kgrid.o
m_get_kpoints_scale.o: get_kpoints_scale.o
m_hamann.o: hamann.o
m_minvec.o: minvec.o
m_object_debug.o: object_debug.o
m_radfft.o: radfft.o
m_timer_tree.o: timer_tree.o
m_timestamp.o: timestamp.o
mtprng.o: m_uuid.o
printmat.o: printmatrix.o
t_spin.o: m_spin.o
timer_options.o: timer.o
