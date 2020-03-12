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

VPATH=$(TOPDIR)/Util/Unfolding/Src:$(TOPDIR)/Src
TOPSRC=$(TOPDIR)/Src

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

XC = $(GRIDXC_LIBS)
INCFLAGS += $(GRIDXC_INCFLAGS)

SYSOBJ=$(SYS).o

SIESTA_OBJS= precision.o m_io.o alloc.o memory_log.o memory.o pxf.o \
             moreParallelSubs.o m_mpi_utils.o atmfuncs.o atm_types.o \
             radial.o reclat.o idiag.o minvec.o \
             parallel.o io.o files.o cellsubs.o sorting.o \
	     timestamp.o m_wallclock.o \
             diag.o diag_option.o bloch_unfold.o \
	     debugmpi.o find_kgrid.o units.o posix_calls.o xml.o \
             interpolation.o spher_harm.o basis_types.o atmparams.o basis_io.o \
             atom_options.o chemical.o \
             get_kpoints_scale.o siesta_geom.o \
             radfft.o bessph.o m_fft_gpfa.o \
             broadcast_basis.o dftu_specs.o basis_specs.o \
             atom.o hamann.o m_filter.o \
             m_cite.o \
             coor.o chkdim.o \
             zmatrix.o m_cell.o  periodic_table.o volcel.o digcel.o \
             arw.o redcel.o broadcast_fdf_struct.o \
             sys.o parallelsubs.o spatial.o \
             domain_decom.o printmatrix.o schecomm.o \
             pspltm1.o mmio.o qsort.o \
             sparse_matrices.o

SIESTA_OBJS += m_walltime.o timer_tree.o m_timer.o timer.o
SIESTA_OBJS += m_spin.o m_vee_integrals.o


# Add specific class objects
# class Objects
SIESTA_OBJS += class_OrbitalDistribution.o \
	class_Sparsity.o \
	class_Data1D.o class_Data2D.o \
	class_SpData1D.o class_SpData2D.o \
        class_Fstack_Pair_Geometry_SpData2D.o \
        class_Pair_Geometry_SpData2D.o \
        class_Geometry.o \
        m_uuid.o object_debug.o

# As the classes are based on inputs we force the dependency to
# alloc
class_Sparsity.o class_Data1D.o class_Data2D.o: alloc.o

DEP_OBJS = $(SIESTA_OBJS) $(TS_OBJS) $(SYSOBJ) $(EXTRA)

# To separate unfold objects from siesta/transiesta
UNFOLD_OBJS  = array.o hsx_m.o
UNFOLD_OBJS += unfold.o 

# Handlers and resolvers of last resort
UNFOLD_OBJS += handlers.o
#

#######################

$(DEP_OBJS): $(COMP_LIBS) 

unfold: $(MPI_INTERFACE) $(COMP_LIBS) $(FDF_LIBS) $(NCPS) $(PSOP) $(DEP_OBJS) $(UNFOLD_OBJS) 
	$(FC) $(FFLAGS) $(LDFLAGS) -o unfold \
              $(DEP_OBJS) $(UNFOLD_OBJS) $(MPI_INTERFACE) \
              $(XC) $(NCPS) $(PSOP) $(PSML_LIBS) $(XMLF90_LIBS) $(FDF_LIBS) $(COMP_LIBS) $(LIBS)

.PHONY: clean
clean: 
	@echo "==> Cleaning object, library, and executable files"
	-rm -f unfold *.o *.a *.mod

PROGS:= unfold
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin

# Dependencies
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o \
		$(addprefix $(TOPSRC)/,$(DEP_OBJS:.o=.f) $(DEP_OBJS:.o=.f90)) \
		$(addprefix $(TOPSRC)/,$(DEP_OBJS:.o=.F) $(DEP_OBJS:.o=.F90)) \
		$(TOPSRC)/class*.T90 \
		$(UNFOLD_OBJS:.o=.f90) $(UNFOLD_OBJS:.o=.F90) \
		$(UNFOLD_OBJS:.o=.f) $(UNFOLD_OBJS:.o=.F) || true

	@sed -i -e 's/\.T90\.o:/.o:/g' Makefile

# DO NOT DELETE THIS LINE - used by make depend
