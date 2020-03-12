# ---
# Copyright (C) 1996-2021       The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---

# Makefile for tbtrans

.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

.PHONY: default what dep version clean-tbt clean

default: what tbtrans 

TOPDIR=.
MAIN_OBJDIR=.

VPATH=$(TOPDIR)/Util/TS/TBtrans:$(TOPDIR)/Src
TOPSRC=$(TOPDIR)/Src


ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
include $(MAIN_OBJDIR)/check_for_build_mk.mk

INCFLAGS += $(FDF_INCFLAGS)
INCFLAGS += $(MPI_WRAPPERS_INCFLAGS)

# Set the default non-mpi routines
FC_DEFAULT := $(FC)
FC_SERIAL  ?= $(FC_DEFAULT)

# We need to at least add the TBTRANS flag
FPPFLAGS += -DTBTRANS

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

SYSOBJ=$(SYS).o

SIESTA_OBJS= precision.o m_io.o alloc.o memory_log.o memory.o pxf.o \
             moreParallelSubs.o m_mpi_utils.o \
             reclat.o idiag.o minvec.o \
             parallel.o io.o files.o cellsubs.o sorting.o \
             kpoint_convert.o timestamp.o m_wallclock.o spin_subs.o m_spin.o \
             diag.o diag_option.o bloch_unfold.o \
             debugmpi.o find_kgrid.o units.o m_cite.o posix_calls.o \
             runinfo_m.o cli_m.o sys.o

# Siesta timings
SIESTA_OBJS += m_walltime.o timer_tree.o m_timer.o timer.o byte_count.o

# Siesta kpoint-scale
SIESTA_OBJS += get_kpoints_scale.o siesta_geom.o

# Add specific class objects
# class Objects
SIESTA_OBJS += class_OrbitalDistribution.o \
	class_Sparsity.o \
	class_Data1D.o class_Data2D.o \
	class_SpData1D.o class_SpData2D.o \
	class_TriMat.o \
        m_uuid.o object_debug.o

# As the classes are based on inputs we force the dependency to
# alloc
class_Sparsity.o class_Data1D.o class_Data2D.o class_TriMat.o: alloc.o

# Some of NPA modules
SIESTA_OBJS += m_char.o cross.o m_os.o
SIESTA_OBJS += intrinsic_missing.o geom_helper.o io_sparse.o m_sparse.o
SIESTA_OBJS += m_handle_sparse.o m_iodm.o
SIESTA_OBJS += m_integrate.o m_interpolate.o m_mat_invert.o
SIESTA_OBJS += m_iterator.o fdf_extra.o
SIESTA_OBJS += create_Sparsity_SC.o create_Sparsity_Union.o
SIESTA_OBJS += m_region.o m_sparsity_handling.o m_pivot_array.o
SIESTA_OBJS += m_pivot.o m_pivot_methods.o
SIESTA_OBJS += ncdf_io.o

# Geometry modules
SIESTA_OBJS += m_geom_aux.o m_geom_objects.o
SIESTA_OBJS += m_geom_box.o m_geom_coord.o m_geom_square.o m_geom_plane.o

# Broadcast fdf data structure
SIESTA_OBJS += broadcast_fdf_struct.o

# The needed transiesta objects
TS_OBJS  = ts_electrode.o m_ts_method.o m_ts_chem_pot.o
TS_OBJS += m_ts_electrode.o m_ts_elec_se.o m_ts_gf.o
TS_OBJS += m_ts_cctype.o m_ts_sparse_helper.o
TS_OBJS += m_ts_io.o m_ts_iodm.o m_ts_io_ctype.o m_ts_debug.o
TS_OBJS += m_ts_io_contour.o m_gauss_quad.o

# These could be left out, however, there are dependencies
TS_OBJS += m_ts_contour_eq.o m_ts_contour_neq.o
TS_OBJS += m_gauss_fermi_inf.o m_gauss_fermi_30.o m_gauss_fermi_28.o
TS_OBJS += m_gauss_fermi_26.o m_gauss_fermi_24.o m_gauss_fermi_22.o
TS_OBJS += m_gauss_fermi_20.o m_gauss_fermi_19.o m_gauss_fermi_18.o
TS_OBJS += m_gauss_fermi_17.o m_ts_aux.o

# The block-tri-diagonal matrix routines
TS_OBJS += m_ts_tri_scat.o m_trimat_invert.o m_ts_trimat_invert.o
TS_OBJS += m_ts_tri_common.o m_ts_rgn2trimat.o

# For creating the global-unit-cell sparsity pattern
TS_OBJS += m_ts_sparse.o m_ts_tdir.o m_ts_pivot.o

DEP_OBJS = $(SIESTA_OBJS) $(TS_OBJS) $(SYSOBJ) $(EXTRA)

# To separate tbtrans objects from siesta/transiesta
TBTOBJS  = version.o m_verbosity.o

# Output stuff
TBTOBJS += tbt_init_output.o

# Add generic tbtrans routines
TBTOBJS += tbt_init.o tbt_end.o tbt_reinit_m.o m_tbt_hs.o m_tbt_regions.o
TBTOBJS += m_tbt_options.o m_tbt_tri_init.o
TBTOBJS += m_tbt_contour.o m_tbt_kpoint.o m_tbt_gf.o
TBTOBJS += sorted_search.o

# Different solution methods
TBTOBJS += m_tbtrans.o m_tbt_trik.o m_tbt_tri_scat.o
TBTOBJS += m_tbt_proj.o m_tbt_diag.o
TBTOBJS += m_tbt_delta.o m_tbt_dH.o m_tbt_dSE.o

# The saving methods
TBTOBJS += m_tbt_sigma_save.o m_tbt_save.o

# The k-region files
TBTOBJS += m_tbt_kregions.o m_tbt_sparse_helper.o

# The main program
TBTOBJS += tbtrans.o 

# Handlers and resolvers of last resort
TBTOBJS += handlers.o

#
# Extract information of the compiler version used.
# This is very basic but for intel, gcc, pgi and open64
# this at least contain the version string of the compiler.
# Possibly later the used compiler version may be extracted.
COMPILER_VERSION := $(shell $(FC) --version | sed '/^[ ]*$$/d;q')
ifeq ($(COMPILER_VERSION),)
  # Just try and capture the case were it can not figure out anything
  COMPILER_VERSION := Unknown version
endif

# The version program
.PHONY: version.o FORCE
SIESTA_VERSION_FILE:=../../../SIESTA.version
version.o: version
version: version.F90 posix_calls.o $(SIESTA_VERSION_FILE)
	@echo
	@echo "==> Incorporating information about present compilation (compiler and flags)"
	@sed "s'SIESTA_ARCH'${SIESTA_ARCH}'g;s'FFLAGS'${FC} ${FFLAGS}'g;\
		s^COMPILER_VERSION^$(COMPILER_VERSION)^;\
		s'TBTRANS_VERSION'$$(cat ${<D}/$(SIESTA_VERSION_FILE))'g;\
		s'LIBS'$(COMP_LIBS) $(LIBS)'g;\
		s'FPPFLAGS'$(FPPFLAGS)'g" $<  > tmp.F90
	@awk '{if (length>80) { printf "%s&\n",substr($$0,1,77); cur=78; \
		while(length-cur>76){  \
		printf "&%s&\n",substr($$0,cur,76) ; cur=cur+76;\
		} printf "&%s\n",substr($$0,cur)} \
		else { print $$0 }}' tmp.F90 > compinfo.F90
	@($(MAKE) "FPPFLAGS=$(FPPFLAGS)" compinfo.o)
	@mv compinfo.o version.o
	@echo
$(SIESTA_VERSION_FILE): FORCE
	@(cd $(TOPDIR) && ./SIESTA_vgen.sh)


# Interfaces to libraries

#######################
#  TS-TBT libs:
#
# Libraries used by transiesta and tbtrans
# They are currently supplied by this distribution
# All rights reserved for the author of these libraries.
# Author: Nick R. Papior
# (FDICT, NCDF: Handled by build.mk through COMP_LIBS)
#
INCFLAGS += $(FDICT_INCFLAGS) $(NCDF_INCFLAGS)
#######################

# This will setup the correct compilation sequence

$(DEP_OBJS): $(COMP_LIBS) $(FDICT_LIBS)
$(TBTOBJS): $(FDF_LIBS) $(FDICT_LIBS)

.PHONY: check-tbtrans
check-tbtrans:
	@if [ -f .phtrans ]; then \
              echo "** You were building phtrans in this directory." ;\
              echo "** Please type 'make clean' first " ;\
              false ;\
	fi
	@touch .tbtrans
tbtrans: check-tbtrans
.PHONY: tb tbt
tb: tbtrans
tbt: tbtrans


.PHONY: check-phtrans
check-phtrans:
	@if [ -f .tbtrans ]; then \
              echo "** You were building tbtrans in this directory." ;\
              echo "** Please type 'make clean' first " ;\
              false ;\
	fi
	@touch .phtrans
phtrans: check-phtrans
.PHONY: ph pht
ph: phtrans
pht: phtrans



# We append the -DTBTRANS as some special routines
# which could take advantage of this preprocessing flag
tbtrans: version $(MPI_INTERFACE) $(FDF_LIBS) $(COMP_LIBS) $(TBT_COMP_LIBS) $(FDICT_LIBS) $(DEP_OBJS) $(TBTOBJS) 
	$(FC) $(FFLAGS) $(LDFLAGS) -o tbtrans \
              $(DEP_OBJS) $(TBTOBJS) $(MPI_INTERFACE) \
              $(FDF_LIBS) $(COMP_LIBS) $(TBT_COMP_LIBS) $(FDICT_LIBS) $(LIBS)

# We append the -DTBTRANS as some special routines
# which could take advantage of this preprocessing flag
# Tbtrans has TBT_PHONON flags to compile special capabilities
phtrans: FPPFLAGS += -DTBT_PHONON
phtrans: version $(MPI_INTERFACE) $(FDF_LIBS) $(COMP_LIBS) $(TBT_COMP_LIBS) $(FDICT_LIBS) $(DEP_OBJS) $(TBTOBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o phtrans \
              $(DEP_OBJS) $(TBTOBJS) $(MPI_INTERFACE) \
              $(FDF_LIBS) $(COMP_LIBS) $(TBT_COMP_LIBS) $(FDICT_LIBS) $(LIBS)

.PHONY: clean-tbt
clean-tbt:
	@echo "==> Cleaning tb/ph-trans only object, and executable files"
	-rm -f m_ts_chem_pot.o m_ts_elec_se.o m_ts_electype.o m_ts_gf.o
	-rm -f tbtrans phtrans $(TBTOBJS)
	-rm -f .phtrans .tbtrans

.PHONY: clean
clean: 
	@echo "==> Cleaning object, library, and executable files"
	-rm -f tbtrans phtrans *.o *.a *.mod
	-rm -f .phtrans .tbtrans

PROGS:= tbtrans
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin

# Dependencies
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o \
		$(addprefix $(TOPSRC)/,$(DEP_OBJS:.o=.f) $(DEP_OBJS:.o=.f90)) \
		$(addprefix $(TOPSRC)/,$(DEP_OBJS:.o=.F) $(DEP_OBJS:.o=.F90)) \
		$(TOPSRC)/class*.T90 \
		$(TBTOBJS:.o=.f90) $(TBTOBJS:.o=.F90) \
		$(TBTOBJS:.o=.f) $(TBTOBJS:.o=.F) || true

	@sed -i -e 's/\.T90\.o:/.o:/g' Makefile


# DO NOT DELETE THIS LINE - used by make depend
