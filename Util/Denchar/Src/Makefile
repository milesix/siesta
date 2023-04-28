# ---
# Copyright (C) 1996-2021       The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for stand-alone Denchar
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

.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90

TOPDIR=.
MAIN_OBJDIR=.

VPATH=$(TOPDIR)/Util/Denchar/Src:$(TOPDIR)/Src

default: denchar

override WITH_MPI=
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
include $(MAIN_OBJDIR)/check_for_build_mk.mk

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC := $(FC_SERIAL)         # Make it non-recursive
#

INCFLAGS:=$(NCPS_INCFLAGS) $(FDF_INCFLAGS) $(PSML_INCFLAGS) $(INCFLAGS)
INCFLAGS += $(UNITS_INCFLAGS)

DENCHAR_UNDEF  = $(DEFS_PREFIX)-UMPI
DENCHAR_UNDEF += $(DEFS_PREFIX)-UCDF
DENCHAR_UNDEF += $(DEFS_PREFIX)-USIESTA__METIS
DENCHAR_UNDEF += $(DEFS_PREFIX)-UON_DOMAIN_DECOMP

DEFS := $(DEFS) $(DENCHAR_UNDEF)
FPPFLAGS := $(FPPFLAGS) $(DENCHAR_UNDEF)

# Uncomment the following line for debugging support
#FFLAGS=$(FFLAGS_DEBUG)

TOP_SRCS=  precision.F bessph.f chkdim.f dismin.f dot.f m_io.f io.f \
           alloc.F90 files.f parallel.F printmatrix.F schecomm.F \
           class_Geometry.F90 class_Data2D.F90 class_OrbitalDistribution.F90 \
           class_Sparsity.F90 class_SpData2D.F90 class_Pair_Geometry_SpData2D.F90 \
           class_Fstack_Pair_Geometry_SpData2D.F90 class_Data1D.F90 class_SpData1D.F90 \
           sparse_matrices.F90 domain_decom.F spatial.F parallelsubs.F iodm.F memory_log.F90 \
           memory.F pxf.F90 m_fft_gpfa.F radfft.f volcel.f debugmpi.F listsc.f \
           atmparams.f atom_options.F90 interpolation.f90 xml.f radial.f atm_types.f spher_harm.f \
           atmfuncs.f  basis_types.f  chemical.f  basis_io.F  qsort.F  mmio.F  pspltm1.F \
           m_getopts.f90 m_uuid.f90 object_debug.F90 class_Pair_Data1D.F90 class_Fstack_Pair_Data1D.F90

LOCAL_SRCS=  local_sys.f  atompla.f  colinear.f  denchar.f  length.f \
             local_die.F  local_reinit.f  matvect.f  neighb.f  planed.f \
             ranger.f  readpla.f  redata_denchar.f  rhoofr.f  timer_local.f \
             wavofr.f  wrout.f  alloc_handlers.f90


SIESTA_OBJS =  $(call change_extensions, $(TOP_SRCS))
LOCAL_OBJS= $(call change_extensions, $(LOCAL_SRCS))

ALL_OBJS_DENCHAR= $(SIESTA_OBJS) $(LOCAL_OBJS)

denchar: $(NCPS) $(UNITS) $(ALL_OBJS_DENCHAR)
	$(FC) -o denchar \
	       $(LDFLAGS) $(ALL_OBJS_DENCHAR) $(FDF_LIBS) $(NCPS) $(UNITS) $(PSML_LIBS) \
               $(XMLF90_LIBS) $(METIS_LIB)

clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f denchar *.o  *.a
	rm -f *.mod
	rm -f _tmp_deps deps.list _deps.list* protomake*

PROGS:= denchar
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin

# Dependencies
ABS_TOP_SRCS:= $(addprefix $(TOPDIR)/Src/,$(TOP_SRCS))
ABS_LOCAL_SRCS:= $(addprefix $(TOPDIR)/Util/Denchar/Src/,$(LOCAL_SRCS))
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o $(ABS_TOP_SRCS) $(ABS_LOCAL_SRCS) || true

# DO NOT DELETE THIS LINE - used by make depend
atm_types.o: precision.o radial.o
atmfuncs.o: atm_types.o local_sys.o precision.o radial.o spher_harm.o
atom_options.o: local_sys.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: local_sys.o precision.o radial.o xml.o
basis_types.o: alloc.o atmparams.o local_sys.o precision.o
bessph.o: local_sys.o precision.o
chemical.o: local_sys.o parallel.o precision.o
chkdim.o: local_sys.o
class_Fstack_Pair_Data1D.o: class_Pair_Data1D.o class_Pair_Data1D.o
class_Fstack_Pair_Geometry_SpData2D.o: class_Pair_Geometry_SpData2D.o
class_Geometry.o: alloc.o
class_Pair_Data1D.o: class_Data1D.o class_Data1D.o
class_Pair_Geometry_SpData2D.o: class_Geometry.o class_SpData2D.o
class_SpData1D.o: class_Data1D.o class_Data1D.o class_Data1D.o class_Data1D.o
class_SpData1D.o: class_Data1D.o
class_SpData2D.o: class_Data2D.o class_Data2D.o class_Data2D.o
class_Sparsity.o: alloc.o
debugmpi.o: parallel.o
domain_decom.o: alloc.o local_sys.o parallel.o precision.o printmatrix.o
domain_decom.o: schecomm.o sparse_matrices.o
io.o: m_io.o
iodm.o: alloc.o files.o local_sys.o parallel.o parallelsubs.o precision.o
listsc.o: alloc.o
m_io.o: local_sys.o
memory.o: memory_log.o
memory_log.o: m_io.o parallel.o precision.o
parallelsubs.o: domain_decom.o local_sys.o parallel.o precision.o spatial.o
printmatrix.o: alloc.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
radial.o: alloc.o interpolation.o precision.o xml.o
schecomm.o: alloc.o
sparse_matrices.o: alloc.o class_Fstack_Pair_Geometry_SpData2D.o
sparse_matrices.o: class_OrbitalDistribution.o class_SpData1D.o
sparse_matrices.o: class_SpData2D.o class_SpData2D.o class_Sparsity.o
sparse_matrices.o: precision.o
spatial.o: precision.o
spher_harm.o: alloc.o local_sys.o precision.o
xml.o: precision.o
atompla.o: local_sys.o precision.o
colinear.o: precision.o
denchar.o: basis_io.o listsc.o m_getopts.o parallel.o precision.o
length.o: precision.o
local_reinit.o: files.o
matvect.o: precision.o
neighb.o: precision.o
planed.o: precision.o
ranger.o: precision.o
readpla.o: precision.o
redata_denchar.o: precision.o
rhoofr.o: atmfuncs.o chemical.o listsc.o planed.o precision.o
wavofr.o: atmfuncs.o chemical.o planed.o precision.o
wrout.o: precision.o
class_ddata1d.o: class_Data1D.o
class_ddata2d.o: class_Data2D.o
class_dspdata1d.o: class_SpData1D.o
class_dspdata2d.o: class_SpData2D.o
class_fstack_pair_ddata1d.o: class_Fstack_Pair_Data1D.o
class_fstack_pair_gdata1d.o: class_Fstack_Pair_Data1D.o
class_fstack_pair_geometry_dspdata2d.o: class_Fstack_Pair_Geometry_SpData2D.o
class_fstack_pair_sdata1d.o: class_Fstack_Pair_Data1D.o
class_gdata1d.o: class_Data1D.o
class_gspdata1d.o: class_SpData1D.o
class_idata1d.o: class_Data1D.o
class_idata2d.o: class_Data2D.o
class_ispdata1d.o: class_SpData1D.o
class_ispdata2d.o: class_SpData2D.o
class_ldata1d.o: class_Data1D.o
class_lspdata1d.o: class_SpData1D.o
class_pair_ddata1d.o: class_Pair_Data1D.o
class_pair_gdata1d.o: class_Pair_Data1D.o
class_pair_geometry_dspdata2d.o: class_Pair_Geometry_SpData2D.o
class_pair_sdata1d.o: class_Pair_Data1D.o
class_sdata1d.o: class_Data1D.o
class_sspdata1d.o: class_SpData1D.o
class_zdata1d.o: class_Data1D.o
class_zdata2d.o: class_Data2D.o
class_zspdata1d.o: class_SpData1D.o
class_zspdata2d.o: class_SpData2D.o
gpfa_core_dp.o: m_fft_gpfa.o
gpfa_core_sp.o: m_fft_gpfa.o
gpfa_fft.o: m_fft_gpfa.o
listsc_module.o: listsc.o
m_bessph.o: bessph.o
m_object_debug.o: object_debug.o
m_radfft.o: radfft.o
mtprng.o: m_uuid.o
printmat.o: printmatrix.o
sys.o: local_sys.o
