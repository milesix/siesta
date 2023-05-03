# ---
# Copyright (C) 1996-2021	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for Ordejon-Lorente STM program
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

.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

TOPDIR=.
MAIN_OBJDIR=.

default: what stm

override WITH_MPI=
override WITH_NCDF=

VPATH=$(TOPDIR)/Util/STM/ol-stm/Src:$(TOPDIR)/Util/Grid:$(TOPDIR)/Src

ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
include $(MAIN_OBJDIR)/check_for_build_mk.mk


#------------- NOTE
# Special needs for the use of FFTW
# You must set these symbols properly, either here or 
# (better) in your top arch.make
# 
#FFTW_ROOT ?=/usr/local
FFTW_INCFLAGS ?= -I$(FFTW_ROOT)/include
FFTW_LIBS ?= -L$(FFTW_ROOT)/lib -lfftw3f -lfftw3
#-----------------------------------------

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

STM_UNDEF  = $(DEFS_PREFIX)-UMPI $(DEFS_PREFIX) -UCDF
STM_UNDEF += $(DEFS_PREFIX)-UCDF
STM_UNDEF += $(DEFS_PREFIX)-USIESTA__METIS
STM_UNDEF += $(DEFS_PREFIX)-UON_DOMAIN_DECOMP

#
DEFS := $(DEFS) $(STM_UNDEF)
FPPFLAGS := $(FPPFLAGS) $(STM_UNDEF)

#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)

what:
	@echo
	@echo "Compilation architecture to be used: ${SIESTA_ARCH}"
	@echo "If this is not what you want, create the right"
	@echo "arch.make file using the models in Sys"
	@echo
	@echo "Hit ^C to abort..."
	@if [ -z "$(FFTW_LIBS)" ] ; then \
           echo "You must set the FFTW symbols in your arch.make" ; \
           false ; \
        fi
	@sleep 2

SYSOBJ=$(SYS).o

# Note that machine-specific files are now in top Src directory.

INCFLAGS+=$(NETCDF_INCFLAGS) \
          $(NCPS_INCFLAGS) \
          $(PSML_INCFLAGS) $(FDF_INCFLAGS) 

INCFLAGS+= $(FFTW_INCFLAGS) 
INCFLAGS+= $(UNITS_INCFLAGS) 
#

LOCAL_SRCS=  local_sys.f  fftw3_mymod.f90  extrapolate.f  sts.f  mainstm.F \
  neighb.f  ranger.f  readstm.f  redata.f  reinitstm.f  stm.f  vacpot.f  handlers.f

SIESTA_SRCS=   precision.F  bessph.f  chkdim.f  dismin.f  dot.f \
  m_io.f  io.f  alloc.F90  files.f  parallel.F  printmatrix.F \
  schecomm.F  class_Geometry.F90  class_Data2D.F90 \
  class_OrbitalDistribution.F90  class_Sparsity.F90  class_SpData2D.F90 \
  class_Pair_Geometry_SpData2D.F90  class_Fstack_Pair_Geometry_SpData2D.F90 \
  class_Data1D.F90  class_SpData1D.F90  sparse_matrices.F90 \
  domain_decom.F  spatial.F  parallelsubs.F  iodm.F  memory_log.F90 \
  memory.F  pxf.F90  m_fft_gpfa.F  radfft.f  volcel.f  listsc.f \
  atmparams.f  atom_options.F90  interpolation.f90  xml.f  radial.f \
  atm_types.f  spher_harm.f  atmfuncs.f  basis_types.f  chemical.f \
  basis_io.F  qsort.F  mmio.F  pspltm1.F  m_uuid.f90  object_debug.F90 \
  class_Pair_Data1D.F90  class_Fstack_Pair_Data1D.F90

GRID_SRCS= m_gridfunc.F90

LOCAL_OBJS = $(call change_extensions, $(LOCAL_SRCS))
SIESTA_OBJS = $(call change_extensions, $(SIESTA_SRCS))
GRID_OBJS = $(call change_extensions, $(GRID_SRCS))

ALL_OBJS_STM= $(SIESTA_OBJS) $(LOCAL_OBJS) $(GRID_OBJS)

# m_gridfunc is compiled through the Util/Grid VPATH component
vacpot.o: m_gridfunc.o

stm: $(NCPS) $(ALL_OBJS_STM) $(UNITS)
	$(FC) -o stm \
	       $(LDFLAGS) $(ALL_OBJS_STM) $(UNITS) $(FDF_LIBS) $(NCPS) \
               $(PSML_LIBS) $(XMLF90_LIBS) $(FFTW_LIBS) 
clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f stm *.o  *.a
	rm -f *.mod

PROGS:= stm
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin
	cp -p stm $(SIESTA_INSTALL_DIRECTORY)/bin/wfs2ldos


ABS_TOP_SRCS:= $(addprefix $(TOPDIR)/Src/,$(SIESTA_SRCS))
ABS_LOCAL_SRCS:= $(addprefix $(TOPDIR)/Util/STM/ol-stm/Src/,$(LOCAL_SRCS))
ABS_GRID_SRCS:= $(addprefix $(TOPDIR)/Util/Grid/,$(GRID_SRCS))

.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o \
            $(ABS_TOP_SRCS) $(ABS_LOCAL_SRCS) $(ABS_GRID_SRCS) || true


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
extrapolate.o: fftw3_mymod.o precision.o
handlers.o: memory_log.o
mainstm.o: basis_io.o listsc.o precision.o sts.o
ranger.o: precision.o
stm.o: atmfuncs.o chemical.o precision.o
sts.o: atmfuncs.o chemical.o precision.o
vacpot.o: m_gridfunc.o precision.o
broadening.o: sts.o
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
