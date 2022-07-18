#
# This file (build.mk) takes care of the low-level details for
# building.  It needs to be included after the users' arch.make file,
# either with an explicit separate include, or maybe using the
# (uncommented) lines *at the end* of the arch.make file:
#
#SELF_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
#include $(SELF_DIR)build.mk
#------------------------------------------------------------
#
ifdef BUILD_MK_H__
$(info multiple inclusion of build.mk...)
else
  BUILD_MK_H__=1

#
SIESTA_ARCH=build-mk-scheme

# Options and machine specific settings (paths, libraries, etc)
# might be:
#
# 1. Set here, perhaps inherited from environmental variables
#    (See examples below)
#
# 2. Set from a 'arch.make' file that includes this file.
#    This is the recommended operation mode.
#
#===========================================================
# Symbols to request particular features or options
# To turn on, set '=1'.
#--------------
# These are mandatory for PSML-enabled Versions
WITH_PSML=1
WITH_GRIDXC=1
#-------------

# Other options supported:
#
# WITH_MPI=
#
# WITH_ELPA=
# WITH_EXTERNAL_ELPA=    (synonymous with the above)

# WITH_POST_2020_ELPA=
# WITH_POST_2020_EXTERNAL_ELPA=
# GPU_TYPE=NVIDIA_GPU  # or AMD_GPU, or INTEL_GPU (for >=2021 external ELPA)

# WITH_FLOOK=

# WITH_NETCDF=
# WITH_EXPLICIT_NETCDF_SYMBOLS=
# WITH_NCDF=
# WITH_NCDF_PARALLEL=

# WITH_GRID_SP=
# WITH_LEGACY_GRIDXC_INSTALL=

# WITH_EXTRA_FPPFLAGS=

#===========================================================
# Symbols for locating the appropriate libraries
# (Either explicit or through shell variables, perhaps
#  set by a module system)
#--------------------------------------------------------

#XMLF90_ROOT=
#PSML_ROOT=
#GRIDXC_ROOT=
#
#ELPA_ROOT=
#ELPA_INCLUDE_DIRECTORY=
#FLOOK_ROOT=
#
#NETCDF_ROOT=$(NETCDF_HOME)
#NETCDF_LIBS=-lnetcdff -lnetcdf -L/opt/hdf5/lib -lhdf5_hl -lhdf5 -lcurl -lz
#NETCDF_INCFLAGS=-I/usr/include
#
#SCALAPACK_LIBS=-lscalapack
#LAPACK_LIBS=-llapack -lblas
#FFTW_ROOT=/apps/FFTW/3.3.8/GCC/OPENMPI/
# Needed for PEXSI (ELSI) support
#LIBS_CPLUS=-lstdc++ -lmpi_cxx

#===========================================================
# Compiler names and flags
#
# FC_PARALLEL=mpif90
# FC_SERIAL=gfortran
# FFLAGS = -O2
# IPO_FLAG = -ipo  # (keep it separate from FFLAGS)
# FFLAGS_DEBUG= -g -O0
#
# FPP = $(FC_SERIAL) -E -P -x c
# RANLIB=echo

#===========================================================
# Possible section on specific recipes for troublesome files, using
# a lower optimization level.
#
#atom.o: atom.F
#	$(FC) -c $(FFLAGS_DEBUG) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) $< 
#state_analysis.o: 
#create_Sparsity_SC.o:

# Note that simply using target-specific variables, such as:
#atom.o: FFLAGS=$(FFLAGS_DEBUG)
# would compile *all* dependencies of atom.o with that setting...

#===========================================================
# Support for idiosynchratic compilers
#
# In case your compiler does not understand the special meaning of 
# the .F and .F90 extensions ("files in need of preprocessing"), you
# will need to use an explicit preprocessing step.
#WITH_EXPLICIT_FPP = 1

# Explicit Fortran preprocessor. Typically this is sufficient to be the
# compiler with options '-E -P -x c'.
#FPP = $(FC) -E -P -x c

# Some compilers (notably IBM's) are not happy with the standard syntax for
# definition of preprocessor symbols (-DSOME_SYMBOL), and they need a prefix
# (i.e. -WF,-DSOME_SYMBOL). This is used in some utility makefiles. Typically
# this need not be defined.
#DEFS_PREFIX = -WF,

# This enables specific preprocessing options for certain source files.
# (example for the IBM compilers)
#FPPFLAGS_fixed_f = -qsuffix=cpp=F -qfixed
#FPPFLAGS_free_f90 = -qsuffix=cpp=F90 -qfree=F90

#===========================================================
#--------------------------------------------------------
# Nothing should need to be changed below
#--------------------------------------------------------

FC_ASIS=$(FC_SERIAL)

# These are for initialization of variables added to below
FPPFLAGS= $(DEFS_PREFIX)-DF2003 
LIBS=
COMP_LIBS=

# ---- ELPA configuration -----------

# For backwards compatibility, we allow both WITH_ELPA and WITH_EXTERNAL_ELPA to enable
# the native ELPA interface.

ifeq ($(WITH_ELPA),1)
   WITH_EXTERNAL_ELPA=1
endif
ifeq ($(WITH_POST_2020_ELPA),1)
   WITH_POST_2020_EXTERNAL_ELPA=1
endif

ifeq ($(WITH_EXTERNAL_ELPA),1)
   ifndef ELPA_ROOT	
     $(error you need to define ELPA_ROOT in your arch.make)
   endif
   ifndef ELPA_INCLUDE_DIRECTORY
     # It cannot be generated directly from ELPA_ROOT...
     $(error you need to define ELPA_INCLUDE_DIRECTORY in your arch.make)
   endif

   FPPFLAGS_ELPA=$(DEFS_PREFIX)-DSIESTA__ELPA
   ifeq ($(WITH_POST_2020_EXTERNAL_ELPA),1)
     ifndef GPU_TYPE
       $(info NVIDIA_GPU used for post-2020 ELPA kernel interface compilation)
       $(info Set GPU_TYPE if you need another kind)
       $(info or disregard if you do not have gpus)
       GPU_TYPE=NVIDIA_GPU
     endif
     FPPFLAGS_ELPA+=$(DEFS_PREFIX)-DELPA_2STAGE_REAL_GPU=ELPA_2STAGE_REAL_$(GPU_TYPE)
     FPPFLAGS_ELPA+=$(DEFS_PREFIX)-DELPA_2STAGE_COMPLEX_GPU=ELPA_2STAGE_COMPLEX_$(GPU_TYPE)
   endif

   ELPA_INCFLAGS= -I$(ELPA_INCLUDE_DIRECTORY)
   INCFLAGS += $(ELPA_INCFLAGS)
   FPPFLAGS += $(FPPFLAGS_ELPA)
   ELPA_LIB = -L$(ELPA_ROOT)/lib -lelpa
   LIBS +=$(ELPA_LIB) 
endif
# ---- end of ELPA configuration -----------


ifeq ($(WITH_NETCDF),1)

   ifeq ($(WITH_EXPLICIT_NETCDF_SYMBOLS),1)

     ifndef NETCDF_INCFLAGS
      $(error you need to define NETCDF_INCFLAGS in your arch.make)
     endif
     ifndef NETCDF_LIBS
      $(error you need to define NETCDF_LIBS in your arch.make)
     endif

   else

     ifndef NETCDF_ROOT
       $(error you need to define NETCDF_ROOT in your arch.make)
     endif

     NETCDF_INCFLAGS = -I$(NETCDF_ROOT)/include
     NETCDF_LIBS = -L$(NETCDF_ROOT)/lib -lnetcdff
   endif

   FPPFLAGS_CDF = $(DEFS_PREFIX)-DCDF
   FPPFLAGS += $(FPPFLAGS_CDF) 
   INCFLAGS += $(NETCDF_INCFLAGS)
   LIBS += $(NETCDF_LIBS)

   ifeq ($(WITH_NCDF),1)
     FPPFLAGS += $(DEFS_PREFIX)-DNCDF $(DEFS_PREFIX)-DNCDF_4
     COMP_LIBS += libncdf.a libfdict.a

     ifeq ($(WITH_NCDF_PARALLEL),1)
       FPPFLAGS += $(DEFS_PREFIX)-DNCDF_PARALLEL
     endif

   endif
endif

ifeq ($(WITH_FLOOK),1)
 ifndef FLOOK_ROOT
   $(error you need to define FLOOK_ROOT in your arch.make)
 endif
 FLOOK_INCFLAGS=-I$(FLOOK_ROOT)/include
 INCFLAGS += $(FLOOK_INCFLAGS)
 FLOOK_LIBS= -L$(FLOOK_ROOT)/lib -lflookall -ldl
 FPPFLAGS_FLOOK = $(DEFS_PREFIX)-DSIESTA__FLOOK
 FPPFLAGS += $(FPPFLAGS_FLOOK) 
 LIBS += $(FLOOK_LIBS)
 COMP_LIBS += libfdict.a
endif

ifeq ($(WITH_MPI),1)
 FC=$(FC_PARALLEL)
 MPI_INTERFACE=$(MPI_WRAPPERS)
 MPI_INCLUDE=.      # Note . for no-op
 FPPFLAGS_MPI = $(DEFS_PREFIX)-DMPI $(DEFS_PREFIX)-DMPI_TIMING
 LIBS += $(SCALAPACK_LIBS)
 FPPFLAGS += $(FPPFLAGS_MPI) 
else
 FC = $(FC_SERIAL)
endif

ifeq ($(WITH_BUILTIN_LAPACK),1)
  ifeq ($(WITH_MPI),1)
    $(error You should not use the built-in LAPACK routines with Scalapack)
  else
    COMP_LIBS += $(BUILTIN_LAPACK) $(BUILTIN_BLAS)
    FPPFLAGS += -DSIESTA__DIAG_2STAGE
    FPPFLAGS += -DSIESTA__MRRR
  endif
else
  LIBS += $(LAPACK_LIBS)
endif

# ------------- libGridXC configuration -----------

ifeq ($(WITH_GRID_SP),1)
  GRIDXC_CONFIG_PREFIX=sp
  FPPFLAGS_GRID= $(DEFS_PREFIX)-DGRID_SP
else
  GRIDXC_CONFIG_PREFIX=dp
endif
ifeq ($(WITH_MPI),1)
  GRIDXC_CONFIG_PREFIX:=$(GRIDXC_CONFIG_PREFIX)_mpi
endif
FPPFLAGS += $(FPPFLAGS_GRID) 
# -------------------------------------------------


SYS=nag

# These lines make use of a custom mechanism to generate library lists and
# include-file management. The mechanism is not implemented in all libraries.
#---------------------------------------------
ifeq ($(WITH_PSML),1)
 include $(XMLF90_ROOT)/share/org.siesta-project/xmlf90.mk
 include $(PSML_ROOT)/share/org.siesta-project/psml.mk
endif

# A legacy libGridXC installation will have dual 'serial' and 'mpi' subdirectories,
# whereas a modern one, generated with the 'multiconfig' option,  will have split
# include directories but a flat lib directory. The details are still handled by
# appropriate .mk files in the installation directories.
#
# The multiconfig option appeared in 0.9.X, but the legacy compilation option is
# still allowed. For single-precision support with the 'legacy' option, you need to
# make sure that your installation is 'single'...
#
ifeq ($(WITH_GRIDXC),1)
  ifeq ($(WITH_LEGACY_GRIDXC_INSTALL),1)
    include $(GRIDXC_ROOT)/gridxc.mk
  else
    include $(GRIDXC_ROOT)/share/org.siesta-project/gridxc_$(GRIDXC_CONFIG_PREFIX).mk
  endif
endif
#
EXTRA_FPPFLAGS:= $(foreach flag,$(WITH_EXTRA_FPPFLAGS),$(DEFS_PREFIX)-D$(flag))
FPPFLAGS+=$(EXTRA_FPPFLAGS)
#
# Truly external libraries can be phony targets
#
.PHONY: $(PSML_LIBS)
$(PSML_LIBS):
	@echo ">>> PSML_LIBS: $(PSML_LIBS)"
#
.PHONY: $(XMLF90_LIBS)
$(XMLF90_LIBS):
	@echo ">>> XMLF90_LIBS: $(XMLF90_LIBS)"
#
.PHONY: $(GRIDXC_LIBS)
$(GRIDXC_LIBS):
	@echo ">>> GRIDXC_LIBS: $(GRIDXC_LIBS)"
#
# Built-in libraries
#
# Each of these snippets will define XXX_INCFLAGS and XXX_LIBS (or simply XXX, for
# backwards compatibility) for the "internal" libraries. Their use elsewhere is thus
# very similar to that of external libraries.
#
.PHONY: DO_FDF
FDF_LIBS=$(MAIN_OBJDIR)/Src/fdf/libfdf.a
FDF_INCFLAGS=-I$(MAIN_OBJDIR)/Src/fdf
$(FDF_LIBS): DO_FDF
DO_FDF:
	@echo "+++ Compiling internal FDF library"
	(cd $(MAIN_OBJDIR)/Src/fdf ; $(MAKE) -j 1 FFLAGS="$(FFLAGS:$(IPO_FLAG)=)" module)
#--------------------
.PHONY: DO_NCPS
NCPS=$(MAIN_OBJDIR)/Src/ncps/src/libncps.a
NCPS_INCFLAGS=-I$(MAIN_OBJDIR)/Src/ncps/src
$(NCPS): DO_NCPS
DO_NCPS:
	@echo "+++ Compiling internal ncps library"
	(cd $(MAIN_OBJDIR)/Src/ncps/src ; $(MAKE) -j 1 FFLAGS="$(FFLAGS:$(IPO_FLAG)=)" module)
#--------------------
.PHONY: DO_PSOP
PSOP=$(MAIN_OBJDIR)/Src/psoplib/src/libpsop.a
PSOP_INCFLAGS=-I$(MAIN_OBJDIR)/Src/psoplib/src
$(PSOP): DO_PSOP
DO_PSOP:
	@echo "+++ Compiling internal psoplib library"
	(cd $(MAIN_OBJDIR)/Src/psoplib/src ; $(MAKE) -j 1 FFLAGS="$(FFLAGS:$(IPO_FLAG)=)" module)
#--------------------
.PHONY: DO_MS
MS=$(MAIN_OBJDIR)/Src/MatrixSwitch/src/libMatrixSwitch.a
MS_INCFLAGS=-I$(MAIN_OBJDIR)/Src/MatrixSwitch/src
$(MS): DO_MS
DO_MS:
	@echo "+++ Compiling internal MatrixSwitch library"
	(cd $(MAIN_OBJDIR)/Src/MatrixSwitch/src ; \
         $(MAKE) -j 1 FFLAGS="$(FFLAGS:$(IPO_FLAG)=)" module)
#--------------------
.PHONY: DO_BUILTIN_LAPACK
BUILTIN_LAPACK=$(MAIN_OBJDIR)/Src/Libs/libsiestaLAPACK.a
$(BUILTIN_LAPACK): DO_BUILTIN_LAPACK
DO_BUILTIN_LAPACK:
	@echo "+++ Compiling internal LAPACK library"
	(cd $(MAIN_OBJDIR)/Src/Libs ; $(MAKE) -j 1 FFLAGS="$(FFLAGS:$(IPO_FLAG)=)" libsiestaLAPACK.a)
#--------------------
.PHONY: DO_BUILTIN_BLAS
BUILTIN_BLAS=$(MAIN_OBJDIR)/Src/Libs/libsiestaBLAS.a
$(BUILTIN_BLAS): DO_BUILTIN_BLAS
DO_BUILTIN_BLAS:
	@echo "+++ Compiling internal BLAS library"
	(cd $(MAIN_OBJDIR)/Src/Libs ; $(MAKE) -j 1 FFLAGS="$(FFLAGS:$(IPO_FLAG)=)" libsiestaBLAS.a)
#--------------------
.PHONY: DO_MPI_WRAPPERS
MPI_WRAPPERS=$(MAIN_OBJDIR)/Src/MPI/libmpi_f90.a
MPI_WRAPPERS_INCFLAGS=-I$(MAIN_OBJDIR)/Src/MPI
$(MPI_WRAPPERS): DO_MPI_WRAPPERS
DO_MPI_WRAPPERS:
	@echo "+++ Compiling MPI wrappers library"
	(cd $(MAIN_OBJDIR)/Src/MPI ; $(MAKE) -j 1 FFLAGS="$(FFLAGS:$(IPO_FLAG)=)" )


# Define default compilation methods
.c.o:
	$(CC) -c $(CFLAGS) $(INCFLAGS) $(CPPFLAGS) $< 
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F)  $< 
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_free_F90) $< 
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FFLAGS_fixed_f)  $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FFLAGS_free_f90)  $<

endif
