#
# Makefile for SIESTA-QMMM 
#
.SUFFIXES: .f .F .o .a  .f90 .F90
#
#
VPATH:=$(shell pwd)/../../Src
OBJDIR=Obj

SIESTA_OBJS_DIR=../../$(OBJDIR)
ARCH_MAKE=../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)

FFLAGS:= $(FFLAGS) -DQMMM_BSC 

INCFLAGS:= -I$(SIESTA_OBJS_DIR) $(INCFLAGS)

LIBS:= libsiesta.a $(LIBS)

SIESTA-QMMM_PATH=$(PWD)
##VPATH=$(SIESTA_OBJS_DIR)

#LD_LIBRARY_PATH:=$(SIESTA_OBJS_DIR):$(LD_LIBRARY_PATH)
#export LD_LIBRARY_PATH
#
default: what libsiesta.a siesta-qmmm
#
######include siesta-qmmm.make
#
what:
	@echo
	@echo "Compilation architecture to be used: ${SIESTA_ARCH}"
	@echo "If this is not what you want, create the right"
	@echo "arch.make file using the models in Sys"
	@echo
	@echo "Hit ^C to abort..."
	@sleep 1

# Note that machine-specific files are now in top Src directory.
OBJS =  m_qmmm_fdf.o siesta_qmmm_options.o qmmm_fixed.o functions.o   \
        qmmm_mneighb.o pbc.o ewald_coeff.o ioxvconstr.o dgesv.o \
	m_fill_auxrho.o qmmm_reinit.o qmmm_timer.o \
        linkatoms.o fsiesta.o \
	pcpot.o centerall.o readall.o subconstr.o \
	subfree.o writeall.o qmmm_lst_blk.o graphite.o solv_assign.o \
	mmforce.o ljef.o solv_ene_fce.o   \
	assignation.o qmmm_mulliken.o elecfield.o qmmm_ioxv.o siesta_qmmm.o \

# FoX build targets:
FoX_configured=$(SIESTA_OBJS_DIR)/FoX/.config
FoX_built=$(SIESTA_OBJS_DIR)/FoX/.FoX
#
# This is how we pick up modules and libraries for FoX:
FoX_FCFLAGS=`$(SIESTA_OBJS_DIR)/FoX/FoX-config --fcflags`
FoX_LIBS=`$(SIESTA_OBJS_DIR)/FoX/FoX-config --libs --wcml`

# And add them to global compiler flags:
INCFLAGS:=$(INCFLAGS) $(FoX_FCFLAGS)
#
# First, it needs to be configured. This may have been done
# by the main SIESTA configure, but in case not:
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
#
#--------------------------------------------------------------
$(FoX_configured):
	(cd $(SIESTA_OBJS_DIR)/FoX; touch arch.make ; \
         CONFIGURE="$(VPATH)/FoX/configure"; \
         $$CONFIGURE VPATH="$(VPATH)/FoX" \
         FC="$(FC_SERIAL)" FCFLAGS="$(FFLAGS)" \
         --enable-wcml $(DUMMY_FOX) || false )
#
# Note we have to make sure to use the same compiler as SIESTA,
# and pick up all the same FFLAGS, and also remember to maybe ask for a dummy library.
# Note also that the "false" clause will stop the 'make' process if the configuration fails.
#
# Then we want to make FoX itself. Like so:
$(FoX_built): $(FoX_configured)
	(cd $(SIESTA_OBJS_DIR)/FoX; $(MAKE))
#
#--------------------------------------------------------------
XMLPARSER=libxmlparser.a
$(XMLPARSER): 
	(cd xmlparser ; $(MAKE) "VPATH=$(VPATH)/xmlparser" module)
#
#--------------------------------------------------------------
XC=libSiestaXC.a
$(XC): 
	(cd SiestaXC ; $(MAKE) \
                    "VPATH_ROOT=$(VPATH)"  \
	            "INCFLAGS=-I../MPI" \
                    "FFLAGS=$(FFLAGS)"  module )
#--------------------------------------------------------------
FDF=libfdf.a
$(FDF): 
#
#
	(cd $(SIESTA_OBJS_DIR)/fdf ; $(MAKE) "VPATH=$(VPATH)/fdf" module)
ALL_OBJS=$(OBJS) $(SYSOBJ)
#
# Interfaces to libraries
#
libmpi_f90.a: 
	@(cd $(SIESTA_OBJS_DIR)/MPI ; $(MAKE) "VPATH=$(VPATH)/MPI")
#
# Libraries that might need to be compiled
#
libblas.a: $(SIESTA_OBJS_DIR)/Libs/blas.f
	@echo "==> Compiling libblas.a in Libs..."
	@(cd $(SIESTA_OBJS_DIR)/Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" libblas.a)
liblapack.a: $(SIESTA_OBJS_DIR)/Libs/lapack.f  Libs/extra_lapack.f
	@echo "==> Compiling liblapack.a in Libs..."
	@(cd $(SIESTA_OBJS_DIR)/Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" liblapack.a)
linalg.a: $(SIESTA_OBJS_DIR)/Libs/blas.f $(SIESTA_OBJS_DIR)Libs/lapack.f
	@echo "==> Compiling linalg.a in Libs..."
	@(cd $(SIESTA_OBJS_DIR)/Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" linalg.a)
dc_lapack.a:  $(SIESTA_OBJS_DIR)/Libs/dc_lapack.f
	@echo "==> Compiling dc_lapack.a in Libs..."
	@(cd $(SIESTA_OBJS_DIR)/Libs ; $(MAKE) "VPATH=$(VPATH)/Libs" dc_lapack.a)

version: version.F90
	@echo
	@echo "==> Incorporating information about present compilation (compiler and flags)"
	@sed "s'SIESTA_ARCH'${SIESTA_ARCH}'g;s'FFLAGS'${FC} ${FFLAGS}'g;s'SIESTA_VERSION'\
                                       $$(cat ${<D}/../version.info)'g" $<  > compinfo.F90
	@($(MAKE) compinfo.o)
	@mv compinfo.o version.o
	#@rm -f compinfo.F90
	@echo

siesta-qmmm: version $(OBJS)
	cd $(SIESTA-QMMM_PATH)
	@if [ -f $(SIESTA_OBJS_DIR)/libfdf.a ]; then ln -sf $(SIESTA_OBJS_DIR)/libfdf.a .; fi
	@if [ -f $(SIESTA_OBJS_DIR)/libxmlparser.a ]; then ln -sf $(SIESTA_OBJS_DIR)/libxmlparser.a .; fi
	@if [ -f $(SIESTA_OBJS_DIR)/libmpi_f90.a ]; then ln -sf $(SIESTA_OBJS_DIR)/libmpi_f90.a .; fi
	@if [ -f $(SIESTA_OBJS_DIR)/dc_lapack.a ]; then ln -sf $(SIESTA_OBJS_DIR)/dc_lapack.a .; fi
	@if [ -f $(SIESTA_OBJS_DIR)/liblapack.a ]; then ln -sf $(SIESTA_OBJS_DIR)/liblapack.a .; fi
	@if [ -f $(SIESTA_OBJS_DIR)/libblas.a ]; then ln -sf $(SIESTA_OBJS_DIR)/libblas.a .; fi
	@if [ -f $(SIESTA_OBJS_DIR)/libSiestaXC.a ]; then ln -sf $(SIESTA_OBJS_DIR)/libSiestaXC.a .; fi
	@if [ -f $(SIESTA_OBJS_DIR)/siesta ]; then ln -sf $(SIESTA_OBJS_DIR)/siesta .; fi
	@if [ -f ./libSiestaXC.a ]; then \
	echo "$(FC) -o siesta-qmmm \
               $(LDFLAGS) $(ALL_OBJS) $(LIBS)  $(FDF) $(WXML) $(XMLPARSER) \
               $(XC) \
               $(MPI_INTERFACE) $(COMP_LIBS) $(FoX_LIBS)"; \
	$(FC) -o siesta-qmmm \
               $(LDFLAGS) $(ALL_OBJS) $(LIBS)  $(FDF) $(WXML) $(XMLPARSER) \
               $(XC) \
               $(MPI_INTERFACE) $(COMP_LIBS) $(FoX_LIBS); \
	else \
	echo "$(FC) -o siesta-qmmm \
               $(LDFLAGS) $(ALL_OBJS) $(LIBS)  $(FDF) $(WXML) $(XMLPARSER) \
               $(MPI_INTERFACE) $(COMP_LIBS) $(FoX_LIBS)"; \
	$(FC) -o siesta-qmmm \
               $(LDFLAGS) $(ALL_OBJS) $(LIBS)  $(FDF) $(WXML) $(XMLPARSER) \
               $(MPI_INTERFACE) $(COMP_LIBS) $(FoX_LIBS); \
	fi

siesta:
	@( cd $(SIESTA_OBJS_DIR); $(MAKE) siesta) 

libsiesta.a:
#	@( cd $(SIESTA_OBJS_DIR); $(MAKE) siesta)
	cd $(SIESTA-QMMM_PATH)
	rm -f libsiesta.a
	ar rvs libsiesta.a $(SIESTA_OBJS_DIR)/*.o
	ar d libsiesta.a siesta.o
	ranlib libsiesta.a

clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f siesta-qmmm siesta *.o  *.a *.pcl *.pc
	rm -f *.mod
	rm -f libsiesta.a

