#
# Makefile for MatrixSwitch in Siesta style
#
.SUFFIXES:
.SUFFIXES: .o .a .f90 .F90

module: libMatrixSwitch.a

ARCH_MAKE_DEFAULT=../../arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)

# Add -DSLAP for MPI compilations only
#
# Use the override command in case we still want
# to pass FPPFLAGS in the command line. (Now it
# is superfluous since we include the arch.make file)
#

ifeq ($(WITH_MPI),1)
  override FPPFLAGS += -DHAVE_LAPACK -DHAVE_SCALAPACK -DCONV -DHAVE_MPI
else
  override FPPFLAGS += -DHAVE_LAPACK -DCONV
endif
#
#
ALLOBJ = MatrixSwitch_ops.o MatrixSwitch_mm_multiply.o MatrixSwitch_m_add.o \
         MatrixSwitch_m_set.o MatrixSwitch_m_copy.o MatrixSwitch_m_register.o \
         MatrixSwitch.o MatrixSwitch_wrapper_params.o MatrixSwitch_wrapper.o

libMatrixSwitch.a: $(ALLOBJ)
	@echo "Current FPPFLAGS: " $(FPPFLAGS)
	$(AR) $(ARFLAGS_EXTRA) cru $@ $(ALLOBJ)
	-$(RANLIB) $@

clean : 
	rm -f *.a *.o *.mod
#
MatrixSwitch.o: MatrixSwitch_ops.o MatrixSwitch_mm_multiply.o MatrixSwitch_m_add.o
MatrixSwitch.o: MatrixSwitch_m_set.o MatrixSwitch_m_copy.o MatrixSwitch_m_register.o
