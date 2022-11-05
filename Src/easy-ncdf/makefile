#
#  Makefile for NCDF
#  
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#

TOPDIR=.
MAIN_OBJDIR=.

default: module
module: libncdf.a

include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

VPATH=$(TOPDIR)/Src/easy-ncdf
#
OBJS=netcdf_ncdf.o
$(OBJS): $(FDICT_LIBS)
INCFLAGS += $(FDICT_INCFLAGS)

libncdf.a:  $(OBJS)
	$(AR) $(ARFLAGS_EXTRA) cru libncdf.a $(OBJS)
	-$(RANLIB) libncdf.a

#
# Cleaning
#
clean:
	rm -f *.o *.mod libncdf.a

install: libncdf.a
	cp -p libncdf.a $(SIESTA_INSTALL_DIRECTORY)/lib
