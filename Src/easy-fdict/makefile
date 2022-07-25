#
#  Makefile for FDICT
#  
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#

TOPDIR=.
MAIN_OBJDIR=.

default: module
module: libfdict.a

include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

VPATH=$(TOPDIR)/Src/easy-fdict
FC:= $(FC_SERIAL)

#
OBJS=variable.o dictionary.o
libfdict.a: $(OBJS)
	$(AR) $(ARFLAGS_EXTRA) cru libfdict.a $(OBJS)
	-$(RANLIB) libfdict.a

dictionary.o: variable.o
#
# Cleaning
#
clean:
	rm -f *.o *.mod libfdict.a

install: libfdict.a
	cp -p libfdict.a $(SIESTA_INSTALL_DIRECTORY)/lib

