#
#  Makefile for FDF library and FDF example.
#  
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#

TOPDIR=.
MAIN_OBJDIR=.

# Library module libfdf.a
default: module
module: libfdf.a

include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

VPATH=$(TOPDIR)/Src/fdf

#
OBJS=iso_fortran_env.o fdf.o io_fdf.o parse.o utils.o prec.o 
libfdf.a: $(OBJS)
	$(AR) $(ARFLAGS_EXTRA) cru libfdf.a $(OBJS)
	-$(RANLIB) libfdf.a

fdf.o: io_fdf.o parse.o prec.o utils.o
io_fdf.o: utils.o prec.o
parse.o: utils.o prec.o
utils.o: prec.o

#
# Tests for FDF library 
#
sample: sample.F90 libfdf.a
	$(FC) $(LDFLAGS) -o $@ $^ 

#
# Cleaning
#
clean:
	rm -f *.o *.mod libfdf.a
	rm -f sample 
	rm -f fdf.debug* sample.out* sample.fdf.*
