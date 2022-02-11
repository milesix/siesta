#
#  Makefile for FDF library and FDF example.
#  
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
# This makefile can also be used "remotely", so we allow
# for an external specification of the (relative) location 
# of the arch.make file.
#
ARCH_MAKE_DEFAULT=../arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
#
# Include copying operations in rule to make sure that
# they are always performed.
#
# Library module libfdf.a
default: module
module: libfdf.a
	cp libfdf.a ..
	@cp -p *.*d ..     
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
