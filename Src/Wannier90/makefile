#
#  Makefile for Wannier90 library (subroutines directly borrowed from Wannier90)
#  
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
.PHONY: check_flags
#
# This makefile can also be used "remotely", so we allow
# for an external specification of the (relative) location 
# of the arch.make file.
#
#
ARCH_MAKE_DEFAULT=../arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
#
# Add items to the INCFLAGS as set by the arch.make file
#
INCFLAGS:=-I$(VPATH) $(INCFLAGS)    # For VPATH operation
#
# This is needed to pick up top-level module files, such as mpi_siesta
#
LOCAL_INCFLAGS=-I../
INCFLAGS:=$(INCFLAGS) $(LOCAL_INCFLAGS)
#
# Include copying operations in rule to make sure that
# they are always performed.
#
# Library module libwann90.a
default: module
module: libwann90.a
	cp libwann90.a ..
	@cp -p *.*d ..     
#
libwann90.a: parameters.o kmesh_proj.o utility.o io_proj.o constants.o 
	$(AR) $(ARFLAGS_EXTRA) cru libwann90.a \
                parameters.o kmesh_proj.o utility.o io_proj.o constants.o 
	-$(RANLIB) libwann90.a

check_flags:
	@echo "In Wannier90, INCFLAGS is: $(INCFLAGS)"
#
io_proj.o: constants.o
kmesh_proj.o: io_proj.o utility.o parameters.o constants.o
utility.o: constants.o io_proj.o
parameters.o: io_proj.o constants.o utility.o

#
# Cleaning
#
clean:
	rm -f *.o *.mod libwann90.a
