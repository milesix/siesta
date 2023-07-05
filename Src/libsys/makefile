.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
default:  libsys.a
module:  libsys.a

TOPDIR=.
MAIN_OBJDIR=.
#
include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

VPATH=$(TOPDIR)/Src/libsys

#
SRCS= sys.F90 external_entries.f90
# Generate local .o file list
OBJFILES:=$(call change_extensions, $(SRCS))

libsys.a: $(OBJFILES)
	$(AR) $(ARFLAGS_EXTRA) cru libsys.a $(OBJFILES)
	-$(RANLIB) libsys.a
#
install:
	cp -p libsys.a $(SIESTA_INSTALL_DIRECTORY)/lib
#
# Generate sources list with proper placement for sfmakedepend

LOCAL_SRCS:= $(addprefix $(TOPDIR)/Src/libsys/,$(SRCS))
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o $(LOCAL_SRCS) 
#
clean:
	rm -f *.o *.*d *.a
#
# DO NOT DELETE THIS LINE - used by make depend
external_entries.o: sys.o

