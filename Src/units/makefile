.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
default:  libunits.a
module:  libunits.a

TOPDIR=.
MAIN_OBJDIR=.
#
include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

VPATH=$(TOPDIR)/Src/units

#
SRCS=  units.F90  units_legacy.f90  units_codata2018.f90  units_common.f90
# Generate local .o file list
OBJFILES:=$(call change_extensions, $(SRCS))

libunits.a: $(OBJFILES)
	$(AR) $(ARFLAGS_EXTRA) cru libunits.a $(OBJFILES)
	-$(RANLIB) libunits.a
#
install:
	cp -p libunits.a $(SIESTA_INSTALL_DIRECTORY)/lib
#
# Generate sources list with proper placement for sfmakedepend

LOCAL_SRCS:= $(addprefix $(TOPDIR)/Src/units/,$(SRCS))
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o $(LOCAL_SRCS) 
#
clean:
	rm -f *.o *.*d *.a
#
# DO NOT DELETE THIS LINE - used by make depend
units.o: units_codata2018.o units_legacy.o
units_codata2018.o: units_common.o
units_legacy.o: units_common.o
units_codata2018_m.o: units_codata2018.o
units_common_m.o: units_common.o
units_legacy_m.o: units_legacy.o
