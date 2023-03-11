.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
#
default:  libpsop.a
module: libpsop.a
#
TOPDIR=.
MAIN_OBJDIR=.
#
include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

VPATH=$(TOPDIR)/Src/psoplib/src

OBJFILES= m_psop.o m_kbgen.o m_localgen.o schrodinger.o arw.o \
          psop_params.o 
#
INCFLAGS=
#
libpsop.a: $(OBJFILES)
	$(AR) $(ARFLAGS_EXTRA) cru libpsop.a $(OBJFILES) 
	-$(RANLIB) libpsop.a
#
dep:
	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90 
#
clean:
	rm -f *.o *.*d *.a
#
# DO NOT DELETE THIS LINE - used by make depend
m_kbgen.o: psop_params.o schrodinger.o
m_psop.o: m_kbgen.o m_localgen.o psop_params.o
schrodinger.o: psop_params.o
schrodinger_m.o: schrodinger.o
