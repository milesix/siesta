.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
.PHONY: libxc-compat
#
default:  libncps.a
module:  libncps.a

TOPDIR=.
MAIN_OBJDIR=.
#
include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

VPATH=$(TOPDIR)/Src/ncps/src

#
libxc-compat:
	@echo "Compiling the libxc-compat library..."
	(cd libxc-compat ; $(MAKE) )
#
OBJFILES= m_ncps.o m_ncps_froyen_reader.o m_ncps_translators.o \
          m_ncps_froyen_ps_t.o m_ncps_reader.o	m_ncps_writers.o \
          interpolation.o m_ncps_psml_plugins.o m_ncps_utils.o
#
INCFLAGS:= -I libxc-compat $(PSML_INCFLAGS)
#
libncps.a: libxc-compat $(OBJFILES)
	$(AR) $(ARFLAGS_EXTRA) cru libncps.a $(OBJFILES) libxc-compat/*.o
	-$(RANLIB) libncps.a
#
install:
	cp -p libncps.a $(SIESTA_INSTALL_DIRECTORY)/lib
#
dep:
	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90 
#
clean:
	rm -f *.o *.*d *.a
	(cd libxc-compat ; $(MAKE) \
                        ARCH_MAKE="$(ARCH_MAKE)" clean)
#
# DO NOT DELETE THIS LINE - used by make depend
m_libxc_sxc_translation.o: m_libxc_list.o
m_ncps.o: m_ncps_froyen_ps_t.o m_ncps_reader.o m_ncps_translators.o
m_ncps.o: m_ncps_writers.o m_ncps_utils.o
m_ncps_froyen_reader.o: interpolation.o m_ncps_froyen_ps_t.o
m_ncps_reader.o: m_ncps_froyen_ps_t.o m_ncps_froyen_reader.o
m_ncps_reader.o: m_ncps_translators.o m_ncps_writers.o
m_ncps_translators.o: m_ncps_froyen_ps_t.o
m_ncps_translators.o: m_ncps_psml_plugins.o
m_ncps_writers.o: m_ncps_froyen_ps_t.o
m_ncps_utils.o: m_ncps_froyen_ps_t.o
