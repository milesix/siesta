#
# Top-level makefile for Siesta plus auxiliary programs (Util, Pseudo)
# 
# Simply typing "make" will build the Siesta executable.
# The whole set of auxiliary programs can be built by "make utils".
#
# Installation: 'make install_siesta' or 'make install_utils'
#
# Two auxiliary targets, "extlibs" and "create_install_directory" are
# called implicitly by the targets above. They can also be called
# directly to carry out the pre-compilation of the dependencies needed
# and to set up the installation directory, respectively. Once these
# two tasks are done, individual auxiliary programs can be built and
# installed independently by going to the relevant directory and
# typing "make" or "make install".
#

default: siesta

TOPDIR=.
MAIN_OBJDIR=.
include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

UTILS = Util/Unfolding/Src \
	Util/Macroave/Src \
        Util/VCA \
	Util/Grid \
	Util/WFS \
        Util/Eig2DOS \
	Util/COOP \
        Util/Denchar/Src \
	Util/DensityMatrix \
        Util/TS/TBtrans \
	Util/Optical \
	Util/Optimizer \
        Util/Contrib/APostnikov \
        Util/Contour \
	Util/JobList/Src \
	Util/pdosxml \
	Util/Sockets \
        Pseudo/converters/psml2psf \
        Pseudo/vnl-operator

# These pending...
UTILS_PENDING =	Util/HSX \
	Util/Bands \
	Util/Gen-basis \
	Util/Grimme \
	Util/Helpers \
	Util/TS/tscontour \
	Util/TS/tshs2tshs \
	Util/TS/ts2ts \
	Util/Vibra/Src \
	Util/ON \
	Util/STM/ol-stm/Src \
	Util/STM/simple-stm \
	Util/SpPivot \
        Util/SiestaSubroutine/SimpleTest/Src \
	Util/SiestaSubroutine/FmixMD/Src \
	Util/SiestaSubroutine/ProtoNEB/Src \
	Util/Projections

.PHONY: utils clean_utils install_utils siesta $(UTILS)
.PHONY: install_siesta create_install_directory

extlibs: $(EXTLIBS)

siesta: $(EXTLIBS) 
	(cd Src; $(MAKE) siesta )

clean_siesta:
	(cd Src; $(MAKE) clean)
#
utils:  MODE= 
utils:  $(EXTLIBS) $(UTILS)

clean_utils: MODE=clean
clean_utils: $(UTILS)

install_utils: create_install_directory
install_utils: MODE=install
install_utils: $(EXTLIBS) $(UTILS)

$(UTILS):
	@echo "This revision cannot build all the programs in Util"
	@echo "using the makefile-based build system".
	@echo "Only TBtrans and Unfold are enabled..."
	@echo "Use the CMake build system if you want a full installation."
	@echo
	@echo "Hit ^C to abort..."
	@sleep 1
	$(MAKE) -C $@ $(MODE)

install_siesta: create_install_directory 
	(cd Src; $(MAKE) install)

create_install_directory:
	mkdir -p $(SIESTA_INSTALL_DIRECTORY)/bin
	mkdir -p $(SIESTA_INSTALL_DIRECTORY)/lib
	mkdir -p $(SIESTA_INSTALL_DIRECTORY)/include
	mkdir -p $(SIESTA_INSTALL_DIRECTORY)/share

install: install_siesta install_utils

clean: clean_siesta clean_utils
	rm -rf local_install

veryclean: clean clean_extlibs

