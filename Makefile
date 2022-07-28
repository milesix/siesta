default: siesta

TOPDIR=.
MAIN_OBJDIR=.
include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

SIESTA_INSTALL_DIRECTORY?=$(MAIN_OBJDIR)/local_install

UTILS = Util/Denchar/Src \
        Util/VCA \
        Util/Eig2DOS \
	Util/Grid \
	Util/WFS \
	Util/HSX \
	Util/DensityMatrix \
	Util/Bands \
	Util/Gen-basis \
	Util/Grimme \
	Util/Helpers \
	Util/COOP \
	Util/Sockets \
	Util/TS/TBtrans \
	Util/TS/tscontour \
	Util/TS/tshs2tshs \
	Util/TS/ts2ts \
	Util/Vibra/Src \
	Util/Macroave/Src \
	Util/STM/ol-stm/Src \
	Util/STM/simple-stm \
	Util/SpPivot \
	Util/Unfolding/Src \
        Pseudo/converters/psml2psf \
        Pseudo/vnl-operator \
        Util/SiestaSubroutine/SimpleTest/Src 

.PHONY: utils clean_utils install_utils siesta $(UTILS)
.PHONY: install_siesta create_install_directory


siesta: $(EXTLIBS) 
	(cd Src; $(MAKE) $(EXTLIBS_SPECS)  siesta )

clean_siesta:
	(cd Src; $(MAKE) clean)
#
utils:  MODE= 
utils:  $(EXTLIBS) $(UTILS)

clean_utils: MODE=clean
clean_utils: $(UTILS)

install_utils: create_install_directory
install_utils: MODE=SIESTA_INSTALL_DIRECTORY="$(SIESTA_INSTALL_DIRECTORY)" install
install_utils: $(EXTLIBS) $(UTILS)

$(UTILS):
	$(MAKE) -C $@ $(MODE) $(EXTLIBS_SPECS)

install_siesta: create_install_directory 
	(cd Src; $(MAKE) SIESTA_INSTALL_DIRECTORY="$(SIESTA_INSTALL_DIRECTORY)" $(EXTLIBS_SPECS) install)

create_install_directory:
	mkdir -p $(SIESTA_INSTALL_DIRECTORY)/bin
	mkdir -p $(SIESTA_INSTALL_DIRECTORY)/lib
	mkdir -p $(SIESTA_INSTALL_DIRECTORY)/include
	mkdir -p $(SIESTA_INSTALL_DIRECTORY)/share

install: install_siesta install_utils

clean: clean_siesta clean_utils
	rm -rf local_install

veryclean: clean clean_extlibs

