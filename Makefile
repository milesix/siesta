default: siesta

TOPDIR=/Users/ag/code/GITLAB/Siesta/garalb
MAIN_OBJDIR=/Users/ag/code/GITLAB/Siesta/garalb/_util
include $(MAIN_OBJDIR)/arch.make
include $(MAIN_OBJDIR)/check_for_build_mk.mk

SIESTA_INSTALL_DIRECTORY?=$(MAIN_OBJDIR)/local_install

UTILS = Util/Denchar/Src \
        Util/VCA \
        Util/Eig2DOS \
	Util/Grid \
	Util/Gen-basis \
	Util/COOP \
	Util/TS/TBtrans \
	Util/TS/tscontour \
	Util/TS/tshs2tshs \
	Util/TS/ts2ts \
	Util/Vibra/Src \
	Util/Macroave/Src \
	Util/STM/ol-stm/Src \
	Util/Unfolding/Src \
	Pseudo/converters/psml2psf \
        Pseudo/vnl-operator

.PHONY: utils clean_utils install_utils siesta $(UTILS)
.PHONY: install_siesta create_install_directory

siesta:
	(cd Src; $(MAKE))
clean_siesta:
	(cd Src; $(MAKE) clean)
#
utils:  MODE=
utils:  $(UTILS)

clean_utils: MODE=clean
clean_utils: $(UTILS)

install_utils: create_install_directory
install_utils: MODE=SIESTA_INSTALL_DIRECTORY="$(SIESTA_INSTALL_DIRECTORY)" install
install_utils: $(UTILS)

$(UTILS):
	$(MAKE) -C $@ $(MODE)

install_siesta: create_install_directory 
	(cd Src; $(MAKE) SIESTA_INSTALL_DIRECTORY="$(SIESTA_INSTALL_DIRECTORY)" install)

create_install_directory:
	mkdir -p $(SIESTA_INSTALL_DIRECTORY)/bin

install: install_siesta install_utils
clean: clean_siesta clean_utils
