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
	Util/Unfolding/Src \
        Pseudo/converters/psml2psf \
        Pseudo/vnl-operator \
        Util/SiestaSubroutine/SimpleTest/Src 

.PHONY: utils clean_utils install_utils siesta $(UTILS)
.PHONY: install_siesta create_install_directory

#-----------------------------------------
CMAKE_BUILD_DIR_xmlf90 := $(MAIN_OBJDIR)/__extlib_xmlf90
CMAKE_SOURCE_DIR_xmlf90 := $(TOPDIR)/ExtLibs/xmlf90

$(CMAKE_BUILD_DIR_xmlf90)/Makefile: $(CMAKE_SOURCE_DIR_xmlf90)/CMakeLists.txt
	cmake -S $(<D) -B $(@D) -DCMAKE_INSTALL_PREFIX=$(MAIN_OBJDIR)/ExtLibs_installs

.PHONY: $(CMAKE_BUILD_DIR_xmlf90)/libxmlf90.a  # to allow CMake's make check the build
$(CMAKE_BUILD_DIR_xmlf90)/libxmlf90.a: $(CMAKE_BUILD_DIR_xmlf90)/Makefile
	cmake --build $(CMAKE_BUILD_DIR_xmlf90) 
	cmake --install $(CMAKE_BUILD_DIR_xmlf90)

xmlf90: $(CMAKE_BUILD_DIR_xmlf90)/libxmlf90.a

#-----------------------------------------
CMAKE_BUILD_DIR_psml := $(MAIN_OBJDIR)/__extlib_psml
CMAKE_SOURCE_DIR_psml := $(TOPDIR)/ExtLibs/libpsml

$(CMAKE_BUILD_DIR_psml)/Makefile: $(CMAKE_SOURCE_DIR_psml)/CMakeLists.txt
	cmake -S $(<D) -B $(@D) \
             -DCMAKE_INSTALL_PREFIX=$(MAIN_OBJDIR)/ExtLibs_installs \
             -DCMAKE_PREFIX_PATH=$(MAIN_OBJDIR)/ExtLibs_installs

.PHONY: $(CMAKE_BUILD_DIR_psml)/libpsml.a  # to allow CMake's make check the build
$(CMAKE_BUILD_DIR_psml)/libpsml.a: $(CMAKE_BUILD_DIR_psml)/Makefile
	cmake --build $(CMAKE_BUILD_DIR_psml) 
	cmake --install $(CMAKE_BUILD_DIR_psml)

psml: xmlf90 $(CMAKE_BUILD_DIR_psml)/libpsml.a
#-----------------------------------------
CMAKE_BUILD_DIR_gridxc := $(MAIN_OBJDIR)/__extlib_gridxc
CMAKE_SOURCE_DIR_gridxc := $(TOPDIR)/ExtLibs/libgridxc

ifeq ($(WITH_MPI), 1)
  MPI_FLAG=ON
else
  MPI_FLAG=OFF
endif

$(CMAKE_BUILD_DIR_gridxc)/Makefile: $(CMAKE_SOURCE_DIR_gridxc)/CMakeLists.txt
	cmake -S $(<D) -B $(@D) \
             -DCMAKE_INSTALL_PREFIX=$(MAIN_OBJDIR)/ExtLibs_installs \
             -DCMAKE_PREFIX_PATH=$(MAIN_OBJDIR)/ExtLibs_installs \
             -DWITH_MPI=$(MPI_FLAG)

.PHONY: $(CMAKE_BUILD_DIR_gridxc)/libgridxc.a  # to allow CMake's make check the build
$(CMAKE_BUILD_DIR_gridxc)/libgridxc.a: $(CMAKE_BUILD_DIR_gridxc)/Makefile
	cmake --build $(CMAKE_BUILD_DIR_gridxc) 
	cmake --install $(CMAKE_BUILD_DIR_gridxc)

gridxc: $(CMAKE_BUILD_DIR_gridxc)/libgridxc.a
#-----------------------------------------

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
	mkdir -p $(SIESTA_INSTALL_DIRECTORY)/lib
	mkdir -p $(SIESTA_INSTALL_DIRECTORY)/include
	mkdir -p $(SIESTA_INSTALL_DIRECTORY)/share

install: install_siesta install_utils

clean: clean_siesta clean_utils
	rm -rf local_install
