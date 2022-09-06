#
# Recipes for compilation of external libraries
#
#
# Credit for mechanism for integration of CMake rules in
# makefile-based framework:
#
# https://stackoverflow.com/questions/62218250/wrapping-cmake-build-with-makefile
#
# Source code for the libraries is handled through git submodules
# deployed in directory ExtLibs
#
#-----------------------------------------
EXTLIBS_INSTALL_PREFIX=$(MAIN_OBJDIR)/ExtLibs_installs
## EXTLIBS_TOOLCHAIN_FILE=$(MAIN_OBJDIR)/extlibs_toolchain.cmake
EXTLIBS_VERSION_CHECK_FILE=$(TOPDIR)/Config/mk-build/extlibs_cmake_version_check.cmake

CMAKE_BUILD_DIR_xmlf90 := $(MAIN_OBJDIR)/__extlib_xmlf90
CMAKE_SOURCE_DIR_xmlf90 := $(TOPDIR)/ExtLibs/xmlf90

check_cmake_version:
	@cmake -P $(EXTLIBS_VERSION_CHECK_FILE)

$(CMAKE_BUILD_DIR_xmlf90)/Makefile: $(CMAKE_SOURCE_DIR_xmlf90)/CMakeLists.txt
	cmake -S $(<D) -B $(@D) \
             -DCMAKE_INSTALL_PREFIX=$(EXTLIBS_INSTALL_PREFIX) \
             -DCMAKE_Fortran_COMPILER=$(FC_SERIAL) \
             -DCMAKE_INSTALL_LIBDIR=$(LIBPREFIX)

.PHONY: $(CMAKE_BUILD_DIR_xmlf90)/libxmlf90.a  # to allow CMake's make check the build
$(CMAKE_BUILD_DIR_xmlf90)/libxmlf90.a: $(CMAKE_BUILD_DIR_xmlf90)/Makefile
	cmake --build $(CMAKE_BUILD_DIR_xmlf90) 
	cmake --install $(CMAKE_BUILD_DIR_xmlf90)

xmlf90: check_cmake_version
xmlf90: $(CMAKE_BUILD_DIR_xmlf90)/libxmlf90.a

#-----------------------------------------
CMAKE_BUILD_DIR_psml := $(MAIN_OBJDIR)/__extlib_psml
CMAKE_SOURCE_DIR_psml := $(TOPDIR)/ExtLibs/libpsml

$(CMAKE_BUILD_DIR_psml)/Makefile: $(CMAKE_SOURCE_DIR_psml)/CMakeLists.txt
	cmake -S $(<D) -B $(@D) \
             -DCMAKE_INSTALL_PREFIX=$(EXTLIBS_INSTALL_PREFIX) \
             -DCMAKE_Fortran_COMPILER=$(FC_SERIAL) \
             -DCMAKE_PREFIX_PATH=$(EXTLIBS_INSTALL_PREFIX)  \
             -DCMAKE_INSTALL_LIBDIR=$(LIBPREFIX)

.PHONY: $(CMAKE_BUILD_DIR_psml)/libpsml.a  # to allow CMake's make check the build
$(CMAKE_BUILD_DIR_psml)/libpsml.a: $(CMAKE_BUILD_DIR_psml)/Makefile
	cmake --build $(CMAKE_BUILD_DIR_psml) 
	cmake --install $(CMAKE_BUILD_DIR_psml)

psml: check_cmake_version
psml: xmlf90 $(CMAKE_BUILD_DIR_psml)/libpsml.a
#-----------------------------------------
CMAKE_BUILD_DIR_gridxc := $(MAIN_OBJDIR)/__extlib_gridxc
CMAKE_SOURCE_DIR_gridxc := $(TOPDIR)/ExtLibs/libgridxc

ifeq ($(WITH_MPI), 1)
  MPI_FLAG=ON
  FORTRAN_COMPILER=$(FC_PARALLEL)
else
  MPI_FLAG=OFF
  FORTRAN_COMPILER=$(FC_SERIAL)
endif
ifeq ($(WITH_LIBXC), 1)
  LIBXC_FLAG=ON
else
  LIBXC_FLAG=OFF
endif
ifeq ($(WITH_GRID_SP), 1)
  SP_FLAG=ON
else
  SP_FLAG=OFF
endif

$(CMAKE_BUILD_DIR_gridxc)/Makefile: $(CMAKE_SOURCE_DIR_gridxc)/CMakeLists.txt
	cmake -S $(<D) -B $(@D) \
             -DCMAKE_INSTALL_PREFIX=$(EXTLIBS_INSTALL_PREFIX) \
             -DCMAKE_Fortran_COMPILER=$(FORTRAN_COMPILER) \
             -DCMAKE_INSTALL_LIBDIR=$(LIBPREFIX) \
             -DCMAKE_PREFIX_PATH="$(EXTLIBS_INSTALL_PREFIX);$(LIBXC_ROOT)" \
             -DWITH_MPI=$(MPI_FLAG) \
             -DWITH_LIBXC=$(LIBXC_FLAG) \
             -DWITH_GRID_SP=$(SP_FLAG)

.PHONY: $(CMAKE_BUILD_DIR_gridxc)/libgridxc.a  # to allow CMake's make check the build
$(CMAKE_BUILD_DIR_gridxc)/libgridxc.a: $(CMAKE_BUILD_DIR_gridxc)/Makefile
	cmake --build $(CMAKE_BUILD_DIR_gridxc) 
	cmake --install $(CMAKE_BUILD_DIR_gridxc)

gridxc: check_cmake_version
gridxc: $(CMAKE_BUILD_DIR_gridxc)/libgridxc.a
#-----------------------------------------

clean_extlibs:
	-rm -rf $(MAIN_OBJDIR)/__extlib_*
	-rm -rf $(EXTLIBS_INSTALL_PREFIX)
