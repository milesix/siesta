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
# deployed in directory External
#
#-----------------------------------------
EXTLIBS_INSTALL_PREFIX=$(MAIN_OBJDIR)/External_installs
## EXTLIBS_TOOLCHAIN_FILE=$(MAIN_OBJDIR)/extlibs_toolchain.cmake
EXTLIBS_VERSION_CHECK_FILE=$(TOPDIR)/Config/mk-build/extlibs_cmake_version_check.cmake

CMAKE_BUILD_DIR_xmlf90 := $(MAIN_OBJDIR)/__extlib_xmlf90
CMAKE_SOURCE_DIR_xmlf90 := $(TOPDIR)/External/xmlf90

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
CMAKE_BUILD_DIR_libpsml := $(MAIN_OBJDIR)/__extlib_libpsml
CMAKE_SOURCE_DIR_libpsml := $(TOPDIR)/External/libpsml

$(CMAKE_BUILD_DIR_libpsml)/Makefile: $(CMAKE_SOURCE_DIR_libpsml)/CMakeLists.txt
	cmake -S $(<D) -B $(@D) \
             -DCMAKE_INSTALL_PREFIX=$(EXTLIBS_INSTALL_PREFIX) \
             -DCMAKE_Fortran_COMPILER=$(FC_SERIAL) \
             -DCMAKE_PREFIX_PATH=$(EXTLIBS_INSTALL_PREFIX)  \
             -DCMAKE_INSTALL_LIBDIR=$(LIBPREFIX)

.PHONY: $(CMAKE_BUILD_DIR_libpsml)/libpsml.a  # to allow CMake's make check the build
$(CMAKE_BUILD_DIR_libpsml)/libpsml.a: $(CMAKE_BUILD_DIR_libpsml)/Makefile
	cmake --build $(CMAKE_BUILD_DIR_libpsml) 
	cmake --install $(CMAKE_BUILD_DIR_libpsml)

libpsml: check_cmake_version
libpsml: xmlf90 $(CMAKE_BUILD_DIR_libpsml)/libpsml.a
#-----------------------------------------
CMAKE_BUILD_DIR_libgridxc := $(MAIN_OBJDIR)/__extlib_libgridxc
CMAKE_SOURCE_DIR_libgridxc := $(TOPDIR)/External/libgridxc

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

$(CMAKE_BUILD_DIR_libgridxc)/Makefile: $(CMAKE_SOURCE_DIR_libgridxc)/CMakeLists.txt
	cmake -S $(<D) -B $(@D) \
             -DCMAKE_INSTALL_PREFIX=$(EXTLIBS_INSTALL_PREFIX) \
             -DCMAKE_Fortran_COMPILER=$(FORTRAN_COMPILER) \
             -DCMAKE_INSTALL_LIBDIR=$(LIBPREFIX) \
             -DCMAKE_PREFIX_PATH="$(EXTLIBS_INSTALL_PREFIX);$(LIBXC_ROOT)" \
             -DWITH_MPI=$(MPI_FLAG) \
             -DWITH_LIBXC=$(LIBXC_FLAG) \
             -DWITH_GRID_SP=$(SP_FLAG)

.PHONY: $(CMAKE_BUILD_DIR_libgridxc)/libgridxc.a  # to allow CMake's make check the build
$(CMAKE_BUILD_DIR_libgridxc)/libgridxc.a: $(CMAKE_BUILD_DIR_libgridxc)/Makefile
	cmake --build $(CMAKE_BUILD_DIR_libgridxc) 
	cmake --install $(CMAKE_BUILD_DIR_libgridxc)

libgridxc: check_cmake_version
libgridxc: $(CMAKE_BUILD_DIR_libgridxc)/libgridxc.a
#-----------------------------------------
CMAKE_BUILD_DIR_libwannier90 := $(MAIN_OBJDIR)/__extlib_libwannier90
CMAKE_SOURCE_DIR_libwannier90 := $(MAIN_OBJDIR)/wannier90-3.1.0

ifeq ($(WITH_MPI), 1)
  MPI_FLAG=ON
  FORTRAN_COMPILER=$(FC_PARALLEL)
else
  MPI_FLAG=OFF
  FORTRAN_COMPILER=$(FC_SERIAL)
endif

$(CMAKE_SOURCE_DIR_libwannier90)/CMakeLists.txt:
	tar xzf $(WANNIER90_PACKAGE) -C $(MAIN_OBJDIR) wannier90-3.1.0/src \
                                                       wannier90-3.1.0/LICENSE
	(cd $(CMAKE_SOURCE_DIR_libwannier90); \
            patch -p1 -i $(TOPDIR)/External/Wannier/Patches/3.1.0.patch)

$(CMAKE_BUILD_DIR_libwannier90)/Makefile: $(CMAKE_SOURCE_DIR_libwannier90)/CMakeLists.txt
	cmake -S $(<D) -B $(@D) \
             -DCMAKE_INSTALL_PREFIX=$(EXTLIBS_INSTALL_PREFIX) \
             -DCMAKE_Fortran_COMPILER=$(FORTRAN_COMPILER) \
             -DCMAKE_INSTALL_LIBDIR=$(LIBPREFIX) \
             -DLAPACK_LIBRARY=$(LAPACK_LIBS) \
             -DWITH_MPI=$(MPI_FLAG)


.PHONY: $(CMAKE_BUILD_DIR_libwannier90)/libwannier90.a  # to allow CMake's make check the build
$(CMAKE_BUILD_DIR_libwannier90)/libwannier90.a: $(CMAKE_BUILD_DIR_libwannier90)/Makefile
	cmake --build $(CMAKE_BUILD_DIR_libwannier90) 
	cmake --install $(CMAKE_BUILD_DIR_libwannier90)

libwannier90: check_cmake_version
libwannier90: $(CMAKE_BUILD_DIR_libwannier90)/libwannier90.a
#-----------------------------------------

clean_extlibs:
	-rm -rf $(MAIN_OBJDIR)/__extlib_*
	-rm -rf $(MAIN_OBJDIR)/wannier90-3.1.0
	-rm -rf $(EXTLIBS_INSTALL_PREFIX)
