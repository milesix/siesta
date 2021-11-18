
cmake_minimum_required(VERSION 3.14)

project(
  "siesta"
  LANGUAGES Fortran C
  VERSION "4.1.6"
  DESCRIPTION "A great DFT code"
)

option(WITH_MPI "Whether Siesta should support MPI-parallelism" FALSE)
option(WITH_NETCDF "Use netCDF interface" FALSE)
option(WITH_ELPA "Use the ELPA library (direct interface)" FALSE)
option(WITH_LUA "Use the flook library for Lua scripting" FALSE)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# Follow GNU conventions for installing directories
include(GNUInstallDirs)


# General configuration information
add_subdirectory("config")

# The MPI interfaces (or lack thereof) are handled
# in a 'private library' for convenience.
#
if (WITH_MPI)
 add_subdirectory("Src/MPI")
endif()

# Search for flook
if (WITH_LUA)
 if("$ENV{FLOOK_ROOT}" STREQUAL "")
   message(SEND_ERROR "Cannot find FLOOK_ROOT")
 else()
  set(FLOOK_ROOT $ENV{FLOOK_ROOT})
  message("Found FLOOK_ROOT: ${FLOOK_ROOT}")
  set(FLOOK_LIBS "-L${FLOOK_ROOT}/lib -lflookall -ldl")
  set(FLOOK_INCLUDE_DIRS "${FLOOK_ROOT}/include")
  message("FLOOK_LIBS: ${FLOOK_LIBS}")
  message("FLOOK_INCLUDE_DIRS: ${FLOOK_INCLUDE_DIRS}")
 endif()
endif()

#
# -- Search for netCDF using pkg-config data
# -- There should be a fallback
#
if (WITH_NETCDF)
 find_package(PkgConfig REQUIRED)
 pkg_check_modules(NETCDF REQUIRED netcdf-fortran)

 message("Netcdf Libdir: ${NETCDF_LIBDIR}")
 message("Netcdf Libraries: ${NETCDF_LIBRARIES}")
 message("Netcdf Link libraries: ${NETCDF_LINK_LIBRARIES}")
 message("NETCDF_INCLUDEDIR: ${NETCDF_INCLUDEDIR}")
 message("NETCDF_INCLUDE_DIRS: ${NETCDF_INCLUDE_DIRS}")
endif()


#
#
# -- Search for ELPA using pkg-config data
# (This, with further refinemens,
#  should go in a module in the cmake subdir)
# 
# There are several non-standard issues:
#
# 1. The pkgconfig filename is typicall of the form
#
#     elpa-2020.05.001.rc1.pc
#
# instead of simply 'elpa.pc'.
#
# 2. The include directory entry in the pkgconfig file
#    is missing the 'modules' leaf component.
#
# The most expedient solution for (1) is to put a symbolic
# link to the desired version somewhere in the PKG_CONFIG_PATH. For example:
#
#   ln -sf /path/to/original/.pc/file $HOME/lib/pkgconfig/elpa.pc
#
#   export PKG_CONFIG_PATH=$HOME/lib/pkgconfig:${PKG_CONFIG_PATH}
#
# (Note: CMAKE_PREFIX_PATH can also be used, but without the
#  '/lib/pkgconfig' part.
#
# A solution for (2) is encoded below.
#
# 
if (WITH_ELPA)
 find_package(PkgConfig REQUIRED)
 pkg_check_modules(ELPA REQUIRED elpa)

 message("Netcdf Libdir: ${ELPA_LIBDIR}")
 message("Netcdf Libraries: ${ELPA_LIBRARIES}")
 message("Netcdf Link libraries: ${ELPA_LINK_LIBRARIES}")
 message("ELPA_INCLUDEDIR: ${ELPA_INCLUDEDIR}")
 message("ELPA_INCLUDE_DIRS: ${ELPA_INCLUDE_DIRS}")
 #
 # Fix non-standard setting
 #
 set(ELPA_INCFLAGS "${ELPA_INCLUDE_DIRS}/modules")
endif()

# Code organization
# -- Internal library dependencies

add_subdirectory("Src/fdf")
add_subdirectory("Src/wxml")
add_subdirectory("Src/SiestaXC")
add_subdirectory("Src/xmlparser")
add_subdirectory("Src/simple-fdict")

if (WITH_NETCDF)
  add_subdirectory("Src/simple-ncdf")
endif()

# -- Siesta executable
add_subdirectory("Src")

# -- Utilities

add_subdirectory("Util")

# Package license files
install(
  FILES
  "COPYING"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)

# add the testsuite
#enable_testing()
#add_subdirectory(" ... ")
