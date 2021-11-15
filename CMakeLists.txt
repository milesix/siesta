
cmake_minimum_required(VERSION 3.14)

project(
  "siesta"
  LANGUAGES Fortran C
  VERSION "4.1.6"
  DESCRIPTION "A great DFT code"
)

option(WITH_MPI "Whether Siesta should support MPI-parallelism" FALSE)
option(WITH_NETCDF "Use netCDF interface" FALSE)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# Follow GNU conventions for installing directories
include(GNUInstallDirs)


# General configuration information
add_subdirectory("config")

# -- Internal library dependencies
if (WITH_MPI)
 add_subdirectory("Src/MPI")
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
#add_subdirectory("Testers")
