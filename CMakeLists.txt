
cmake_minimum_required(VERSION 3.14)

list(APPEND
   CMAKE_MODULE_PATH
   ${CMAKE_CURRENT_SOURCE_DIR}/config/cmake
  )

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
  include(search_for_flook)
endif()
#
if (WITH_NETCDF)
 include(search_for_netcdf)
endif()

if (WITH_ELPA)
 include(search_for_elpa)
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
