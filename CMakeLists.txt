
cmake_minimum_required(VERSION 3.14)

list(APPEND
   CMAKE_MODULE_PATH
   ${CMAKE_CURRENT_SOURCE_DIR}/config/cmake
  )

include(utils)
ensure_out_of_source_build()

project(
  "siesta"
  LANGUAGES Fortran C
  VERSION "5.1.0"
  DESCRIPTION "A great DFT code"
)

option(WITH_MPI "Whether Siesta should support MPI-parallelism" FALSE)
option(WITH_NETCDF "Use netCDF interface" FALSE)
option(WITH_ELPA "Use the ELPA library (direct interface)" FALSE)
option(WITH_LUA "Use the flook library for Lua scripting" FALSE)
option(WITH_LIBXC "Include libxc support" FALSE)
#
option(WITH_UTILS "Compile all utilities (besides tbtrans)" TRUE)


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

# This might be necessary if we are not compiling libgridxc
# (pending pkg-config support for the libs.required)
#
if(WITH_LIBXC)
  include(search_for_libxc)
endif()  

# Handle external dependencies

option(DOWNLOAD_FALLBACK "Download dependencies if not found" ON)

include(search_for_xmlf90)
include(search_for_libpsml)
include(search_for_gridxc)



# Code organization
# -- Internal library dependencies

add_subdirectory("Src/fdf")
add_subdirectory("Src/simple-fdict")
add_subdirectory("Src/ncps/src")
add_subdirectory("Src/psoplib/src")
add_subdirectory("Src/MatrixSwitch/src")

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
