
cmake_minimum_required(VERSION 3.14)

# CMAKE_MODULE_PATH is locally searched for *.cmake files
# before CMAKE_PREFIX_PATH.
# We will prefer to use CMAKE_MODULE_PATH
list(APPEND
  CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/Config/cmake
  )

# Define project details
project(
  siesta
  # If C++ code is added, we should remember to add it here
  LANGUAGES Fortran C
  VERSION "5.1.0"
  # HOMEPAGE_URL requires >=3.12
  HOMEPAGE_URL "https://siesta-project.org/siesta"
  DESCRIPTION "A first-principles materials simulation code using DFT."
  )

# Default build-type is Release (Debug for now)
option(WITH_LIBXC "Include libxc support" FALSE)
option(WITH_UTILS "Compile all utilities (besides tbtrans)" TRUE)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

include(flags)

# First setup the environment
# Define library specific content
set(SIESTA_AUTHOR  "Siesta group")
set(SIESTA_LICENSE "GPLv3")

# Follow GNU conventions for installing directories
# This will use the bin/ include/ etc. folder structure
include(GNUInstallDirs)


# Include Siesta specific utility functions
include(SiestaUtils)
siesta_util_ensure_out_of_source_build()


# Top-level handlers for compilers
# Ensure all sub-projects gets added a flag for the GCC fortran compiler (for GCC>=10)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # Append -fallow-argument-mismatch for specific versions
  if ( CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10.0 )
    message(STATUS "Using GNU compilers with a version >= 10.0 requires -fallow-argument-mismatch, forcefully added")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif()
endif()


# We will search for packages and base defaults on these
# Use Pkg-Config if available
find_package(PkgConfig)



# Find MPI
find_package(MPI OPTIONAL_COMPONENTS Fortran C)
if( NOT MPI_FOUND )
  message(STATUS "MPI could not be found by the CMake configuration files. Will default to not build Siesta with parallel support.")
  message(WARNING "To enable MPI support redefine MPI_Fortran and MPI_C variables and set -DWITH_MPI=on.")
endif()
option(WITH_MPI "Build Siesta with MPI support" ${MPI_FOUND})


# Find OpenMP
find_package(OpenMP COMPONENTS Fortran C)
option(WITH_OPENMP "Build Siesta with OpenMP support" ${OpenMP_FOUND})
if( WITH_OPENMP )
  if( OpenMP_FOUND )
    # Check that the OpenMP compiler has the omp_lib.mod file
    if(NOT OpenMP_Fortran_HAVE_OMPLIB_MODULE)
      message(ERROR_FATAL "OpenMP does not have the 'omp_lib' module which is required for Siesta OpenMP usage")
    endif()
  else()
    message(ERROR_FATAL "Requesting OpenMP support but could not find any OpenMP flags (default FindOpenMP search)")
  endif()
endif()



# Find NetCDF, some fallbacks are provided for pkg-config
if( (NOT DEFINED NetCDF_ROOT) AND (NOT DEFINED NetCDF_PATH) )
  if( PKG_CONFIG_FOUND )
    # Define NetCDF_PATH if not defined
    # Use pkg-config to retrieve the path
    pkg_get_variable(NetCDF_PATH netcdf prefix)
  else()
    if( DEFINED NETCDF_ROOT )
      set(NetCDF_ROOT ${NETCDF_ROOT})
    elseif( DEFINED NETCDF_PATH )
      set(NetCDF_PATH ${NETCDF_PATH})
    endif()
  endif()
endif()

# Require at least 4.0.0 version (we need >=4)
find_package(NetCDF 4.0.0 OPTIONAL_COMPONENTS Fortran C)
option(WITH_NETCDF "Build Siesta with NetCDF support" ${NetCDF_FOUND})

option(WITH_NCDF "Use the NCDF library" ${WITH_NETCDF})
option(WITH_NCDF_PARALLEL "Use NCDF_PARALLEL" FALSE)



# TODO
option(WITH_ELPA "Use the ELPA library (direct interface)" FALSE)
option(WITH_LUA "Use the flook library for Lua scripting" FALSE)
option(WITH_UTILS "Compile all utilities (besides tbtrans)" TRUE)


# Search for flook
if (WITH_LUA)
  include(search_for_flook)
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

if (WITH_UTILS)
  # For STM/ol-stm
  include(search_for_fftw)
endif()


# For now we add each directory manually.
# This follows the code structure since each
# folder contains a library.
# Each added sub-directory should thus be a directory

if( WITH_MPI )
  # optional mpi-siesta library
  add_subdirectory("Src/MPI")
  # This is the safest to have a simple unified approach
  # I.e. namespaces for modules
  add_library(Siesta::MPI ALIAS mpi_siesta)
endif()

# FDF
add_subdirectory("Src/fdf")
add_library(Siesta::libfdf ALIAS fdf)

add_subdirectory("Src/easy-fdict")
add_subdirectory("Src/ncps/src/libxc-compat")
add_subdirectory("Src/ncps/src")
add_subdirectory("Src/psoplib/src")
add_subdirectory("Src/MatrixSwitch/src")
if (WITH_NCDF)
  message("---- Using the NCDF library")
  add_subdirectory("Src/easy-ncdf")
endif()

# Siesta sources
add_subdirectory("Src")

# -- Utilities
add_subdirectory("Util")
add_subdirectory("Pseudo/vnl-operator")

# Package license files
install(
  FILES
  "COPYING"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)

# add the testsuite
#enable_testing()
add_subdirectory("Src/ncps/examples")

