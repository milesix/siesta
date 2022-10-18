
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

set(SIESTA_AUTHOR  "Siesta group")
set(SIESTA_LICENSE "GPLv3")

# Include Siesta specific utility functions
include(SiestaUtils)
siesta_util_ensure_out_of_source_build()

# Default build-type 
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

# -- Options ----------------------------------------
#
option(WITH_MPI "Build Siesta with MPI support" TRUE)

option(WITH_NETCDF "Build Siesta with NetCDF support" TRUE)
option(WITH_NCDF "Use the NCDF library" ${WITH_NETCDF})
option(WITH_NCDF_PARALLEL "Use NCDF_PARALLEL" FALSE)

option(WITH_FFTW "FFTW support (for STM/ol-stm utility)" TRUE)
option(WITH_LIBXC "Include libxc support" FALSE)
option(WITH_GRID_SP "Use single-precision for grid magnitudes" FALSE)

option(WITH_OPENMP "Build Siesta and TBTrans with OpenMP support" FALSE)

option(WITH_ELPA "Use the ELPA library (direct interface)" FALSE)
option(WITH_LUA "Use the flook library for Lua scripting" FALSE)
#--------------------------------------

# Compiler flags
include(flags)

# Follow GNU conventions for installing directories
# This will use the bin/ include/ etc. folder structure
include(GNUInstallDirs)

# We will search for packages and base defaults on these
# Use Pkg-Config if available
find_package(PkgConfig)


if (WITH_MPI)
  if (NOT TARGET MPI::MPI_Fortran)
   # Find MPI
   find_package(MPI OPTIONAL_COMPONENTS Fortran C)
   if( NOT MPI_Fortran_FOUND )
     # To change the setting to not build Siesta with parallel support,
     # use the following line
     #set(WITH_MPI OFF CACHE BOOL "Could not find MPI" FORCE)
     message(STATUS "MPI could not be found by CMake")
     message(STATUS "Try setting MPI_Fortran_COMPILER to an MPI wrapper, or the MPI_HOME env variable")
     message(STATUS "If unsuccessful, set -DWITH_MPI=OFF")
     message(FATAL_ERROR "MPI could not be found by CMake")
   endif()
  endif()
endif()

#---------------------------- OpenMP

if( WITH_OPENMP )

  if(NOT TARGET OpenMP::OpenMP_Fortran)
    find_package(OpenMP COMPONENTS Fortran C)

    if( OpenMP_FOUND )

      if(NOT OpenMP_Fortran_HAVE_OMPLIB_MODULE)
             message(ERROR_FATAL
	             "OpenMP does not have the 'omp_lib' module used by Siesta OpenMP")
      endif()
      
    else()
       message(ERROR_FATAL
               "OpenMP support requested but could not find any OpenMP flags")
    endif()
  endif()
  
else()

  # Create dummy empty targets
  if(NOT TARGET OpenMP::OpenMP_Fortran)
    add_library(OpenMP::OpenMP_Fortran INTERFACE IMPORTED)
  endif()
  if(NOT TARGET OpenMP::OpenMP_C)
    add_library(OpenMP::OpenMP_C INTERFACE IMPORTED)
  endif()

endif()


if (WITH_NETCDF)
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

endif()

#
# Look for BLAS and LAPACK
#
find_package(CustomBlas REQUIRED)
find_package(CustomLapack REQUIRED)
if (WITH_MPI)
  find_package(CustomScalapack REQUIRED)
endif()

include(check_for_linalg_features)

# Handle external dependencies
# Use proper order for xmlf90 and libpsml dependencies

if (NOT TARGET xmlf90::xmlf90)
  find_package(xmlf90)
endif()

if (NOT TARGET libpsml::libpsml)
  find_package(libpsml)
endif()

if(WITH_LIBXC)
  if (NOT TARGET libxc::XC_Fortran)
     include(search_for_libxc)
  endif()
endif()  

if (NOT TARGET libgridxc::libgridxc)
  include(search_for_gridxc)
endif()

# Search for flook
if (WITH_LUA)
  include(search_for_flook)
endif()

if (WITH_ELPA)
 include(search_for_elpa)
endif()

if (WITH_FFTW)
  # For STM/ol-stm only for now
  find_package(FFTW QUIET COMPONENTS DOUBLE_LIB)
  if (FFTW_DOUBLE_LIB_FOUND)
     message(STATUS "Found FFTW Double")
  endif()
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
add_subdirectory("Pseudo")

# Package license files
install(
  FILES
  "COPYING"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)

# add the testsuite

enable_testing()
add_subdirectory(Examples)
if(WITH_MPI)
  add_subdirectory("Util/MPI_test")
endif(WITH_MPI)

add_subdirectory("Src/ncps/examples")

