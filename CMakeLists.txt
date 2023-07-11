cmake_minimum_required(VERSION 3.17)

# CMAKE_MODULE_PATH is locally searched for *.cmake files
# before CMAKE_PREFIX_PATH.
# We will prefer to use CMAKE_MODULE_PATH
list(APPEND
  CMAKE_MODULE_PATH
  # add the generic configuration cmake directory
  "${CMAKE_CURRENT_SOURCE_DIR}/Config/cmake"
  # add the external modules directory
  "${CMAKE_CURRENT_SOURCE_DIR}/Config/cmake/Modules"
  )

# Define project details
project(
  siesta
  LANGUAGES Fortran C
  HOMEPAGE_URL "https://siesta-project.org/siesta"
  DESCRIPTION "A first-principles materials simulation code using DFT."
  )

set(SIESTA_AUTHOR  "Siesta group")
set(SIESTA_LICENSE "GPLv3")

# Include Siesta specific utility functions
include(SiestaUtils)
siesta_util_ensure_out_of_source_build()

# Compiler flags
include(SiestaFlags)

# -- Options ----------------------------------------
option(WITH_GRID_SP "Use single-precision for grid magnitudes" FALSE)

option(BUILD_DOCS "Create Doxygen-based documentation (WIP)" FALSE)
option(WITH_PEXSI "Use the PEXSI library (direct interface)" FALSE)
#--------------------------------------


# Follow GNU conventions for installing directories
# This will use the bin/ include/ etc. folder structure
include(GNUInstallDirs)

# Use Pkg-Config if available
# Many of our sub-inclusions requires this one, so we will catch errors here
# To not confuse users!
find_package(PkgConfig REQUIRED)

# Look for blas, if needed, it is not required if LAPACK ships BLAS
find_package(CustomBlas)

# Look for LAPACK
find_package(CustomLapack REQUIRED)

if(NOT TARGET MPI::MPI_Fortran)
 # Find MPI
 find_package(MPI OPTIONAL_COMPONENTS Fortran C)
endif()

if(MPI_Fortran_FOUND)
  find_package(CustomScalapack)
  if( NOT SCALAPACK_FOUND )
    message(WARNING
      "MPI is found, but ScaLAPACK library cannot be found (or compiled against).\n"
      "If parallel support is required please supply the ScaLAPACK library with appropriate "
      "flags:\n"
      " -DSCALAPACK_LIBRARY=<lib>")
    set(MPI_Fortran_FOUND FALSE)
  endif()
endif()

option(WITH_MPI "Build Siesta with MPI support" ${MPI_Fortran_FOUND})
if( WITH_MPI )

  if(NOT MPI_Fortran_FOUND)
    message(WARNING
      "MPI could not be found by CMake"
      "Try setting MPI_Fortran_COMPILER to an MPI wrapper, or the MPI_HOME env variable"
      "If unsuccessful, set -DWITH_MPI=OFF")
    message(FATAL_ERROR "MPI was requested, but could not be found by CMake")
  endif()
  if( NOT SCALAPACK_FOUND )
    message(WARNING
      "ScaLAPACK could not be found by CMake"
      "Try adding SCALAPACK_LIBRARIES with appropriate flags"
      "If ScaLAPACK is not present MPI needs to be turned off:\n"
      "  -DWITH_MPI=OFF"
      )
    message(FATAL_ERROR "MPI was requested, but ScaLAPACK could not be found by CMake")
  endif()

endif()

# Now we can check for linear algebra stuff
include(SiestaCheckLinalg)
include(SiestaCheckLinalgFeatures)

#---------------------------- OpenMP
# OpenMP *must* be user-requested since it brings additional
# overhead. We should also, for consistency check whether the LAPACK
# and/or BLAS libraries support OpenMP in some form.
# There is little gain when using OpenMP and BLAS libraries are not threaded.
option(WITH_OPENMP "Build with OpenMP support (very few utilities will benefit)" FALSE)
if( WITH_OPENMP )

  if(NOT TARGET OpenMP::OpenMP_Fortran)
    find_package(OpenMP COMPONENTS Fortran C)

    if( OpenMP_FOUND )

      if(NOT OpenMP_Fortran_HAVE_OMPLIB_MODULE)
        message(ERROR_FATAL "OpenMP does not have the 'omp_lib' module used by Siesta OpenMP")
      endif()
      
    else()
      message(ERROR_FATAL "OpenMP support requested but could not find any OpenMP flags")
    endif()
  endif()
  
endif()


# First, a fault check for case-sensitive arguments
# Instead of relying on multiple variants, we should just
# kill the build.
# The problem is that NetCDF *may* require additional
# variables in correct case format. And if we allow upper-case
# for the _ROOT variable, we should do so for all of the variables
# which turns out to be a mess... (we can do a loop if we have everything,
# but for now lets not do this.
foreach(name IN ITEMS "PATH" "ROOT")
  if( (DEFINED NETCDF_${name}) AND (NOT DEFINED NetCDF_${name}) )
    message(WARNING
      "Searching for NetCDF requires the correct case of the search."
      "The build script detected that NETCDF_${name} is defined whereas "
      "NetCDF_${name} is not defined.\n"
      "The search for the NetCDF library requires that you use the "
      "NetCDF_* format for variables. E.g. do this:\n"
      " -DNetCDF_${name}=${NETCDF_${name}}\n"
      "on the command-line instead.")
    message(FATAL_ERROR "Wrong case-notation of NetCDF search variables.")
  endif()
endforeach()
find_package(NetCDF 4.0.0 COMPONENTS Fortran C)


# TODO catch WITH_NetCDF
option(WITH_NETCDF "Build Siesta with NetCDF support" ${NetCDF_Fortran_FOUND})
if( WITH_NETCDF )
  # One should not be able to request NCDF without NetCDF
  # either they should be fault-checked, or WITH_NETCDF enabled if WITH_NCDF is true
  option(WITH_NCDF "Use the NCDF library" ${WITH_NETCDF})
  # Only if netcdf-supports it will we add it.
  option(WITH_NCDF_PARALLEL "Use NCDF_PARALLEL" ${NetCDF_PARALLEL})
else()
  set(WITH_NCDF FALSE CACHE BOOL "Use the NCDF library" FORCE)
  set(WITH_NCDF_PARALLEL FALSE CACHE BOOL "Use parallel NCDF library" FORCE)
endif()
if (WITH_NETCDF AND (NOT NetCDF_Fortran_FOUND))
  message(WARNING
    "NetCDF_Fortran is requested but could not be found by CMake.\n"
    "Try setting variables these variables for searching:\n"
    " - NetCDF_ROOT or NetCDF_PATH\n"
    " - NetCDF_INCLUDE_DIR or NetCDF_Fortran_INCLUDE_DIR\n"
    " - NetCDF_INCLUDE_DIRS or NetCDF_Fortran_INCLUDE_DIRS\n"
    "to appropriate variables for discovering NetCDF fortran libraries."
    "Another possibility would be to add the NetCDF cmake package path to\n"
    "  CMAKE_PREFIX_PATH\n"
    "for automatic discovery."
    "If still unsuccessful, please disable NetCDF support with:\n"
    "  -DWITH_NETCDF=OFF\n")
  message(FATAL_ERROR "NetCDF could not be found by CMake")
endif()


# Handle external dependencies
# Use proper order for xmlf90 and libpsml dependencies
if (NOT TARGET libfdf::libfdf)
  find_package(Customlibfdf REQUIRED)
endif()

if (NOT TARGET xmlf90::xmlf90)
  find_package(Customxmlf90 REQUIRED)
endif()

if (NOT TARGET libpsml::libpsml)
  find_package(Customlibpsml REQUIRED)
endif()
# ensure libpsml links with xmlf90
target_link_libraries(libpsml::libpsml INTERFACE xmlf90::xmlf90)

if(NOT LIBPSML_USES_PROCEDURE_POINTER )

  # Work around bug in shared library linking
  # This should only temporarily be a hack since future
  # psml versions will always rely on procedure pointers
  add_library(${PROJECT_NAME}-psml_wrappers
    STATIC
    "${CMAKE_CURRENT_SOURCE_DIR}/Src/psml_wrappers.F90"
  )

  # One cannot use OBJECT libraries as their objects are not
  # passed down to subsequent dependent usages. I.e. that would cause
  # name-collisions for static libraries.
  target_link_libraries(libpsml::libpsml
    INTERFACE
    ${PROJECT_NAME}-psml_wrappers
  )

endif()

if (NOT TARGET Libxc::xc_Fortran)
  find_package(CustomLibxc)
  if (TARGET Libxc::xc_Fortran)
    message(VERBOSE "Defining new target libxc::XC_Fortran...")
    add_library(libxc::XC_Fortran INTERFACE IMPORTED)
    target_link_libraries(libxc::XC_Fortran INTERFACE Libxc::xc_Fortran)
  endif()
endif()
option(WITH_LIBXC "Include libxc support" ${LIBXC_Fortran_FOUND})
if( WITH_LIBXC AND (NOT LIBXC_Fortran_FOUND) )
  message(WARNING
    "Libxc support has been requested ${WITH_LIBXC} : ${LIBXC_Fortran_FOUND}, but the fortran libxc library cannot be found.")
  message(FATAL_ERROR "Libxc could not be found by CMake")
endif()

if (NOT TARGET libgridxc::libgridxc)
  find_package(CustomLibGridxc REQUIRED)
endif()
if(NOT LIBGRIDXC_USES_PROCEDURE_POINTER)

  add_library(${PROJECT_NAME}-gridxc_wrappers
    STATIC
    "${CMAKE_CURRENT_SOURCE_DIR}/Src/gridxc_wrappers.F90"
  )

  target_link_libraries(libgridxc::libgridxc
    INTERFACE
    ${PROJECT_NAME}-gridxc_wrappers
  )

endif()

option(WITH_DFTD3  "Support for DFT-D3 corrections" TRUE)
if(WITH_DFTD3)
  add_subdirectory("External/DFTD3")
  if(NOT TARGET s-dftd3::s-dftd3)
    message(WARNING "Cannot configure DFT-D3 module. Check sources in External/DFTD3")
    set(WITH_DFTD3  FALSE CACHE BOOL "Support for DFT-D3 corrections" FORCE)
  endif()
endif()

# Search for DFTD3
# TODO use the utility for locating it, but then revert here
#### siesta_add_subdirectory_option(
####   DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/External/DFTD3/"
####   HELP "Include support for DFT-D3 corrections"
####   NITEMS 2
####   OPTION WITH_DFTD3
####   )
# TODO add a fallback for failure in DFTD3

# Search for flook
siesta_add_subdirectory_option(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/External/Lua-Engine/"
  HELP "Include support for Lua engine through flook (custom MD + SCF interaction)"
  NITEMS 2
  OPTION WITH_FLOOK
  )


if( WITH_MPI )
  # no need to search for ELPA in non-mpi environments
  find_package(CustomElpa)
  option(WITH_ELPA "Include ELPA support" ${ELPA_FOUND})
endif()

if (WITH_PEXSI)

   # We want to support PEXSI-2.1, which has some issues
   # with target installation. Hence we provide a "Custom" finder.
   # If we are dealing with PEXSI-2.0, a proper cmake package is
   # available

   if (PEXSI_2_0)
     find_package(PEXSI CONFIG REQUIRED)
     enable_language(CXX)
     # Avoid to link everything associated with PEXSI with C++
     set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 0)
   else()
     find_package(CustomPEXSI REQUIRED)
     enable_language(CXX)
   endif()

   # Some of the following settings might apply only to 2.0,
   # since 2.1 discovery is through library archives and not
   # through a target that might have some properties set

   if (WITH_MPI)
    if (NOT TARGET MPI::MPI_CXX)
      # Find MPI
      find_package(MPI OPTIONAL_COMPONENTS CXX)
      if( NOT MPI_CXX_FOUND )
        message(FATAL_ERROR "MPI_CXX, needed by PEXSI, could not be found by CMake")
      endif()
    endif()
   endif()

endif()

# For STM/ol-stm only for now
find_package(FFTW QUIET COMPONENTS DOUBLE_LIB)
option(WITH_FFTW "FFTW support (currently only STM/ol-stm utility)" ${FFTW_DOUBLE_LIB_FOUND})


if( WITH_MPI )
  # optional mpi-siesta library
  add_subdirectory("Src/MPI")
  # This is the safest to have a simple unified approach
  # I.e. namespaces for modules
  # To be done for all
  add_library(Siesta::MPI ALIAS mpi_siesta)
endif()

add_subdirectory("Src/libsys")
add_subdirectory("Src/easy-fdict")
add_subdirectory("Src/ncps/src/libxc-compat")
add_subdirectory("Src/ncps/src")
add_subdirectory("Src/psoplib/src")
add_subdirectory("Src/MatrixSwitch/src")
if (WITH_NCDF)
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

# Docs
# Enable BUILD_DOCS to create Doxygen-based documentation (wip)
add_subdirectory(Docs)

# add a few very basic tests for now

enable_testing()

include(SiestaTestMPIConfig)

if(WITH_MPI)
  add_subdirectory("Util/MPI_test")
endif()

add_subdirectory(Tests)

include(SiestaPrintBuildInfo)

#add_subdirectory("Src/ncps/examples")
