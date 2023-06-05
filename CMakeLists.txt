
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
  LANGUAGES Fortran C CXX
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
  set(CMAKE_BUILD_TYPE "ReleaseWithDebInfo")
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
option(WITH_ELSI "Use the ELSI library interface" FALSE)
option(WITH_LUA "Use the flook library for Lua scripting" FALSE)

option(WITH_DFTD3 "Include support for DFT-D3" FALSE)

# 
option(BUILD_DOCS "Create Doxygen-based documentation (wip)" FALSE)
#--------------------------------------

# Compiler flags
include(flags)

# Follow GNU conventions for installing directories
# This will use the bin/ include/ etc. folder structure
include(GNUInstallDirs)

# Use Pkg-Config if available
find_package(PkgConfig)

if (WITH_MPI)
  if (NOT TARGET MPI::MPI_Fortran)
   # Find MPI
   find_package(MPI OPTIONAL_COMPONENTS Fortran C CXX)
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
  if(NOT TARGET NetCDF::NetCDF_Fortran)
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
    if( NOT NetCDF_Fortran_FOUND )
     message(STATUS "NetCDF_Fortran could not be found by CMake")
     message(STATUS "Try setting variables: ")
     message(STATUS "    NetCDF_ROOT or NetCDF_PATH")
     message(STATUS " or NETCDF_ROOT or NETCDF_PATH")
     message(STATUS " or NetCDF_INCLUDE_DIR or NetCDF_Fortran_INCLUDE_DIR")
     message(STATUS " or NetCDF_INCLUDE_DIRS or NetCDF_Fortran_INCLUDE_DIRS")
     message(STATUS "to appropriate values,")
     message(STATUS "or use an appropriate path in CMAKE_PREFIX_PATH")
     message(STATUS "If unsuccessful, set -DWITH_NETCDF=OFF")
     message(STATUS "   (some Siesta functionalities will not be available)")
     message(FATAL_ERROR "NetCDF could not be found by CMake")
    endif()
  endif()
endif()
#
# Look for LAPACK
#
find_package(CustomLapack REQUIRED)
include(check_lapack)
if (WITH_MPI)
  find_package(CustomScalapack REQUIRED)
endif()

include(check_for_linalg_features)

# Handle external dependencies
# Use proper order for xmlf90 and libpsml dependencies

if (NOT TARGET libfdf::libfdf)
  find_package(libfdf REQUIRED)
endif()

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

if ( WITH_DFTD3 )
 add_subdirectory("External/DFTD3")
endif()

# Search for flook
if (WITH_LUA)
   add_subdirectory("External/Lua-Engine")
endif()

if (WITH_ELPA)
 include(search_for_elpa)
 message(DEBUG "ELPA lib directory (for ELSI linking): ${ELPA_LIBDIR}")
endif()

if (WITH_ELSI AND NOT TARGET elsi::elsi)
  set(ELSI_MIN_VERSION "2.5.0")
  find_package(elsi ${ELSI_MIN_VERSION} REQUIRED)
  message(STATUS "Found ELSI library (CMake package) Version: ${elsi_VERSION}")

  if(NOT TARGET elsi::elpa)
     # ELSI uses an external ELPA library
     # To link it properly, we need to add the right ELPA directory,
     # as the CMake-based ELSI package just links "-lelpa".
     # Alternatively (see https://gitlab.com/elsi_project/elsi_interface/-/issues/57)
     # we could discover ELSI using pkgconfig, as the pkg-config file provides the
     # right directory
     if(NOT ELPA_LIBDIR)
       # We have not searched for ELPA...
       # Since Siesta can use a native interface to ELPA, we just suggest to
       # the user to set -DWITH_ELPA=ON. Alternatively, the user should provide
       # the right value in ELPA_LIBDIR
       message(STATUS "**")
       message(STATUS "The found ELSI library needs an external ELPA library")
       message(STATUS "Either set -DWITH_ELPA=ON")
       message(STATUS "  (recommended, and you get also the native Siesta interface to ELPA)")
       message(STATUS "or set -DELPA_LIBDIR=/path/to/elpa/libdir")
       message(STATUS "  (this might not work with all compilers)")
       message(FATAL_ERROR "**")
     else()
       message(STATUS "The found ELSI library needs an external ELPA library")
       message(STATUS "Adding ELPA directory for ELSI linking: ${ELPA_LIBDIR}")
       target_link_directories(elsi::elsi INTERFACE "${ELPA_LIBDIR}")
       #
       # Note: Some compilers (e.g., maybe the NEC one) migh need access to
       # the full chain of modules. In that case we would also need:
       #     target_include_directories(elsi::elsi INTERFACE ${ELPA_FORTRAN_INC_DIRS})
       # It might be better to just 'include(search_for_elpa)' here
       # and use the variables. 
     endif()
  endif()

  if(TARGET elsi::pexsi)
    message(STATUS "The found ELSI library contains a built-in PEXSI")
    set(ELSI_WITH_PEXSI TRUE)
    enable_language(CXX)
    # Avoid to link everything associated with PEXSI with C++
    set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 0)
  elseif(ELSI_HAS_EXTERNAL_PEXSI)
    message(STATUS "You affirm that the ELSI library links to an external PEXSI")
    find_package(PEXSI REQUIRED)
    set(ELSI_WITH_PEXSI TRUE)
    enable_language(CXX)
    # Avoid to link everything associated with PEXSI with C++
    set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 0)
  else()  
    set(ELSI_WITH_PEXSI FALSE)
  endif()
endif()

if (WITH_FFTW)
  # For STM/ol-stm only for now
  find_package(FFTW QUIET COMPONENTS DOUBLE_LIB)
  if (FFTW_DOUBLE_LIB_FOUND)
     message(STATUS "Found FFTW Double")
  endif()
endif()


if( WITH_MPI )
  # optional mpi-siesta library
  add_subdirectory("Src/MPI")
  # This is the safest to have a simple unified approach
  # I.e. namespaces for modules
  # To be done for all
  add_library(Siesta::MPI ALIAS mpi_siesta)
endif()

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

# Docs
# Enable BUILD_DOCS to create Doxygen-based documentation (wip)
add_subdirectory(Docs)

# add a few very basic tests for now

enable_testing()

include(TestMPIConfig)

if(WITH_MPI)
  add_subdirectory("Util/MPI_test")
endif(WITH_MPI)

add_subdirectory(Tests)


#add_subdirectory("Src/ncps/examples")

