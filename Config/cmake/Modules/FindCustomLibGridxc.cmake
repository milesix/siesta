#
# Copyright (C) 2023 Siesta developers group.
#
#
#  Search for libgridxc library.
#  This is a custom module, as it needs to support several ways of discovery,
#  including two flavors of auto-tools-installations.
#[=======================================================================[.rst:
FindCustomLibGridxc
----------------

This will find and locate the GridXC library.

This is a wrapper around the Findlibgridxc module with additional functionality.

There are a couple of cases it will handle, in the prescribed order.
For all cases it will check for 1) MPI support, 2) grid-precision matching the
requested Siesta precision and 3) libxc support.

1. Use the cmake-package libgridxc, this has all necessary information
2. Search for the libgridxc in multi-config mode using pkg-config
3. Search for the libgridxc in single-config mode using pkg-config
4. If libgridxc directory is found in the External/ folder it may be
   used as a sub-project.
5. It uses the SiestaFindPackage to wrap and thus the options for controlling
  what to search for is applicable here as well.

The target:
   libgridxc::libgridxc
will be usable upon return.

#]=======================================================================]
set(_name libgridxc)

# try and find it using the Siesta_find_package repository
include(SiestaFindPackage)

# Lets just include this always
include(CheckFortranSourceCompiles)
include(CheckFortranSourceRuns)

# For compiling libgridxc we (before procedure pointers) require
# a stub handler for the libraries.
# Here is a listing that has all of these
set(_libgridxc_stub [=[
  subroutine gridxc_timer_start()
  end
  subroutine gridxc_timer_stop()
  end
  subroutine alloc_error_report()
  end
  subroutine alloc_memory_event()
  end
  subroutine die()
  end
  ]=])

# First parse some options as passed to the Siesta compilation
if (WITH_MPI)
  if (WITH_GRID_SP)
    set(pkg_configs "gridxc_multi libgridxc_sp_mpi>=0.10.0")
    #pkg_check_modules(gridxc_multi IMPORTED_TARGET GLOBAL libgridxc_sp_mpi>=0.10.0)
  else()
    set(pkg_configs "gridxc_multi libgridxc_dp_mpi>=0.10.0")
    #pkg_check_modules(gridxc_multi IMPORTED_TARGET GLOBAL libgridxc_dp_mpi>=0.10.0)
  endif()
else()
  if (WITH_GRID_SP)
    set(pkg_configs "gridxc_multi libgridxc_sp>=0.10.0")
    #pkg_check_modules(gridxc_multi IMPORTED_TARGET GLOBAL libgridxc_sp>=0.10.0)
  else()
    set(pkg_configs "gridxc_multi libgridxc_dp>=0.10.0")
    #pkg_check_modules(gridxc_multi IMPORTED_TARGET GLOBAL libgridxc_dp>=0.10.0)
  endif()
endif()

# When dealing with single-installs the name is the same, always.
# Then later down we will check the compatibility.
list(APPEND "gridxc libgridxc>=0.10.0")

# Required flag is from outside
# This should locate some libgridxc and provide the proper target libgridxc::libgridxc
Siesta_find_package(${_name}
  REQUIRED
  PKG_CONFIG ${pkg_configs}
  GIT_REPOSITORY "https://gitlab.com/siesta-project/libraries/libgridxc"
  GIT_TAG "master"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/External/libgridxc"
  )

set(_compat TRUE)

if( "${${_name}_FOUND_METHOD}" STREQUAL "cmake" OR
    "${${_name}_FOUND_METHOD}" STREQUAL "pkgconf")


# Start conversation about compatibility
list(APPEND CMAKE_MESSAGE_INDENT "  ")
message(CHECK_START "Checking for libgridxc compatibility")

# Check libxc support (if requested!)
if (WITH_LIBXC)

  message(CHECK_START "... checking for libxc support")
  list(APPEND CMAKE_MESSAGE_INDENT "  ")

  set(CMAKE_REQUIRED_LIBRARIES libgridxc::libgridxc)
  check_fortran_source_compiles("use gridxc, only: gridxc_setXC_libxc; end;
    ${_libgridxc_stub}"
    GRIDXC_USES_LIBXC SRC_EXT F90)

  # In case the library does not have libxc as a dependency, lets add it
  if (NOT GRIDXC_USES_LIBXC)
    message(VERBOSE " trying by explicitly adding libxc as a dependency")
    target_link_libraries(libgridxc::libgridxc INTERFACE Libxc::xc_Fortran)
    check_fortran_source_compiles("use gridxc, only: gridxc_setXC_libxc; end
      ${_libgridxc_stub}"
      GRIDXC_USES_LIBXC SRC_EXT F90)
  endif()
  unset(CMAKE_REQUIRED_LIBRARIES)

  list(POP_BACK CMAKE_MESSAGE_INDENT)

  if( GRIDXC_USES_LIBXC )
    message(CHECK_PASS "found")
  else()
    message(CHECK_FAIL "not found")
    set(_compat FALSE)
  endif()

endif()

if (WITH_MPI)

  message(CHECK_START "... checking for MPI support")
  list(APPEND CMAKE_MESSAGE_INDENT "  ")

  set(CMAKE_REQUIRED_LIBRARIES libgridxc::libgridxc)
  check_fortran_source_compiles("use gridxc, only: gridxc_init; call gridxc_init(1); end
    ${_libgridxc_stub}"
    GRIDXC_HAS_MPI SRC_EXT F90)
  unset(CMAKE_REQUIRED_LIBRARIES)

  list(POP_BACK CMAKE_MESSAGE_INDENT)

  if( GRIDXC_HAS_MPI )
    message(CHECK_PASS "found")
  else()
    message(CHECK_FAIL "not found")
    set(_compat FALSE)
  endif()
endif()


if (WITH_GRID_SP)
  message(CHECK_START "... checking for single precision")
  list(APPEND CMAKE_MESSAGE_INDENT "  ")

  set(CMAKE_REQUIRED_LIBRARIES libgridxc::libgridxc)
  check_fortran_source_runs("use gridxc, only: grid_p; if (kind(1.0) /= grid_p) ERROR STOP 1; end
    ${_libgridxc_stub}"
    GRIDXC_USES_SP SRC_EXT F90)
  unset(CMAKE_REQUIRED_LIBRARIES)

  list(POP_BACK CMAKE_MESSAGE_INDENT)

  if( GRIDXC_USES_SP )
    message(CHECK_PASS "has support")
  else()
    message(CHECK_FAIL "no support")
    set(_compat FALSE)
  endif()
else()
  message(CHECK_START "... checking for double precision")
  list(APPEND CMAKE_MESSAGE_INDENT "  ")

  set(CMAKE_REQUIRED_LIBRARIES libgridxc::libgridxc)
  check_fortran_source_runs("use gridxc, only: grid_p; if (kind(1.d0) /= grid_p) ERROR STOP 1; end
    ${_libgridxc_stub}"
    GRIDXC_USES_DP SRC_EXT F90)
  unset(CMAKE_REQUIRED_LIBRARIES)

  list(POP_BACK CMAKE_MESSAGE_INDENT)

  if( GRIDXC_USES_DP )
    message(CHECK_PASS "has support")
  else()
    message(CHECK_FAIL "no support")
    set(_compat FALSE)
  endif()

endif()


# Figure out whether psml uses the procedure pointer, or not
set(CMAKE_REQUIRED_LIBRARIES libgridxc::libgridxc)
check_fortran_source_compiles("use gridxc, only: gridxc_set_error_handler; end
  ${_libgridxc_stub}"
  LIBGRIDXC_HAS_ERROR_PROCEDURE_POINTER SRC_EXT F90)
unset(CMAKE_REQUIRED_LIBRARIES)

else()# check for source download

set(LIBGRIDXC_USES_PROCEDURE_POINTER FALSE CACHE BOOL "Whether the sources uses the procedure pointer or not")
mark_as_advanced(LIBGRIDXC_USES_PROCEDURE_POINTER)

endif()

# Final clean-up of the search
if(LIBGRIDXC_FOUND AND _compat)
  if( LIBGRIDXC_HAS_ERROR_PROCEDURE_POINTER )
    set(LIBGRIDXC_USES_PROCEDURE_POINTER TRUE)
  else()
    set(LIBGRIDXC_USES_PROCEDURE_POINTER FALSE)
  endif()
  message(CHECK_PASS "found")
elseif(LIBGRIDXC_FOUND)
  message(CHECK_FAIL "missing compatibility")
  if( CustomLibGridxc_FIND_REQUIRED )
    message(STATUS "Compatibility could not be asserted. Try one of the other search options: libgridxc_FIND_METHOD (not ${${_name}_FOUND_METHOD})")
    message(FATAL_ERROR "Required package libgridxc was not compatible with options")
  endif()
else()
  message(CHECK_FAIL "not found")
  if( CustomLibGridxc_FIND_REQUIRED )
    message(FATAL_ERROR "Required package libgridxc cannot be found")
  endif()
endif()
list(POP_BACK CMAKE_MESSAGE_INDENT)
