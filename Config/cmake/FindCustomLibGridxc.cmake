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
5. If LIBGRIDXC_ALLOW_FETCH is true cmake will try and fetch the remote
   and add it as a sub-project.


The target:
   libgridxc::libgridxc
will be usable upon return.

#]=======================================================================]


set(LIBGRIDXC_ALLOW_FETCH "ON" CACHE BOOL "Allow remote fetching of libgridxc code")

set(_lib "libgridxc")
set(_pkg "LIBGRIDXC")
set(_url "https://gitlab.com/siesta-project/libraries/libgridxc")
set(_tag "master") #TODO revision of version for libgridxc

set(_GOOD_LIBRARY TRUE)

message(CHECK_START "Searching for ${_lib} library")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

macro(found_return)
  unset(_GOOD_LIBRARY)
  list(POP_BACK CMAKE_MESSAGE_INDENT)
  set(LIBGRIDXC_FOUND TRUE)
  message(CHECK_PASS "found with compatibility")
  return()
endmacro()


# -- As a cmake package
find_package(libgridxc CONFIG)
if (libgridxc_FOUND)
  message(STATUS "Found ${_lib} cmake package -- checking compatibility")

  # We must still make sure that it has the right features
  # In the future we might compile different components for the library
  
  if (WITH_LIBXC)
    message(CHECK_START "${_lib}: checking for libxc support")
    list(APPEND CMAKE_MESSAGE_INDENT "  ")
    include(CheckFortranSourceCompiles)
    set(CMAKE_REQUIRED_LIBRARIES libgridxc::libgridxc)
    check_fortran_source_compiles("use gridxc, only: gridxc_setXC_libxc; end"
                               GRIDXC_USES_LIBXC SRC_EXT F90)
    unset(CMAKE_REQUIRED_LIBRARIES)
    list(POP_BACK CMAKE_MESSAGE_INDENT)

    if( GRIDXC_USES_LIBXC )
      message(CHECK_PASS "found")
    else()
      message(CHECK_FAIL "not found")
      set(_GOOD_LIBRARY FALSE)
    endif()
    # If compiled with libxc, the package will have libxc::XC_Fortran as a dependency
    # Hence it is not necessary to add it.
  endif()

  # Conversely, if the library links to libxc, the following checks will fail if WITH_LIBXC is not ON,
  # since we will not have found libxc:XC_Fortran...

  if (WITH_MPI)
    message(CHECK_START "${_lib}: checking for MPI support")
    list(APPEND CMAKE_MESSAGE_INDENT "  ")
    include(CheckFortranSourceCompiles)
    set(CMAKE_REQUIRED_LIBRARIES libgridxc::libgridxc)
    check_fortran_source_compiles("use gridxc, only: gridxc_init; call gridxc_init(1); end"
                               GRIDXC_HAS_MPI SRC_EXT F90)
    unset(CMAKE_REQUIRED_LIBRARIES)
    list(POP_BACK CMAKE_MESSAGE_INDENT)

    if( GRIDXC_HAS_MPI )
      message(CHECK_PASS "found")
    else()
      message(CHECK_FAIL "not found")
      set(_GOOD_LIBRARY FALSE)
    endif()
  endif()

  if (WITH_GRID_SP)
    message(CHECK_START "${_lib}: checking for single precision")
    list(APPEND CMAKE_MESSAGE_INDENT "  ")
    include(CheckFortranSourceRuns)
    set(CMAKE_REQUIRED_LIBRARIES libgridxc::libgridxc)
    check_fortran_source_runs(
        "use gridxc, only: grid_p; if (kind(1.0) /= grid_p) ERROR STOP 1; end"
         GRIDXC_USES_SP SRC_EXT F90)
    unset(CMAKE_REQUIRED_LIBRARIES)
    list(POP_BACK CMAKE_MESSAGE_INDENT)

    if( GRIDXC_USES_SP )
      message(CHECK_PASS "has support")
    else()
      message(CHECK_FAIL "no support")
      set(_GOOD_LIBRARY FALSE)
    endif()
  endif()

  if(_GOOD_LIBRARY)
    # we are in good condition, return
    found_return()
  endif()
  
endif()


# we did not find with cmake package.
# Now we will try with pkg-config + multiconfig

# reset variable check
set(_GOOD_LIBRARY TRUE)

find_package(PkgConfig QUIET)

if (WITH_MPI)
  if (WITH_GRID_SP)
    pkg_check_modules(gridxc_multi IMPORTED_TARGET GLOBAL libgridxc_sp_mpi>=0.10.0)
  else()
    pkg_check_modules(gridxc_multi IMPORTED_TARGET GLOBAL libgridxc_dp_mpi>=0.10.0)
  endif()
else()
  if (WITH_GRID_SP)
    pkg_check_modules(gridxc_multi IMPORTED_TARGET GLOBAL libgridxc_sp>=0.10.0)
  else()
    pkg_check_modules(gridxc_multi IMPORTED_TARGET GLOBAL libgridxc_dp>=0.10.0)
  endif()
endif()  

if(gridxc_multi_FOUND)
  message(STATUS "Found ${_lib} (multi) through pkgconfig -- checking compatibility")
  
  if (WITH_LIBXC)
    message(CHECK_START "${_lib}: checking for libxc support")
    list(APPEND CMAKE_MESSAGE_INDENT "  ")
    include(CheckFortranSourceCompiles)
    set(CMAKE_REQUIRED_LIBRARIES PkgConfig::gridxc_multi)
    check_fortran_source_compiles("use gridxc, only: gridxc_setXC_libxc; end"
                             GRIDXC_USES_LIBXC SRC_EXT F90)
    unset(CMAKE_REQUIRED_LIBRARIES)
    list(POP_BACK CMAKE_MESSAGE_INDENT)

    if (GRIDXC_USES_LIBXC)
      message(CHECK_PASS "found")
      target_link_libraries(PkgConfig::gridxc_multi INTERFACE libxc::XC_Fortran)
    else()      
      message(CHECK_FAIL "not found")
      set(_GOOD_LIBRARY FALSE)
    endif()
  endif()
 
  # quick return, if able
  if(_GOOD_LIBRARY)
    add_library(libgridxc::libgridxc ALIAS PkgConfig::gridxc_multi)
    found_return()
  else()
    unset(_GOOD_LIBRARY)
  endif()
  
endif()


# we did not find the pkg-config multiconfig
# Now we will try pkg-config single-config

set(_GOOD_LIBRARY TRUE)

pkg_check_modules(gridxc IMPORTED_TARGET GLOBAL libgridxc>=0.10.0)
if(gridxc_FOUND)
  message(STATUS "Found ${_lib} (single) through pkgconfig -- checking compatibility")

  if (WITH_MPI)
    message(CHECK_START "${_lib}: checking for MPI support")
    list(APPEND CMAKE_MESSAGE_INDENT "  ")
    include(CheckFortranSourceCompiles)
    set(CMAKE_REQUIRED_LIBRARIES PkgConfig::gridxc)
    check_fortran_source_compiles("use gridxc, only: gridxc_init; call gridxc_init(1); end"
                             GRIDXC_HAS_MPI SRC_EXT F90)
    unset(CMAKE_REQUIRED_LIBRARIES)
    list(POP_BACK CMAKE_MESSAGE_INDENT)

    if( GRIDXC_HAS_MPI )
      message(CHECK_PASS "found")
    else()
      message(CHECK_FAIL "not found")
      set(_GOOD_LIBRARY FALSE)
    endif()
  endif()

  if (WITH_LIBXC)
    message(CHECK_START "${_lib}: checking for libxc support")
    list(APPEND CMAKE_MESSAGE_INDENT "  ")
    include(CheckFortranSourceCompiles)
    set(CMAKE_REQUIRED_LIBRARIES PkgConfig::gridxc)
    check_fortran_source_compiles("use gridxc, only: gridxc_setXC_libxc; end"
                             GRIDXC_USES_LIBXC SRC_EXT F90)
    unset(CMAKE_REQUIRED_LIBRARIES)
    list(POP_BACK CMAKE_MESSAGE_INDENT)

    if (GRIDXC_USES_LIBXC)
      message(CHECK_PASS "found")
      target_link_libraries(PkgConfig::gridxc INTERFACE libxc::XC_Fortran)
    else()      
      message(CHECK_FAIL "not found")
      set(_GOOD_LIBRARY FALSE)
    endif()
  endif()

  if (WITH_GRID_SP)
    message(CHECK_START "${_lib}: checking for single precision")
    list(APPEND CMAKE_MESSAGE_INDENT "  ")
    include(CheckFortranSourceRuns)
    set(CMAKE_REQUIRED_LIBRARIES PkgConfig::gridxc)
    check_fortran_source_runs(
        "use gridxc, only: grid_p; if (kind(1.0) /= grid_p) ERROR STOP 1; end"
         GRIDXC_USES_SP SRC_EXT F90)
    unset(CMAKE_REQUIRED_LIBRARIES)
    list(POP_BACK CMAKE_MESSAGE_INDENT)

    if( GRIDXC_USES_SP )
      message(CHECK_PASS "has support")
    else()
      message(CHECK_FAIL "no support")
      set(_GOOD_LIBRARY FALSE)
    endif()
  endif()
  
  if (_GOOD_LIBRARY)
    add_library(libgridxc::libgridxc ALIAS PkgConfig::gridxc)
    found_return()
  else()
    unset(_GOOD_LIBRARY)
  endif()
endif()


# None of the pkg-config's were found.
# Fallback to sub-projects.

# We will fall-back to shipped libgridxc
set(GRIDXC_SOURCE_DIR "${PROJECT_SOURCE_DIR}/External/libgridxc")
set(GRIDXC_BINARY_DIR "${PROJECT_BINARY_DIR}/External/libgridxc")

if(EXISTS "${GRIDXC_SOURCE_DIR}/CMakeLists.txt")
  message(STATUS "Include ${_lib} from subprojects")
  add_subdirectory(
    "${GRIDXC_SOURCE_DIR}"
    "${GRIDXC_BINARY_DIR}"
  )

  add_library(libgridxc::libgridxc INTERFACE IMPORTED)
  target_link_libraries(libgridxc::libgridxc INTERFACE libgridxc)

  # We need the module directory in the subproject before we finish the configure stage
  if(NOT EXISTS "${GRIDXC_BINARY_DIR}/include")
    make_directory("${GRIDXC_BINARY_DIR}/include")
  endif()

  found_return()
endif()

if (LIBGRIDXC_ALLOW_FETCH)
  # -- Fetch
  message(STATUS "Retrieving ${_lib} from ${_url}")
  include(FetchContent)
  FetchContent_Declare(
      libgridxc
      GIT_REPOSITORY "${_url}"
      GIT_TAG "${_tag}"
  )
  FetchContent_MakeAvailable(libgridxc)

  add_library(libgridxc::libgridxc INTERFACE IMPORTED)
  target_link_libraries(libgridxc::libgridxc INTERFACE libgridxc)

  # We need the module directory in the subproject before we finish the configure stage
  FetchContent_GetProperties(libgridxc BINARY_DIR GRIDXC_BINARY_DIR)
  if(NOT EXISTS "${GRIDXC_BINARY_DIR}/include")
     make_directory("${GRIDXC_BINARY_DIR}/include")
  endif()

  found_return()

else(LIBGRIDXC_ALLOW_FETCH)

  message(WARNING "** Note that fetching is disabled for ${_lib}")
  message(WARNING "** Enable through LIBGRIDXC_ALLOW_FETCH")

endif(LIBGRIDXC_ALLOW_FETCH)

set(LIBGRIDXC_FOUND FALSE)
list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_FAIL "not found")
