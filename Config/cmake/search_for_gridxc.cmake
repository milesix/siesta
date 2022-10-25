#
#  Search for libgridxc library.
#  This is a custom module, as it needs to support several ways of discovery,
#  including two flavors of auto-tools-installations.

set(GRIDXC_ALLOW_FETCH "ON" CACHE BOOL "Allow remote fetching of libgridxc code")

set(_lib "libgridxc")
set(_pkg "LIBGRIDXC")
set(_url "https://gitlab.com/siesta-project/libraries/libgridxc")
set(_tag "cmake-master")

set(BAD_LIBRARY "0")

#  -- As a cmake package
#
find_package(libgridxc QUIET CONFIG)
if (libgridxc_FOUND)
  message(STATUS "Found libgridxc cmake package")

  # We must still make sure that it has the right features
  # In the future we might compile different components for the library
  
    if (WITH_LIBXC)
      include(CheckFortranSourceCompiles)
      set(CMAKE_REQUIRED_LIBRARIES libgridxc::libgridxc)
      check_fortran_source_compiles("use gridxc, only: gridxc_setXC_libxc; end"
                               GRIDXC_USES_LIBXC SRC_EXT F90)
      unset(CMAKE_REQUIRED_LIBRARIES)
      if (NOT GRIDXC_USES_LIBXC)
         set(BAD_LIBRARY "1")
         message(WARNING "Found libgridxc library does not have required libxc support")
      endif()
      # If compiled with libxc, the package will have libxc::XC_Fortran as a dependency
      # Hence it is not necessary to add it.
    endif()

    # Conversely, if the library links to libxc, the following checks will fail if WITH_LIBXC is not ON,
    # since we will not have found libxc:XC_Fortran...

    if (WITH_MPI)
      include(CheckFortranSourceCompiles)
      set(CMAKE_REQUIRED_LIBRARIES libgridxc::libgridxc)
      check_fortran_source_compiles("use gridxc, only: gridxc_init; call gridxc_init(1); end"
                               GRIDXC_HAS_MPI SRC_EXT F90)
      unset(CMAKE_REQUIRED_LIBRARIES)
      if (NOT GRIDXC_HAS_MPI)
         message(WARNING "Found libgridxc library does not have required MPI support")
	 set(BAD_LIBRARY "1")
      endif()
    endif()

    if (WITH_GRID_SP)
      include(CheckFortranSourceRuns)
      set(CMAKE_REQUIRED_LIBRARIES libgridxc::libgridxc)
      check_fortran_source_runs(
          "use gridxc, only: grid_p; if (kind(1.0) /= grid_p) ERROR STOP 1; end"
           GRIDXC_USES_SP SRC_EXT F90)
      unset(CMAKE_REQUIRED_LIBRARIES)
      if (NOT GRIDXC_USES_SP)
        message(WARNING "Found libgridxc library does not have single-precision cellxc")
        set(BAD_LIBRARY "1")
      endif()
    endif()

    if (NOT BAD_LIBRARY)
      return()
    else()
      unset(BAD_LIBRARY)
    endif()
    
endif()

#  -- With pkg-config

#     --- Multiconfig setting

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
    message(STATUS "Found (multiconfig) libgridxc through pkgconfig")
    
    if (WITH_LIBXC)
      include(CheckFortranSourceCompiles)
      set(CMAKE_REQUIRED_LIBRARIES PkgConfig::gridxc_multi)
      check_fortran_source_compiles("use gridxc, only: gridxc_setXC_libxc; end"
                               GRIDXC_USES_LIBXC SRC_EXT F90)
      unset(CMAKE_REQUIRED_LIBRARIES)
      if (GRIDXC_USES_LIBXC)
        target_link_libraries(PkgConfig::gridxc_multi INTERFACE libxc::XC_Fortran)
      else()      
        message(WARNING "Found libgridxc library does not have required libxc support")
        set(BAD_LIBRARY "1")
      endif()
    endif()
    
    if (NOT BAD_LIBRARY)
      add_library(libgridxc::libgridxc ALIAS PkgConfig::gridxc_multi)
      return()
    else()
      unset(BAD_LIBRARY)
    endif()
    
  endif()
  
#     --- Single name

  pkg_check_modules(gridxc IMPORTED_TARGET GLOBAL libgridxc>=0.10.0)
  if(gridxc_FOUND)
    message(STATUS "Found plain libgridxc through pkgconfig")

    if (WITH_MPI)
      include(CheckFortranSourceCompiles)
      set(CMAKE_REQUIRED_LIBRARIES PkgConfig::gridxc)
      check_fortran_source_compiles("use gridxc, only: gridxc_init; call gridxc_init(1); end"
                               GRIDXC_HAS_MPI SRC_EXT F90)
      unset(CMAKE_REQUIRED_LIBRARIES)
      if (NOT GRIDXC_HAS_MPI)
         message(WARNING "Found libgridxc library does not have required MPI support")
	 set(BAD_LIBRARY "1")
      endif()
    endif()

    if (WITH_LIBXC)
      include(CheckFortranSourceCompiles)
      set(CMAKE_REQUIRED_LIBRARIES PkgConfig::gridxc)
      check_fortran_source_compiles("use gridxc, only: gridxc_setXC_libxc; end"
                               GRIDXC_USES_LIBXC SRC_EXT F90)
      unset(CMAKE_REQUIRED_LIBRARIES)
      if (GRIDXC_USES_LIBXC)
        target_link_libraries(PkgConfig::gridxc INTERFACE libxc::XC_Fortran)
      else()      
        message(WARNING "Found libgridxc library does not have required libxc support")
        set(BAD_LIBRARY "1")
      endif()
    endif()

    if (WITH_GRID_SP)
      include(CheckFortranSourceRuns)
      set(CMAKE_REQUIRED_LIBRARIES PkgConfig::gridxc)
      check_fortran_source_runs(
          "use gridxc, only: grid_p; if (kind(1.0) /= grid_p) ERROR STOP 1; end"
           GRIDXC_USES_SP SRC_EXT F90)
      unset(CMAKE_REQUIRED_LIBRARIES)
      if (NOT GRIDXC_USES_SP)
        message(WARNING "Found libgridxc library does not have single-precision cellxc")
        set(BAD_LIBRARY "1")
      endif()
    endif()
    
    if (NOT BAD_LIBRARY)
      add_library(libgridxc::libgridxc ALIAS PkgConfig::gridxc)
      return()
    else()
      unset(BAD_LIBRARY)
    endif()
  endif()

# -- Submodule

  set(GRIDXC_SOURCE_DIR "${PROJECT_SOURCE_DIR}/External/libgridxc")
  set(GRIDXC_BINARY_DIR "${PROJECT_BINARY_DIR}/External/libgridxc")
  if(EXISTS "${GRIDXC_SOURCE_DIR}/CMakeLists.txt")
    message(STATUS "Include libgridxc from subprojects")
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

    return()
  endif()

if (GRIDXC_ALLOW_FETCH)
# -- Fetch

  message(STATUS "Retrieving libgridxc from ${_url}")
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

  return()

else(GRIDXC_ALLOW_FETCH)

  message(STATUS "** Note that fetching is disabled for libgridxc")
  message(STATUS "** Enable through GRIDXC_ALLOW_FETCH")

endif(GRIDXC_ALLOW_FETCH)

message(FATAL_ERROR "Cannot find libgridxc !!")
