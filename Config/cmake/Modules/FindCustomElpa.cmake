#
#
# -- Search for ELPA using pkg-config data
# (This, with further refinemens,
#  should go in a module in the cmake subdir)
# 
# There are several non-standard issues:
#
# 1. The pkgconfig filename might be of the form
#
#     elpa-2020.05.001.rc1.pc
#
# instead of simply 'elpa.pc'.
#
# 2. The include directory entry in the pkgconfig file
#    is missing the 'modules' leaf component.
#
# The most expedient solution for (1) is to put a symbolic
# link to the desired version somewhere in the PKG_CONFIG_PATH. For example:
#
#   ln -sf /path/to/original/.pc/file $HOME/lib/pkgconfig/elpa.pc
#
#   export PKG_CONFIG_PATH=$HOME/lib/pkgconfig:${PKG_CONFIG_PATH}
#
# (Note: CMAKE_PREFIX_PATH can also be used, but without the
#  '/lib/pkgconfig' part.
#
# A solution for (2) is encoded below.
# 
find_package(PkgConfig REQUIRED QUIET)

message(CHECK_START "Searching for ELPA library")
if( WITH_OPENMP )
  pkg_check_modules(ELPA elpa_openmp)
  if( NOT ELPA_FOUND )
    message(WARNING "ELPA could not locate elpa_openmp package, resorting to non-threaded version")
    pkg_check_modules(ELPA elpa)
  endif()
else()
  pkg_check_modules(ELPA elpa)
endif()

if(NOT ELPA_FOUND)
  message(CHECK_FAIL "not found")
  return()
endif()

message(DEBUG "ELPA_LIBDIR: ${ELPA_LIBDIR}")
message(DEBUG "ELPA_LIBRARIES: ${ELPA_LIBRARIES}")
message(DEBUG "ELPA_LINK_LIBRARIES: ${ELPA_LINK_LIBRARIES}")
message(DEBUG "ELPA_INCLUDEDIR: ${ELPA_INCLUDEDIR}")
message(DEBUG "ELPA_INCLUDE_DIRS: ${ELPA_INCLUDE_DIRS}")

# Retrieve compilation flags
pkg_get_variable(ELPA_FCFLAGS elpa fcflags)
message(DEBUG "ELPA_FCFLAGS: ${ELPA_FCFLAGS}")

#
# Fix non-standard setting of Fortran module directory
#
if(ELPA_FCFLAGS)
  message(STATUS "Fixing include directories from FCFLAGS")
  string(REPLACE "-I" "" ELPA_Fortran_INCLUDE_DIRS ${ELPA_FCFLAGS}) 
else()
  message(STATUS "Adding include-directories manually")
  # TODO this should probably be a loop since there may be several include dir's
  set(ELPA_Fortran_INCLUDE_DIRS "${ELPA_INCLUDE_DIRS}/modules")
endif()
 
message(STATUS "ELPA Fortran search path to be used: ${ELPA_Fortran_INCLUDE_DIRS}")

message(CHECK_PASS "found")

# Manual adding the library
add_library(Elpa::elpa INTERFACE IMPORTED)
target_link_libraries(Elpa::elpa INTERFACE ${ELPA_LINK_LIBRARIES})
target_include_directories(Elpa::elpa INTERFACE ${ELPA_Fortran_INCLUDE_DIRS})
