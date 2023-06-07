#
#
# -- Search for ELSI using pkg-config data
#

option(SEARCH_ELSI_PKGCONF "Use ELSI pkg-conf file" TRUE)

if(SEARCH_ELSI_PKGCONF)

 # This seems to work better with external libraries (including CUDA), since somehow all the pieces
 # for the link are specified in the elsi.pc file.

 find_package(PkgConfig REQUIRED)
 pkg_check_modules(ELSI REQUIRED elsi>=2.9.1)
 message(STATUS "Found ELSI library (pkg-conf file) Version: ${elsi_VERSION}")

 message(DEBUG "ELSI Libdir: ${ELSI_LIBDIR}")
 message(DEBUG "ELSI Libraries: ${ELSI_LIBRARIES}")
 message(STATUS "ELSI Link libraries: ${ELSI_LINK_LIBRARIES}")
 message(DEBUG "ELSI_INCLUDEDIR: ${ELSI_INCLUDEDIR}")
 message(STATUS "ELSI_INCLUDE_DIRS: ${ELSI_INCLUDE_DIRS}")

 add_library(elsi::elsi INTERFACE IMPORTED)
 target_link_libraries(elsi::elsi INTERFACE  ${ELSI_LINK_LIBRARIES})
 target_include_directories(elsi::elsi INTERFACE ${ELSI_INCLUDE_DIRS})

#
# We have not yet implemented PEXSI discovery in ELSI through this route.
# If it has PEXSI support, please use -DELSI_WITH_PEXSI=ON in the command line...
# This still does not work...  See https://cmake.cmake.narkive.com/THgbc5It/fortran-program-linked-with-c-library
#   ... maybe I need to use a STATIC qualifier?
if(ELSI_WITH_PEXSI)
  message(STATUS "If you have PEXSI in ELSI you need to use the CMake package")
  message(SEND_ERROR "by using the option -DSEARCH_ELSI_PKGCONF=OFF")
  enable_language(CXX)
  # Avoid to link everything associated with PEXSI with C++
  #set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 0)
  set_target_properties(elsi::elsi PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
  )
endif()

else(SEARCH_ELSI_PKGCONF)

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
  else()
    set(ELSI_WITH_PEXSI FALSE)
  endif()

endif(SEARCH_ELSI_PKGCONF)

