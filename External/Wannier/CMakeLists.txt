if (NOT TARGET libwannier90::libwannier90)

  # Location of wannier90-3.1.0.tar.gz package
  # It can be a local file or a remote URL
  
  if(NOT "$ENV{WANNIER90_PACKAGE}" STREQUAL "")
       set(url "$ENV{WANNIER90_PACKAGE}")
  elseif(EXISTS
     "${CMAKE_SOURCE_DIR}/External/Wannier/libwannier90/LICENSE")
     #
     # We assume that we have somehow put the pristine source in this directory
     # (perhaps using a submodule -- but do not 'init' it by default,
     # as it is huge -- cloned directory is over 170Mb)
     # Better ship tar.gz file around
     #
     set(url "${CMAKE_SOURCE_DIR}/External/Wannier/libwannier90")

  elseif(EXISTS
     "${CMAKE_SOURCE_DIR}/External/Wannier/patched_libwannier90/CMakeLists.txt")

    # Intended for developers
    message(STATUS "Adding External/Wannier/patched_libwannier90 subdirectory...")
    message(STATUS "... containing already-patched sources for w90 wrapper")

    add_subdirectory("${CMAKE_SOURCE_DIR}/External/Wannier/patched_libwannier90")
    add_library("libwannier90::libwannier90" ALIAS "libwannier90-lib")
    RETURN()
    
  else()
    message(STATUS "Please set environment variable WANNIER90_PACKAGE")
    message(STATUS "to point to a pristine (possibly remote) wannier90-3.1.0.tar.gz package,")
    message(STATUS "or populate ${CMAKE_SOURCE_DIR}/External/Wannier/libwannier90/")
    message(STATUS "with a pristine wannier90 v3.1.0 src/ directory plus the LICENSE file,")
    message(STATUS "or populate ${CMAKE_SOURCE_DIR}/External/Wannier/patched_libwannier90/")
    message(STATUS "with the already patched w90 wrapper sources")

    message(FATAL_ERROR "Cannot find wannier90 sources")
  endif()
  
  set(package "libwannier90")
  message(STATUS "Retrieving ${package} from ${url}")
  include(FetchContent)

  # Note that the hash check has no effect for directories
  # It is possible to put, e.g. ' && rm -rf test-suite' after the patch command
  # string to lighten the unpacked wannier90 distribution, but better not.
  
  FetchContent_Declare(libwannier90
    URL ${url}
    URL_HASH SHA256=40651a9832eb93dec20a8360dd535262c261c34e13c41b6755fa6915c936b254
    PATCH_COMMAND patch -p1 -i ${CMAKE_SOURCE_DIR}/External/Wannier/Patches/3.1.0.patch
    UPDATE_DISCONNECTED true
    UPDATE_COMMAND ""
  )

  # Slightly indirect way to avoid re-patching if CMake is invoked again.
  # The call to _Populate, as opposed to _MakeAvailable, seems not to
  # apply the patching step a second time...
  # (set -DFETCHCONTENT_QUIET=OFF to see the internal messages)
  #
  # The variable libwannier90_POPULATED is not stored in the cache,
  # and thus not saved across invocations of CMake. It is only useful
  # when the subproject might be added at several places in a given
  # CMake run.
  #
  FetchContent_GetProperties(libwannier90)

  if(NOT libwannier90_POPULATED)
    # Fetch the content using previously declared details
    message(STATUS "Making-available ${package}...")
    FetchContent_Populate(libwannier90)

    # Set custom variables, policies, etc.
    # ...

    # Bring the populated content into the build
    message(STATUS "Adding subdirectory ${libwannier90_SOURCE_DIR}")
    add_subdirectory("${libwannier90_SOURCE_DIR}" "${libwannier90_BINARY_DIR}")
  endif()

  ###  Direct way, but will involve re-patching, and thus errors, if
  ##  CMake is re-run.
  ##
  ##  FetchContent_MakeAvailable(libwannier90)

  # The libwannier90::libwannier90 target is not defined in the
  # subproject until installation time. Hence the need for the
  # following lines.
  # Note the GLOBAL qualifier, since this is in a subordinate scope
  # I am not sure about the interpretation of 'libwannier90'
  # in the second line.
  
  add_library(libwannier90::libwannier90 INTERFACE IMPORTED GLOBAL)
  target_link_libraries(libwannier90::libwannier90 INTERFACE libwannier90)

  # We need the module directory in the subproject before we finish the configure stage
  FetchContent_GetProperties(${package} BINARY_DIR LIBWANNIER90_BINARY_DIR)
  if(NOT EXISTS "${LIBWANNIER90_BINARY_DIR}/include")
     make_directory("${LIBWANNIER90_BINARY_DIR}/include")
  endif()

endif()
