#
# 
#
set(FLOOK_ROOT "$ENV{FLOOK_ROOT}" CACHE FILEPATH "flook installation path")

#
# Try first a proper installation with a valid pkg-config file,
# augmenting any specified search paths with FLOOK_ROOT

if(NOT "${FLOOK_ROOT}" STREQUAL "" )
  set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" "${FLOOK_ROOT}")
endif()

find_package(PkgConfig REQUIRED)
pkg_check_modules(FLOOK flook>=0.8.1)

if(FLOOK_FOUND)

  message("Flook: Found installation via pkg-config")
  message("FLOOK_LINK_LIBRARIES: ${FLOOK_LINK_LIBRARIES}")
  message("FLOOK_INCLUDE_DIRS: ${FLOOK_INCLUDE_DIRS}")
  add_library(flook::flook INTERFACE IMPORTED GLOBAL)
  target_link_libraries(flook::flook
               INTERFACE "${FLOOK_LINK_LIBRARIES}")
  target_include_directories(flook::flook
               INTERFACE
              "${FLOOK_INCLUDE_DIRS}")

  RETURN()
endif()

#
# Check for an installation without .pc file in FLOOK_ROOT and set libraries
# manually, assuming the recommended recipe has been followed...
#

if(NOT "${FLOOK_ROOT}" STREQUAL "" )

  message("Flook: Using barebones flook installation in ${FLOOK_ROOT}")
  message("Flook: ... assuming basic structure is right")

  set(FLOOK_LIBS "-L${FLOOK_ROOT}/lib -lflookall -ldl")
  set(FLOOK_INCLUDE_DIRS "${FLOOK_ROOT}/include")
  message("FLOOK_LIBS: ${FLOOK_LIBS}")
  message("FLOOK_INCLUDE_DIRS: ${FLOOK_INCLUDE_DIRS}")

  add_library(flook::flook INTERFACE IMPORTED GLOBAL)
  target_link_libraries(flook::flook
               INTERFACE "${FLOOK_LIBS}")
  target_include_directories(flook::flook
               INTERFACE
              "${FLOOK_INCLUDE_DIRS}")


else()

#
# Compile flook
#

include(ExternalProject)

# Compiling the internal Lua package requires tweaking if you are
# using a different platform than `linux`.
#
#    PLATFORM = aix | bsd | c89 | freebsd | generic
#               | linux | macosx | mingw | posix | solaris
#
# where the makefile system will try and guess the correct
# platform. However, if the build fails due to the Lua library, then
# please supply `PLATFORM` with the correct platform in your
# `setup.make` file.
#
# Please add more cases here if needed.

if("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
  set(FLOOK_PLATFORM "macosx")
else()
  set(FLOOK_PLATFORM "linux")
endif()

configure_file(
  setup.make.in
  ${CMAKE_CURRENT_BINARY_DIR}/setup.make
  @ONLY
  )

find_program(MAKE_EXECUTABLE NAMES gmake make mingw32-make REQUIRED)

#
# Possible sources
#
set(url_default "https://github.com/ElectronicStructureLibrary/flook/archive/refs/tags/v0.8.1.tar.gz")
#
if(NOT "$ENV{FLOOK_PACKAGE}" STREQUAL "")
     set(url "$ENV{FLOOK_PACKAGE}")
     message("Flook: will compile on-the-fly from sources: $ENV{FLOOK_PACKAGE}")
elseif(EXISTS
     "${CMAKE_SOURCE_DIR}/External/Lua-Engine/flook/LICENSE")
     #
     # We assume that we have instantiated a submodule
     #
     set(url "${CMAKE_SOURCE_DIR}/External/Lua-Engine/flook")
     message("Flook: will compile on-the-fly from submodule in External/Lua-Engine/flook")
else()

    set(url "${url_default}")
    
    message(STATUS "Flook: We are going to attempt to download the source from:")
    message(STATUS "  ${url_default}")
    message(STATUS "If this is not possible/convenient:")
    message(STATUS "  Please set environment variable FLOOK_PACKAGE")
    message(STATUS "  to point to a pristine (possibly remote) flook-0.8.1.tar.gz package,")
    message(STATUS "  or populate ${CMAKE_SOURCE_DIR}/External/Lua-Engine/flook/")
    message(STATUS "  with the sources for flook. For example, via a submodule:")
    message(STATUS "      git submodule update --init External/Lua-Engine/flook")
endif()

ExternalProject_Add(flook
 URL "${url}"
 UPDATE_DISCONNECTED true
 ##CONFIGURE_HANDLED_BY_BUILD true
 CONFIGURE_COMMAND ""
 BUILD_COMMAND ${MAKE_EXECUTABLE} -j liball
 ## We avoid installing until we know how to specify the prefix...
 ##  ... the following should work but does not.
 ## INSTALL_COMMAND ${MAKE_EXECUTABLE} DESTDIR=<INSTALL_DIR> install
 INSTALL_COMMAND ""
 TEST_COMMAND ""
)

#
# Copy configured setup.make file
# Copy also everything from SOURCE_DIR to BINARY_DIR (this is
# typically done by 'configure' for autotools, but it is not
# done for 'make' projects -- there must be a better way).
# Note special idiom for the 'cp' command. **Check** that it
# works with all shells?
#
ExternalProject_Add_Step(flook
  copy_setup_make
  COMMAND cp ${CMAKE_CURRENT_BINARY_DIR}/setup.make . && cp -a <SOURCE_DIR>/. .
  DEPENDEES download
  DEPENDERS build
  WORKING_DIRECTORY <BINARY_DIR> 
)

add_library(flook::flook INTERFACE IMPORTED GLOBAL)

add_dependencies(flook::flook flook)

#
# We do not install, so the .a and .mod files end up here
#

set(_flook_install_dir "${CMAKE_CURRENT_BINARY_DIR}/flook-prefix/src/flook-build")


SET_TARGET_PROPERTIES(flook::flook PROPERTIES
      INTERFACE_LINK_LIBRARIES
        "${_flook_install_dir}/libflookall${CMAKE_STATIC_LIBRARY_SUFFIX};dl"
      INTERFACE_INCLUDE_DIRECTORIES
        "${_flook_install_dir}"

)

endif(NOT "${FLOOK_ROOT}" STREQUAL "" )
