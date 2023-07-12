#
# Handling of compiler flags
#
# Following ideas and code borrowed from DFTB+ project
#
# The options hard-coded in this file can be overridden by the use of
# command-line variables:
#
#    cmake ......  -DFortran_FLAGS_DEBUG=" -g -O0" -DCMAKE_BUILD_TYPE=Debug
#
#    cmake ......  -DFortran_FLAGS=" -g -O0"  -DCMAKE_BUILD_TYPE=None
#
#  The default flags in Fortran_FLAGS are those passed to the compiler as FFLAGS from
#  the environment.
#
#  Toolchain files, or .cmake files processed with -C  might also be useful for this.
#
include(CheckFortranCompilerFlag)


message(STATUS "Using ${CMAKE_Fortran_COMPILER_ID} compiler")
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel|IntelLLVM)
  set(_toolchain "intel")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  set(_toolchain "gnu")

else()
  set(_toolchain "generic")
endif()

set(SIESTA_TOOLCHAIN "${_toolchain}"
  CACHE STRING
  "Define which default compiler settings to use: none|generic|gnu|intel|mac|marconi100|marenostrum|vega"
  )
set_property(CACHE SIESTA_TOOLCHAIN PROPERTY STRINGS
  "none" "generic" "gnu" "intel" "mac" "marconi100" "marenostrum" "vega")

foreach(_toolchain IN LISTS SIESTA_TOOLCHAIN)
  # Now load the toolchain
  if(EXISTS "${PROJECT_SOURCE_DIR}/Config/cmake/toolchains/${_toolchain}.cmake")
    message(STATUS "Using toolchain: Config/cmake/toolchains/${_toolchain}.cmake")
    include("${PROJECT_SOURCE_DIR}/Config/cmake/toolchains/${_toolchain}.cmake")
  elseif(EXISTS "${_toolchain}.cmake")
    message(STATUS "Using toolchain: ${_toolchain}.cmake")
    include("${_toolchain}.cmake")
  elseif(EXISTS "${_toolchain}")
    message(STATUS "Using toolchain: ${_toolchain}")
    include("${_toolchain}")
  else()
    message(FATAL_ERROR "Unknown toolchain searched: Config/cmake/toolchains/${_toolchain}.cmake, ${_toolchain}.cmake and ${_toolchain}.")
  endif()
endforeach()


# Get all languages, enables easy extensions
get_property(_languages GLOBAL PROPERTY ENABLED_LANGUAGES)

# Default build-type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "
No build type selected. ${PROJECT_NAME} will default to 'Release'.
To override pass -DCMAKE_BUILD_TYPE=<type> in order to configure ${PROJECT_NAME}.
Available options are:
  * -DCMAKE_BUILD_TYPE=Release - For an optimized build with no assertions or debug info.
  * -DCMAKE_BUILD_TYPE=Debug - For an unoptimized build with assertions and debug info.
  * -DCMAKE_BUILD_TYPE=Check - For an unoptimized build with assertions and debug info + code checks.
  * -DCMAKE_BUILD_TYPE=RelWithDebInfo - For an optimized build with no assertions but with debug info.
  * -DCMAKE_BUILD_TYPE=MinSizeRel - For a build optimized for size instead of speed.
")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build" FORCE)
endif()

# For cmake-gui
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
  "Debug" "Check" "Release" "MinRelSize" "RelWithDebInfo")


if(CMAKE_BUILD_TYPE)
  set(_buildtypes ${CMAKE_BUILD_TYPE})
else()
  set(_buildtypes ${CMAKE_CONFIGURATION_TYPES})
endif()

foreach(_buildtype IN LISTS _buildtypes)
  foreach(lang IN LISTS _languages)
    string(TOUPPER "${_buildtype}" _buildtype_upper)
    set(CMAKE_${lang}_FLAGS "${${lang}_FLAGS}")
    set(CMAKE_${lang}_FLAGS_${_buildtype_upper} "${${lang}_FLAGS_${_buildtype_upper}}")
    message(STATUS "Flags for ${lang}-compiler (build type: ${_buildtype}): "
      "${CMAKE_${lang}_FLAGS} ${CMAKE_${lang}_FLAGS_${_buildtype_upper}}")
  endforeach()
endforeach()
unset(_buildtypes)
unset(_buildtype)
unset(_buildtype_upper)
unset(_languages)


if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10.0 )
    message(STATUS "Adding '-fallow-argument-mismatch' for GNU >= 10.0")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif()
endif()

# For tagging in version-info.inc
#
string(TOUPPER ${CMAKE_BUILD_TYPE} _buildtype_upper)
set(Fortran_FLAGS_CURRENT ${CMAKE_Fortran_FLAGS}
                          ${CMAKE_Fortran_FLAGS_${_buildtype_upper}})
