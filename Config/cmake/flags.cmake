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

# Pick up any default flags
set(Fortran_FLAGS ${CMAKE_Fortran_FLAGS} CACHE STRING "Build-type independent flags")



# Just a couple of vendors for now. Others should really be put in a toolchain file.

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    message("Using an Intel compiler")
    
      set(Fortran_FLAGS_RELEASE "-O2 -ip -xHOST -fp-model=strict -prec-div -prec-sqrt" CACHE STRING "Fortran release flags")
      set(Fortran_FLAGS_RELWITHDEBINFO   "-g  ${Fortran_FLAGS_RELEASE} -traceback" CACHE STRING "Fortran rel-with-deb-info flags")
      set(Fortran_FLAGS_CHECK  "-g -O0 -check -traceback" CACHE STRING "Fortran debug flags")
      set(Fortran_FLAGS_DEBUG  "-g -O0 -traceback" CACHE STRING "Fortran debug flags")
    
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    message("Using a GNU compiler")
 
    if ( CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10.0 )
      message(STATUS "  Adding '-fallow-argument-mismatch' for GNU >= 10.0")
      set(Fortran_FLAGS "${Fortran_FLAGS} -fallow-argument-mismatch")
    endif()

    set(Fortran_FLAGS_RELEASE "-Ofast -march=native"
                      ## alt:    set(Fortran_FLAGS_RELEASE ""-O2 -funroll-all-loops")
                CACHE STRING "Fortran release flags")

    set(Fortran_FLAGS_RELWITHDEBINFO   "-g  ${Fortran_FLAGS_RELEASE} -fbacktrace"
                CACHE STRING "Fortran rel-with-debinfo flags")

    set(Fortran_FLAGS_CHECK  "-g -O0 -fbounds-check" CACHE STRING "Fortran debug flags")
    set(Fortran_FLAGS_DEBUG  "-g -O0" CACHE STRING "Fortran debug flags")

endif()

#
#  This is common for all compilers, as a convenience fallback for fuller control
#  Set it explicitly from the command line or use in combination with Fortran_FLAGS.
#
set(Fortran_FLAGS_NONE  " " CACHE STRING "Fortran 'none' flags")

#
# Minimal support for C and CXX
#
set(C_FLAGS ${CMAKE_C_FLAGS} CACHE STRING "Build-type independent flags")
set(C_FLAGS_NONE  " " CACHE STRING "C 'none' flags")
set(C_FLAGS_DEBUG  " -g -O0 " CACHE STRING "C 'Debug' flags")
set(C_FLAGS_CHECK  " -g -O0 " CACHE STRING "C 'Debug' flags")
set(C_FLAGS_RELEASE  " -O2 " CACHE STRING "C 'Release' flags")
set(C_FLAGS_RELWITHDEBINFO  "-g -O2 " CACHE STRING "C 'RelWithDebInfo' flags")

set(CXX_FLAGS ${CMAKE_CXX_FLAGS} CACHE STRING "Build-type independent flags")
set(CXX_FLAGS_NONE  " " CACHE STRING "C++ 'none' flags")
set(CXX_FLAGS_DEBUG  " -g -O0 " CACHE STRING "C++ 'Debug' flags")
set(CXX_FLAGS_CHECK  " -g -O0 " CACHE STRING "C++ 'Debug' flags")
set(CXX_FLAGS_RELEASE  " -O2 " CACHE STRING "C++ 'Release' flags")
set(CXX_FLAGS_RELWITHDEBINFO  "-g -O2 " CACHE STRING "C++ 'RelWithDebInfo' flags")

#
  if(CMAKE_BUILD_TYPE)
    set(_buildtypes ${CMAKE_BUILD_TYPE})
  else()
    set(_buildtypes ${CMAKE_CONFIGURATION_TYPES})
  endif()
  foreach(_buildtype IN LISTS _buildtypes)
    foreach (lang IN ITEMS Fortran C)
      string(TOUPPER "${_buildtype}" _buildtype_upper)
      set(CMAKE_${lang}_FLAGS " ${${lang}_FLAGS}")
      set(CMAKE_${lang}_FLAGS_${_buildtype_upper} " ${${lang}_FLAGS_${_buildtype_upper}}")
      message(STATUS "Flags for ${lang}-compiler (build type: ${_buildtype}): "
        "${CMAKE_${lang}_FLAGS} ${CMAKE_${lang}_FLAGS_${_buildtype_upper}}")
    endforeach()
  endforeach()
  unset(_buildtypes)
  unset(_buildtype)
  unset(_buildtype_upper)

#
# For tagging in version-info.inc
#
string(TOUPPER ${CMAKE_BUILD_TYPE} _buildtype_upper)
set(Fortran_FLAGS_CURRENT ${CMAKE_Fortran_FLAGS}
                          ${CMAKE_Fortran_FLAGS_${_buildtype_upper}})
