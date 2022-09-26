#
#  This can be improved.
#  Toolchain files might also be useful for this
#


if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    message("Using an Intel compiler")
    
    if (NOT CMAKE_Fortran_FLAGS_RELEASE)
      set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -ip -xHOST -fp-model=strict -prec-div -prec-sqrt")
    endif()
    if (NOT CMAKE_Fortran_FLAGS_RELWITHDEBINFO)
      set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO   "-g  ${CMAKE_Fortran_FLAGS_RELEASE} -traceback")
    endif()
    if (NOT CMAKE_Fortran_FLAGS_DEBUG)
      set(CMAKE_Fortran_FLAGS_DEBUG  "-g -O0 -check -traceback")
    endif()
    
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    message("Using a GNU compiler")
 
    if ( CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10.0 )
      message(STATUS "  Adding '-fallow-argument-mismatch' for GNU >= 10.0")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
    endif()

    if (NOT CMAKE_Fortran_FLAGS_RELEASE)
      set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -march=native")
                      ## alt:    set(CMAKE_Fortran_FLAGS_RELEASE ""-O2 -funroll-all-loops")
    endif()
    if (NOT CMAKE_Fortran_FLAGS_RELWITHDEBINFO)
      set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO   "-g  ${CMAKE_Fortran_FLAGS_RELEASE} -traceback")
    endif()
    if (NOT CMAKE_Fortran_FLAGS_DEBUG)
      set(CMAKE_Fortran_FLAGS_DEBUG  "-g -O0 -fbounds-check")
    endif()

endif()


#
# Following snippet borrowed from DFTB+ project
#
  if(CMAKE_BUILD_TYPE)
    set(_buildtypes ${CMAKE_BUILD_TYPE})
  else()
    set(_buildtypes ${CMAKE_CONFIGURATION_TYPES})
  endif()
  foreach(_buildtype IN LISTS _buildtypes)
    foreach (lang IN ITEMS Fortran C)
      string(TOUPPER "${_buildtype}" _buildtype_upper)
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
