# Print details related to the build information
# This should comprise *ALL* components.
# And preferentially in some order to retain information locally.

# Global function for printing stuff
function(print_feature_info)
  set(options REQUIRED)
  set(oneValueArgs OPTION FOUND)
  set(multiValueArgs
    VARIABLES DEPENDENCIES OPTIONAL_DEPENDENCIES
    MSG MSGOFF MSGON
    HEADER FOOTER)
  cmake_parse_arguments(_pi "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # add an empty line
  message(NOTICE "")
  message(NOTICE "-------------------------------------------------------")
  list(APPEND CMAKE_MESSAGE_INDENT "|")

  if(DEFINED _pi_HEADER)
    foreach(line IN LISTS _pi_HEADER)
      message(NOTICE "${line}")
    endforeach()
  endif()
  list(APPEND CMAKE_MESSAGE_INDENT "  ")
  
  if(DEFINED _pi_OPTION)
    if( ${${_pi_OPTION}} )
      message(NOTICE "Feature is turned ON and controlled by '${_pi_OPTION}=${${_pi_OPTION}}'")
    else()
      message(NOTICE "Feature is turned OFF and controlled by '${_pi_OPTION}=${${_pi_OPTION}}'")
    endif()
  elseif(DEFINED _pi_FOUND)
    if( ${${_pi_FOUND}} )
      message(NOTICE "Feature is turned ON")
    else()
      message(NOTICE "Feature is turned OFF")
    endif()
    set("${_pi_OPTION}" "${${_pi_FOUND}}")
  endif()

  if(_pi_REQUIRED AND (NOT ${${_pi_OPTION}}) )
    message(FATAL_ERROR "Logic in library information could not be fulfilled"
      "Neither OPTION or FOUND variable is defined")
  endif()

  # Print out information if on
  if( ${${_pi_OPTION}} )
    if(DEFINED _pi_MSGON)
      foreach(line IN LISTS _pi_MSGON)
        message(NOTICE "${line}")
      endforeach()
    endif()
  elseif(DEFINED _pi_MSGOFF)
    foreach(line IN LISTS _pi_MSGOFF)
      message(NOTICE "${line}")
    endforeach()
  endif()

  if(DEFINED _pi_MSG)
    foreach(line IN LISTS _pi_MSG)
      message(NOTICE "${line}")
    endforeach()
  endif()

  if(DEFINED _pi_DEPENDENCIES)
    message(NOTICE "The following dependencies are required to use this feature:")
    foreach(dep IN LISTS _pi_DEPENDENCIES)
      message(NOTICE " - ${dep}")
    endforeach()
  endif()
  if(DEFINED _pi_OPTIONAL_DEPENDENCIES)
    message(NOTICE "The following dependencies are optional when using this feature:")
    foreach(dep IN LISTS _pi_OPTIONAL_DEPENDENCIES)
      message(NOTICE " - ${dep}")
    endforeach()
  endif()

  if(DEFINED _pi_VARIABLES)

    if( ${${_pi_OPTION}} )

      # We have enabled the feature, so the user has supplied enough information
      message(NOTICE "Variables used to enable feature:")
      list(APPEND CMAKE_MESSAGE_INDENT "  - ")
      foreach(var IN LISTS _pi_VARIABLES)
        if(NOT "${${var}}" STREQUAL "")
          message(NOTICE "${var}=${${var}}")
        endif()
      endforeach()
      list(POP_BACK CMAKE_MESSAGE_INDENT)
      
      message(NOTICE "Empty or undefined variables (possibly not needed to be set!):")
      list(APPEND CMAKE_MESSAGE_INDENT "  - ")
      foreach(var IN LISTS _pi_VARIABLES)
        if("${${var}}" STREQUAL "")
          message(NOTICE "${var}")
        endif()
      endforeach()
      list(POP_BACK CMAKE_MESSAGE_INDENT)

    else( ${${_pi_OPTION}} )
      
      # We have not enabled the feature, so the user may have some missing information
      message(NOTICE "Variables used but the feature is NOT enabled still:")
      list(APPEND CMAKE_MESSAGE_INDENT "  - ")
      foreach(var IN LISTS _pi_VARIABLES)
        if(NOT "${${var}}" STREQUAL "")
          message(NOTICE "${var}=${${var}}")
        endif()
      endforeach()
      list(POP_BACK CMAKE_MESSAGE_INDENT)
      
      message(NOTICE "Empty or undefined variables (possibly some are needed to enable the feature):")
      list(APPEND CMAKE_MESSAGE_INDENT "  - ")
      foreach(var IN LISTS _pi_VARIABLES)
        if("${${var}}" STREQUAL "")
          message(NOTICE "${var}")
        endif()
      endforeach()
      list(POP_BACK CMAKE_MESSAGE_INDENT)

    endif( ${${_pi_OPTION}} )

  endif()
  
  if(DEFINED _pi_FOOTER)
    foreach(line IN LISTS _pi_FOOTER)
      message(NOTICE "${line}")
    endforeach()
  endif()

  list(POP_BACK CMAKE_MESSAGE_INDENT)
  list(POP_BACK CMAKE_MESSAGE_INDENT)
  message(NOTICE "-------------------------------------------------------")
endfunction()

message(NOTICE "")
message(NOTICE
  "Printing out information related to the build about to proceed.\n"
  "Please carefully go through the following lines to assert that the "
  "options are as expected, some default fall-backs may disable features "
  "depending on how variables are passed."
)

print_feature_info(
  HEADER "BLAS support"
  VARIABLES
    BLAS_LIBRARY
    BLAS_LIBRARY_DIR
    BLAS_LIBRARIES
    BLAS_LINKER_FLAG
    BLAS_DETECTION
  FOUND BLAS_FOUND # will generally be TRUE since otherwise the build will crash
  MSG "Required library for fast performance"
  "Recommended libraries are:"
  "  - mkl"
  "  - openblas"
  "  - blis"
  "The NetLib BLAS library is a reference implementation which should be avoided"
  "for performance reasons"
  )

print_feature_info(REQUIRED
  HEADER "LAPACK support"
  VARIABLES
    LAPACK_LIBRARY
    LAPACK_LIBRARY_DIR
    LAPACK_LIBRARIES
    LAPACK_LINKER_FLAG
    LAPACK_DETECTION
  FOUND LAPACK_FOUND # will generally be TRUE since otherwise the build will crash
  MSG "Required library for fast performance"
  "Recommended libraries are:"
  "  - mkl"
  "  - openblas"
  "  - flame"
  "The NetLib LAPACK library is fine to use as long as the linked BLAS library is"
  "NOT the NetLib BLAS library!"
  )

print_feature_info(
  HEADER "MPI (parallel support)"
  OPTION WITH_MPI
  FOUND MPI_Fortran_FOUND
  MSGOFF "Parallel support is highly advised to allow scalable and faster calculations"
  DEPENDENCIES "ScaLAPACK"
)

print_feature_info(
  HEADER "ScaLAPACK library"
  OPTION WITH_MPI
  FOUND SCALAPACK_FOUND
  DEPENDENCIES "MPI"
  VARIABLES
    SCALAPACK_LIBRARY
    SCALAPACK_LIBRARY_DIR
    SCALAPACK_LIBRARIES
    SCALAPACK_LINKER_FLAG
    SCALAPACK_DETECTION
  MSGOFF "Parallel support is highly advised to allow scalable and faster calculations"
  )

print_feature_info(
  HEADER "ELPA library support (faster diagonalizations)"
  OPTION WITH_ELPA
  FOUND ELPA_FOUND
  MSGOFF
    "ELPA provides a significant speedup for diagonalizations, users are "
    "generally advised to add support for this library, if able."
  MSGON
    "ELPA support is controlled in the fdf via:"
    "  Diag.Algorithm ELPA-2stage|ELPA-1stage # former is preferred"
  VARIABLES
    ELPA_LIBDIR
    ELPA_LIBRARIES
    ELPA_LINK_LIBRARIES
    ELPA_INCLUDEDIR
    ELPA_INCLUDE_DIRS
    ELPA_FCFLAGS
  DEPENDENCIES "MPI" "ScaLAPACK"
  )


print_feature_info(
  HEADER "OpenMP (threaded) support"
  MSGOFF
    "Can be used for extremely large systems to reduce memory requirements"
    "at the cost of some overhead since only some parts of the code is parallelized with OpenMP"
  MSGON
    "Carefully analyze your typicals runs for whether this makes sense (performance wise)"
    "It might be that the overhead of using OpenMP is very high in which case it should not be used."
    "Users are encouraged to have both a non-OpenMP AND an OpenMP executable to easily switch on a case-by-case"
  OPTION WITH_OPENMP
  )


  
print_feature_info(
  HEADER "Libxc for reference exchange-correlation functionals"
  MSGOFF
    "Users are advised to use the libxc library for full feature completeness"
    "using the libxc interaction with libgridxc."
    "The XC functionals will be restricted to those directly implemented in Siesta."
  VARIABLES
    LIBXC_Fortran_INTERFACE
    LIBXC_C_LIBDIR
    LIBXC_C_LIBRARIES
    LIBXC_C_LINK_LIBRARIES
    LIBXC_C_INCLUDEDIR
    LIBXC_C_INCLUDE_DIRS
    LIBXC_F03_LIBDIR
    LIBXC_F03_LIBRARIES
    LIBXC_F03_LINK_LIBRARIES
    LIBXC_F03_INCLUDEDIR
    LIBXC_F03_INCLUDE_DIRS
    LIBXC_F90_LIBDIR
    LIBXC_F90_LIBRARIES
    LIBXC_F90_LINK_LIBRARIES
    LIBXC_F90_INCLUDEDIR
    LIBXC_F90_INCLUDE_DIRS
  OPTION WITH_LIBXC
  FOUND LIBXC_Fortran_FOUND
  )

print_feature_info(REQUIRED
  HEADER "LibGridXC support to calculate XC functionals on the grid"
  VARIABLES
    GRIDXC_ALLOW_FETCH
  OPTIONAL_DEPENDENCIES LibXC MPI
  )


print_feature_info(
  HEADER "NetCDF output file support"
  MSGOFF
    "Users are adviced to add support for NetCDF due to many functionalities"
    "relying on NetCDF file outputs."
    "It also allows parallel file outputs which can easily be compressed subsequently"
  VARIABLES
    NetCDF_ROOT
    NetCDF_PATH
    NetCDF_INCLUDE_DIR
    NetCDF_INCLUDE_DIRS
    NetCDF_Fortran_INCLUDE_DIR
    NetCDF_Fortran_INCLUDE_DIRS

  OPTION WITH_NETCDF
  FOUND NetCDF_FOUND 
  )

print_feature_info(
  HEADER "Lua support through flook library"
  MSGON "Interaction with Lua can be enabled by adding"
  "  Lua.Script luafile.lua"
  "For molecular dynamics controlled in Lua, additionally define:"
  "  MD.TypeOfRun Lua"
  VARIABLES
    FLOOK_ROOT
    FLOOK_LINK_LIBRARIES
    FLOOK_INCLUDE_DIRS
  OPTION WITH_FLOOK
  FOUND FLOOK_FOUND
  )

print_feature_info(
  HEADER "DFT-D3 corrections"
  OPTION WITH_DFTD3
  MSGOFF
    "Adding support for DFT-D3 can improve force descriptions"
    "using simple but heuristically determined corrections to forces"
  )




# Empty line
message(NOTICE "")
message(NOTICE "Done with build information")
message(NOTICE "")

