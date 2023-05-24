#
# Toolchain file for
#
# Intel compiler, MKL library
#
# Notes:
#
#  * CMake format: Command line options (e.g. compiler flags) space separated, other kind
#    of lists semicolon separated.
#
#  * Variables containing library search paths are empty by default. The CMAKE_PREFIX_PATH
#    environment variable should be set up correctly, so that CMake can find those libraries
#    automatically. If that is not the case, override those variables to add search paths
#    manually
#


#
# Fortran compiler settings
#
set(Fortran_FLAGS_RELEASE "-O2 -ip -xHost -fp-model=strict -prec-div -prec-sqrt"
  CACHE STRING "Fortran compiler flags for Release build")
  
set(Fortran_FLAGS_MINRELSIZE "-Os -ip -xHost -fp-model=strict -prec-div -prec-sqrt"
    CACHE STRING "Fortran compiler flags for minimum size build")

set(Fortran_FLAGS_RELWITHDEBINFO "-g ${Fortran_FLAGS_RELEASE}"
  CACHE STRING "Fortran compiler flags for Release with debug info build")

set(Fortran_FLAGS_DEBUG "-g -Og -check -traceback"
  CACHE STRING "Fortran compiler flags for Debug build")
  
set(Fortran_FLAGS_CHECK "-g -Og -check -traceback"
  CACHE STRING "Fortran compiler flags for Debug (checking) build")

#
# C compiler settings
#
set(C_FLAGS "${CMAKE_C_FLAGS}"
  CACHE STRING "Build type independent C compiler flags")

set(C_FLAGS_RELEASE "-O2 -ip -xHost"
  CACHE STRING  "C compiler flags for Release build")

set(C_FLAGS_DEBUG "-g -Og -check -traceback"
  CACHE STRING "C compiler flags for Debug build")

