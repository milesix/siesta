#
#  If you have veclibfort, and it works, you could uncomment this line
#
## set(LAPACK_LIBRARY "-lveclibfort" CACHE STRING "lapack library chosen")
#
#  More general settings for use with shell modules
#
set(LAPACK_LIBRARY "$ENV{LAPACK_LIBS}" CACHE STRING "lapack library chosen")
set(SCALAPACK_LIBRARY "$ENV{SCALAPACK_LIBS}" CACHE STRING "scalapack library chosen")
#
# These can be overridden from the command line
#
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "build_type")
set(Fortran_FLAGS_DEBUG  "-g -O0" CACHE STRING "Fortran debug flags")

#  set(LAPACK_DETECTION OFF CACHE BOOL "LAPACK detection")

set(WITH_ELSI "ON" CACHE BOOL "with ELSI")
#
# This is temporarily needed because the discovery of an external PEXSI in the ELSI config file is
# not complete.
#
set(ELSI_HAS_EXTERNAL_PEXSI "ON" CACHE BOOL "the linked ELSI uses an external PEXSI library")
#
# You will need to specify the paths to the installation of ELSI and PEXSI
# in the CMAKE_PREFIX_PATH variable, either here on in the command line.
# Assuming those are in the environment variables $ELSI_ROOT and $PEXSI_ROOT:
#
set(CMAKE_PREFIX_PATH "$ENV{ELSI_ROOT};$ENV{PEXSI_ROOT}" CACHE STRING "paths to ELSI and PEXSI")

