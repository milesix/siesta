set(LAPACK_LIBRARY "-lveclibfort" CACHE STRING "lapack library chosen")
set(SCALAPACK_LIBRARY "$ENV{SCALAPACK_LIBS}" CACHE STRING "scalapack library chosen")
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "build_type")
set(Fortran_FLAGS_DEBUG  "-g -O0" CACHE STRING "Fortran debug flags")

#  set(LAPACK_DETECTION OFF CACHE BOOL "LAPACK detection")

