set(LAPACK_LIBRARY
      "-L/cineca/prod/opt/libraries/lapack/3.9.0/gnu--8.4.0/lib -llapack -L/cineca/prod/opt/libraries/essl/6.2.1/binary/lib64 -lessl"
        CACHE STRING "lapack library chosen")
	
set(SCALAPACK_LIBRARY
     "-L/cineca/prod/opt/libraries/scalapack/2.1.0/spectrum_mpi--10.3.1--binary/lib -lscalapack"
     CACHE STRING "scalapack library chosen")
     
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "build_type")
set(Fortran_FLAGS_DEBUG  "-g -O0" CACHE STRING "Fortran debug flags")


