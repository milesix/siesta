# Run sanity checks on the libraries.
# In particular, check that the return convention is appropriate on MacOS
# (See https://github.com/mcg1969/vecLibFort for more details)
#

include(CheckFortranSourceCompiles)
include(CheckFortranSourceRuns)

set(CMAKE_REQUIRED_QUIET OFF)
message(STATUS
        "Checking that Lapack library links...")
	
set(CMAKE_REQUIRED_LIBRARIES LAPACK::LAPACK)
check_fortran_source_compiles(
"
external :: sgemm, dsysv
call sgemm()
call dsysv()
end
"
LAPACK_LINKS_OK SRC_EXT F90)

unset(CMAKE_REQUIRED_LIBRARIES)

if (NOT LAPACK_LINKS_OK)
 message(STATUS "  ---------------------------------------------")
 message(STATUS "  *** LAPACK library does not link properly")
 message(STATUS "  Please check the library linking string found or used by CMake")
 message(STATUS "  ---------------------------------------------")
 message(FATAL_ERROR "  *** LAPACK library does not link properly")
endif()

set(CMAKE_REQUIRED_LIBRARIES LAPACK::LAPACK)
check_fortran_source_runs(
"
complex :: c
complex, dimension(2) :: a = [ 1, 2 ]
complex :: cdotu
external :: cdotu
c=cdotu(2,a(:),1,a(:),1)
end
"
LAPACK_RETURN_CONVENTION_OK SRC_EXT F90)
unset(CMAKE_REQUIRED_LIBRARIES)

if (NOT LAPACK_RETURN_CONVENTION_OK)
 message(STATUS "  ---------------------------------------------")
 message(STATUS "  *** LAPACK library uses wrong return-value convention!!!")
 message(STATUS "  This is likely to happen on MacOS if the default Accelerate framework is used")
 message(STATUS "  You can install veclibfort (https://github.com/mcg1969/vecLibFort)")
 message(STATUS "  and then set -DLAPACK_LIBRARY=-lveclibfort in your cmake invocation")
 message(STATUS "  You can also set -DBLAS_LIBRARY=-lveclibfort, but it is not strictly needed")
 message(STATUS "    (The key is that -lveclibfort appears first in the link stage,")
 message(STATUS "     so the BLAS library can still point to the raw Accelerate")
 message(STATUS "  ")
 message(STATUS "  Alternatively you can install OpenBLAS and set the variables accordingly")
 message(STATUS "  ---------------------------------------------")
 message(FATAL_ERROR "  *** LAPACK library uses wrong return-value convention!!!")
endif()
