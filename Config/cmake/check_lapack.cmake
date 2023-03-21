# Run sanity checks on the libraries.
# In particular, check that the return convention is appropriate on MacOS
# (See https://github.com/mcg1969/vecLibFort for more details)
#

include(CheckFortranSourceCompiles)
include(CheckFortranSourceRuns)

set(CMAKE_REQUIRED_QUIET OFF)


message(STATUS "Checking that BLAS library works...")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

# Check that sgemm can be found (blas check)
if( TARGET BLAS::BLAS )
  set(CMAKE_REQUIRED_LIBRARIES BLAS::BLAS)
else()
  set(CMAKE_REQUIRED_LIBRARIES LAPACK::LAPACK)
endif()

check_fortran_source_compiles(
"
external :: sgemm
call sgemm()
end
"
  blas_has_sgemm SRC_EXT F90)

if (NOT blas_has_sgemm)
  message(WARNING
    "---------------------------------------------"
    " BLAS library cannot link properly"
    " Please check the library linking string found or used by CMake"
    "---------------------------------------------")
  message(FATAL_ERROR " BLAS library does not link properly")
endif()


check_fortran_source_runs(
"
complex :: c
complex, dimension(2) :: a = [ 1, 2 ]
complex :: cdotu
external :: cdotu
c=cdotu(2,a(:),1,a(:),1)
end
"
blas_cdotu_return_convention SRC_EXT F90)
unset(CMAKE_REQUIRED_LIBRARIES)

if (NOT blas_cdotu_return_convention)
  message(WARNING
    "---------------------------------------------"
    " BLAS library uses wrong return-value convention!"
    " This is likely to happen on MacOS if the default Accelerate framework is used"
    " You can install veclibfort (https://github.com/mcg1969/vecLibFort)"
    " and then set either of:\n"
    "   -DBLAS_LIBRARY=-lveclibfort"
    " or"
    "   -DLAPACK_LIBRARY=-lveclibfort\n"
    " in your cmake invocation."
    " Where the key point is to have the BLAS library linked first."
    " Alternatively you can install OpenBLAS/BLIS and set the variables accordingly."
    "---------------------------------------------")
  message(FATAL_ERROR " BLAS library uses wrong return-value convention!!!")
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)


message(STATUS "Checking that LAPACK library works...")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

# Test for lapack
set(CMAKE_REQUIRED_LIBRARIES LAPACK::LAPACK)
check_fortran_source_compiles(
"
external :: dsysv
call dsysv()
end
"
lapack_has_dsysv SRC_EXT F90)

if (NOT lapack_has_dsysv)
  message(WARNING
    "---------------------------------------------"
    " LAPACK library cannot link properly"
    " Please check the library linking string found or used by CMake"
    "---------------------------------------------")
 message(FATAL_ERROR "  *** LAPACK library does not link properly")
endif()


unset(CMAKE_REQUIRED_LIBRARIES)

list(POP_BACK CMAKE_MESSAGE_INDENT)
