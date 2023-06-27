#
# s-dftd3 depends on a number of libraries
# They can be included as submodules in this directory,
# and this file processes them in order and configures
# those options that can be configured externally.
#
# It would be good to have an option to build only
# the library itself.


set(WITH_TESTS "${WITH_DFTD3_TESTS}" CACHE BOOL "Perform s-dftd3 tests")
# We need to forcefully disable OpenMP, s-dftd3 is defaulting to true for this.
set(WITH_OpenMP "${WITH_OPENMP}")

### Re-use information already found by Siesta itself
### and skip BLAS search in s-dftd3


if (TARGET BLAS::BLAS)
  # Re-use target
  set(BLAS_FOUND "TRUE")
elseif(TARGET LAPACK::LAPACK)
  # Create a BLAS wrapper target
  add_library("BLAS::BLAS" INTERFACE IMPORTED)
  target_link_libraries("BLAS::BLAS" INTERFACE LAPACK::LAPACK)
  set(BLAS_FOUND "TRUE")
endif()

if(NOT TARGET "mctc-lib::mctc-lib")
  find_package("mctc-lib" REQUIRED)
endif()

if(NOT TARGET "mstore::mstore" AND WITH_TESTS)
  find_package("mstore" REQUIRED)
endif()

if(NOT TARGET "test-drive::test-drive")
  find_package("test-drive" REQUIRED)
endif()

if(NOT TARGET "toml-f::toml-f")
  find_package("toml-f" REQUIRED)
endif()

if (NOT TARGET s-dftd3::s-dftd3)
  find_package(s-dftd3 REQUIRED)
#  find_package(Custom-s-dftd3 REQUIRED)
endif()
#
# This is necessary in case a cmake package is found
# for s-dftd3, as imported targets are not global by
# default, and we are running this in a "subdirectory"
# scope. If s-dftd3 is compiled ("source" method), this would not be
# necessary, since we manually add the library as IMPORTED GLOBAL
# in SiestaFindPackage.
#
set_target_properties(s-dftd3::s-dftd3
                      PROPERTIES
                      IMPORTED_GLOBAL TRUE)

