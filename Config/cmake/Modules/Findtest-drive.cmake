include(SiestaFindPackage)

Siesta_find_package(test-drive
  REQUIRED
  GIT_REPOSITORY "https://github.com/fortran-lang/test-drive"
  GIT_TAG "v0.4.0"
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/External/DFTD3/test-drive
  )
