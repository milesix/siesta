include(SiestaFindPackage)

Siesta_find_package(test-drive
  URL "https://github.com/fortran-lang/test-drive"
  REVISION "v0.4.0"
  SUBMODULE_DIR ${PROJECT_SOURCE_DIR}/External/DFTD3/test-drive
  )
