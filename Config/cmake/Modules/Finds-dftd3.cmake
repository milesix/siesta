include(SiestaFindPackage)

Siesta_find_package(s-dftd3
  URL "https://github.com/dftd3/simple-dftd3"
  REVISION "v0.7.0"
  SUBMODULE_DIR ${PROJECT_SOURCE_DIR}/External/DFTD3/s-dftd3
  )
