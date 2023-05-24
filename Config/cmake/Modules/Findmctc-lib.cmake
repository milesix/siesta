include(SiestaFindPackage)

Siesta_find_package(mctc-lib
  URL "https://github.com/grimme-lab/mctc-lib"
  REVISION "v0.3.1"
  SUBMODULE_DIR ${PROJECT_SOURCE_DIR}/External/DFTD3/mctc-lib
  )
