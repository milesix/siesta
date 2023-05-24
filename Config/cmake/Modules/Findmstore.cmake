include(SiestaFindPackage)

Siesta_find_package(mstore
  URL "https://github.com/grimme-lab/mstore"
  REVISION "HEAD"
  SUBMODULE_DIR ${PROJECT_SOURCE_DIR}/External/DFTD3/mstore
  )
