include(SiestaFindPackage)

Siesta_find_package(mstore
  GIT_REPO "https://github.com/grimme-lab/mstore"
  GIT_TAG "HEAD"
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/External/DFTD3/mstore
  )
