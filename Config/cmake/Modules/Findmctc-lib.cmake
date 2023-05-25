include(SiestaFindPackage)

Siesta_find_package(mctc-lib
  GIT_REPOSITORY "https://github.com/grimme-lab/mctc-lib"
  GIT_TAG "v0.3.1"
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/External/DFTD3/mctc-lib
  )
