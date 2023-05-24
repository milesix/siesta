include(SiestaFindPackage)

Siesta_find_package(libfdf
  URL "https://gitlab.com/siesta-project/libraries/libfdf"
  REVISION "master"
  SUBMODULE_DIR ${PROJECT_SOURCE_DIR}/External/libfdf
  )

