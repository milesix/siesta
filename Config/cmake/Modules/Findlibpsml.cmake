include(SiestaFindPackage)

Siesta_find_package(libpsml
  URL "https://gitlab.com/siesta-project/libraries/libpsml"
  REVISION "master"
  SUBMODULE_DIR ${PROJECT_SOURCE_DIR}/External/libpsml
  )
