include(SiestaFindPackage)

Siesta_find_package(libpsml
  GIT_REPO "https://gitlab.com/siesta-project/libraries/libpsml"
  GIT_TAG "master"
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/External/libpsml
  )
