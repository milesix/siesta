include(SiestaFindPackage)

Siesta_find_package(xmlf90
  GIT_REPO "https://gitlab.com/siesta-project/libraries/xmlf90"
  GIT_TAG "master"
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/External/xmlf90
  )
