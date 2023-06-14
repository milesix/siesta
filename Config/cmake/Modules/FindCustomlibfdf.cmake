include(SiestaFindPackage)

Siesta_find_package(libfdf
  REQUIRED
  GIT_REPOSITORY "https://gitlab.com/siesta-project/libraries/libfdf"
  GIT_TAG "master"
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/External/libfdf
  )

