include(SiestaFindPackage)

Siesta_find_package(libfdf
  REQUIRED
  MIN_VERSION 0.5.0
  GIT_REPOSITORY "https://gitlab.com/siesta-project/libraries/libfdf"
  GIT_TAG "0.5.0"
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/External/libfdf
  )

