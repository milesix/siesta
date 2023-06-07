include(SiestaFindPackage)

Siesta_find_package(libpsml
  GIT_REPOSITORY "https://gitlab.com/siesta-project/libraries/libpsml"
  GIT_TAG "master"
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/External/libpsml
  )

if( "${${_name}_FOUND_METHOD}" STREQUAL "cmake" OR
    "${${_name}_FOUND_METHOD}" STREQUAL "pkgconf")

include(CheckFortranSourceCompiles)

# Figure out whether psml uses the procedure pointer, or not
set(CMAKE_REQUIRED_LIBRARIES libpsml::libpsml)
check_fortran_source_compiles("use m_psml, only: ps_set_error_handler; end"
  LIBPSML_HAS_ERROR_PROCEDURE_POINTER SRC_EXT F90)
unset(CMAKE_REQUIRED_LIBRARIES)

if( LIBPSML_HAS_ERROR_PROCEDURE_POINTER )
  set(LIBPSML_USES_PROCEDURE_POINTER TRUE)
else()
  set(LIBPSML_USES_PROCEDURE_POINTER FALSE)
endif()

else()

set(LIBPSML_USES_PROCEDURE_POINTER TRUE CACHE BOOL "Whether the sources uses the procedure pointer or not")
mark_as_advanced(LIBPSML_USES_PROCEDURE_POINTER)

endif()
