add_library(
  ${PROJECT_NAME}-libncps

  interpolation.f90
  m_ncps.f90
  m_ncps_froyen_ps_t.f90
  m_ncps_froyen_reader.f
  m_ncps_psml_plugins.f90
  m_ncps_reader.f
  m_ncps_translators.f90
  m_ncps_utils.f90
  m_ncps_writers.f90
  )

### target_include_directories(
###   ${PROJECT_NAME}-libncps
###   PRIVATE
###   ${CMAKE_CURRENT_BINARY_DIR}
### )
target_link_libraries(
  ${PROJECT_NAME}-libncps
  PRIVATE
  ${PROJECT_NAME}-libxc-trans
  libpsml
)
target_include_directories(
  ${PROJECT_NAME}-libncps
  INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}
)


