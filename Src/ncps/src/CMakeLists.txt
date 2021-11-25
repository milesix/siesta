add_library(
  ${PROJECT_NAME}-libncps

  interpolation.f90
  libxc-compat/m_libxc_compat.f90
  libxc-compat/m_libxc_list.F90
  libxc-compat/m_libxc_sxc_translation.F90
  libxc-compat/m_siestaxc_list.F90
  m_ncps.f90
  m_ncps_froyen_ps_t.f90
  m_ncps_froyen_reader.f
  m_ncps_psml_plugins.f90
  m_ncps_reader.f
  m_ncps_translators.f90
  m_ncps_writers.f90

  )

target_include_directories(
  ${PROJECT_NAME}-libncps
  INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}
)
target_link_libraries(
  ${PROJECT_NAME}-libncps
  libpsml
)


