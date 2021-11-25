add_library(
  ${PROJECT_NAME}-libpsop

  arw.f
  m_kbgen.F
  m_localgen.f
  m_psop.f90
  psop_params.f
  schrodinger.F

)

target_include_directories(
  ${PROJECT_NAME}-libpsop
  INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}
)


