set(top_srcdir "${CMAKE_SOURCE_DIR}/Src")

add_library(psop_top_objs OBJECT
    "${top_srcdir}/m_getopts.f90"
    "${top_srcdir}/m_uuid.f90"
    "${top_srcdir}/periodic_table.f"
)

add_executable(
  psop

  check_grid.f90
  dpnint.f90
  handlers.f
  local_die.F
  m_kb.f90
  psml_die.f
  psop.f90
  psop_options.f90
)

# For semicore_info_froyen...
 target_link_libraries(psop_top_objs
  ${PROJECT_NAME}-libncps
 )
 
 target_link_libraries(
  psop
  
  psop_top_objs
  libgridxc::libgridxc
  libpsml::libpsml
  xmlf90::xmlf90
  ${PROJECT_NAME}-libncps
  ${PROJECT_NAME}-libpsop
  ${PROJECT_NAME}-libxc-trans
 )

install(
  TARGETS psop
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

