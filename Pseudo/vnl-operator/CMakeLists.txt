set(top_srcdir "${CMAKE_SOURCE_DIR}/Src")

add_library(psop_top_objs OBJECT
    "${top_srcdir}/m_getopts.f90"
    "${top_srcdir}/m_uuid.f90"
)

add_executable(
  psop

  check_grid.f90
  dpnint.f90
  handlers.f
  local_die.F
  m_kb.f90
  m_libxc_list.F90
  m_libxc_sxc_translation.F90
  m_siestaxc_list.F90
  psml_die.f
  psop.f90
  psop_options.f90
  semicore_info_froyen.f90
)

#
 target_link_libraries(
  psop
  
  psop_top_objs
  gridxc
  libpsml
  xmlf90
  ${PROJECT_NAME}-libncps
  ${PROJECT_NAME}-libpsop

 )

#
# This is needed with Ninja.
# In the compilation logs there is a '-I../Pseudo/vnl-operator'
# instead of '-IPseudo/vnl-operator' if I comment this out.
# Hence, it cannot find the modules provided by its own sources.
# The reason is that the 'ncps' library provides a (different) set
# of {libxc,siestaxc} compatibility wrappers.

target_include_directories(
  psop
  PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
  )

install(
  TARGETS psop
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

