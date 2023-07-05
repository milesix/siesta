set(top_srcdir "${CMAKE_SOURCE_DIR}/Src")

add_library(psop_top_objs OBJECT
    "${top_srcdir}/alloc.F90"
    "${top_srcdir}/m_getopts.f90"
    "${top_srcdir}/m_uuid.f90"
    "${top_srcdir}/m_walltime.f90"
    "${top_srcdir}/m_io.f"
    "${top_srcdir}/parallel.F"
    "${top_srcdir}/precision.F"
    "${top_srcdir}/periodic_table.f"
    "${top_srcdir}/timer_tree.f90"
    "${top_srcdir}/m_timer.F90"
    "${top_srcdir}/timer.F90"
    "${top_srcdir}/moreParallelSubs.F90"
)

add_executable(
  psop

  check_grid.f90
  dpnint.f90
  m_kb.f90
  psop.f90
  psop_options.f90
)

# For semicore_info_froyen...
target_link_libraries(psop_top_objs
  PRIVATE
  ${PROJECT_NAME}-libncps
  ${PROJECT_NAME}-libsys
 )
 
target_link_libraries(
  psop
  PRIVATE
  ${PROJECT_NAME}-libsys
  ${PROJECT_NAME}-libncps
  ${PROJECT_NAME}-libpsop
  ${PROJECT_NAME}-libxc-trans
  psop_top_objs
  libpsml::libpsml
  libgridxc::libgridxc
 )

install(
  TARGETS psop
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

