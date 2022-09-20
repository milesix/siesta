set(top_srcdir "${CMAKE_SOURCE_DIR}/Src")
set(grid_srcdir "${CMAKE_SOURCE_DIR}/Util/Grid")


add_executable( plstm 
  plstm.f90
  ${top_srcdir}/m_getopts.f90
  ${grid_srcdir}/m_gridfunc.F90
)

add_executable( plsts
  plsts.f
  ${grid_srcdir}/m_gridfunc.F90
)

install(
  TARGETS plstm plsts
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

