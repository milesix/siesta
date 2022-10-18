set(top_srcdir "${CMAKE_SOURCE_DIR}/Src")
set(grid_srcdir "${CMAKE_SOURCE_DIR}/Util/Grid")

# To avoid "multiple rules" in Ninja
# Maybe turn into global-utility targets

add_library(stm_gridfunc   ${grid_srcdir}/m_gridfunc.F90)

add_executable( plstm 
  plstm.f90
  ${top_srcdir}/m_getopts.f90
)

add_executable( plsts
  plsts.f
)

target_link_libraries(plstm PRIVATE stm_gridfunc)
target_link_libraries(plsts PRIVATE stm_gridfunc)

install(
  TARGETS plstm plsts
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

