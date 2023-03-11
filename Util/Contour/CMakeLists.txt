add_executable( grid1d grid1d.f )
add_executable( grid2d grid2d.f )

install(
  TARGETS grid1d grid2d
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

