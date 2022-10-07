add_executable( orbmol_proj orbmol_proj.f90)

install(
  TARGETS orbmol_proj
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

