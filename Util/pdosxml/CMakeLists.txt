add_executable( pdosxml
     m_orbital_chooser.f90 m_pdos.f90 pdosxml.f90)

target_link_libraries(pdosxml xmlf90::xmlf90)

install(
  TARGETS pdosxml
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

