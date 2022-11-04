add_executable( optical optical.f )
add_executable( optical_input optical_input.f )

install(
  TARGETS optical optical_input
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

