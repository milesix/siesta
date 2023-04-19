add_executable( optical optical_cond.f )
add_executable( optical_input optical_input.f )
target_link_libraries(optical_input
  PRIVATE ${PROJECT_NAME}-libunits)

install(
  TARGETS optical optical_input
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

