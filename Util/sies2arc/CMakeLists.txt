#
add_executable(sies2arc
   sies2arc.f linepro.f wtof.f cell.f uncell.f
)

install(
  TARGETS sies2arc
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

