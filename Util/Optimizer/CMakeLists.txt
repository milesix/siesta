add_executable( simplex
    io.f
    parse.f
    minimizer.f90
    vars_module.f90
    precision.f90
    amoeba.f
    system.f90
    sys.f
    simplex.f90
)

add_executable( swarm
    io.f
    parse.f
    minimizer.f90
    precision.f90
    vars_module.f90
    system.f90
    sys.f
    swarm.f90
)

target_link_libraries(swarm OpenMP::OpenMP_Fortran)

install(
  TARGETS simplex swarm
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

