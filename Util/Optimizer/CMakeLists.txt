add_library(aux_optim
    io.f
    parse.f
    minimizer.f90
    vars_module.f90
    precision.f90
    system.f90
    sys.f
)    

add_executable( simplex
    amoeba.f
    simplex.f90
)

add_executable( swarm
    swarm.f90
)

target_link_libraries(swarm PRIVATE OpenMP::OpenMP_Fortran)
target_link_libraries(swarm PRIVATE aux_optim)
target_link_libraries(simplex PRIVATE aux_optim)

install(
  TARGETS simplex swarm
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

