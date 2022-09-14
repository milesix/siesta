
add_executable(pi3
  # Sources
  pi3.F
  handlers.f
  )

target_link_libraries(pi3 PRIVATE mpi_siesta)
add_test(NAME pi3_Runs COMMAND pi3)
add_test(NAME pi3_Runs_inMPI
         COMMAND
	 ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
         ${MPIEXEC_PREFLAGS} pi3 ${MPIEXEC_POSTFLAGS})


add_executable(newcomm
  # Sources
  newcomm.F
  handlers.f
  )

target_link_libraries(newcomm PRIVATE mpi_siesta)
add_test(NAME newcomm_Runs_inMPI
         COMMAND
	 ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
         ${MPIEXEC_PREFLAGS} newcomm ${MPIEXEC_POSTFLAGS})

add_executable(blacs_prb
  # Sources
  blacs_prb.f90
  )

target_link_libraries(blacs_prb
  PRIVATE
  "$ENV{SCALAPACK_LIBS}"
  MPI::MPI_Fortran
)
add_test(NAME blacs_prb_Runs_inMPI
         COMMAND
	 ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
         ${MPIEXEC_PREFLAGS} blacs_prb ${MPIEXEC_POSTFLAGS})



