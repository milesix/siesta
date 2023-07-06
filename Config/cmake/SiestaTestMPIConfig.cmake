if (WITH_MPI)
  set(max_num_ranks "${MPIEXEC_MAX_NUMPROCS}")
  message(STATUS "Maximum number of usable MPI ranks: ${max_num_ranks}")
  if(NOT DEFINED TEST_MPI_PROCS)
    if (${max_num_ranks} GREATER "4")
       set(num_ranks 4)
    else()
       set(num_ranks ${max_num_ranks})
    endif()
  else()
    if ("${TEST_MPI_PROCS}" GREATER "${max_num_ranks}")
       set(num_ranks ${max_num_ranks})
    else()
       set(num_ranks "${TEST_MPI_PROCS}")
    endif()
  endif()
  message(STATUS "Number of MPI ranks in tests: ${num_ranks}")
  set(COMMAND_PREAMBLE
         "${MPIEXEC_EXECUTABLE}
          ${MPIEXEC_NUMPROC_FLAG} ${num_ranks}
          ${MPIEXEC_PREFLAGS}"
     )
  set(COMMAND_POSTAMBLE "${MPIEXEC_POSTFLAGS}")
  set(TEST_LABEL_MPI "_mpi")     
else(WITH_MPI)

 set(COMMAND_PREAMBLE)
 set(COMMAND_POSTAMBLE)
 set(TEST_LABEL_MPI)

endif(WITH_MPI)
