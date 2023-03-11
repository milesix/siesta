#
# This function uses 'num_ranks' as set in a previous
# call to TestMPIConfig.
#

# Note TARGET_FILE idiom for generality. This also covers the case in
# which the current directory is not in the PATH ('./prog' could work
# too). Otherwise mpiexec would not find it.

function(setup_test prog)

   set(mpiexec_args
          ${MPIEXEC_NUMPROC_FLAG} ${num_ranks}
          ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:${prog}>
          ${MPIEXEC_POSTFLAGS})
	  
   add_test( NAME    ${prog}_mpi_np_${num_ranks}
             COMMAND ${MPIEXEC_EXECUTABLE} ${mpiexec_args} )

endfunction()

add_executable(pi3
  pi3.F
  handlers.f
  )

target_link_libraries(pi3 PRIVATE mpi_siesta)
#
# To check "singleton" MPI runs (i.e., that we can run
# a parallel program without 'mpirun'). This is not
# possible in some systems (e.g. IBM's MPI)
#
add_test(NAME pi3_Runs_Singleton COMMAND pi3)

# Normal MPI run with ${num_ranks}
setup_test(pi3)

#
add_executable(newcomm
  newcomm.F
  handlers.f
  )

target_link_libraries(newcomm PRIVATE mpi_siesta)
setup_test(newcomm)

#
add_executable(blacs_prb  blacs_prb.f90 )
target_link_libraries(blacs_prb
  PRIVATE
  Scalapack::Scalapack
  MPI::MPI_Fortran
  LAPACK::LAPACK
)
setup_test(blacs_prb)

#
  add_executable(pblas_prb  pblas_prb.f90 )
  target_link_libraries(pblas_prb
    PRIVATE
    Scalapack::Scalapack
    MPI::MPI_Fortran
    LAPACK::LAPACK
  )

  file(COPY "${CMAKE_CURRENT_LIST_DIR}/pblas.dat" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

#
# Run the test automatically only if we have enough ranks available.
#
if("${num_ranks}" GREATER_EQUAL 4)
  setup_test(pblas_prb)
endif()


