#
# FIXME: The names in this directory are very confusing
#        Only a subset of the executables are built
#
set(top_srcdir "${CMAKE_SOURCE_DIR}/Src")

if (WITH_MPI)
#
# Phonons: Uses the MPI dispatch
#
add_executable(phonons
    phonons.f90
    handlers.f

    ${top_srcdir}/fsiesta_mpi.F90
)

target_link_libraries(phonons ${PROJECT_NAME}-libsiesta)

install(
  TARGETS phonons
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

# Simple MPI-aware driver

add_executable(mpi_driver
    simple_mpi_parallel.f90
    handlers.f

    ${top_srcdir}/fsiesta_mpi.F90
)

target_link_libraries(mpi_driver ${PROJECT_NAME}-libsiesta)

install(
  TARGETS mpi_driver
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

#
# Sockets interface to parallel-version of Siesta
#
add_executable(sockets_parallel
    simple_parallel.f90

    ${top_srcdir}/fsiesta_sockets.F90
    ${top_srcdir}/fsockets.f90
    ${top_srcdir}/sockets.c
    ${top_srcdir}/posix_calls.f90
)

install(
  TARGETS sockets_parallel
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

# Sockets interface to parallel-version of Siesta
#
add_executable(pipes_parallel
    simple_parallel.f90

    ${top_srcdir}/fsiesta_pipes.F90
    ${top_srcdir}/pxf.F90
    ${top_srcdir}/posix_calls.f90
)

install(
  TARGETS pipes_parallel
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

else(WITH_MPI)

#
add_executable(sockets_serial
    simple_serial.f90

    ${top_srcdir}/fsiesta_sockets.F90
    ${top_srcdir}/fsockets.f90
    ${top_srcdir}/sockets.c
    ${top_srcdir}/posix_calls.f90
)

install(
  TARGETS sockets_serial
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

#
add_executable(pipes_serial
    simple_serial.f90

    ${top_srcdir}/fsiesta_pipes.F90
    ${top_srcdir}/pxf.F90
    ${top_srcdir}/posix_calls.f90
)

install(
  TARGETS pipes_serial
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

endif(WITH_MPI)




  
