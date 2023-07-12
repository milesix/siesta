#
# FIXME: The names in this directory are very confusing
#        Only a subset of the executables are built
#
set(top_srcdir "${PROJECT_SOURCE_DIR}/Src")


#
# Create auxiliary libraries, with their own module directories,
# to satisfy Ninja. This is also good practice in general
# 
# -- Sockets
add_library(sockets-siesta-dispatch 
    ${top_srcdir}/fsiesta_sockets.F90
    ${top_srcdir}/fsockets.f90
    ${top_srcdir}/sockets.c
    ${top_srcdir}/posix_calls.f90
)
set_target_properties(sockets-siesta-dispatch
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/sockets-include")
target_include_directories(sockets-siesta-dispatch
  INTERFACE
  "${CMAKE_CURRENT_BINARY_DIR}/sockets-include")
#
# -- Pipes
add_library(pipes-siesta-dispatch 
    ${top_srcdir}/fsiesta_pipes.F90
    ${top_srcdir}/pxf.F90
    ${top_srcdir}/posix_calls.f90
)
set_target_properties(pipes-siesta-dispatch
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pipes-include")
target_include_directories(pipes-siesta-dispatch
  INTERFACE
  "${CMAKE_CURRENT_BINARY_DIR}/pipes-include")

if (WITH_MPI)

# -- MPI dispatch -- uses 'libSiestaForces'
add_library(mpi-siesta-dispatch  ${top_srcdir}/fsiesta_mpi.F90)
target_link_libraries(mpi-siesta-dispatch PUBLIC ${PROJECT_NAME}-libsiesta)

set_target_properties(mpi-siesta-dispatch
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mpi-include")
target_include_directories(mpi-siesta-dispatch
  INTERFACE
  "${CMAKE_CURRENT_BINARY_DIR}/mpi-include")

# ---- 
#
# Phonons: Uses the MPI dispatch
#
add_executable(phonons phonons.f90)
target_link_libraries(phonons PRIVATE mpi-siesta-dispatch)

install(
  TARGETS phonons
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

# Simple MPI-aware driver

add_executable(mpi_driver
    simple_mpi_parallel.f90
)
target_link_libraries(mpi_driver PRIVATE mpi-siesta-dispatch)

install(
  TARGETS mpi_driver
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

#
# Sockets interface to parallel-version of Siesta
#
add_executable(sockets_parallel  simple_parallel.f90)
target_link_libraries(sockets_parallel PRIVATE sockets-siesta-dispatch)

install(
  TARGETS sockets_parallel
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

# Pipes interface to parallel-version of Siesta
#
add_executable(pipes_parallel
    simple_parallel.f90
)
target_link_libraries(pipes_parallel PRIVATE pipes-siesta-dispatch)

install(
  TARGETS pipes_parallel
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

else(WITH_MPI)

#
add_executable(sockets_serial  simple_serial.f90)
target_link_libraries(sockets_serial PRIVATE sockets-siesta-dispatch)

install(
  TARGETS sockets_serial
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

#
add_executable(pipes_serial simple_serial.f90)
target_link_libraries(pipes_serial PRIVATE pipes-siesta-dispatch)

install(
  TARGETS pipes_serial
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

endif(WITH_MPI)




  
