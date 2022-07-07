set(fdf_sources
  fdf.F90
  io_fdf.F90
  parse.F90
  prec.F90
  utils.F90
  )

# fdf library without MPI support (serial)
# Plain, serial version
add_library(fdf-serial ${fdf_sources})

# Use a known top-level place to store the modules, since
# we have to have two separate versions to link to different
# targets in the suite (see below)

set_target_properties(
  fdf-serial
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod-fdf-serial"
)

# For stuff that depends on libfdf, add these include directories
target_include_directories(fdf-serial INTERFACE "${CMAKE_BINARY_DIR}/mod-fdf-serial")

# In case the user also has MPI enabled, also build this one
# In principle we should have a stricter CMakeLists.txt file
# that encompass fdf as a separate project. This would make it
# easier
if( WITH_MPI )

  # fdf library with MPI support
  # TODO missing features such as -DCLUSTER and -DBLOCKING
  add_library(fdf-mpi ${fdf_sources})

  set_target_properties(
   fdf-mpi
   PROPERTIES
   Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod-fdf-mpi"
)
  # For stuff that depends on libfdf, add these include directories
  target_include_directories(fdf-mpi INTERFACE "${CMAKE_BINARY_DIR}/mod-fdf-mpi" )

  # Add stuff for MPI
  # This enables correct compilation and also correct flags for
  # MPI builds.
  # Also, it ensures that the include directories are correctly setup
  target_compile_options(fdf-mpi PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
  target_include_directories(fdf-mpi PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
  target_link_directories(fdf-mpi PUBLIC ${MPI_Fortran_LIBRARY_DIRS})
  target_link_libraries(fdf-mpi PUBLIC ${MPI_Fortran_LIBRARIES})

  # TODO this should be adapted when fdf discards the mpi_siesta dependency
  target_link_libraries(fdf-mpi PUBLIC Siesta::MPI)

  # Add -DMPI flag for the local compiled files
  # This need not be the same flag in other libraries or codes
  target_compile_definitions(fdf-mpi PRIVATE -DMPI)

endif()

