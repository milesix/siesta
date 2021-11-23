#
set(fdf_sources
  fdf.F90
  io_fdf.F90
  parse.F90
  prec.F90
  utils.F90
)

#
# Plain, serial version
#
add_library(
  ${PROJECT_NAME}-libfdf
  ${fdf_sources}  
)

# Use a known top-level place to store the modules, since
# we have to have two separate versions to link to different
# targets in the suite (see below)

set_target_properties(
  ${PROJECT_NAME}-libfdf
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod-fdf"
)

#–----------------------------
# We now prepare a MPI version. This is needed since fdf contains
# some MPI code for broadcasting the fdf data object.
# This MPI code is removed in the stand-alone library, and relegated
# to an extra routine that clients can use on their side.
#

if (WITH_MPI)
 add_library(
  ${PROJECT_NAME}-libfdf-mpi
  ${fdf_sources}  
 )

 target_compile_definitions(
  ${PROJECT_NAME}-libfdf-mpi
  PUBLIC
  MPI
 )

 target_link_libraries(
  ${PROJECT_NAME}-libfdf-mpi
  PRIVATE
  ${PROJECT_NAME}-mpi-int
 )

#
# Put the module in a different place
#
set_target_properties(
  ${PROJECT_NAME}-libfdf-mpi
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod-fdf-mpi"
 )
endif()

