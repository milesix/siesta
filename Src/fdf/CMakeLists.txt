set(fdf_sources
  fdf.F90
  io_fdf.F90
  parse.F90
  prec.F90
  utils.F90
)

add_library(
  ${PROJECT_NAME}-libfdf
  ${fdf_sources}  
)

# Use a known top-level place to store the modules, since
# we have to have two separate versions (see below)

set_target_properties(
  ${PROJECT_NAME}-libfdf
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod-fdf"
)

#–----------------------------
# We now prepare a MPI version. This is needed since fdf contains
# some MPI code (it is removed in the library version)

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

set_target_properties(
  ${PROJECT_NAME}-libfdf-mpi
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod-fdf-mpi"
 )
endif()

