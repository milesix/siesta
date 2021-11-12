add_library(
  ${PROJECT_NAME}-libfdf
   
  fdf.F90
  io_fdf.F90
  parse.F90
  prec.F90
  utils.F90
)

set_target_properties(
  ${PROJECT_NAME}-libfdf
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod-fdf"
)

if (WITH_MPI)
 add_library(
  ${PROJECT_NAME}-libfdf-mpi
   
  fdf.F90
  io_fdf.F90
  parse.F90
  prec.F90
  utils.F90
  "${CMAKE_CURRENT_SOURCE_DIR}/../MPI/mpi_siesta.F90"
 )

  target_compile_definitions(
  ${PROJECT_NAME}-libfdf-mpi
  PUBLIC
  MPI
  NO_MPI_INTERFACES
  )

set_target_properties(
  ${PROJECT_NAME}-libfdf-mpi
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod-fdf-mpi"
 )
endif()

