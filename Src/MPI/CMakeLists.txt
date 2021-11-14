add_library(
  ${PROJECT_NAME}-mpi-int
   
  mpi_siesta.F90
)

# For now, without explicit interfaces


target_include_directories(
  ${PROJECT_NAME}-mpi-int
  INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}
)  

target_compile_definitions(
    ${PROJECT_NAME}-mpi-int
    PRIVATE
    NO_MPI_INTERFACES
)
