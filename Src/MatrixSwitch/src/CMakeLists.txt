#
set(ms_sources

  MatrixSwitch.F90
  MatrixSwitch_m_add.F90
  MatrixSwitch_m_copy.F90
  MatrixSwitch_m_register.F90
  MatrixSwitch_m_set.F90
  MatrixSwitch_mm_multiply.F90
  MatrixSwitch_ops.F90
  MatrixSwitch_wrapper.F90
  MatrixSwitch_wrapper_params.F90
)

 add_library(
  ${PROJECT_NAME}-libms
  ${ms_sources}  
 )

if (WITH_MPI)
 # specific preprocessor variable
 target_compile_definitions(
  ${PROJECT_NAME}-libms
  PUBLIC
  HAVE_MPI
  HAVE_SCALAPACK
 )

# Note that this library does not use
# mpi_siesta...
# target_link_libraries(
#  ${PROJECT_NAME}-libms
#  PRIVATE
#  mpi_siesta
# )

  # Hence, add ground stuff for MPI
  # This enables correct compilation and also correct flags for
  # MPI builds.
  # Also, it ensures that the include directories are correctly setup

 target_compile_options(${PROJECT_NAME}-libms PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
 target_include_directories(${PROJECT_NAME}-libms PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
 target_link_directories(${PROJECT_NAME}-libms PUBLIC ${MPI_Fortran_LIBRARY_DIRS})
 target_link_libraries(${PROJECT_NAME}-libms PUBLIC  ${MPI_Fortran_LIBRARIES})
  

endif()

target_include_directories(
  ${PROJECT_NAME}-libms
  INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}
)

