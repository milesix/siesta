add_library(
  ${PROJECT_NAME}-libncdf
   
  netcdf_ncdf.f90
  )

target_include_directories(
  ${PROJECT_NAME}-libncdf
  INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}
)


