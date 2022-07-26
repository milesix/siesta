#
# Use a pre-generated source file (see top of file)
#
add_library(
  ${PROJECT_NAME}-libncdf
   
  netcdf_ncdf.F90
  )

target_link_libraries(
   ${PROJECT_NAME}-libncdf
   # 'PUBLIC' to bring in fdict
   PRIVATE
   ${PROJECT_NAME}-libfdict
   NetCDF::NetCDF_Fortran
)

target_include_directories(
  ${PROJECT_NAME}-libncdf
  INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}
)

if(WITH_NETCDF_PARALLEL)
  target_compile_definitions(
  "${PROJECT_NAME}-libncdf"
  PRIVATE
  NCDF_PARALLEL
  )
endif()


