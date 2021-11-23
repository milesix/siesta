#
# I got this "pre-processed" file from the ncdf distribution, but I do
# not remember under which conditions.
#
add_library(
  ${PROJECT_NAME}-libncdf
   
  netcdf_ncdf.f90
  )

target_link_libraries(
   ${PROJECT_NAME}-libncdf
   PRIVATE
   ${PROJECT_NAME}-libfdict
)

target_include_directories(
  ${PROJECT_NAME}-libncdf
  INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}
)

# The conditional should be superfluous
if(WITH_NETCDF)
  target_compile_definitions(
  "${PROJECT_NAME}-libncdf"
  PRIVATE
  CDF
  NCDF
  NCDF_4
  )
  target_include_directories(
    ${PROJECT_NAME}-libncdf
    PRIVATE
    ${NETCDF_INCLUDEDIR}
    ${NETCDF_INCLUDE_DIRS}
)
#  target_link_libraries(
#     ${PROJECT_NAME}-siesta
#     PRIVATE
#    ${NETCDF_LINK_LIBRARIES}
# )
endif()


