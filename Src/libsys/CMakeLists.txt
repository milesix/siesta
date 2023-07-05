#
#
add_library(
  ${PROJECT_NAME}-libsys
   
  sys.F90
  external_entries.f90
  )

target_include_directories(
  ${PROJECT_NAME}-libsys
  INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}
)


