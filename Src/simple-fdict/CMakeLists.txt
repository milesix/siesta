add_library(
  ${PROJECT_NAME}-libfdict
   
  variable.f90
  dictionary.f90
  )

target_include_directories(
  ${PROJECT_NAME}-libfdict
  INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}
)


