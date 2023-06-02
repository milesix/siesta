# This internal library uses a "pre-processed" file from the fdict
# distribution, containing the needed variable types.
  
add_library(
  ${PROJECT_NAME}-libfdict
   
  variable.F90
  dictionary.f90
  )

target_include_directories(
  ${PROJECT_NAME}-libfdict
  INTERFACE
  "${CMAKE_CURRENT_BINARY_DIR}"
)


