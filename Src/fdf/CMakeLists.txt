set(fdf_sources
  fdf.F90
  io_fdf.F90
  parse.F90
  prec.F90
  utils.F90
  )

add_library(fdf ${fdf_sources})

# For stuff that depends on libfdf, add these include directories
target_include_directories(fdf INTERFACE ${CMAKE_CURRENT_BINARY_DIR})

