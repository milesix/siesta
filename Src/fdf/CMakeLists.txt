set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND srcs
  "${dir}/fdf.F90"
  "${dir}/io_fdf.F90"
  "${dir}/parse.F90"
  "${dir}/prec.F90"
  "${dir}/utils.F90"
)

set(srcs "${srcs}" PARENT_SCOPE)
