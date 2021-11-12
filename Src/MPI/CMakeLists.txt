set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND srcs
  "${dir}/mpi_siesta.F90"
)

set(srcs "${srcs}" PARENT_SCOPE)

