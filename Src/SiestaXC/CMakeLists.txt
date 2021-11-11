set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND srcs


  "${dir}/am05.f90"
  "${dir}/array.F90"
  "${dir}/atomxc.F90"
  "${dir}/cellxc.F90"
  "${dir}/debugxc.f90"
  "${dir}/fft3d.F90"
  "${dir}/fftr.F90"
  "${dir}/ggaxc.f"
  "${dir}/ldaxc.f"
  "${dir}/mesh1d.f90"
  "${dir}/mesh3d.F90"
  "${dir}/siestaxc.F90"
  "${dir}/vdwxc.F90"
  "${dir}/vv_vdwxc.F90"
  "${dir}/xcmod.f90"

)

set(srcs "${srcs}" PARENT_SCOPE)
