set(dir "${CMAKE_CURRENT_SOURCE_DIR}")
set(sources

  am05.f90
  array.F90
  atomxc.F90
  cellxc.F90
  debugxc.f90
  fft3d.F90
  fftr.F90
  ggaxc.f
  ldaxc.f
  mesh1d.f90
  mesh3d.F90
  siestaxc.F90
  vdwxc.F90
  vv_vdwxc.F90
  xcmod.f90
  sys.F
  io.f
)
list(
    APPEND
    sources

    "${dir}/../alloc.F90"
    "${dir}/../bessph.f"
    "${dir}/../cellsubs.f"
    "${dir}/../debugmpi.F"
    "${dir}/../chkgmx.f"
    "${dir}/../fft1d.F"
    "${dir}/../interpolation.f90"
    "${dir}/../minvec.f"
    "${dir}/../m_io.f"
    "${dir}/../moreParallelSubs.F90"
    "${dir}/../class_OrbitalDistribution.F90"
    "${dir}/../m_timer.F90"
    "${dir}/../timer_tree.F90"
    "${dir}/../m_walltime.f90"
    "${dir}/../parallel.F"
    "${dir}/../precision.F"
    "${dir}/../pxf.F90"
    "${dir}/../radfft.f"
    "${dir}/../m_fft_gpfa.F"
    "${dir}/../sorting.f"
    "${dir}/../timer.F90"

)


add_library(
  ${PROJECT_NAME}-libSiestaXC
  ${sources}
)

target_include_directories(
  ${PROJECT_NAME}-libSiestaXC
  INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}
)  

if (WITH_MPI)

  target_compile_definitions(
    ${PROJECT_NAME}-libSiestaXC
    PRIVATE
    MPI
  )

 target_link_libraries(
  ${PROJECT_NAME}-libSiestaXC
  PRIVATE
  ${PROJECT_NAME}-mpi-int
  )
  
endif()

