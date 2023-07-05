#
set(top_src_dir "${CMAKE_SOURCE_DIR}/Src")

add_executable(pvtsp

   ${top_src_dir}/precision.F
   ${top_src_dir}/parallel.F
   ${top_src_dir}/alloc.F90
   ${top_src_dir}/m_io.f
   ${top_src_dir}/reclat.f
   ${top_src_dir}/io.f
   ${top_src_dir}/memory_log.F90
   ${top_src_dir}/memory.F
   ${top_src_dir}/m_uuid.f90
   ${top_src_dir}/object_debug.F90
   ${top_src_dir}/class_Sparsity.F90
   ${top_src_dir}/class_OrbitalDistribution.F90
   ${top_src_dir}/class_Data1D.F90
   ${top_src_dir}/class_Data2D.F90
   ${top_src_dir}/class_SpData1D.F90
   ${top_src_dir}/class_SpData2D.F90
   ${top_src_dir}/intrinsic_missing.F90
   ${top_src_dir}/geom_helper.f90
   ${top_src_dir}/io_sparse.F90
   ${top_src_dir}/m_os.F90
   ${top_src_dir}/create_Sparsity_SC.F90
   ${top_src_dir}/m_sparse.F90
   ${top_src_dir}/m_ts_io.F90
   ${top_src_dir}/m_region.F90
   ${top_src_dir}/m_pivot_methods.F90
   ${top_src_dir}/m_pivot.F90
   ${top_src_dir}/m_interpolate.F90
   ${top_src_dir}/m_sparsity_handling.F90

   pvtsp.F90
   local_timer.f90
)

target_link_libraries(pvtsp PRIVATE ${PROJECT_NAME}-libsys)

install(
  TARGETS pvtsp
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

