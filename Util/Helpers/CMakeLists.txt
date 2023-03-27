#
set(top_src_dir "${CMAKE_SOURCE_DIR}/Src")

add_executable(get_chem_labels
   ${top_src_dir}/parallel.F
   ${top_src_dir}/chemical.f
   ${top_src_dir}/atom_options.F90
   ${top_src_dir}/m_io.f
   ${top_src_dir}/io.f
   ${top_src_dir}/alloc.F90
   ${top_src_dir}/memory.F
   ${top_src_dir}/memory_log.F90
   ${top_src_dir}/periodic_table.f
   ${top_src_dir}/precision.F
   ${top_src_dir}/pxf.F90
   ${top_src_dir}/units.f90
   get_chem_labels.f90
   local_sys.f
   handlers.f
)

target_link_libraries(get_chem_labels libfdf::libfdf)


install(
  TARGETS get_chem_labels
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

