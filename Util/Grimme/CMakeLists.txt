#
set(top_src_dir "${CMAKE_SOURCE_DIR}/Src")

add_executable(fdf2grimme
   ${top_src_dir}/parallel.F
   ${top_src_dir}/chemical.f
   ${top_src_dir}/atom_options.F90
   ${top_src_dir}/m_io.f
   ${top_src_dir}/io.f
   ${top_src_dir}/periodic_table.f
   ${top_src_dir}/precision.F
   grimme.f90
   local_sys.f90
   handlers.f
)

target_link_libraries(fdf2grimme
  libfdf::libfdf
  ${PROJECT_NAME}-libunits)

install(
  TARGETS fdf2grimme
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

