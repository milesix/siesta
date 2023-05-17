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
   ${top_src_dir}/units.f90
   grimme.f90
   local_timer.f90
)

target_link_libraries(fdf2grimme
                      PRIVATE
		      libfdf::libfdf
      		      ${PROJECT_NAME}-libsys
)

install(
  TARGETS fdf2grimme
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

