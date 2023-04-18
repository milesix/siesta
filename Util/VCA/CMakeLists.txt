#
# There are two targets in this directory. The 'OBJECT library' idiom
# was needed to avoid ninja complaining about "two paths for the same
# .mod file"
#
set(top_src_dir "${CMAKE_SOURCE_DIR}/Src")

add_library(vca_objs OBJECT
   ${top_src_dir}/interpolation.f90
   ${top_src_dir}/m_io.f
   ${top_src_dir}/io.f
   ${top_src_dir}/periodic_table.f
   ${top_src_dir}/precision.F
   local_sys.f
   handlers.f
)
target_link_libraries(vca_objs PUBLIC
			       ${PROJECT_NAME}-libncps)

add_executable( mixps mixps.f )
add_executable( fractional fractional.f)


# These will pick up the ncps dependency from vca_objs
target_link_libraries(mixps      vca_objs )
target_link_libraries(fractional vca_objs )
    
install(
  TARGETS mixps fractional
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

