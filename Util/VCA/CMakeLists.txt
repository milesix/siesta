#
# There are two targets in this directory. The 'OBJECT library' idiom
# was needed to avoid ninja complaining about "two paths for the same
# .mod file"
#
set(top_src_dir "${PROJECT_SOURCE_DIR}/Src")

add_library(vca_objs OBJECT
   ${top_src_dir}/interpolation.f90
   ${top_src_dir}/m_io.f
   ${top_src_dir}/io.f
   ${top_src_dir}/periodic_table.f
   ${top_src_dir}/precision.F
)
target_link_libraries(vca_objs
                      PRIVATE
		      ${PROJECT_NAME}-libsys)

add_executable( mixps mixps.f )
add_executable( fractional fractional.f)


target_link_libraries(mixps    PRIVATE
                      vca_objs
                      ${PROJECT_NAME}-libncps
		      ${PROJECT_NAME}-libsys)
target_link_libraries(fractional PRIVATE
                      vca_objs
                      ${PROJECT_NAME}-libncps
		      ${PROJECT_NAME}-libsys)
    
install(
  TARGETS mixps fractional
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

