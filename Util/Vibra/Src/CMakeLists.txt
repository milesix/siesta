#
# There are two targets in this directory. The 'OBJECT library' idiom
# was needed to avoid ninja complaining about "two paths for the same
# .mod file"
#
set(top_src_dir "${CMAKE_SOURCE_DIR}/Src")

add_library(vibra_objs OBJECT
   ${top_src_dir}/io.f
   ${top_src_dir}/m_io.f
   ${top_src_dir}/precision.F
   ${top_src_dir}/reclat.f
   local_sys.f
   recoor.f
)

add_executable(
    fcbuild
    fcbuild.f
)

add_executable(
   vibra
   hermdp.F
   klines.f
   outbands.f
   vibra.F
)

target_link_libraries(
    vibra_objs
    PRIVATE
    libfdf::libfdf
    )
    
target_link_libraries(fcbuild
    vibra_objs
    libfdf::libfdf
    ${PROJECT_NAME}-libunits
    )
    
target_link_libraries(vibra
    vibra_objs
    libfdf::libfdf
    ${PROJECT_NAME}-libunits
    LAPACK::LAPACK
    )
    
install(
  TARGETS fcbuild vibra
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

