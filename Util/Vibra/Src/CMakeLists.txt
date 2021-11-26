set(top_src_dir "${CMAKE_SOURCE_DIR}/Src")

add_executable(
   fcbuild
   fcbuild.f
   local_sys.F
   recoor.f
   ${top_src_dir}/io.f
   ${top_src_dir}/m_io.f
   ${top_src_dir}/precision.F
   ${top_src_dir}/pxf.F90
   ${top_src_dir}/reclat.f
   ${top_src_dir}/units.f90

)
add_executable(
   vibra
   hermdp.F
   klines.f
   local_sys.F
   outbands.f
   recoor.f
   vibra.F
   ${top_src_dir}/io.f
   ${top_src_dir}/m_io.f
   ${top_src_dir}/precision.F
   ${top_src_dir}/pxf.F90
   ${top_src_dir}/reclat.f
   ${top_src_dir}/units.f90

)

 target_link_libraries(fcbuild
    ${PROJECT_NAME}-libfdf
    )
 target_link_libraries(vibra
    ${PROJECT_NAME}-libfdf
    )
 target_link_libraries(vibra
    $ENV{LAPACK_LIBS}
    )
    
 target_include_directories(
     fcbuild
     PUBLIC
    "${CMAKE_BINARY_DIR}/mod-fdf")
   
 target_include_directories(
     vibra
     PUBLIC
    "${CMAKE_BINARY_DIR}/mod-fdf")
   
install(
  TARGETS fcbuild vibra
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

