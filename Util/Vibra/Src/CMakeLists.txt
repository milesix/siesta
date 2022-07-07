#
# There are two targets in this directory. The 'OBJECT library' idiom
# was needed to avoid ninja complaining about "two paths for the same
# .mod file"
#
# Since the top objects are compiled separately, their 'sys' module
# must be the top one. Here only external 'die' et al are bundled in
# the 'sys_externals.f90' file.
#
set(top_src_dir "${CMAKE_SOURCE_DIR}/Src")

add_library(vibra_top_objs OBJECT
   ${top_src_dir}/io.f
   ${top_src_dir}/m_io.f
   ${top_src_dir}/precision.F
   ${top_src_dir}/pxf.F90
   ${top_src_dir}/reclat.f
   ${top_src_dir}/sys.F
   ${top_src_dir}/units.f90
)

add_executable(
    fcbuild
    fcbuild.f
    recoor.f
    sys_externals.f90
)

add_executable(
   vibra
   hermdp.F
   klines.f
   outbands.f
   recoor.f
   sys_externals.f90
   vibra.F
)

target_link_libraries(fcbuild
    vibra_top_objs
    fdf-serial
    )
target_link_libraries(vibra
    vibra_top_objs
    fdf-serial
    $ENV{LAPACK_LIBS}
    )
    
install(
  TARGETS fcbuild vibra
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

