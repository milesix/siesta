#
# There are two targets in this directory. The 'OBJECT library' idiom
# was needed to avoid ninja complaining about "two paths for the same
# .mod file"
#
set(top_src_dir "${CMAKE_SOURCE_DIR}/Src")

add_library(bands_top_objs OBJECT
   ${top_src_dir}/m_getopts.f90
)

add_executable(
    gnubands
    gnubands.f90
)

add_executable(
   eigfat2plot
   eigfat2plot.f90
)

target_link_libraries(gnubands
    bands_top_objs
    )
target_link_libraries(eigfat2plot
    bands_top_objs
    )
    
install(
  TARGETS gnubands eigfat2plot
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

