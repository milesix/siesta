#
set(top_src_dir "${CMAKE_SOURCE_DIR}/Src")

add_library(sockets_objs OBJECT
   ${top_src_dir}/fsockets.f90
   ${top_src_dir}/sockets.c
)

add_executable( f2fmaster f2fmaster.f90 )
add_executable( f2fslave  f2fslave.f90 )

target_link_libraries(f2fmaster sockets_objs )
target_link_libraries(f2fslave sockets_objs )
    
install(
  TARGETS f2fmaster f2fslave
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

