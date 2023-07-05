set(top_src_dir "${PROJECT_SOURCE_DIR}/Src")
add_executable(
   Eig2DOS
   Eig2DOS.f90
   ${top_src_dir}/m_getopts.f90

)

install(
  TARGETS Eig2DOS
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

