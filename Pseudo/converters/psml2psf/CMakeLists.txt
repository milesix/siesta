#
set(top_src_dir "${CMAKE_SOURCE_DIR}/Src")

add_executable(psml2psf
   ${top_src_dir}/m_getopts.f90

   psml2psf.f90
   handlers.f
)

target_link_libraries(psml2psf
                      ${PROJECT_NAME}-libncps
		      libpsml::libpsml)

install(
  TARGETS psml2psf
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

