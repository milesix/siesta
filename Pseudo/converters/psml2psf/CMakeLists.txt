#
set(top_src_dir "${PROJECT_SOURCE_DIR}/Src")

add_executable(psml2psf
   ${top_src_dir}/m_getopts.f90

   psml2psf.f90
)

target_link_libraries(psml2psf
                      ${PROJECT_NAME}-libncps
                      ${PROJECT_NAME}-libsys
		      libpsml::libpsml)

install(
  TARGETS psml2psf
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

