#
set(top_src_dir "${CMAKE_SOURCE_DIR}/Src")

#-----------------------------------------
add_executable(readwfx
   ${top_src_dir}/m_getopts.f90
   readwfx.f90
)

add_executable(readwf readwf.f)  
add_executable(wfs2wfsx  wfs2wfsx.f)
add_executable(wfsx2wfs  wfsx2wfs.f)

install(TARGETS
        readwf
        readwfx
        wfs2wfsx
        wfsx2wfs
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )

if (WITH_NETCDF)

  add_executable(wfsnc2wfsx wfsnc2wfsx.F90)

  target_link_libraries(wfsnc2wfsx PRIVATE
    NetCDF::NetCDF_Fortran
    ${PROJECT_NAME}-libunits
    )
  target_compile_definitions(wfsnc2wfsx  PRIVATE CDF )
		     
  install(TARGETS wfsnc2wfsx RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )

endif(WITH_NETCDF)
