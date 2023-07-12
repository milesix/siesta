add_library(objs_getopts OBJECT
            ${PROJECT_SOURCE_DIR}/Src/m_getopts.f90 )
	    
add_executable(dmfilter dmfilter.f90)
add_executable(dmbs2dm dmbs2dm.F90 )
target_link_libraries(dmbs2dm objs_getopts)
add_executable(dmUnblock dmUnblock.F90)
target_link_libraries(dmUnblock objs_getopts)
add_executable(dm_noncol_sign_flip4 dm_noncol_sign_flip4.f90)

install(
  TARGETS dmfilter dmbs2dm dmUnblock dm_noncol_sign_flip4
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

if(WITH_NETCDF)

 add_executable(dm2cdf dm2cdf.F90)
 add_executable(cdf2dm cdf2dm.F90)

 install(
  TARGETS dm2cdf cdf2dm
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
  
endif(WITH_NETCDF)


