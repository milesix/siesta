if (WITH_NETCDF)

  add_executable( lwf2cdf lwf2cdf.F90 )
  install(TARGETS lwf2cdf  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

endif()


