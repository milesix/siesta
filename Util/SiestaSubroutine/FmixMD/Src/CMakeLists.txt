#
set(top_srcdir "${PROJECT_SOURCE_DIR}/Src")

# Driver for mixed (coarse/fine)  MD using pipes
#
add_executable(fmixmd-driver

    driver.f90
    fmixmd.f90
    sample.f90
    ${top_srcdir}/fsiesta_pipes.F90
    ${top_srcdir}/posix_calls.f90
    ${top_srcdir}/pxf.F90
)

install(
  TARGETS fmixmd-driver
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
