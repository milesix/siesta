set(top_srcdir "${CMAKE_SOURCE_DIR}/Src" )

set(sources

  defs_basis.f90
  defs_common.f90
  io.f
  iorho.f
  macroave.f
  thetaft.f
  surpla.f
  volcel.f
  hdr_io.f90
  local_die.f
  dpnint1.f90

)

list(
  APPEND
  sources

  ${top_srcdir}/interpolation.f90
  ${top_srcdir}/m_fft_gpfa.F

)

add_executable(macroave "${sources}")

TARGET_LINK_LIBRARIES(
    macroave
    PRIVATE
    fdf-serial
)

install(
  TARGETS macroave
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )


