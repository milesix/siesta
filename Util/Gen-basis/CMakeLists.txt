set(top_srcdir "${CMAKE_SOURCE_DIR}/Src")

set(sources

  gen-basis.F
  handlers.f
  
)
list(
    APPEND
    sources


    "${top_srcdir}/alloc.F90"
    "${top_srcdir}/arw.f"
    "${top_srcdir}/atm_types.f"
    "${top_srcdir}/atmparams.f"
    "${top_srcdir}/atom.F"
    "${top_srcdir}/atom_options.F90"
    "${top_srcdir}/basis_io.F"
    "${top_srcdir}/basis_specs.f"
    "${top_srcdir}/basis_types.f"
    "${top_srcdir}/bessph.f"
    "${top_srcdir}/cellsubs.f"
    "${top_srcdir}/chemical.f"
    "${top_srcdir}/dftu_specs.f"
    "${top_srcdir}/dot.f"
    "${top_srcdir}/files.f"
    "${top_srcdir}/get_kpoints_scale.f90"
    "${top_srcdir}/hamann.f90"
    "${top_srcdir}/interpolation.f90"
    "${top_srcdir}/io.f"
    "${top_srcdir}/minvec.f"
    "${top_srcdir}/m_cite.F90"
    "${top_srcdir}/m_fft_gpfa.F"
    "${top_srcdir}/m_filter.f90"
    "${top_srcdir}/m_io.f"
    "${top_srcdir}/m_spin.F90"
    "${top_srcdir}/m_walltime.f90"
    "${top_srcdir}/memory.F"
    "${top_srcdir}/memory_log.F90"
    "${top_srcdir}/nag.f"
    "${top_srcdir}/parallel.F"
    "${top_srcdir}/parsing.f"
    "${top_srcdir}/periodic_table.f"
    "${top_srcdir}/precision.F"
    "${top_srcdir}/pxf.F90"
    "${top_srcdir}/radfft.f"
    "${top_srcdir}/radial.f"
    "${top_srcdir}/reclat.f"
    "${top_srcdir}/read_xc_info.F"
    "${top_srcdir}/siesta_geom.F90"
    "${top_srcdir}/sorting.f"
    "${top_srcdir}/sys.F"
    "${top_srcdir}/units.f90"
    "${top_srcdir}/xml.f"
     ${top_srcdir}/m_spin.F90
     ${top_srcdir}/m_vee_integrals.F90 

)

add_executable(
   gen-basis
   ${sources}
)

install(
  TARGETS gen-basis
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

#
# lapack is needed for m_filter...
#
 target_link_libraries(
  gen-basis
  PRIVATE
  fdf
  gridxc
  libpsml::libpsml
  ${PROJECT_NAME}-libncps
  ${PROJECT_NAME}-libpsop
  "$ENV{LAPACK_LIBS}"
 )


