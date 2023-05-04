set(top "${CMAKE_SOURCE_DIR}/Src")

set(gen-basis-sources

  gen-basis.F
  handlers.f
  
)
list(
    APPEND
    gen-basis-sources

    ${top}/alloc.F90
    ${top}/arw.f
    ${top}/atm_types.f
    ${top}/atmparams.f
    ${top}/atom.F
    ${top}/atom_options.F90
    ${top}/basis_io.F
    ${top}/basis_specs.f
    ${top}/basis_types.f
    ${top}/bessph.f
    ${top}/cellsubs.f
    ${top}/chemical.f
    ${top}/dftu_specs.f
    ${top}/dot.f
    ${top}/files.f
    ${top}/get_kpoints_scale.f90
    ${top}/hamann.f90
    ${top}/interpolation.f90
    ${top}/io.f
    ${top}/minvec.f
    ${top}/m_cite.F90
    ${top}/m_fft_gpfa.F
    ${top}/m_filter.f90
    ${top}/m_io.f
    ${top}/m_spin.F90
    ${top}/m_walltime.f90
    ${top}/memory.F
    ${top}/memory_log.F90
    ${top}/nag.f
    ${top}/parallel.F
    ${top}/parsing.f
    ${top}/periodic_table.f
    ${top}/precision.F
    ${top}/pxf.F90
    ${top}/radfft.f
    ${top}/radial.f
    ${top}/reclat.f
    ${top}/read_xc_info.F
    ${top}/siesta_geom.F90
    ${top}/sorting.f
    ${top}/sys.F
    ${top}/xml.f
    ${top}/m_spin.F90
    ${top}/m_vee_integrals.F90 

)

add_executable(
   gen-basis
   ${gen-basis-sources}
)

#------

set(ioncat-sources

  ioncat.f
  handlers.f
  
)

list(
   APPEND
   ioncat-sources

   ${top}/m_getopts.f90 ${top}/alloc.F90 ${top}/arw.f ${top}/atm_types.f ${top}/atmparams.f ${top}/atom.F
   ${top}/atmfuncs.f  ${top}/atom_options.F90  ${top}/basis_io.F   ${top}/basis_specs.f
   ${top}/basis_types.f   ${top}/bessph.f 
   ${top}/cellsubs.f  ${top}/chkdim.f ${top}/chemical.f   ${top}/dftu_specs.f   ${top}/dot.f
   ${top}/files.f   ${top}/get_kpoints_scale.f90
   ${top}/hamann.f90   ${top}/interpolation.f90   ${top}/io.f   ${top}/minvec.f   ${top}/m_cite.F90
   ${top}/m_fft_gpfa.F ${top}/m_filter.f90   ${top}/m_io.f ${top}/m_walltime.f90   ${top}/memory.F
   ${top}/memory_log.F90   ${top}/nag.f  ${top}/parallel.F   ${top}/parsing.f   ${top}/periodic_table.f
   ${top}/precision.F   ${top}/pxf.F90  ${top}/radfft.f   ${top}/radial.f   ${top}/reclat.f
   ${top}/read_xc_info.F   ${top}/siesta_geom.F90   ${top}/sorting.f 
   ${top}/sys.F   ${top}/xml.f   ${top}/m_spin.F90  ${top}/posix_calls.f90
   ${top}/m_vee_integrals.F90 ${top}/version.F90 ${top}/spher_harm.f
    
)

add_executable(
   ioncat
   ${ioncat-sources}
)


install(
  TARGETS gen-basis ioncat
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

#
# lapack is needed for m_filter...
#
 target_link_libraries(
  gen-basis
  PRIVATE
  libfdf::libfdf
  libgridxc::libgridxc
  libpsml::libpsml
  ${PROJECT_NAME}-libncps
  ${PROJECT_NAME}-libpsop
  ${PROJECT_NAME}-libunits
  LAPACK::LAPACK
 )

 target_link_libraries(
  ioncat
  PRIVATE
  libfdf::libfdf
  libgridxc::libgridxc
  libpsml::libpsml
  ${PROJECT_NAME}-libncps
  ${PROJECT_NAME}-libpsop
  ${PROJECT_NAME}-libunits
  LAPACK::LAPACK
 )


