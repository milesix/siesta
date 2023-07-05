set(top_srcdir "${CMAKE_SOURCE_DIR}/Src")

set(sources

  gen-basis.F
  
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
    "${top_srcdir}/domain_decom.F"
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
    "${top_srcdir}/moreParallelSubs.F90"
    "${top_srcdir}/nag.f"
    "${top_srcdir}/parallel.F"
    "${top_srcdir}/parallelsubs.F"
    "${top_srcdir}/parsing.f"
    "${top_srcdir}/periodic_table.f"
    "${top_srcdir}/precision.F"
    "${top_srcdir}/printmatrix.F"
    "${top_srcdir}/pxf.F90"
    "${top_srcdir}/radfft.f"
    "${top_srcdir}/radial.f"
    "${top_srcdir}/reclat.f"
    "${top_srcdir}/read_xc_info.F"
    "${top_srcdir}/schecomm.F"
    "${top_srcdir}/siesta_geom.F90"
    "${top_srcdir}/spatial.F"
    "${top_srcdir}/sparse_matrices.F90"
    "${top_srcdir}/sorting.f"
    "${top_srcdir}/timer_tree.f90"
    "${top_srcdir}/m_timer.F90"
    "${top_srcdir}/timer.F90"
    "${top_srcdir}/xml.f"
     ${top_srcdir}/m_spin.F90
     ${top_srcdir}/m_vee_integrals.F90 

  ${top_srcdir}/pspltm1.F
  ${top_srcdir}/qsort.F
  ${top_srcdir}/mmio.F
  ${top_srcdir}/class_Geometry.F90
  ${top_srcdir}/class_Data2D.F90
  ${top_srcdir}/class_OrbitalDistribution.F90
  ${top_srcdir}/class_Sparsity.F90
  ${top_srcdir}/class_SpData2D.F90
  ${top_srcdir}/class_Data1D.F90
  ${top_srcdir}/class_SpData1D.F90
  ${top_srcdir}/class_Pair_Geometry_SpData2D.F90
  ${top_srcdir}/class_Fstack_Pair_Geometry_SpData2D.F90
)

add_executable(
   gen-basis
   ${sources}
)

target_link_libraries(
  gen-basis
  PRIVATE
  ${PROJECT_NAME}-libsys
  ${PROJECT_NAME}-libncps
  ${PROJECT_NAME}-libpsop
  ${PROJECT_NAME}-libunits
  libfdf::libfdf
  libpsml::libpsml
  libgridxc::libgridxc
  LAPACK::LAPACK
 )

install(
  TARGETS gen-basis
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

