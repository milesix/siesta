set(top_srcdir "${CMAKE_SOURCE_DIR}/Src")
set(grid_srcdir "${CMAKE_SOURCE_DIR}/Util/Grid")

set(sources

  local_sys.f
  fftw3_mymod.f90
  extrapolate.f
  sts.f
  mainstm.F
  neighb.f
  ranger.f
  readstm.f
  redata.f
  reinitstm.f
  stm.f
  vacpot.f
  handlers.f
  
)

list(
    APPEND
    sources

  ${grid_srcdir}/m_gridfunc.F90

 )

list(
    APPEND
    sources

  ${top_srcdir}/precision.F
  ${top_srcdir}/bessph.f
  ${top_srcdir}/chkdim.f
  ${top_srcdir}/dismin.f
  ${top_srcdir}/dot.f
  ${top_srcdir}/m_io.f
  ${top_srcdir}/io.f
  ${top_srcdir}/alloc.F90
  ${top_srcdir}/files.f
  ${top_srcdir}/parallel.F
  ${top_srcdir}/printmatrix.F
  ${top_srcdir}/schecomm.F
  ${top_srcdir}/class_Geometry.F90
  ${top_srcdir}/class_Data2D.F90
  ${top_srcdir}/class_OrbitalDistribution.F90
  ${top_srcdir}/class_Sparsity.F90
  ${top_srcdir}/class_SpData2D.F90
  ${top_srcdir}/class_Pair_Geometry_SpData2D.F90
  ${top_srcdir}/class_Fstack_Pair_Geometry_SpData2D.F90
  ${top_srcdir}/class_Data1D.F90
  ${top_srcdir}/class_SpData1D.F90
  ${top_srcdir}/sparse_matrices.F90
  ${top_srcdir}/domain_decom.F
  ${top_srcdir}/spatial.F
  ${top_srcdir}/parallelsubs.F
  ${top_srcdir}/iodm.F
  ${top_srcdir}/memory_log.F90
  ${top_srcdir}/memory.F
  ${top_srcdir}/pxf.F90
  ${top_srcdir}/m_fft_gpfa.F
  ${top_srcdir}/radfft.f
  ${top_srcdir}/volcel.f
  ${top_srcdir}/listsc.f
  ${top_srcdir}/atmparams.f
  ${top_srcdir}/atom_options.F90
  ${top_srcdir}/interpolation.f90
  ${top_srcdir}/xml.f
  ${top_srcdir}/radial.f
  ${top_srcdir}/atm_types.f
  ${top_srcdir}/spher_harm.f
  ${top_srcdir}/atmfuncs.f
  ${top_srcdir}/basis_types.f
  ${top_srcdir}/chemical.f
  ${top_srcdir}/basis_io.F
  ${top_srcdir}/qsort.F
  ${top_srcdir}/mmio.F
  ${top_srcdir}/pspltm1.F
  ${top_srcdir}/units.f90
  ${top_srcdir}/m_uuid.f90
  ${top_srcdir}/object_debug.F90
  ${top_srcdir}/class_Pair_Data1D.F90
  ${top_srcdir}/class_Fstack_Pair_Data1D.F90

)

add_executable(
   ol-stm
   ${sources}
)

install(
  TARGETS ol-stm
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

 target_link_libraries(
  ol-stm
  PRIVATE
  fdf
  libpsml::libpsml
  ${PROJECT_NAME}-libncps
  FFTW::Double
 )

