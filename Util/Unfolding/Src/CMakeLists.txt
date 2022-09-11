set(top_srcdir "${CMAKE_SOURCE_DIR}/Src" )

set(sources

  array.F90
  hsx_m.f90
  unfold.F90
  handlers.f
)

list(
  APPEND
  sources

  ${top_srcdir}/precision.F
  ${top_srcdir}/broadcast_fdf_struct.F90
  ${top_srcdir}/sys.F
  ${top_srcdir}/m_io.f
  ${top_srcdir}/alloc.F90
  ${top_srcdir}/parallel.F
  ${top_srcdir}/memory_log.F90
  ${top_srcdir}/memory.F
  ${top_srcdir}/pxf.F90
  ${top_srcdir}/moreParallelSubs.F90
  ${top_srcdir}/m_mpi_utils.F
  ${top_srcdir}/interpolation.f90
  ${top_srcdir}/xml.f
  ${top_srcdir}/radial.f
  ${top_srcdir}/atm_types.f
  ${top_srcdir}/spher_harm.f
  ${top_srcdir}/atmfuncs.f
  ${top_srcdir}/reclat.f
  ${top_srcdir}/idiag.f
  ${top_srcdir}/cellsubs.f
  ${top_srcdir}/sorting.f
  ${top_srcdir}/minvec.f
  ${top_srcdir}/io.f
  ${top_srcdir}/files.f
  ${top_srcdir}/timestamp.f90
  ${top_srcdir}/m_walltime.f90
  ${top_srcdir}/m_wallclock.f90
  ${top_srcdir}/diag_option.F90
  ${top_srcdir}/diag.F90
  ${top_srcdir}/units.f90
  ${top_srcdir}/bloch_unfold.F90
  ${top_srcdir}/debugmpi.F
  ${top_srcdir}/find_kgrid.F
  ${top_srcdir}/posix_calls.f90
  ${top_srcdir}/atmparams.f
  ${top_srcdir}/basis_types.f
  ${top_srcdir}/atom_options.F90
  ${top_srcdir}/chemical.f
  ${top_srcdir}/basis_io.F
  ${top_srcdir}/siesta_geom.F90
  ${top_srcdir}/get_kpoints_scale.f90
  ${top_srcdir}/bessph.f
  ${top_srcdir}/m_fft_gpfa.F
  ${top_srcdir}/radfft.f
  ${top_srcdir}/periodic_table.f
  ${top_srcdir}/basis_specs.f
  ${top_srcdir}/hamann.f90
  ${top_srcdir}/m_filter.f90
  ${top_srcdir}/atom.F
  ${top_srcdir}/m_cite.F90
  ${top_srcdir}/dftu_specs.f
  ${top_srcdir}/broadcast_basis.F
  ${top_srcdir}/m_cell.f
  ${top_srcdir}/zmatrix.F
  ${top_srcdir}/coor.F
  ${top_srcdir}/chkdim.f
  ${top_srcdir}/volcel.f
  ${top_srcdir}/digcel.f
  ${top_srcdir}/arw.f
  ${top_srcdir}/redcel.F
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
  ${top_srcdir}/pspltm1.F
  ${top_srcdir}/mmio.F
  ${top_srcdir}/qsort.F
  ${top_srcdir}/timer_tree.f90
  ${top_srcdir}/m_timer.F90
  ${top_srcdir}/timer.F90
  ${top_srcdir}/m_uuid.f90
  ${top_srcdir}/object_debug.F90
  ${top_srcdir}/nag.f
  ${top_srcdir}/m_spin.F90
  ${top_srcdir}/m_vee_integrals.F90 
)

add_executable(unfold "${sources}")

if(WITH_MPI)
  target_compile_definitions(
  unfold
  PUBLIC
  MPI
  )
# MPI interfaces
  target_link_libraries(
    unfold
    PRIVATE
    mpi_siesta
  )
endif()


TARGET_LINK_LIBRARIES(
    unfold
    PRIVATE
    fdf
)

# External libraries.


if(WITH_MPI)
  target_link_libraries(
    unfold
    PRIVATE
    "$ENV{SCALAPACK_LIBS}"
  )
endif()
  
TARGET_LINK_LIBRARIES(
  unfold
  PRIVATE
  "$ENV{LAPACK_LIBS}")

if(WITH_NETCDF)

  TARGET_LINK_LIBRARIES(
    unfold
    PRIVATE
    NetCDF::NetCDF_Fortran
)
  target_compile_definitions(
    unfold
    PUBLIC
    CDF
    NCDF
    NCDF_4
  )

endif()

 target_link_libraries(
  unfold
  PRIVATE
  libgridxc::libgridxc
  libpsml::libpsml
  ${PROJECT_NAME}-libncps
  ${PROJECT_NAME}-libpsop
 )

install(
  TARGETS unfold
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
