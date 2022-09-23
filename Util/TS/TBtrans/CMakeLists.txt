set(top_srcdir "${CMAKE_SOURCE_DIR}/Src" )

set(sources

 handlers.F
 m_tbt_contour.F90
 m_tbt_dH.F90
 m_tbt_dSE.F90
 m_tbt_delta.F90
 m_tbt_diag.F90
 m_tbt_gf.F90
 m_tbt_hs.F90
 m_tbt_kpoint.F90
 m_tbt_kregions.F90
 m_tbt_options.F90
 m_tbt_proj.F90
 m_tbt_regions.F90
 m_tbt_save.F90
 m_tbt_sigma_save.F90
 m_tbt_sparse_helper.F90
 m_tbt_tri_init.F90
 m_tbt_tri_scat.F90
 m_tbt_trik.F90
 m_tbtrans.F90
 m_verbosity.f90
 simple-version.F90
 sorted_search.F90
 tbt_end.F90
 tbt_init.F90
 tbt_init_output.f90
 tbt_reinit_m.F90
 tbtrans.F90

)
 
list(
  APPEND
  sources

  "${top_srcdir}/alloc.F90"
  "${top_srcdir}/basic_func.inc"
  "${top_srcdir}/basic_type.inc"
  "${top_srcdir}/bloch_unfold.F90"
   ${top_srcdir}/broadcast_fdf_struct.F90
  "${top_srcdir}/byte_count.F90"
  "${top_srcdir}/cellsubs.f"
  "${top_srcdir}/class_Data1D.F90"
  "${top_srcdir}/class_Data1D.T90"
  "${top_srcdir}/class_Data2D.F90"
  "${top_srcdir}/class_Data2D.T90"
  "${top_srcdir}/class_Distribution.F90"
  "${top_srcdir}/class_OrbitalDistribution.F90"
  "${top_srcdir}/class_SpData1D.F90"
  "${top_srcdir}/class_SpData1D.T90"
  "${top_srcdir}/class_SpData2D.F90"
  "${top_srcdir}/class_SpData2D.T90"
  "${top_srcdir}/class_Sparsity.F90"
  "${top_srcdir}/class_TriMat.F90"
  "${top_srcdir}/class_TriMat.T90"
  "${top_srcdir}/cli_m.f90"
  "${top_srcdir}/create_Sparsity_SC.F90"
  "${top_srcdir}/create_Sparsity_Union.F90"
  "${top_srcdir}/cross.f"
  "${top_srcdir}/debugmpi.F"
  "${top_srcdir}/densematrix.f90"
  "${top_srcdir}/diag.F90"
  "${top_srcdir}/diag_option.F90"
  "${top_srcdir}/fdf_extra.F90"
  "${top_srcdir}/files.f"
  "${top_srcdir}/find_kgrid.F"
  "${top_srcdir}/geom_helper.f90"
  "${top_srcdir}/get_kpoints_scale.f90"
  "${top_srcdir}/idiag.f"
  "${top_srcdir}/intrinsic_missing.F90"
  "${top_srcdir}/io.f"
  "${top_srcdir}/kpoint_convert.f90"
  "${top_srcdir}/m_char.f90"
  "${top_srcdir}/m_cite.F90"
  "${top_srcdir}/m_gauss_fermi_17.f90"
  "${top_srcdir}/m_gauss_fermi_18.f90"
  "${top_srcdir}/m_gauss_fermi_19.f90"
  "${top_srcdir}/m_gauss_fermi_20.f90"
  "${top_srcdir}/m_gauss_fermi_22.f90"
  "${top_srcdir}/m_gauss_fermi_24.f90"
  "${top_srcdir}/m_gauss_fermi_26.f90"
  "${top_srcdir}/m_gauss_fermi_28.f90"
  "${top_srcdir}/m_gauss_fermi_30.f90"
  "${top_srcdir}/m_gauss_fermi_inf.f90"
  "${top_srcdir}/m_gauss_quad.f90"
  "${top_srcdir}/m_geom_aux.f90"
  "${top_srcdir}/m_geom_box.f90"
  "${top_srcdir}/m_geom_coord.f90"
  "${top_srcdir}/m_geom_objects.f90"
  "${top_srcdir}/m_geom_plane.f90"
  "${top_srcdir}/m_geom_square.f90"
  "${top_srcdir}/m_handle_sparse.F90"
  "${top_srcdir}/m_integrate.f90"
  "${top_srcdir}/m_interpolate.F90"
  "${top_srcdir}/m_io.f"
  "${top_srcdir}/io_sparse.F90"
  "${top_srcdir}/m_iodm.F90"
  "${top_srcdir}/m_iterator.f90"
  "${top_srcdir}/m_mat_invert.F90"
  "${top_srcdir}/m_mpi_utils.F"
  "${top_srcdir}/m_os.F90"
  "${top_srcdir}/m_pivot.F90"
  "${top_srcdir}/m_pivot_array.f90"
  "${top_srcdir}/m_pivot_methods.F90"
  "${top_srcdir}/m_region.F90"
  "${top_srcdir}/m_sparse.F90"
  "${top_srcdir}/m_sparsity_handling.F90"
  "${top_srcdir}/m_spin.F90"
  "${top_srcdir}/m_timer.F90"
  "${top_srcdir}/m_trimat_invert.F90"
  "${top_srcdir}/m_ts_aux.F90"
  "${top_srcdir}/m_ts_cctype.f90"
  "${top_srcdir}/m_ts_chem_pot.F90"
  "${top_srcdir}/m_ts_contour.f90"
  "${top_srcdir}/m_ts_contour_eq.f90"
  "${top_srcdir}/m_ts_contour_neq.f90"
  "${top_srcdir}/m_ts_debug.F90"
  "${top_srcdir}/m_ts_elec_se.F90"
  "${top_srcdir}/m_ts_electrode.F90"
  "${top_srcdir}/m_ts_method.f90"
  "${top_srcdir}/m_ts_gf.F90"
  "${top_srcdir}/m_ts_io.F90"
  "${top_srcdir}/m_ts_io_contour.f90"
  "${top_srcdir}/m_ts_io_ctype.f90"
  "${top_srcdir}/m_ts_iodm.F90"
  "${top_srcdir}/m_ts_pivot.F90"
  "${top_srcdir}/m_ts_rgn2trimat.F90"
  "${top_srcdir}/m_ts_sparse.F90"
  "${top_srcdir}/m_ts_sparse_helper.F90"
  "${top_srcdir}/m_ts_tdir.f90"
  "${top_srcdir}/m_ts_tri_common.F90"
  "${top_srcdir}/m_ts_tri_scat.F90"
  "${top_srcdir}/m_ts_trimat_invert.F90"
  "${top_srcdir}/m_uuid.f90"
  "${top_srcdir}/m_wallclock.f90"
  "${top_srcdir}/m_walltime.f90"
  "${top_srcdir}/memory.F"
  "${top_srcdir}/memory_log.F90"
  "${top_srcdir}/minvec.f"
  "${top_srcdir}/moreParallelSubs.F90"
  "${top_srcdir}/ncdf_io.F90"
  "${top_srcdir}/object_debug.F90"
  "${top_srcdir}/parallel.F"
  "${top_srcdir}/posix_calls.f90"
  "${top_srcdir}/precision.F"
  "${top_srcdir}/pxf.F90"
  "${top_srcdir}/reclat.f"
  "${top_srcdir}/runinfo_m.F90"
  "${top_srcdir}/siesta_geom.F90"
  "${top_srcdir}/sorting.f"
  "${top_srcdir}/sys.F"
  "${top_srcdir}/spin_subs.F90"
  "${top_srcdir}/timer.F90"
  "${top_srcdir}/timer_tree.f90"
  "${top_srcdir}/timestamp.f90"
  "${top_srcdir}/ts_electrode.F90"
  "${top_srcdir}/units.f90"
)

add_executable(tbtrans "${sources}")

# -- Version information

# 
configure_file(
  version-info-template.inc
  ${CMAKE_CURRENT_BINARY_DIR}/generated/version-info.inc
  @ONLY
  )

target_include_directories(
  tbtrans
  PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/generated
  )


# Add -DTBTRANS
target_compile_definitions(
  tbtrans
  PRIVATE
  TBTRANS
)

if(WITH_MPI)
  target_compile_definitions(
  tbtrans
  PRIVATE
  MPI
  )
# MPI interfaces
  target_link_libraries(
    tbtrans
    PRIVATE
    mpi_siesta
  )
endif()

TARGET_LINK_LIBRARIES(
    tbtrans
    PRIVATE
    fdf
)

# External libraries.

# Due to some design decisions in the code,
# TBtrans depends on FDICT even if not
# using NCDF

    target_link_libraries(
      tbtrans
      PRIVATE
      ${PROJECT_NAME}-libfdict
      )

if(WITH_MPI)
  target_link_libraries(
    tbtrans
    PRIVATE
    Scalapack::Scalapack
  )
endif()
  
TARGET_LINK_LIBRARIES(
  tbtrans
  PRIVATE
  LAPACK::LAPACK)

if(WITH_NETCDF)

  TARGET_LINK_LIBRARIES(
    tbtrans
    PRIVATE
    NetCDF::NetCDF_Fortran
  )
  target_compile_definitions(
    tbtrans
    PRIVATE
    CDF
  )

  if(WITH_NCDF)
    target_compile_definitions(
    tbtrans
    PRIVATE 
    NCDF
    NCDF_4
    )
#
    if(WITH_NCDF_PARALLEL)
      target_compile_definitions(
      tbtrans
      PRIVATE
      NCDF_PARALLEL
      )
    endif()

    target_link_libraries(
      tbtrans
      PRIVATE
      ${PROJECT_NAME}-libncdf
      )
    target_link_libraries(
      tbtrans
      PRIVATE
      ${PROJECT_NAME}-libfdict
      )
  endif()
endif()

install(
  TARGETS tbtrans
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

     
