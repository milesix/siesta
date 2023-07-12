set(top_srcdir "${PROJECT_SOURCE_DIR}/Src")

set(sources

   m_ts_io_version.f90
   tshs2tshs.f90
)

list(
    APPEND
    sources

   ${top_srcdir}/precision.F
   ${top_srcdir}/parallel.F
   ${top_srcdir}/alloc.F90
   ${top_srcdir}/units.f90
   ${top_srcdir}/m_io.f
   ${top_srcdir}/memory_log.F90
   ${top_srcdir}/io.f
   ${top_srcdir}/class_OrbitalDistribution.F90
   ${top_srcdir}/class_Data2D.F90
   ${top_srcdir}/class_Sparsity.F90
   ${top_srcdir}/class_SpData2D.F90
   ${top_srcdir}/intrinsic_missing.F90
   ${top_srcdir}/geom_helper.f90
   ${top_srcdir}/m_sparse.F90
   ${top_srcdir}/reclat.f
   ${top_srcdir}/class_Data1D.F90
   ${top_srcdir}/class_SpData1D.F90
   ${top_srcdir}/io_sparse.F90
   ${top_srcdir}/ncdf_io.F90
   ${top_srcdir}/m_os.F90
   ${top_srcdir}/m_ts_io.F90
   ${top_srcdir}/m_uuid.f90
   ${top_srcdir}/object_debug.F90
   ${top_srcdir}/nag.f
)

add_executable(
   tshs2tshs
   ${sources}
)

install(
  TARGETS tshs2tshs
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

#
 target_link_libraries(
  tshs2tshs
  PRIVATE
  ${PROJECT_NAME}-libsys
  $<$<BOOL:${WITH_NETCDF}>:NetCDF::NetCDF_Fortran>
  "$<$<BOOL:${WITH_NCDF}>:${PROJECT_NAME}-libncdf>"
  "$<$<BOOL:${WITH_NCDF}>:${PROJECT_NAME}-libfdict>"
 )

target_compile_definitions(
  tshs2tshs
  PRIVATE
  "$<$<BOOL:${WITH_NETCDF}>:CDF>"
  "$<$<BOOL:${WITH_NCDF}>:NCDF>"
  "$<$<BOOL:${WITH_NCDF}>:NCDF_4>"
  "$<$<BOOL:${WITH_NCDF_PARALLEL}>:NCDF_PARALLEL>"
)

