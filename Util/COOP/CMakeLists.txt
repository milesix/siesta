set(top_srcdir "${CMAKE_SOURCE_DIR}/Src")

set (sources_top

        ${top_srcdir}/precision.F
        ${top_srcdir}/alloc.F90
        ${top_srcdir}/parallel.F
        ${top_srcdir}/sys.F
        ${top_srcdir}/m_io.f
        ${top_srcdir}/pxf.F90
        ${top_srcdir}/moreParallelSubs.F90
        ${top_srcdir}/m_getopts.f90
)

set(sources_fat
	${sources_top}
        io.f
        subs.f90
        units.f90
        main_vars.f90
        io_hs.f90
        orbital_set.f90
        read_curves.f90
        handlers.f
        fat.f90
)
set(sources_mprop
	${sources_top}
        io.f
        subs.f90
        units.f90
        main_vars.f90
        io_hs.f90
        orbital_set.f90
        read_curves.f90
        handlers.f
        mprop.f90
)
set(sources_spin_texture
	${sources_top}

        io.f
        subs.f90
        units.f90
        main_vars.f90
        io_hs.f90
        handlers.f
        spin_texture.f90
)
set(sources_dm_creator
	${sources_top}

        io.f
        subs.f90
        units.f90
        main_vars.f90
        io_hs.f90
        iodm_netcdf.F90
        write_dm.f
        handlers.f
        dm_creator.F90
)

add_executable(
   fat
   ${sources_fat}
)
add_executable(
   mprop
   ${sources_mprop}
)
add_executable(
   spin_texture
   ${sources_spin_texture}
)
add_executable(
   dm_creator
   ${sources_dm_creator}
)

if(WITH_NETCDF)
 target_link_libraries(dm_creator
   PRIVATE
   NetCDF::NetCDF_Fortran
 )
  target_compile_definitions(
    dm_creator
    PRIVATE
    CDF
  )
endif()

install(
  TARGETS fat mprop spin_texture dm_creator
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
