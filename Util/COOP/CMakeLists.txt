#
# There are several targets in this directory. The 'OBJECT library' idiom
# was needed to avoid ninja complaining about "two paths for the same
# .mod file". Common files are grouped in a pseudo-library object.
#
# This is a major annoyance.
#
set(top_srcdir "${PROJECT_SOURCE_DIR}/Src")

add_library(coop_common_objs OBJECT

        ${top_srcdir}/alloc.F90
        ${top_srcdir}/parallel.F
        ${top_srcdir}/m_io.f
        ${top_srcdir}/moreParallelSubs.F90
        ${top_srcdir}/m_getopts.f90
        ${top_srcdir}/precision.F
        io.f
        subs.f90
        units.f90
        main_vars.f90
        io_hs.f90
        orbital_set.f90
        read_curves.f90
	local_timer.f90
)

target_link_libraries(
   coop_common_objs
   PRIVATE
   ${PROJECT_NAME}-libsys
   )

set(sources_fat
        fat.f90
)
set(sources_mprop
        mprop.f90
)
set(sources_spin_texture
        spin_texture.f90
)
set(sources_dm_creator
        iodm_netcdf.F90
        write_dm.f
        dm_creator.F90
)

add_executable(
   fat
   ${sources_fat}
)
target_link_libraries(fat
   coop_common_objs
   ${PROJECT_NAME}-libsys
   )
   
add_executable(
   mprop
   ${sources_mprop}
)
target_link_libraries(mprop
   PRIVATE
   coop_common_objs
   ${PROJECT_NAME}-libsys
   )

add_executable(
   spin_texture
   ${sources_spin_texture}
)
target_link_libraries(spin_texture
   PRIVATE
   coop_common_objs
   ${PROJECT_NAME}-libsys
   )

add_executable(
   dm_creator
   ${sources_dm_creator}
)
target_link_libraries(dm_creator
   PRIVATE
   ${PROJECT_NAME}-libsys
   coop_common_objs
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
