#
# There are several targets in this directory. The 'OBJECT library' idiom
# was needed to avoid ninja complaining about "two paths for the same
# .mod file"
#
set(top_src_dir "${CMAKE_SOURCE_DIR}/Src")

# To avoid "multiple rules" in Ninja
# Maybe turn into global-utility targets

add_library(aux_opts     ${top_src_dir}/m_getopts.f90)
add_library(aux_gridfunc  m_gridfunc.F90)

#-----------------------------------------
add_executable(g2c_ng
   ${top_src_dir}/reclat.f
   m_struct.f90
   local_die.f90
   g2c_ng.f
)

# For future enabling of NetCDF in m_gridfunc
target_link_libraries(g2c_ng PRIVATE $<$<BOOL:${WITH_NETCDF}>:NetCDF::NetCDF_Fortran>)
target_link_libraries(g2c_ng PRIVATE aux_opts aux_gridfunc)
target_compile_definitions(g2c_ng  PRIVATE $<$<BOOL:${WITH_NETCDF}>:CDF> )
		     
install(TARGETS g2c_ng RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
#---------------------------------------------------------------------

add_executable(v_info v_info.f90)

# For future enabling of NetCDF in m_gridfunc
target_link_libraries(v_info PRIVATE $<$<BOOL:${WITH_NETCDF}>:NetCDF::NetCDF_Fortran>)
target_link_libraries(v_info PRIVATE aux_gridfunc)
target_compile_definitions(v_info  PRIVATE $<$<BOOL:${WITH_NETCDF}>:CDF> )
		     
install(TARGETS v_info RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )

#---------------------------------------------------------------------
add_executable(grid2cube grid2cube.f)
add_executable(grid2val grid2val.F)
add_executable(grid_rotate grid_rotate.F90)
add_executable(grid_supercell grid_supercell.F90)
install(TARGETS grid2cube grid2val RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
install(TARGETS grid_rotate grid_supercell RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
#---------------------------------------------------------------------

if (WITH_NETCDF)
#-----------------------------------------
add_executable(cdf2grid  cdf2grid.F90)

target_link_libraries(cdf2grid PRIVATE NetCDF::NetCDF_Fortran )
target_compile_definitions(cdf2grid  PRIVATE CDF )
		     
install(TARGETS cdf2grid RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
#---------------------------------------------------------------------
add_executable(grid2cdf  grid2cdf.F90)

target_link_libraries(grid2cdf PRIVATE NetCDF::NetCDF_Fortran )
target_compile_definitions(grid2cdf  PRIVATE CDF )
		     
install(TARGETS grid2cdf RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
#-----------------------------------------
add_executable(cdf2xsf  cdf2xsf.F90)

target_link_libraries(cdf2xsf PRIVATE NetCDF::NetCDF_Fortran )
target_compile_definitions(cdf2xsf  PRIVATE CDF )
		     
install(TARGETS cdf2xsf RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
#-----------------------------------------

# To avoid "multiple rules" in Ninja
# (could add OBJECT qualifier to avoid the lib appearing)
add_library(aux_m_grid  m_grid.F90)
target_compile_definitions(aux_m_grid  PRIVATE CDF )
target_link_libraries(aux_m_grid PUBLIC NetCDF::NetCDF_Fortran )



add_executable(cdf_diff cdf_diff.F90)

target_link_libraries(cdf_diff PRIVATE aux_m_grid NetCDF::NetCDF_Fortran )
target_compile_definitions(cdf_diff  PRIVATE CDF )
		     
install(TARGETS cdf_diff RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
#-----------------------------------------
add_executable(cdf_get_cell cdf_get_cell.F90)

target_link_libraries(cdf_get_cell PRIVATE aux_m_grid NetCDF::NetCDF_Fortran )
target_compile_definitions(cdf_get_cell  PRIVATE CDF )
		     
install(TARGETS cdf_get_cell RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
#-----------------------------------------

# To avoid "multiple rules" in Ninja
add_library(aux_fft  ${top_src_dir}/m_fft_gpfa.F)

add_executable(cdf_laplacian
   ${top_src_dir}/reclat.f
   fft3d.f
   local_die.f90
   cdf_laplacian.F90
)

target_link_libraries(cdf_laplacian PRIVATE aux_fft aux_opts)
target_link_libraries(cdf_laplacian PRIVATE NetCDF::NetCDF_Fortran )
target_compile_definitions(cdf_laplacian  PRIVATE CDF )
		     
install(TARGETS cdf_laplacian RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )

#-----------------------------------------

add_executable(cdf_fft
   ${top_src_dir}/reclat.f
   fft3d.f
   local_die.f90
   cdf_fft.F90
)

target_link_libraries(cdf_fft PRIVATE aux_fft)
target_link_libraries(cdf_fft PRIVATE NetCDF::NetCDF_Fortran )
target_compile_definitions(cdf_fft  PRIVATE CDF )
		     
install(TARGETS cdf_fft RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )

endif(WITH_NETCDF)
