# -- Search for netCDF using pkg-config data
# -- There should be a fallback
#
 find_package(PkgConfig REQUIRED)
 pkg_check_modules(NETCDF REQUIRED netcdf-fortran)

 message(STATUS "Netcdf Libdir: ${NETCDF_LIBDIR}")
 message(STATUS "Netcdf Libraries: ${NETCDF_LIBRARIES}")
 message(STATUS "Netcdf Link libraries: ${NETCDF_LINK_LIBRARIES}")
 message(STATUS "NETCDF_INCLUDEDIR: ${NETCDF_INCLUDEDIR}")
 message(STATUS "NETCDF_INCLUDE_DIRS: ${NETCDF_INCLUDE_DIRS}")

