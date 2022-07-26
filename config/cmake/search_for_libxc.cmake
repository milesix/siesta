# -- Search for libxc using pkg-config data
#
 find_package(PkgConfig REQUIRED)
 
 pkg_check_modules(LIBXC_F90 REQUIRED libxcf90)

 message(STATUS "Libxc_f90 Libdir: ${LIBXC_F90_LIBDIR}")
 message(STATUS "Libxc_f90 Libraries: ${LIBXC_F90_LIBRARIES}")
 message(STATUS "Libxc_f90 Link libraries: ${LIBXC_F90_LINK_LIBRARIES}")
 message(STATUS "LIBXC_F90_INCLUDEDIR: ${LIBXC_F90_INCLUDEDIR}")
 message(STATUS "LIBXC_F90_INCLUDE_DIRS: ${LIBXC_F90_INCLUDE_DIRS}")

 pkg_check_modules(LIBXC_C REQUIRED libxc)
 message(STATUS "Libxc_c Libdir: ${LIBXC_C_LIBDIR}")
 message(STATUS "Libxc_c Libraries: ${LIBXC_C_LIBRARIES}")
 message(STATUS "Libxc_c Link libraries: ${LIBXC_C_LINK_LIBRARIES}")
 message(STATUS "LIBXC_C_INCLUDEDIR: ${LIBXC_C_INCLUDEDIR}")
 message(STATUS "LIBXC_C_INCLUDE_DIRS: ${LIBXC_C_INCLUDE_DIRS}")

 pkg_check_modules(LIBXC_F03 libxcf03)

 if(LIBXC_F03_FOUND)
   message(STATUS "Found F2003 interface for libxc. Not used yet")
   message(STATUS "Libxc_f03 Libdir: ${LIBXC_F03_LIBDIR}")
   message(STATUS "Libxc_f03 Libraries: ${LIBXC_F03_LIBRARIES}")
   message(STATUS "Libxc_f03 Link libraries: ${LIBXC_F03_LINK_LIBRARIES}")
   message(STATUS "LIBXC_F03_INCLUDEDIR: ${LIBXC_F03_INCLUDEDIR}")
   message(STATUS "LIBXC_F03_INCLUDE_DIRS: ${LIBXC_F03_INCLUDE_DIRS}")
 endif()
 
 set(LIBXC_LINK_LIBRARIES
     ${LIBXC_F90_LINK_LIBRARIES}
     ${LIBXC_C_LINK_LIBRARIES}
 )
 
 set(LIBXC_INCLUDE_DIRS
     ${LIBXC_F90_INCLUDE_DIRS}
     ${LIBXC_C_INCLUDE_DIRS}
 )
  