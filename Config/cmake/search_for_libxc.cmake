# -- Search for libxc using pkg-config data
#
 find_package(PkgConfig REQUIRED)

 
 pkg_check_modules(LIBXC_F90 REQUIRED libxcf90)

 message(DEBUG "Libxc_f90 Libdir: ${LIBXC_F90_LIBDIR}")
 message(DEBUG "Libxc_f90 Libraries: ${LIBXC_F90_LIBRARIES}")
 message(DEBUG "Libxc_f90 Link libraries: ${LIBXC_F90_LINK_LIBRARIES}")
 message(DEBUG "LIBXC_F90_INCLUDEDIR: ${LIBXC_F90_INCLUDEDIR}")
 message(DEBUG "LIBXC_F90_INCLUDE_DIRS: ${LIBXC_F90_INCLUDE_DIRS}")

        add_library(libxc::xcf90 INTERFACE IMPORTED)
        target_link_libraries(
          libxc::xcf90
          INTERFACE
          "${LIBXC_F90_LINK_LIBRARIES}"
        )
        target_include_directories(
          libxc::xcf90
          INTERFACE
          "${LIBXC_F90_INCLUDE_DIRS}"
        )

 pkg_check_modules(LIBXC_C REQUIRED libxc)

 message(DEBUG "Libxc_c Libdir: ${LIBXC_C_LIBDIR}")
 message(DEBUG "Libxc_c Libraries: ${LIBXC_C_LIBRARIES}")
 message(DEBUG "Libxc_c Link libraries: ${LIBXC_C_LINK_LIBRARIES}")
 message(DEBUG "LIBXC_C_INCLUDEDIR: ${LIBXC_C_INCLUDEDIR}")
 message(DEBUG "LIBXC_C_INCLUDE_DIRS: ${LIBXC_C_INCLUDE_DIRS}")

        add_library(libxc::xc INTERFACE IMPORTED)
        target_link_libraries(
          libxc::xc
          INTERFACE
          "${LIBXC_C_LINK_LIBRARIES}"
        )
        target_include_directories(
          libxc::xc
          INTERFACE
          "${LIBXC_C_INCLUDE_DIRS}"
        )

pkg_check_modules(LIBXC_F03 libxcf03)

 if(LIBXC_F03_FOUND)

        add_library(libxc::xcf03 INTERFACE IMPORTED)
        target_link_libraries(
          libxc::xcf03
          INTERFACE
          "${LIBXC_F03_LINK_LIBRARIES}"
        )
        target_include_directories(
          libxc::xcf03
          INTERFACE
          "${LIBXC_F03_INCLUDE_DIRS}"
        )

   message(DEBUG "Libxc_f03 Libdir: ${LIBXC_F03_LIBDIR}")
   message(DEBUG "Libxc_f03 Libraries: ${LIBXC_F03_LIBRARIES}")
   message(DEBUG "Libxc_f03 Link libraries: ${LIBXC_F03_LINK_LIBRARIES}")
   message(DEBUG "LIBXC_F03_INCLUDEDIR: ${LIBXC_F03_INCLUDEDIR}")
   message(DEBUG "LIBXC_F03_INCLUDE_DIRS: ${LIBXC_F03_INCLUDE_DIRS}")
 
   message(STATUS "   ---> Found F2003 interface for libxc. Using it as default")
   add_library(libxc::XC_Fortran INTERFACE IMPORTED)
   target_link_libraries(libxc::XC_Fortran INTERFACE libxc::xcf03 libxc::xc)

 else()
 
  message(STATUS "   ---> Did not find F2003 interface for libxc")
  message(STATUS "   ---> Found F90 interface for libxc. Using it...")

  add_library(libxc::XC_Fortran INTERFACE IMPORTED)
  target_link_libraries(libxc::XC_Fortran INTERFACE libxc::xcf90 libxc::xc)
  
 endif()
