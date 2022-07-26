# -- Search for FFTW using pkg-config data
#
 find_package(PkgConfig QUIET)
 
 pkg_check_modules(FFTW_FORTRAN IMPORTED_TARGET GLOBAL fftw3f>=3.0)

 message(STATUS "FFTW_Fortran Libdir: ${FFTW_FORTRAN_LIBDIR}")
 message(STATUS "FFTW_Fortran Libraries: ${FFTW_FORTRAN_LIBRARIES}")
 message(STATUS "FFTW_Fortran Link libraries: ${FFTW_FORTRAN_LINK_LIBRARIES}")
 message(STATUS "FFTW_FORTRAN_INCLUDEDIR: ${FFTW_FORTRAN_INCLUDEDIR}")
 message(STATUS "FFTW_FORTRAN_INCLUDE_DIRS: ${FFTW_FORTRAN_INCLUDE_DIRS}")

 pkg_check_modules(FFTW_C IMPORTED_TARGET GLOBAL fftw3>=3.0)
 message(STATUS "fftw_c Libdir: ${FFTW_C_LIBDIR}")
 message(STATUS "fftw_c Libraries: ${FFTW_C_LIBRARIES}")
 message(STATUS "fftw_c Link libraries: ${FFTW_C_LINK_LIBRARIES}")
 message(STATUS "FFTW_C_INCLUDEDIR: ${FFTW_C_INCLUDEDIR}")
 message(STATUS "FFTW_C_INCLUDE_DIRS: ${FFTW_C_INCLUDE_DIRS}")

set(FFTW_LINK_LIBRARIES
     ${FFTW_FORTRAN_LINK_LIBRARIES}
     ${FFTW_C_LINK_LIBRARIES}
 )
 
 set(FFTW_INCLUDE_DIRS
     ${FFTW_FORTRAN_INCLUDE_DIRS}
     ${FFTW_C_INCLUDE_DIRS}
 )


 