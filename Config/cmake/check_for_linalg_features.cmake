include(CheckFortranSourceCompiles)

set(CMAKE_REQUIRED_LIBRARIES LAPACK::LAPACK)
check_fortran_source_compiles("external :: zheevr, dsyevr; call zheevr(); call dsyevr(); end"
                               LAPACK_HAS_MRRR SRC_EXT F90)

check_fortran_source_compiles("external :: dsyevr_2stage, zheevr_2stage; call zheevr_2stage(); call dsyevr_2stage(); end"
                               LAPACK_HAS_MRRR_2STAGE SRC_EXT F90)

check_fortran_source_compiles("external :: zheevx_2stage; call zheevx_2stage(); end"
                               LAPACK_HAS_ZVX_2STAGE SRC_EXT F90)
check_fortran_source_compiles("external :: zheevd_2stage; call zheevd_2stage(); end"
                               LAPACK_HAS_ZVD_2STAGE SRC_EXT F90)
check_fortran_source_compiles("external :: zheev_2stage; call zheev_2stage(); end"
                               LAPACK_HAS_ZV_2STAGE SRC_EXT F90)

check_fortran_source_compiles("external :: dsyevx_2stage; call dsyevx_2stage(); end"
                               LAPACK_HAS_DVX_2STAGE SRC_EXT F90)
check_fortran_source_compiles("external :: dsyevd_2stage; call dsyevd_2stage(); end"
                               LAPACK_HAS_DVD_2STAGE SRC_EXT F90)
check_fortran_source_compiles("external :: dsyev_2stage; call dsyev_2stage(); end"
                               LAPACK_HAS_DV_2STAGE SRC_EXT F90)

unset(CMAKE_REQUIRED_LIBRARIES)

if(WITH_MPI)

  set(CMAKE_REQUIRED_LIBRARIES Scalapack::Scalapack LAPACK::LAPACK)
  check_fortran_source_compiles("external :: pdsyevr, pzheevr; call pdsyevr(); call pzheevr(); end"
                               SCALAPACK_HAS_MRRR SRC_EXT F90)
  unset(CMAKE_REQUIRED_LIBRARIES)

  set(HAS_MRRR ${LAPACK_HAS_MRRR} AND ${SCALAPACK_HAS_MRRR} CACHE BOOL "MRRR routines available")

else()

  set(HAS_MRRR ${LAPACK_HAS_MRRR} CACHE BOOL "MRRR routines available")

endif(WITH_MPI)

