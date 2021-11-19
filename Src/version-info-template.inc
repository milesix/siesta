! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module version_info

implicit none

character(len=*), parameter :: version_str =  "@PROJECT_VERSION@"
character(len=*), parameter :: compiler_version = "@CMAKE_Fortran_COMPILER@"

character(len=*), parameter :: siesta_arch= &
"SIESTA_ARCH"
character(len=*), parameter :: fflags= &
"FFLAGS"
character(len=*), parameter :: fppflags= &
"FPPFLAGS"
character(len=*), parameter :: libs= &
"LIBS"


private

public :: version_str
public :: siesta_arch, fflags, fppflags, libs
public :: compiler_version

end module version_info
!================================================================

subroutine prversion

! Simple routine to print the version string. Could be extended to
! provide more information, if needed.

! Use free format in file to make more room for long option strings...

use version_info
implicit none

write(6,'(2a)') 'Siesta Version  : ', trim(version_str)
write(6,'(2a)') 'Architecture    : ', trim(siesta_arch)
write(6,'(2a)') 'Compiler version: ', trim(compiler_version)
write(6,'(2a)') 'Compiler flags  : ', trim(fflags)
write(6,'(2a)') 'PP flags        : ', trim(fppflags)
write(6,'(2a)') 'Libraries       : ', trim(libs)

#ifdef MPI
write(6,'(a)') 'PARALLEL version'
#else
write(6,'(a)') 'SERIAL version'
#endif

!$OMP parallel
!$OMP master
!$ write(*,'(a)') 'THREADED version'
#ifdef _OPENMP
!$ write(*,'(a,i0)') '* OpenMP version ', _OPENMP
#endif
!$OMP end master
!$OMP end parallel

#ifdef USE_GEMM3M
write(6,'(a)') 'GEMM3M support'
#endif
#ifdef CDF
write(6,'(a)') 'NetCDF support'
#endif
#ifdef NCDF_4
write(6,'(a)') 'NetCDF-4 support'
#ifdef NCDF_PARALLEL
write(6,'(a)') 'NetCDF-4 MPI-IO support'
#endif
#endif
#if defined(ON_DOMAIN_DECOMP) || defined(SIESTA__METIS)
write(6,'(a)') 'METIS ordering support'
#endif
#ifdef SIESTA__FLOOK
write(6,'(a)') 'Lua support'
#endif
#ifdef TRANSIESTA
write(6,'(a)') '******************************************************'
write(6,'(a)') 'transiesta executable is deprecated, please use siesta'
write(6,'(a)') '******************************************************'
#endif

end subroutine prversion
!----------------------------------------------------------
