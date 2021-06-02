!
! Copyright (C) 1996-2021       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module version_info

implicit none

character(len=*), parameter :: version_str =  &
"SIESTA_VERSION"
character(len=*), parameter :: compiler_version = &
"COMPILER_VERSION"
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

public :: prversion

!================================================================

CONTAINS

subroutine prversion

! Simple routine to print the version string. Could be extended to
! provide more information, if needed.

! Use free format in file to make more room for long option strings...

implicit none

logical :: has_parallel
!$ integer :: omp_version
!$ character(len=:), allocatable :: omp_name

write(6,'(2a)') 'Siesta Version  : ', trim(adjustl(version_str))
write(6,'(2a)') 'Architecture    : ', trim(adjustl(siesta_arch))
write(6,'(2a)') 'Compiler version: ', trim(adjustl(compiler_version))
write(6,'(2a)') 'Compiler flags  : ', trim(adjustl(fflags))
write(6,'(2a)') 'PP flags        : ', trim(adjustl(fppflags))
write(6,'(2a)') 'Libraries       : ', trim(adjustl(libs))

write(6,'(a)',ADVANCE='NO') 'Parallelisations: '
#ifdef MPI
has_parallel = .true.
write(6,'(a)',ADVANCE='NO') 'MPI'
#else
has_parallel = .false.
#endif

!$ if (has_parallel) write(6,'(a)', ADVANCE='NO') '; and '
!$ write(6,'(a)',ADVANCE='NO') 'OpenMP threads'
!$ has_parallel = .true.
if (.not. has_parallel) write(6,'(a)',ADVANCE='NO') 'none'
write(6,'(a)') '.'
#ifdef _OPENMP
!$ omp_version = _OPENMP
!$ select case (omp_version)
!$    case (202011)
!$       omp_name = 'OpenMP 5.1'
!$    case (201811)
!$       omp_name = 'OpenMP 5.0'
!$    case (201611)
!$       ! jme52: Many versions of ifort report this value
!$       ! (I don't think this should be a valid value),
!$       ! despite not even (always?) providing full OpenMP 4.0.
!$       omp_name = 'OpenMP 5.0 Preview 1 non-normative Technical Report'
!$    case (201511)
!$       omp_name = 'OpenMP 4.5'
!$    case (201307)
!$       omp_name = 'OpenMP 4.0'
!$    case (201107)
!$       omp_name = 'OpenMP 3.1'
!$    case (200805)
!$       omp_name = 'OpenMP 3.0'
!$    case (200505)
!$       omp_name = 'OpenMP 2.5'
!$    case (200011)
!$       omp_name = 'OpenMP 2.0'
!$    ! Earlier versions of OpenMP (1.x) did not specify
!$    ! the value of _OPENMP
!$    case default
!$       omp_name = 'unknown OpenMP version'
!$ end select
!$ write(6,'(a,i0,a)') '* OpenMP version: ', omp_version, ' ('//omp_name//').'
#endif

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

end module version_info
