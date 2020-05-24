! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
program matel_test

! The PAO and KB information can optionally be read from ASCII files
! (those produced by a standard run of Siesta or Base, but with the extension 
! renamed to '.input' instead of '.dump'), or from NetCDF files (if NetCDF
! is available). Note that there must be files for *all* species.
!
      use precision

      use m_ion_io
      use atm_types, only: species, nspecies
      use m_matel_registry, only: show_pool
      use matel_mod, only: init_matel, new_matel
      
      implicit none

      
      character(len=20), allocatable, dimension(:) :: species_label
      integer :: io
      real(dp) :: val, grad(3)
      
      nspecies = 1
      allocate(species_label(1))
      species_label(1) = "Si"
      
      write(6,'(/a)') 'Reading PAOs and KBs from ascii files...'
      call read_basis_ascii(nspecies,species_label)

      write(*,*) species(1)%label
      
      call register_rfs()
      call show_pool()
      call init_matel()

      io = species(1)%orb_gindex(1)

      call new_MATEL('S', io, io, [ 1.0_dp, 0.0_dp, 0.0_dp ], val, grad)
      print "(a,f10.4,2x,3f10.4)", "<3s|  | 3s(1.0)> : ", val, grad
      call new_MATEL('X', io, io, [ 1.0_dp, 0.0_dp, 0.0_dp ], val, grad)
      print "(a,f10.4,2x,3f10.4)", "<3s| X | 3s(1.0)> : ", val, grad
      call new_MATEL('Y', io, io, [ 1.0_dp, 0.0_dp, 0.0_dp ], val, grad)
      print "(a,f10.4,2x,3f10.4)", "<3s| Y | 3s(1.0)> : ", val, grad
      call new_MATEL('Z', io, io, [ 1.0_dp, 0.0_dp, 0.0_dp ], val, grad)
      print "(a,f10.4,2x,3f10.4)", "<3s| Z | 3s(1.0)> : ", val, grad
    end program  matel_test




