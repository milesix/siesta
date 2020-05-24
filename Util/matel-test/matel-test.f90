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
      
      implicit none

      character(len=20), allocatable, dimension(:) :: species_label
      
      nspecies = 1
      allocate(species_label(1))
      species_label(1) = "Si"
      
      write(6,'(/a)') 'Reading PAOs and KBs from ascii files...'
      call read_basis_ascii(nspecies,species_label)

      write(*,*) species(1)%label
      
!      call register_rfs()

    end program  matel_test




