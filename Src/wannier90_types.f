! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors
!
      module wannier90_types
!
!=======================================================================
!
!     This module defines data structures to handle the specification
!     of the Lowdin orthonormalization. This specification
!     is read from the fdf file by routine 'read_lowdin_specs'.

      use precision, only: dp

      implicit none


!
! Variables related with the atomic structure of the system
!

      real(dp) :: reclatvec_wannier(3,3) ! Reciprocal lattice vectors.
                                         !   Cartesian coordinates.
                                         !   Readed in Angstroms^-1 
                                         !   and transformed
                                         !   to Bohr^-1 internally
                                         !   First  index: component
                                         !   Second index: vector
                                         !   This is consistent with the 
                                         !     reciprocal lattice read in Siesta
                                         !     but it has changed with respect
                                         !     the first version implemented 
                                         !     by R. Korytar


!
! Variables related with the k-point list for which the overlap
! matrices Mmn between a k-point and its neighbor will be computed
!
      integer  :: numkpoints_wannier ! Total number of k-points
                                     !   for which the overlap of the
                                     !   periodic part of the wavefunct
                                     !   with a neighbour k-point will
                                     !   be computed
      real(dp), pointer :: kpointsfrac_wannier(:,:)
                                     ! List of k points relative
                                     !   to the reciprocal lattice vectors.
                                     !   First  index: component
                                     !   Second index: k-point index in the list
!
! Variables related with the number of bands considered for wannierization
!
      integer :: numbands_wannier(2)    ! Number of bands for wannierization
                                        !    before excluding bands
      integer :: numincbands_wannier(2) ! Number of included bands in the calc.
                                        !   of the overlap and projection 
                                        !   matrices after excluding the bands
      integer :: nincbands_loc_wannier  ! Number of bands for wannierizatio
                                        !   after excluding bands
                                        !   in the local node
      integer :: blocksizeincbands_wannier
                                        ! Maximum number of bands
                                        !   considered for wannierization
                                        !   per node
      logical, pointer :: isexcluded_wannier(:) 
                                        ! Masks excluded bands

      end module wannier90_types
