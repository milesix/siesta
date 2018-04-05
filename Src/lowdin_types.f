! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors
!
      module lowdin_types
!
!=======================================================================
!
!     This module defines data structures to handle the specification
!     of the Lowdin orthonormalization. This specification
!     is read from the fdf file by routine 'read_lowdin_specs'.

      use precision, only: dp

      implicit none

      type, public ::  lowdin_manifold_t
          integer                :: initial_band
          integer                :: final_band
          integer                :: number_of_bands
          integer, pointer       :: orbital_indices(:)
          logical, pointer       :: isexcluded(:)  ! List of bands excluded for
                                                   !   Lowdin orthonormaliza
          integer, pointer       :: isincluded(:)  ! List of bands included for
                                                   !   Lowdin orthonormaliza
          logical, pointer       :: orbexcluded(:) ! List of orbitals excluded 
                                                   !   for Lowdin orthonormaliza
          integer, pointer       :: orb_in_manifold(:) ! List of orbitals excluded 
                                                   !   for Lowdin orthonormaliza
          integer                :: blocksizeincbands_lowdin 
                                                   ! Maximum number of bands
                                                   !    considered for 
                                                   !    Lowdin orthonormalizat

          integer                :: nincbands_loc_lowdin     
                                                   ! Number of included bands 
                                                   !   in the calc.
                                                   !   of the overlap and 
                                                   !   projection matrices.
                                                   !   in the local Node
      end type lowdin_manifold_t

      type(lowdin_manifold_t), public,
     .     allocatable, save, target     :: manifold_bands_lowdin(:)

      complex(dp), pointer, save ::   overlaptilde(:,:)
      complex(dp), pointer, save ::   coeffs_k(:,:)
      complex(dp), pointer, save ::   overlap_sq(:,:)
      complex(dp), pointer, save ::   phitilde(:,:)
      complex(dp), pointer, save ::   invsqrtover(:,:)
      complex(dp), pointer, save ::   coeffshatphi(:,:)

      
 

      end module lowdin_types

