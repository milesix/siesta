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
      use files,     only : label_length         ! Number of characters in slabel
      use trialorbitalclass

      implicit none

      type, public ::  lowdin_manifold_t
          character(label_length+3)  :: seedname_lowdin  
                                         ! Name of the file where the Wannier90
                                         !   code, when used as a postprocessing
                                         !   tool, reads or dumps the
                                         !   information.

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
          integer, pointer       :: orb_in_manifold(:) 
                                                   ! List of orbitals excluded 
                                                   !   for Lowdin orthonormaliza
          integer                :: blocksizeincbands_lowdin 
                                                   ! Maximum number of bands
                                                   !    considered for 
                                                   !    Lowdin orthonormalizat
                                                   !    per node
          integer                :: numbands_lowdin
          integer                :: nincbands_loc_lowdin     
                                                   ! Number of included bands 
                                                   !   in the calc.
                                                   !   of the overlap and 
                                                   !   projection matrices
                                                   !   in the local node
          type(trialorbital), allocatable  :: proj_lowdin(:)
      end type lowdin_manifold_t

      type(lowdin_manifold_t), public,
     .     allocatable, save, target     :: manifold_bands_lowdin(:)


      integer :: numkpoints_lowdin   ! Total number of k-points
                                     !   used in the Lowdin normalization

      real(dp), pointer :: kpointsfrac_lowdin(:,:)
                                     ! List of k points relative
                                     !   to the reciprocal lattice vectors.
                                     !   First  index: component
                                     !   Second index: k-point index in the list
!
! Variables related with the neighbours of the k-points
!
      integer           :: nncount_lowdin
                             ! The number of nearest neighbours belonging to
                             !   each k-point of the Monkhorst-Pack mesh
      integer, pointer  :: nnlist_lowdin(:,:)
                             ! nnlist(ikp,inn) is the index of the
                             !   inn-neighbour of ikp-point
                             !   in the Monkhorst-Pack grid folded to the
                             !   first Brillouin zone
      integer, pointer  :: nnfolding_lowdin(:,:,:)
                             ! nnfolding(i,ikp,inn) is the i-component
                             !   of the reciprocal lattice vector
                             !   (in reduced units) that brings
                             !   the inn-neighbour specified in nnlist
                             !   (which is in the first BZ) to the
                             !   actual \vec{k} + \vec{b} that we need.
                             !   In reciprocal lattice units. 
      real(dp) :: latvec_lowdin(3,3)        !< Lattice vectors
      real(dp) :: reclatvec_lowdin(3,3)     !< Reciprocal lattice vectors
                                     !!  Cartesian coordinates in Bohr^-1
                                     !!  First  index: component
                                     !!  Second index: vector
      
      real(dp), pointer :: bvectorsfrac_lowdin(:,:)
                                         !! The vectors b that connect
                                         !!   each mesh-point k
                                         !!   to its nearest neighbours
      complex(dp), pointer, save ::   overlaptilde(:,:)
      complex(dp), pointer, save ::   coeffs_k(:,:)
      complex(dp), pointer, save ::   overlap_sq(:,:)
      complex(dp), pointer, save ::   phitilde(:,:)
      complex(dp), pointer, save ::   invsqrtover(:,:)
      complex(dp), pointer, save ::   coeffshatphi(:,:)

      integer :: kmeshlowdin(3)      ! Number of divisions along the three
                                     !   reciprocal lattice vectors that will be
                                     !   used in the Lowdinn projections



      end module lowdin_types

