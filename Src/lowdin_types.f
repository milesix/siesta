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
!     of the Wannier90 transformation within SIESTA. This specification
!     is read from the fdf file by routine 'read_lowdin_specs'.

      use precision, only: dp
      use files,     only: label_length         ! Number of characters in slabel
      use trialorbitalclass

      implicit none

      type, public ::  lowdin_manifold_t
          character(label_length+3)  :: seedname_lowdin  
                                       ! Seed of the name of the file 
                                       !   that will be used to dump the 
                                       !   information of the matrix elements
                                       !   produced by the WANNIER90 subroutines
          integer                :: initial_band
          integer                :: final_band
          integer                :: number_of_bands
                                       ! Number of bands passed to Wannier90
          integer                :: numbands_lowdin
                                       ! Number of bands that will be 
                                       !   orthonormalized
          integer                :: num_iter       ! Number of iterations for th
                                                   !   minimization of \Omega
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
                                                   ! Blocking size for orbital
                                                   !   distribution across 
                                                   !   the nodes
          integer                :: nincbands_loc_lowdin     
                                                   ! Number of included bands 
                                                   !   in the calc.
                                                   !   of the overlap and 
                                                   !   projection matrices
                                                   !   in the local node
          real(dp)               :: dis_win_min   !< Bottom of the outer energy
                                                  !!   window for band 
                                                  !!   disentanglement
                                                  !!   Units: energy
                                                  !!   Read from the fdf input
                                                  !!   in the 
                                                  !!   %block WannierProjections
          real(dp)               :: dis_win_max   !< Top of the outer energy
                                                  !!   window for band 
                                                  !!   disentanglement
                                                  !!   Units: energy
                                                  !!   Read from the fdf input
                                                  !!   in the 
                                                  !!   %block WannierProjections
          real(dp)               :: dis_froz_min  !< Bottom of the inner 
                                                  !!   (frozen) energy window 
                                                  !!   for band disentanglement
                                                  !!   Units: energy
                                                  !!   Read from the fdf input
                                                  !!   in the 
                                                  !!   %block WannierProjections
          real(dp)               :: dis_froz_max  !< Top of the inner 
                                                  !!   (frozen) energy window 
                                                  !!   for band disentanglement
                                                  !!   Units: energy
                                                  !!   Read from the fdf input
                                                  !!   in the 
                                                  !!   %block WannierProjections
          logical                :: disentanglement !< Is the disentanglement 
                                                  !!   procedure required for
                                                  !!   this manifold?
          logical                :: wannier_plot  !! Plot the Wannier functions?
          integer                :: wannier_plot_supercell(3)
                                                  !! Size of the supercell for 
                                                  !!   plotting the WF
          logical                :: fermi_surface_plot  
                                                  !! Plot the Fermi surface?
          logical                :: write_hr      !! Write the Hamiltonian in 
                                                  !!   the WF basis?
          logical                :: write_tb      !! Write lattice vectors, 
                                                  !!   Hamiltonian, and position
                                                  !!   operator in WF basis
          type(trialorbital), allocatable  :: proj_lowdin(:)
      end type lowdin_manifold_t

      type(lowdin_manifold_t), public,
     .     allocatable, save, target     :: manifold_bands_lowdin(:)


      integer :: numkpoints_w90_in   ! Total number of k-points
                                     !   used in Wannier90

      real(dp), pointer :: kpointsfrac_w90_in(:,:)
                                     ! List of k points relative
                                     !   to the reciprocal lattice vectors.
                                     !   First  index: component
                                     !   Second index: k-point index in the list
!
! Variables related with the neighbours of the k-points
!
      integer           :: nncount_w90_in
                             ! The number of nearest neighbours belonging to
                             !   each k-point of the Monkhorst-Pack mesh
      integer, pointer  :: nnlist_w90_in(:,:)
                             ! nnlist(ikp,inn) is the index of the
                             !   inn-neighbour of ikp-point
                             !   in the Monkhorst-Pack grid folded to the
                             !   first Brillouin zone
      integer, pointer  :: nnfolding_w90_in(:,:,:)
                             ! nnfolding(i,ikp,inn) is the i-component
                             !   of the reciprocal lattice vector
                             !   (in reduced units) that brings
                             !   the inn-neighbour specified in nnlist
                             !   (which is in the first BZ) to the
                             !   actual \vec{k} + \vec{b} that we need.
                             !   In reciprocal lattice units. 
      real(dp) :: latvec_w90_in(3,3)        !< Lattice vectors
      real(dp) :: reclatvec_w90_in(3,3)     !< Reciprocal lattice vectors
                                     !!  Cartesian coordinates in Bohr^-1
                                     !!  First  index: component
                                     !!  Second index: vector
      
      real(dp), pointer :: bvectorsfrac_w90_in(:,:)
                                         !! The vectors b that connect
                                         !!   each mesh-point k
                                         !!   to its nearest neighbours

      integer :: kmesh_w90_in(3)     ! Number of divisions along the three
                                     !   reciprocal lattice vectors that will be
                                     !   used in Wannier90



      end module lowdin_types

