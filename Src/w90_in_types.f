! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors

!> \brief General purpose of the w90_in_siesta types module:
!! defines data structures to handle the specification
!!     of the Wannier90 transformation within SIESTA. 
!!
!! This specification
!! is read from the fdf file by routine 'read_w90_in_specs'.
!!
!! Here we define the derived type w90_in_manifold_t, where all the
!! relevant information for a manifold of bands is stored.
!!
!! Then, we also define some relevant information regarding the k-point
!! sampling, the atomic structure, and eventually,
!! the chemical potential that can be applied to a particular
!! Wannier function

      module w90_in_siesta_types
!
!

      use precision, only: dp          ! Double precision
      use files,     only: label_length! Number of characters in slabel
      use trialorbitalclass            ! Determines all the variables
                                       !   that define the localized
                                       !   functions used as first guess
                                       !   during the minimization of the
                                       !   spreading

      implicit none

      type, public ::  w90_in_manifold_t
          character(label_length+3)  :: seedname_w90_in
                                       ! Seed of the name of the file 
                                       !   that will be used to dump the 
                                       !   information of the matrix elements
                                       !   produced by the WANNIER90 subroutines
          integer                :: initial_band
                                       !   Initial band of the manifold
                                       !     to be wannierized
          integer                :: final_band
                                       !   Final band of the manifold
                                       !     to be wannierized
          integer                :: number_of_bands
                                       ! Number of bands within the manifold
                                       !     passed to WANNIER90
          integer                :: numbands_w90_in
                                       ! Number of bands that will be 
                                       !     orthonormalized.
                                       !     It might be different from the 
                                       !     number of bands in the manifold.
                                       !     In this case, a disentanglement
                                       !     procedure is required
          integer                :: num_iter       
                                       ! Number of iterations for th
                                       !     minimization of \Omega
          integer, pointer       :: orbital_indices(:)
                                       ! Indices of the orbitals that will be
                                       !     used as localized functions in the
                                       !     initial guess of the minimization
                                       !     of the spreading
          logical, pointer       :: isexcluded(:)  
                                       ! List of bands excluded for
                                       !     Wannier transformation
          integer, pointer       :: isincluded(:)  
                                       ! List of bands included for
                                       !     Wannier transformation
          logical, pointer       :: orbexcluded(:) 
                                       ! List of orbitals excluded 
                                       !     for Wannier transformation
          integer, pointer       :: orb_in_manifold(:) 
                                       ! Sequential index of the orbitals that
                                       !     will be used in the Wannier 
                                       !     transformation
          integer                :: blocksizeincbands_w90_in
                                       ! Blocking size for orbital
                                       !     distribution across the nodes
          integer                :: nincbands_loc_w90_in
                                       ! Number of included bands 
                                       !     in the calculation 
                                       !     of the overlap and 
                                       !     projection matrices
                                       !     in the local node
          real(dp)               :: dis_win_min   
                                       ! Bottom of the outer energy
                                       !     window for band disentanglement
                                       !     Units: energy
                                       !     Read from the fdf input
                                       !     in the 
                                       !     %block WannierProjections
          real(dp)               :: dis_win_max 
                                       ! Top of the outer energy
                                       !     window for band disentanglement
                                       !     Units: energy
                                       !     Read from the fdf input
                                       !     in the 
                                       !     %block WannierProjections
          real(dp)               :: dis_froz_min  
                                       ! Bottom of the inner (frozen) energy
                                       !     window for band disentanglement
                                       !     Units: energy
                                       !     Read from the fdf input
                                       !     in the 
                                       !     %block WannierProjections
          real(dp)               :: dis_froz_max  
                                       ! Top of the inner (frozen) energy
                                       !     window for band disentanglement
                                       !     Units: energy
                                       !     Read from the fdf input
                                       !     in the 
                                       !     %block WannierProjections
          logical                :: disentanglement  
                                       ! Is the disentanglement 
                                       !   procedure required for this manifold?
          logical                :: wannier_plot  
                                       ! Plot the Wannier functions?
          integer                :: wannier_plot_supercell(3)
                                       ! Size of the supercell for 
                                       !   plotting the WF
          logical                :: fermi_surface_plot  
                                       ! Plot the Fermi surface?
          logical                :: write_hr      
                                       ! Write the Hamiltonian in 
                                       !   the WF basis?
          logical                :: write_tb      
                                       ! Write lattice vectors, 
                                       !   Hamiltonian, and position
                                       !   operator in WF basis
          type(trialorbital), allocatable  :: proj_w90_in(:)
                                       ! Initial guess for the localized 
                                       !   functions
      end type w90_in_manifold_t

      type(w90_in_manifold_t), public,
     .     allocatable, save, target     :: manifold_bands_w90_in(:)

!
! Variables related with the neighbours of the k-points
!

      integer :: numkpoints_w90_in     ! Total number of k-points
                                       !   used in Wannier90

      real(dp), pointer :: kpointsfrac_w90_in(:,:)
                                       ! List of k points relative
                                       !   to the reciprocal lattice vectors.
                                       !   First  index: component
                                       !   Second index: k-point index in 
                                       !      the list
      integer           :: nncount_w90_in
                                       ! The number of nearest neighbours
                                       !      belonging to each k-point of the 
                                       !      Monkhorst-Pack mesh
      integer, pointer  :: nnlist_w90_in(:,:)
                                       ! nnlist(ikp,inn) is the index of the
                                       !   inn-neighbour of ikp-point
                                       !   in the Monkhorst-Pack grid folded 
                                       !   to the first Brillouin zone
      integer, pointer  :: nnfolding_w90_in(:,:,:)
                                       ! nnfolding(i,ikp,inn) is the i-component
                                       !   of the reciprocal lattice vector
                                       !   (in reduced units) that brings
                                       !   the inn-neighbour specified in nnlist
                                       !   (which is in the first BZ) to the
                                       !   actual \vec{k} + \vec{b} that we need
      real(dp), pointer :: bvectorsfrac_w90_in(:,:)
                                       ! The vectors b that connect
                                       !   each mesh-point k
                                       !   to its nearest neighbours
                                       !   In reciprocal lattice units. 
      integer           :: kmesh_w90_in(3)    
                                       ! Number of divisions along the three
                                       !   reciprocal lattice vectors that
                                       !   will be used in Wannier90

!
! Variables related with the atomic structure
!
      real(dp) :: latvec_w90_in(3,3)   ! Lattice vectors
      real(dp) :: reclatvec_w90_in(3,3)! Reciprocal lattice vectors
                                       !  Cartesian coordinates in Bohr^-1
                                       !  First  index: component
                                       !  Second index: vector

!
! Variables related with the chemical potential
!
      real(dp), pointer :: chempotwann_val(:)
                                       ! Chemical potential 
                                       !   applied to shift the energy of a 
                                       !   the matrix elements in real space
                                       !   associated with a given 
                                       !   Wannier function
      integer           :: num_proj_local
                                       ! Number of projections that will be 
                                       !   handled in the local node

      end module w90_in_siesta_types

