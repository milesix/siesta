!
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!> \brief General purpose of the switch_local_projection module:
!!
!! The interface between Siesta and Wannier90 shares many similarities
!! with the Lowdin projection algorithm.
!! Here, depending on the kind of projection to localized orbitals 
!! that will be carried out,
!! the relevant matrices and arrays for the projections will be populated.

module m_switch_local_projection

  use precision,      only: dp                !< Real double precision type
  use siesta_options, only: w90_processing    !< Will we call the interface with
                                              !!    Wannier90
  use siesta_options, only: lowdin_processing !< Will we call the interface with
                                              !!   the Lowdin orthonormalization
  use atomlist,       only: no_u              !< Number of orbitals in unit cell
  use siesta_geom,    only: ucell             !< Unit cell lattice vectors

  implicit none

  integer           :: numkpoints        !< Total number of k-points
                                         !!   in the chosen method of projection
  real(dp), pointer :: kpointsfrac(:,:)  !< List of k points relative
                                         !!   to the reciprocal lattice vectors.
                                         !!   First  index: components.
                                         !!   Second index: k-point index in 
                                         !!   the list
  integer           :: numbands(2)       !< Number of bands for wannierization
                                         !!   before excluding bands
  integer           :: numincbands(2)    !< Number of included bands in the 
                                         !!   calculation of the chosen 
                                         !!   projection method after excluding
                                         !!   bands
  logical, pointer :: isexcluded(:)      !< Masks excluded bands
  integer          :: blocksizeincbands  !< Maximum number of bands
                                         !!   considered for wannierization 
                                         !!   per node
  integer          :: nincbands_loc      !! Number of included bands in the 
                                         !!   calculation of the overlap and 
                                         !!   projection matrices.
                                         !!   in the local Node
  real(dp)         :: reclatvec(3,3)     !< Reciprocal lattice vectors
                                         !!  Cartesian coordinates in Bohr^-1 
                                         !!  First  index: component
                                         !!  Second index: vector


! Routines
  public :: switch_local_projection

  CONTAINS

  subroutine switch_local_projection( index_manifold )
    use lowdin_types, only: manifold_bands_lowdin !< Derived type where all 
                                                  !!   the details of the 
                                                  !!   band manifolds for 
                                                  !!   Lowdin orthonomalization
                                                  !!   are stored
    use lowdin_types, only: numkpoints_lowdin     !< Number of k-points in 
                                                  !!   the Monkhorst-Pack grid
                                                  !!   that will be used in
                                                  !!   the Lowdin 
                                                  !!   orthonormalization
    use lowdin_types, only: kpointsfrac_lowdin    !< List of k-points in 
                                                  !!   the Monkhorst-Pack grid
                                                  !!   that will be used in
                                                  !!   the Lowdin 
                                                  !!   orthonormalization
                                                  !!   (in fractional units, 
                                                  !!   i.e. relative to the
                                                  !!   reciprocal space 
                                                  !!   lattice vectors)
    use wannier90_types, only: numkpoints_wannier !< Number of k-points in 
                                                  !!   the Monkhorst-Pack grid
                                                  !!   for which the overlap of
                                                  !!   the periodic part of the
                                                  !!   wavefunct with a 
                                                  !!   k-point will be computed
    use wannier90_types, only: kpointsfrac_wannier!< List of k-points in 
                                                  !!   the Monkhorst-Pack grid
                                                  !!   relative to the 
                                                  !!   reciprocal lattice 
                                                  !!   vectors.
                                                  !!   First  index: components
                                                  !!   Second index: k-point 
                                                  !!   index in the list
    use wannier90_types, only: numbands_wannier   !< Number of bands for 
                                                  !!   wannierization before
                                                  !!   excluding bands
    use wannier90_types, only: numincbands_wannier!< Number of included bands
                                                  !!   in the calculation of
                                                  !!   the overlap and projec.
                                                  !!   matrices after excluding
                                                  !!   bands
    use wannier90_types, only: nincbands_loc_wannier
                                                  !< Number of bands in the 
                                                  !!   local node for 
                                                  !!   wannierization after
                                                  !!   after excluding bands
    use wannier90_types, only: blocksizeincbands_wannier
                                                  !< Maximum number of bands per
                                                  !!   node considered for 
                                                  !!   for wannierization
    use wannier90_types, only: isexcluded_wannier !< Masks excluded bands for
                                                  !!   Wannier
    use wannier90_types, only: reclatvec_wannier  !< Reciprocal lattice vectors
    use alloc,           only: re_alloc           !< Reallocation routines

    integer, intent(in) :: index_manifold         !< Index of the band manifold
                                                  !!   in the Lowdin 
                                                  !!   orthonormalization

!
! Internal variables
!
    integer :: ik             ! Counter for loop on ik points
    integer :: iband          ! Counter for loop on bands
    integer :: iorb           ! Counter for loop on atomic orbitals
    integer :: ivec           ! Counter for loop on vectors

    if( lowdin_processing ) then
      write(6,'(/a)')                                                     & 
 &      'switch_local_projection: Populating the relevant matrices for '  
      write(6,'(a)')                                                      & 
 &      'switch_local_projection: the Lowdin orthonormalization'
      write(6,'(a,i5)')                                                   & 
 &      'switch_local_projection: band manifold = ', index_manifold
      numkpoints = numkpoints_lowdin
!     Initialize the list of k-points
      nullify( kpointsfrac )
      call re_alloc( kpointsfrac, 1, 3, 1, numkpoints,  &
 &                   name='kpointsfrac', routine='switch_local_projection')
      do ik = 1, numkpoints
        kpointsfrac(:,ik) = kpointsfrac_lowdin(:,ik)
      enddo

      numbands(1)   = manifold_bands_lowdin(index_manifold)%numbands_lowdin
      numincbands(1)= manifold_bands_lowdin(index_manifold)%number_of_bands
      nincbands_loc = manifold_bands_lowdin(index_manifold)%nincbands_loc_lowdin
      blocksizeincbands =  &
 &       manifold_bands_lowdin(index_manifold)%blocksizeincbands_lowdin

!     Initialize the list of excluded bands
      nullify( isexcluded )
      call re_alloc( isexcluded, 1, no_u, name='isexcluded', &
 &                   routine='switch_local_projection' )
      isexcluded = manifold_bands_lowdin(index_manifold)%isexcluded

!     Reciprocal lattice vectors
      call reclat( ucell, reclatvec, 1 )

      goto 100
    endif

    if( w90_processing ) then
      write(6,'(/a)')                                                     & 
 &      'switch_local_projection: Populating the relevant matrices for '
      write(6,'(a)')                                                      & 
 &      'switch_local_projection: the interface with Wannier 90'
      numkpoints = numkpoints_wannier

!     Initialize the list of k-points
      nullify( kpointsfrac )
      call re_alloc( kpointsfrac, 1, 3, 1, numkpoints,  &
 &                   name='kpointsfrac', routine='switch_local_projection')
      do ik = 1, numkpoints
        kpointsfrac(:,ik) = kpointsfrac_wannier(:,ik)
      enddo

      numbands          = numbands_wannier
      numincbands       = numincbands_wannier
      nincbands_loc     = nincbands_loc_wannier
      blocksizeincbands = blocksizeincbands_wannier

!     Initialize the list of excluded bands
      nullify( isexcluded )
      call re_alloc( isexcluded, 1, no_u, name='isexcluded', &
 &                   routine='switch_local_projection' )
      isexcluded = isexcluded_wannier

!     Reciprocal lattice vector
      reclatvec = reclatvec_wannier
    
    endif

100 continue

!!   For debugging
!    write(6,'(/a,i5)')                                         & 
! &    'switch_local_projection: index_manifold    = ', index_manifold
!    write(6,'(a,i5)')                                          & 
! &    'switch_local_projection: numkpoints        = ', numkpoints
!!    do ik = 1, numkpoints
!!      write(6,'(a,i5,3f12.5)')                                &
!! &      'switch_local_projection: ik, kpointsfrac = ',        &
!! &      ik, kpointsfrac(:,ik) 
!!    enddo
!    write(6,'(a,2i5)')                                         &
! &    'switch_local_projection: numbands          = ',         &
! &      numbands
!    write(6,'(a,2i5)')                                         &
! &    'switch_local_projection: numincbands       = ',         &
! &     numincbands
!    write(6,'(a,2i5)')                                         &
! &    'switch_local_projection: nincbands_loc     = ',         &
! &     nincbands_loc
!    write(6,'(a,2i5)')                                         &
! &    'switch_local_projection: blocksizeincbands = ',         &
! &     blocksizeincbands
!    do iband = 1, no_u
!      write(6,'(a,i5,l5)')                                     &
! &      'switch_local_projection: iband, isexcluded = ',       &
! &      iband, isexcluded(iband)
!    enddo
!    do ivec = 1, 3
!      write(6,'(a,i5,3f12.5)')                                 &
! &      'switch_local_projection: reciprocal lattice= ',       &
! &      ivec, reclatvec(:,ivec)
!    enddo
!!   End debugging

  end subroutine switch_local_projection

end module m_switch_local_projection
