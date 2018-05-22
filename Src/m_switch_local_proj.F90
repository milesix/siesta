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
  use files,     only : label_length         ! Number of characters in slabel
  use trialorbitalclass

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
  real(dp)         :: latvec(3,3) 
  real(dp)         :: reclatvec(3,3)     !< Reciprocal lattice vectors
                                         !!  Cartesian coordinates in Bohr^-1 
                                         !!  First  index: component
                                         !!  Second index: vector
!
! Variables related with the k-point list for which the overlap
! matrices Mmn between a k-point and its neighbor will be computed
!
  integer          :: nncount            !< The number of nearest
                                         !!   neighbours belonging to
                                         !!   each k-point of the 
                                         !!   Monkhorst-Pack mesh
  integer, pointer :: nnlist(:,:)        !< nnlist(ikp,inn) is the index of the
                                         !!   inn-neighbour of ikp-point
                                         !!   in the Monkhorst-Pack grid 
                                         !!   folded to the first Brillouin zone
  integer, pointer :: nnfolding(:,:,:)   !< nnfolding(i,ikp,inn) is the
                                         !!   i-component of the reciprocal 
                                         !!   lattice vector
                                         !!   (in reduced units) that brings
                                         !!   the inn-neighbour specified in 
                                         !!   nnlist (which is in the first BZ)
                                         !!   to the actual \vec{k} + \vec{b} 
                                         !!   that we need.
                                         !!   In reciprocal lattice units.
  real(dp), pointer :: bvectorsfrac(:,:)
                                         !! The vectors b that connect
                                         !!   each mesh-point k
                                         !!   to its nearest neighbours


!
! Variables related with the coefficients of the wavefunctions and
! eigenvalues at the Wannier90 k-point mesh
!
  complex(dp), pointer :: coeffs(:,:,:) => null() ! Coefficients of the wavefunctions.
                                         !   First  index: orbital
                                         !   Second index: band
                                         !   Third  index: k-point
  real(dp),    pointer :: eo(:,:) => null()        ! Eigenvalues of the Hamiltonian
                                         !   at the numkpoints introduced in
                                         !   kpointsfrac
                                         !   First  index: band index
                                         !   Second index: k-point index

!
! Variables related with the projections with trial functions,
! initial approximations to the MLWF
!
  integer  :: numproj        ! Total number of projection centers,
                             !   equal to the number of MLWF

  type(trialorbital), target, allocatable  :: projections(:)


!
! Output matrices
!

  complex(dp), pointer :: Mmnkb(:,:,:,:) => null()  
                                         ! Matrix of the overlaps of
                                         !   periodic parts of Bloch waves.
                                         !   <u_{ik}|u_{jk+b}>
                                         !   The first two indices refer to
                                         !   the number of occupied bands
                                         !   (indices m and n in standard
                                         !   notation, see for instance,
                                         !   Eq. (27) of the paper by
                                         !   Marzari et al., RMP 84, 1419 (2012)
                                         !   The third index refer to the kpoint
                                         !   The fourth index refer to the neig
  complex(dp), pointer :: Amnmat(:,:,:) => null() 
                                         ! Projections of a trial function
                                         !   with a Bloch orbital
                                         !   <\psi_{m k}|g_n>

!
! Variables related with the input/output files
!
  character(label_length+3)  :: seedname ! Name of the file where the Wannier90
                                         !   code, when used as a postprocessing
                                         !   tool, reads or dumps the
                                         !   information.



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
     use lowdin_types, only: nncount_lowdin       !< The number of nearest
                                                  !!   neighbours belonging to
                             !   each k-point of the Monkhorst-Pack mesh
     use lowdin_types, only: nnlist_lowdin
                             ! nnlist(ikp,inn) is the index of the
                             !   inn-neighbour of ikp-point
                             !   in the Monkhorst-Pack grid folded to the
                             !   first Brillouin zone
     use lowdin_types, only: nnfolding_lowdin
                             ! nnfolding(i,ikp,inn) is the i-component
                             !   of the reciprocal lattice vector
                             !   (in reduced units) that brings
                             !   the inn-neighbour specified in nnlist
                             !   (which is in the first BZ) to the
                             !   actual \vec{k} + \vec{b} that we need.
                             !   In reciprocal lattice units.
     use lowdin_types, only: bvectorsfrac_lowdin
     use lowdin_types, only: latvec_lowdin
 
    use wannier90_types, only: numbands_wannier   !< Number of bands for 
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
     use wannier90_types, only: nncount_wannier   !< The number of nearest
                                                  !!   neighbours belonging to
                             !   each k-point of the Monkhorst-Pack mesh
     use wannier90_types, only: nnlist_wannier
                             ! nnlist(ikp,inn) is the index of the
                             !   inn-neighbour of ikp-point
                             !   in the Monkhorst-Pack grid folded to the
                             !   first Brillouin zone
     use wannier90_types, only: nnfolding_wannier
                             ! nnfolding(i,ikp,inn) is the i-component
                             !   of the reciprocal lattice vector
                             !   (in reduced units) that brings
                             !   the inn-neighbour specified in nnlist
                             !   (which is in the first BZ) to the
                             !   actual \vec{k} + \vec{b} that we need.
                             !   In reciprocal lattice units.
 
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
    use wannier90_types, only: latvec_wannier     !< Reciprocal lattice vectors
    use wannier90_types, only: reclatvec_wannier  !< Reciprocal lattice vectors
    use wannier90_types, only: bvectorsfrac_wannier
    use wannier90_types, only: numproj_wannier    !< Number of projectors
    use wannier90_types, only: seedname_wannier
    use wannier90_types, only: projections_wannier

    use alloc,           only: re_alloc           !< Reallocation routines

    integer, intent(in) :: index_manifold         !< Index of the band manifold
                                                  !!   in the Lowdin 
                                                  !!   orthonormalization

!
! Internal variables
!
    integer :: ik             ! Counter for loop on ik points
    integer :: nn             ! 
    integer :: i              ! 
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

!     Initialize the list of neighbour k-points
      nullify( nnlist        )
      nullify( nnfolding     )

!     Broadcast information regarding the number of k-points neighbours
!     and allocate in all nodes the corresponding arrays containing information
!     about k-point neighbours

      nncount = nncount_lowdin

      call re_alloc( nnlist, 1, numkpoints, 1, nncount,           &
 &                   name = "nnlist", routine = "read_nnkp" )
      call re_alloc( nnfolding, 1, 3, 1, numkpoints, 1, nncount,  &
 &                   name = "nnfolding", routine = "read_nnkp" )
      nnlist     = nnlist_lowdin
      nnfolding  = nnfolding_lowdin

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

      latvec = ucell

!     Initialize number of projectors
      numproj = numincbands(1) 
      if( allocated(projections) ) deallocate( projections )
      allocate(projections(numproj))
      projections = manifold_bands_lowdin(index_manifold)%proj_lowdin

!     Reciprocal lattice vectors
      call reclat( ucell, reclatvec, 1 )


      nullify( bvectorsfrac )
      call re_alloc( bvectorsfrac, 1, 3, 1, nncount,    &
                     name="bvectorsfrac", routine = "chosing_b_vectors")
      bvectorsfrac = bvectorsfrac_lowdin

      seedname = manifold_bands_lowdin(index_manifold)%seedname_lowdin
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

!     Initialize the list of neighbour k-points
      nullify( nnlist        )
      nullify( nnfolding     )

!     Broadcast information regarding the number of k-points neighbours
!     and allocate in all nodes the corresponding arrays containing information
!     about k-point neighbours

      nncount = nncount_wannier

      call re_alloc( nnlist, 1, numkpoints, 1, nncount,           &
 &                   name = "nnlist", routine = "read_nnkp" )
      call re_alloc( nnfolding, 1, 3, 1, numkpoints, 1, nncount,  &
 &                   name = "nnfolding", routine = "read_nnkp" )
      nnlist     = nnlist_wannier
      nnfolding  = nnfolding_wannier

      numbands          = numbands_wannier
      numincbands       = numincbands_wannier
      nincbands_loc     = nincbands_loc_wannier
      blocksizeincbands = blocksizeincbands_wannier

!     Initialize the list of excluded bands
      nullify( isexcluded )
      call re_alloc( isexcluded, 1, no_u, name='isexcluded', &
 &                   routine='switch_local_projection' )
      isexcluded = isexcluded_wannier

      latvec    = latvec_wannier
!     Reciprocal lattice vector
      reclatvec = reclatvec_wannier

      nullify( bvectorsfrac )
      call re_alloc( bvectorsfrac, 1, 3, 1, nncount,    &
                     name="bvectorsfrac", routine = "chosing_b_vectors")
      bvectorsfrac = bvectorsfrac_wannier

      numproj = numproj_wannier
      if( allocated(projections) ) deallocate( projections )
      allocate(projections(numproj))
      projections = projections_wannier

      seedname = seedname_wannier
    
    endif

100 continue

!!   For debugging
!    write(6,'(/a,i5)')                                         & 
! &    'switch_local_projection: index_manifold    = ', index_manifold
!    write(6,'(a,i5)')                                          & 
! &    'switch_local_projection: numkpoints        = ', numkpoints
!    do ik = 1, numkpoints
!      write(6,'(a,i5,3f12.5)')                                &
! &      'switch_local_projection: ik, kpointsfrac = ',        &
! &      ik, kpointsfrac(:,ik) 
!    enddo
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
! &      'switch_local_projection: lattice vectors   = ',       &
! &      ivec, latvec(:,ivec)
!    enddo
!    do ivec = 1, 3
!      write(6,'(a,i5,3f12.5)')                                 &
! &      'switch_local_projection: reciprocal lattice= ',       &
! &      ivec, reclatvec(:,ivec)
!    enddo
!    write(6,'(a)') 'begin nnkpts'
!    write(6,'(i4)') nncount
!    do ik = 1,numkpoints
!      do nn = 1, nncount
!        write(6,'(2i6,3x,3i4)') &
! &        ik,nnlist(ik,nn),(nnfolding(i,ik,nn),i=1,3)
!      end do
!    end do
!    write(6,'(a)') 'begin bvectorsfrac'
!    do nn = 1, nncount
!      write(6,'(i6,3x,3f12.5)') &
! &      nn,(bvectorsfrac(nn,i),i=1,3)
!    end do
!    write(6,'(a/)') 'end bvectorsfrac'
!      
!!   End debugging

  end subroutine switch_local_projection

end module m_switch_local_projection
