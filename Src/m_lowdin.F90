! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!> \brief General purpose of the Lowdin module 
!! 
!! In this module we perform a Lowdin orthonormalization of the Bloch functions
!! corresponding to a given manifold of bands
!! We follow the recipe given in \cite Iniguez-04. 
!!
!! Example: band structure of SrTiO3
!! \image html Bands_STO.png
module m_lowdin


  use precision,      only: dp                     !< Real double precision type
  use siesta_options, only: n_lowdin_manifolds     !< Number of bands manifolds 
                                                   !!  that will be considered
                                                   !!  for Lowdin transformation
  use parallel,       only: Node, Nodes, IOnode, BlockSize
  use m_spin,         only: nspin                  ! Number of spin components
  use m_switch_local_projection, only: coeffs
  use lowdin_types,   only: lowdin_manifold_t    
  use lowdin_types,   only: manifold_bands_lowdin  ! Variable where the initial
                                                   !   and final band of each
                                                   !   manifold are stored
  use lowdin_types,   only: coeffshatphi           !
  use lowdin_types,   only: invsqrtover            !
  use lowdin_types,   only: overlap_sq             !
  use lowdin_types,   only: overlaptilde           !
  use lowdin_types,   only: phitilde               !
  use lowdin_types,   only: numkpoints_lowdin
  use lowdin_types,   only: kpointsfrac_lowdin     !
  use lowdin_types,   only: nncount_lowdin
  use lowdin_types,   only: nnlist_lowdin
  use lowdin_types,   only: nnfolding_lowdin
  use lowdin_types,   only: bvectorsfrac_lowdin

  use sparse_matrices,    only: maxnh        ! Maximum number of orbitals
                                             !   interacting
                                             ! NOTE: While running in parallel,
                                             !   maxnh changes from one core to
                                             !   the other
  use sparse_matrices,    only: numh         ! Number of nonzero element of each
                                             !   row of the hamiltonian matrix
  use sparse_matrices,    only: listh        ! Nonzero hamiltonian-matrix elemen
  use sparse_matrices,    only: listhptr     ! Pointer to start of each row
                                             !   of the hamiltonian matrix
  use sparse_matrices,    only: xijo         ! Vectors between orbital centers
  use sparse_matrices,    only: tight_binding_param ! Parameters of the
                                             !   tight-binding Hamiltonian
  use sparse_matrices,    only: H            ! Real space Hamiltonian
  use sparse_matrices,    only: S            ! Real space Overlap
  use atomlist,           only: indxuo       ! Index of equivalent orbital in
                                             !   the unit cell
  use atomlist,           only: no_u         ! Number of orbitals in unit cell
                                             ! NOTE: When running in parallel,
                                             !   this is core independent
  use sys,                only : die         ! Termination routine
  use w90_kmesh,          only : kmesh_get
  use w90_kmesh,          only : kmesh_write
  use w90_parameters,     only : nntot
  use w90_parameters,     only : nnlist
  use w90_parameters,     only : neigh
  use w90_parameters,     only : nncell
  use fdf

!
! Allocation/Deallocation routines
!
  use alloc,              only: re_alloc     ! Reallocation routines
  use alloc,              only: de_alloc     ! Deallocation routines

  use parallelsubs,       only: LocalToGlobalOrb

  implicit none

! Routines
  public :: read_lowdin_specs
  public :: check_normalization
  public :: define_phitilde
  public :: define_overlap_phitilde
  public :: define_invsqover_phitilde
  public :: define_coeffshatphi
  public :: compute_tight_binding_param
  public :: compute_position
  public :: write_tight_binding_param
  public :: set_excluded_bands_lowdin
  public :: setup_Lowdin 
  public :: allocate_matrices_Lowdin
  public :: deallocate_matrices_Lowdin
  public :: define_trial_orbitals

  private

  CONTAINS

! subroutine read_lowdin_specs         : Subroutine that reads all the
!                                        info in the fdf file related with the
!                                        Lowdin orthogonalization

  subroutine read_lowdin_specs( )
!
!   Processes the information in an fdf file
!   regarding the bands that will enter into the orthonormalization procedure
!   and the atomic orbitals that will be orthonormalized.
!
    use fdf
    use m_cite,   only: add_citation
    use parallel, only: IONode        ! Node for Input/output
    use siesta_geom,  only: ucell      ! Unit cell lattice vectors
    use lowdin_types, only: reclatvec_lowdin
    use lowdin_types, only: kmeshlowdin
    use w90_parameters, only: mp_grid
    use w90_parameters, only: gamma_only
    use w90_parameters, only: num_kpts
    use w90_parameters, only: recip_lattice
    use w90_parameters, only: real_lattice
    use w90_parameters, only: kpt_latt
    use w90_parameters, only: param_read
    use m_digest_nnkp,  only: chosing_b_vectors 
                                       ! Subroutine that computes the b
                                         ! vectors that connect each mesh
                                         ! k-point to its nearest neighbours.
    use files,         only: slabel        ! Short system label,
                                         !   used to generate file names

!   Internal variables
    integer :: index_manifold      ! Counter for the number of manifolds
    integer :: iorb                ! Counter for the number of atomic orbitals
    integer :: ik                  ! Counter for loop on k-points
    integer :: ikx                 ! Counter for loop on k-points along 1-direc
    integer :: iky                 ! Counter for loop on k-points along 2-direc
    integer :: ikz                 ! Counter for loop on k-points along 3-direc
    integer :: nkp                 ! Counter for loop on k-points along 3-direc
    integer :: nn                  ! Counter for loop on k-points along 3-direc
    integer :: i                   ! Counter for loop on k-points along 3-direc

    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline

!   Allocate the pointer where the initial and final band of every 
!   manifold will be stored
    allocate(manifold_bands_lowdin(n_lowdin_manifolds))
    
!   Read the LowdinProjections block
    if (.not. fdf_block('LowdinProjections',bfdf)) RETURN

!   Add citation
    if ( IONode ) then
      call add_citation("arXiv:cond-mat/0407677")
    end if

    do while(fdf_bline(bfdf, pline))     !! over band manifolds to be 
                                         !! orthonormalized
      if (.not. fdf_bmatch(pline,'iii'))        &    ! We expect that each line
                                                     !   contains three integers
                                                     !   That is the meaning of
                                                     !   iii
 &      call die('Wrong format in LowdinProjections')
        
!     Assign the initial and final band of each manifold
      index_manifold = fdf_bintegers(pline,1)
      manifold_bands_lowdin(index_manifold)%initial_band=fdf_bintegers(pline,2)
      manifold_bands_lowdin(index_manifold)%final_band  =fdf_bintegers(pline,3)
      manifold_bands_lowdin(index_manifold)%number_of_bands=      &
 &        ( manifold_bands_lowdin(index_manifold)%final_band   -  &
 &          manifold_bands_lowdin(index_manifold)%initial_band )  + 1
      manifold_bands_lowdin(index_manifold)%numbands_lowdin =     &
 &        manifold_bands_lowdin(index_manifold)%final_band 

      nullify( manifold_bands_lowdin(index_manifold)%orbital_indices )
      call re_alloc( manifold_bands_lowdin(index_manifold)%orbital_indices,  &
 &        1, manifold_bands_lowdin(index_manifold)%number_of_bands )

      do iorb = 1, manifold_bands_lowdin(index_manifold)%number_of_bands
        manifold_bands_lowdin(index_manifold)%orbital_indices(iorb) =        &
 &         fdf_bintegers(pline,3+iorb)
      enddo

      write(manifold_bands_lowdin(index_manifold)%seedname_lowdin,"(a,'.',i1.1)") &
 &       trim(slabel), index_manifold

!!     For debugging
!         write(6,*)manifold_bands_lowdin(index_manifold)%seedname_lowdin
!      write(6,'(a,4i5)') &
! &      'index manifold, initial band, final band, number of bands = ',      &
! &      index_manifold,                                                      &
! &      manifold_bands_lowdin(index_manifold)%initial_band,                  &
! &      manifold_bands_lowdin(index_manifold)%final_band,                    &
! &      manifold_bands_lowdin(index_manifold)%number_of_bands
!     write(6,'(a,2i5)')'read_lowdin_specs: Number of bands for Lowdin = ',   &
! &      index_manifold, manifold_bands_lowdin(index_manifold)%numbands_lowdin
!
!      do iorb = 1, manifold_bands_lowdin(index_manifold)%number_of_bands
!        write(6,*) manifold_bands_lowdin(index_manifold)%orbital_indices(iorb)
!      enddo
!!     End debugging
    enddo ! end loop over band manifolds

!   Read the data to generate the grid in reciprocal space that will be used
!   for the Lowdin Projections
    if (.not. fdf_block('kMeshforLowdin',bfdf)) RETURN

    do while(fdf_bline(bfdf, pline))     
      if (.not. fdf_bmatch(pline,'iii'))        &   ! We expect that each line
                                                    !   contains three integers
                                                    !   That is the meaning of
                                                    !   iii
                                                    ! The first integer is the
                                                    !   number of divisions 
                                                    !   the first reciprocal 
                                                    !   lattice vector and so on
 &      call die('Wrong format in kMeshforLowdin')
      kmeshlowdin(1) = fdf_bintegers(pline,1)
      kmeshlowdin(2) = fdf_bintegers(pline,2)
      kmeshlowdin(3) = fdf_bintegers(pline,3)
    enddo 

!     Define the total number of k-points used in the Lowdin projection
      numkpoints_lowdin = kmeshlowdin(1) * kmeshlowdin(2) * kmeshlowdin(3)

!     Compute and store the components of the k-points in fractional units
      nullify( kpointsfrac_lowdin )
      call re_alloc( kpointsfrac_lowdin, 1, 3, 1, numkpoints_lowdin,  &
 &                   name='kpointsfrac_lowdin', routine='read_lowdin_specs')

      ik = 0
      do ikx = 0, kmeshlowdin(1) - 1
        do iky = 0, kmeshlowdin(2) - 1
          do ikz = 0, kmeshlowdin(3) - 1
            ik = ik + 1
            kpointsfrac_lowdin(1,ik) = (ikx*1.0_dp)/kmeshlowdin(1)
            kpointsfrac_lowdin(2,ik) = (iky*1.0_dp)/kmeshlowdin(2)
            kpointsfrac_lowdin(3,ik) = (ikz*1.0_dp)/kmeshlowdin(3)
          enddo 
        enddo 
      enddo 

!!     For debugging
!      write(6,'(a,a,3i5)')'read_lowdin_specs: Number of subdivisions ',  &  
! &      'of the reciprocal vectors for Lowdin = ', kmeshlowdin(:)
!      write(6,'(a,a,3i5)')'read_lowdin_specs: Number of k-points used ', & 
! &      'in the Lowdin projection = ', numkpoints_lowdin
!      do ik = 1, numkpoints_lowdin
!        write(6,'(a,a,i5,3f12.5)')'read_lowdin_specs: k-points in ',     & 
! &        'fractional units:', ik, kpointsfrac_lowdin(:,ik)
!      enddo
!!     End debugging

!     Reciprocal lattice vectors
      real_lattice  = ucell * 0.529177_dp
      call reclat( real_lattice, reclatvec_lowdin, 1 )
      mp_grid       = kmeshlowdin
      num_kpts      = numkpoints_lowdin
      recip_lattice = reclatvec_lowdin 
      allocate ( kpt_latt(3,num_kpts) )
      kpt_latt      = kpointsfrac_lowdin
      gamma_only    = .false.

      call param_read
      call kmesh_get( )


      nncount_lowdin = nntot
!     Initialize the list of neighbour k-points
      nullify( nnlist_lowdin    )
      nullify( nnfolding_lowdin )

!     Broadcast information regarding the number of k-points neighbours
!     and allocate in all nodes the corresponding arrays containing information
!     about k-point neighbours

      call re_alloc( nnlist_lowdin, 1, numkpoints_lowdin, 1, nncount_lowdin,   &
 &                   name = "nnlist_lowdin", routine = "read_lowdin_specs" )
      call re_alloc( nnfolding_lowdin, 1, 3, 1, numkpoints_lowdin,             &
 &                   1, nncount_lowdin, name = "nnfolding_lowdin",             &
 &                   routine = "read_lowdin_specs" )
      nnlist_lowdin     = nnlist
      nnfolding_lowdin  = nncell

!     Compute the vectors that connect each mesh k-point 
!     to its nearest neighbours
      call chosing_b_vectors( kpointsfrac_lowdin, nncount_lowdin,  &
 &                          nnlist_lowdin, nnfolding_lowdin,       &
 &                          bvectorsfrac_lowdin )


!!     For debugging
!      write(6,'(a)') 'begin nnkpts'
!      write(6,'(i4)') nntot
!      do nkp=1,num_kpts
!         do nn=1,nntot
!            write(6,'(2i6,3x,3i4)') &
!               nkp,nnlist(nkp,nn),(nncell(i,nkp,nn),i=1,3)
!         end do
!      end do
!      write(6,'(a/)') 'end nnkpts'
!!     End debugging

  endsubroutine read_lowdin_specs

  subroutine setup_Lowdin

    use parallel,  only : Node                 ! This process node
    use parallel,  only : Nodes                ! Total number of processor nodes
    use parallel,  only : BlockSize            ! Total number of processor nodes

#ifdef MPI
    use parallelsubs,       only : set_blocksizedefault
!
!   Subroutine to order the indices of the different bands after
!   excluding some of them for wannierization
!
    use m_orderbands,       only: order_index
#endif

    integer :: index_manifold
    integer :: numincbands_tmp
    integer :: blocksizeincbands_tmp 
    integer :: nincbands_loc_tmp

#ifdef MPI
    integer, external :: numroc
#endif

    do index_manifold = 1, n_lowdin_manifolds
       call set_excluded_bands_lowdin( index_manifold )
       numincbands_tmp = manifold_bands_lowdin(index_manifold)%number_of_bands

!      Allocate memory related with the coefficients of the wavefunctions
#ifdef MPI
!      Find the number of included bands for Wannierization that will be stored
!      per node. Use a block-cyclic distribution of nincbands over Nodes.
!
       call set_blocksizedefault( Nodes, numincbands_tmp, 
 &                                blocksizeincbands_tmp )

!       write(6,'(a,3i5)')' diagonalizeHk: Node, Blocksize = ', &
! &                               Node, blocksizeincbands_tmp

       nincbands_loc_tmp = numroc( numincbands_tmp,
                               blocksizeincbands_tmp, Node, 0, Nodes )
#else
       nincbands_loc_tmp = numincbands_tmp
#endif

       manifold_bands_lowdin(index_manifold)%nincbands_loc_lowdin = &
 &        nincbands_loc_tmp
       manifold_bands_lowdin(index_manifold)%blocksizeincbands_lowdin = &
 &        blocksizeincbands_tmp

#ifdef MPI
!      Set up the arrays that control the indices of the bands to be
!      considered after excluding some of them for wannierization
!      This is done once and for all the k-points
       call order_index( no_l, no_u, numincbands_tmp )
#endif

    enddo  

    nullify( tight_binding_param )
    call re_alloc( tight_binding_param,                                 &
 &                 1, n_lowdin_manifolds,                               &
 &                 1, maxnh,                                            &
 &                 1, nspin,                                            &
 &                 name='tight_binding_param',                          &
 &                 routine='setup_Lowdin',                              &
 &                 shrink=.false., copy=.false.)
    tight_binding_param = 0.0_dp

  endsubroutine setup_Lowdin

  subroutine allocate_matrices_Lowdin( index_manifold )
    integer, intent(in) :: index_manifold

    integer number_of_orbitals_to_project
    integer number_of_bands_in_manifold_local

    number_of_orbitals_to_project =                                      &
 &        manifold_bands_lowdin(index_manifold)%number_of_bands
    number_of_bands_in_manifold_local =                                  &
 &        manifold_bands_lowdin(index_manifold)%nincbands_loc_lowdin

!   Allocate the different matrices required for Lowdin orthonormalization

    nullify( overlaptilde )
    call re_alloc( overlaptilde,                                            &
 &                 1, number_of_orbitals_to_project,                        &
 &                 1, number_of_orbitals_to_project,                        &
 &                 name='overlaptilde', routine='allocate_matrices_Lowdin', &
 &                 shrink=.false., copy=.false.)
    overlaptilde = cmplx(0.0_dp, 0.0_dp, kind=dp)

    nullify( invsqrtover )
    call re_alloc( invsqrtover,                                             &
 &                 1, number_of_orbitals_to_project,                        &
 &                 1, number_of_orbitals_to_project,                        &
 &                 name='invsqrtover', routine='allocate_matrices_Lowdin',  &
 &                 shrink=.false., copy=.false.)
    invsqrtover = cmplx(0.0_dp, 0.0_dp, kind=dp)

    nullify( overlap_sq )
    call re_alloc( overlap_sq,                                              &
 &                 1, no_u,                                                 &
 &                 1, no_u,                                                 &
 &                 name='overlap_sq', routine='allocate_matrices_Lowdin',   &
 &                 shrink=.false., copy=.false.)
    overlap_sq = cmplx(0.0_dp, 0.0_dp, kind=dp)

    nullify( coeffshatphi )
    call re_alloc( coeffshatphi,                                            &
 &                 1, number_of_bands_in_manifold_local,                    &
 &                 1, number_of_orbitals_to_project,                        &
 &                 name='coeffshatphi', routine='allocate_matrices_Lowdin', &
 &                 shrink=.false., copy=.false.)
    coeffshatphi = cmplx(0.0_dp, 0.0_dp, kind=dp)

  endsubroutine allocate_matrices_Lowdin

  subroutine deallocate_matrices_Lowdin
    call de_alloc( overlaptilde,                             &
 &                 name='overlaptilde',                      & 
 &                 routine='deallocate_matrices_Lowdin' )
    call de_alloc( overlap_sq,                               &
 &                 name='overlap_sq',                        & 
 &                 routine='deallocate_matrices_Lowdin' )
    call de_alloc( invsqrtover,                              &
 &                 name='invsqrtover',                       & 
 &                 routine='deallocate_matrices_Lowdin' )
    call de_alloc( coeffshatphi,                             &
 &                 name='coeffshatphi',                      & 
 &                 routine='deallocate_matrices_Lowdin' )
  endsubroutine deallocate_matrices_Lowdin



  subroutine set_excluded_bands_lowdin( index_manifold )
    
    integer, intent(in) :: index_manifold

    integer :: iband
    integer :: iorb
    integer :: iterex
    integer :: index
    integer :: index_orb_included

    nullify( manifold_bands_lowdin(index_manifold)%isexcluded )
    call re_alloc( manifold_bands_lowdin(index_manifold)%isexcluded,      &
 &                 1, no_u,                                               &
 &                 name='isexcluded',                                     &
 &                 routine='set_excluded_bands_lowdin' )

    nullify( manifold_bands_lowdin(index_manifold)%orbexcluded )
    call re_alloc( manifold_bands_lowdin(index_manifold)%orbexcluded,     &
 &                 1, no_u,                                               &
 &                 name='orbexcluded',                                    &
 &                 routine='set_excluded_bands_lowdin' )

    nullify( manifold_bands_lowdin(index_manifold)%orb_in_manifold )
    call re_alloc( manifold_bands_lowdin(index_manifold)%orb_in_manifold, &
 &                 1, no_u,                                               &
 &                 name='orb_in_manifold',                                &
 &                 routine='set_excluded_bands_lowdin' )


!   By default, all the bands are excluded from the calculation
    manifold_bands_lowdin(index_manifold)%isexcluded(:) = .true.

    do iband = 1, no_u
!     Exclude the corresponding bands from the computation
      if( (iband .ge. manifold_bands_lowdin(index_manifold)%initial_band) &
 &        .and.                                                           &
 &        (iband .le. manifold_bands_lowdin(index_manifold)%final_band) ) then
        manifold_bands_lowdin(index_manifold)%isexcluded( iband ) = .false.
      endif
    enddo

!!   For debugging
!    do iterex = 1, no_u
!     write(6,'(a,i5,l5)')                                               &
! &     'excluded_bands_lowind (iband, excluded)',                       &
! &     iterex, manifold_bands_lowdin(index_manifold)%isexcluded(iterex)
!    enddo
!!   End debugging

!   By default, all the orbitals are excluded from the calculation
    manifold_bands_lowdin(index_manifold)%orbexcluded(:) = .true.
    do iorb = 1, manifold_bands_lowdin(index_manifold)%number_of_bands
      index_orb_included =                                          & 
 &      manifold_bands_lowdin(index_manifold)%orbital_indices(iorb)
      manifold_bands_lowdin(index_manifold)%orbexcluded(index_orb_included) = &
 &      .false. 
    enddo

    index = 0
    manifold_bands_lowdin(index_manifold)%orb_in_manifold(:) = 0
    do iorb = 1, manifold_bands_lowdin(index_manifold)%number_of_bands
      index = index + 1
      index_orb_included = manifold_bands_lowdin(index_manifold)%orbital_indices(iorb)
      manifold_bands_lowdin(index_manifold)%orb_in_manifold(index_orb_included) = index
    enddo

!!   For debugging
!    do iterex = 1, no_u
!     write(6,'(a,2i5,l5,i5)')                                            &
! &     'orbitals excluded: manifold, orb_unit cell, excluded, orb_mani', &
! &     index_manifold, iterex,                                           &
! &     manifold_bands_lowdin(index_manifold)%orbexcluded(iterex),        &
! &     manifold_bands_lowdin(index_manifold)%orb_in_manifold(iterex)    
!    enddo
!!   End debugging

  endsubroutine set_excluded_bands_lowdin


!> \f{eqnarray*}{
!!    \langle \psi_{j} (\vec{k}) \vert \psi_{i} (\vec{k}) \rangle & = &
!!    \sum_{\mu \nu} c_{j \nu} (\vec{k})
!!    \langle \phi_{\nu} (\vec{k}) \vert \phi_{\mu} (\vec{k}) \rangle
!!    c_{\mu i} (\vec{k})   \\
!!    & = & \sum_{\mu \nu} c_{\nu j}^{\ast} (\vec{k})
!!    S_{\nu \mu} (\vec{k}) c_{\mu i} (\vec{k})  \\
!!    & = & \delta_{ij},
!! \f}
  subroutine check_normalization( ik, no_u, overlap_sq )

    integer, intent(in)       :: ik
    integer, intent(in)       :: no_u
    complex(dp), intent(in)   :: overlap_sq(no_u,no_u)

    complex(dp), pointer, save ::   normalization(:,:)

!
!   Internal variables
!
    integer iuo      ! Counter for loop on atomic orbitals
    integer juo      ! Counter for loop on atomic orbitals

    nullify( normalization )
    call re_alloc(normalization, 1, no_u, 1, no_u,                     &
 &                name='normalization', routine='check_normalization', &
 &                shrink=.false., copy=.false.)
    normalization = cmplx(0.0_dp, 0.0_dp, kind=dp)

!   Check the normalization of the eigenfunctions
!   that come out of the diagonalization
    normalization = matmul( overlap_sq, coeffs(:,:,ik) )
    normalization = matmul( transpose(conjg(coeffs(:,:,ik))),normalization)

!   For debugging
    write(6,'(a)') ' Normalization: '
    do iuo = 1, no_u
      do juo = 1, no_u
        write(6,'(2i5,2f12.5)') iuo, juo, normalization(iuo,juo)
      enddo
    enddo
!   End debugging
   
    call de_alloc( normalization,                &
 &                 name='normalization',         & 
 &                 routine='check_normalization' )

  endsubroutine check_normalization

  subroutine define_phitilde( ik, no_u, overlap_sq, phitilde )

    integer, intent(in)       :: ik
    integer, intent(in)       :: no_u
    complex(dp), intent(in)   :: overlap_sq(no_u,no_u)
    complex(dp), intent(out)  :: phitilde(no_u,no_u)

!
!   Internal variables
!
    integer jband    ! Counter for loop on bands
    integer mu       ! Counter for loop on atomic orbitals
    integer nu       ! Counter for loop on atomic orbitals
    integer lambda   ! Counter for loop on atomic orbitals

    complex(dp), pointer, save ::   aux(:,:) 

    nullify( aux )
    call re_alloc(aux, 1, no_u, 1, no_u,                  &
 &                name='aux', routine='define_phitilde',  &
 &                shrink=.false., copy=.false.)
    aux = cmplx(0.0_dp, 0.0_dp, kind=dp)

    do mu  = 1, no_u
      do jband = 1, no_u
        do nu = 1, no_u
          aux(jband,mu) = aux(jband,mu) +                    &
 &          conjg(coeffs(nu,jband,ik)) * overlap_sq(nu,mu)
        enddo 
      enddo
    enddo

    do mu  = 1, no_u
      do lambda = 1, no_u
        do jband = 1, no_u
          phitilde(lambda,mu) = phitilde(lambda,mu) +        &
 &          coeffs(lambda,jband,ik) * aux(jband,mu)
        enddo 
      enddo
    enddo

!
!   If we want to perform the previous operation directly using the
!   matrix multiplication internal routines.
!   This is only for testing, since while using this matmul function:
!   - We sum over all the bands (and not only over a chosen set)
!   - It is not well adapted for multiplication
!   In this case, if the band index runs over all the bands, 
!   phitilde has to be equal to delta_(\lambda,\mu)
!
!   phitilde = cmplx(0.0_dp,0.0_dp,kind=dp)
!   phitilde = matmul( transpose(conjg(coeffs(:,:,ik))), overlap_sq )
!   phitilde = matmul( coeffs(:,:,ik), phitilde )

!!   For debugging
!    write(6,'(a)')'Phitilde = '
!    do lambda = 1, no_u
!      do mu  = 1, no_u
!        write(6,'(2i5,2f12.5)') lambda, mu, phitilde(lambda,mu)
!      enddo
!    enddo
!!   End debugging

    call de_alloc( aux,                          &
 &                 name='aux',                   & 
 &                 routine='define_phitilde' )


  endsubroutine define_phitilde

!
!
!

  subroutine define_overlap_phitilde( index_manifold, ik )

    integer, intent(in)       :: index_manifold
    integer, intent(in)       :: ik

!
!   Internal variables
!
    integer jband    ! Counter for loop on bands
    integer mu       ! Counter for loop on atomic orbitals
    integer mu_index ! Counter for loop on atomic orbitals
    integer nu_index ! Counter for loop on atomic orbitals
    integer nu       ! Counter for loop on atomic orbitals
    integer lambda   ! Counter for loop on atomic orbitals
    integer rho      ! Counter for loop on atomic orbitals
    integer number_of_orbitals_to_project
    integer number_of_bands_in_manifold

    complex(dp), pointer, save ::   aux(:,:), aux2(:,:)

    number_of_orbitals_to_project =                                      &
 &        manifold_bands_lowdin(index_manifold)%number_of_bands
    number_of_bands_in_manifold =                                        &
 &        manifold_bands_lowdin(index_manifold)%number_of_bands

    nullify( aux, aux2 )
    call re_alloc( aux,                                                  &
 &                 1, number_of_bands_in_manifold,                       &
 &                 1, number_of_orbitals_to_project,                     &
 &                 name='aux', routine='define_overlap_phitilde',        &
 &                 shrink=.false., copy=.false. )
    aux = cmplx(0.0_dp, 0.0_dp, kind=dp)

    call re_alloc( aux2,                                                 &
 &                 1, no_u,                                              &
 &                 1, number_of_orbitals_to_project,                     &
 &                 name='aux2', routine='define_overlap_phitilde',       &
 &                 shrink=.false., copy=.false. )
    aux2 = cmplx(0.0_dp, 0.0_dp, kind=dp)

    overlaptilde = cmplx(0.0_dp,0.0_dp,kind=dp)

!!   For debugging
!    do mu  = 1, number_of_orbitals_to_project
!      write(6,'(a,3i5)') ' index_manifold, mu, orbital = ',        &
! &      index_manifold, mu,                                        &
! &      manifold_bands_lowdin(index_manifold)%orbital_indices(mu)
!    enddo
!
!    do jband = 1, number_of_bands_in_manifold
!      write(6,'(a,3i5)') ' index_manifold, jband = ',              &
! &      index_manifold, jband,                                     &
! &      manifold_bands_lowdin(index_manifold)%initial_band+ (jband - 1)
!    enddo
!!   End debugging

    do mu  = 1, number_of_orbitals_to_project
      mu_index = manifold_bands_lowdin(index_manifold)%orbital_indices(mu)
      do jband = 1, number_of_bands_in_manifold
        do lambda = 1, no_u
          aux(jband,mu) = aux(jband,mu) +                    &
 &          conjg(coeffs(lambda,jband,ik)) * overlap_sq(lambda,mu_index)
        enddo 
      enddo
    enddo

    do mu   = 1, number_of_orbitals_to_project
      do rho = 1, no_u
        do jband = 1, number_of_bands_in_manifold
          aux2(rho,mu) = aux2(rho,mu) +                     &
 &          coeffs(rho,jband,ik) * aux(jband,mu)
        enddo 
      enddo
    enddo

    do mu  = 1, number_of_orbitals_to_project
      do nu = 1, number_of_orbitals_to_project
        nu_index = manifold_bands_lowdin(index_manifold)%orbital_indices(nu)
        do rho = 1, no_u
          overlaptilde(nu,mu) = overlaptilde(nu,mu) +       &
 &          overlap_sq(nu_index,rho) * aux2(rho,mu)
        enddo 
      enddo
    enddo

!!   For debugging
!    write(6,'(a)')'Overlap phitilde = '
!    do nu = 1, number_of_orbitals_to_project
!      do mu  = 1, number_of_orbitals_to_project
!        write(6,'(2i5,2f12.5)') nu, mu, overlaptilde(nu,mu)
!      enddo
!    enddo
!!   End debugging

    call de_alloc( aux,                                &
 &                 name='aux',                         & 
 &                 routine='define_overlap_phitilde' )

    call de_alloc( aux2,                                &
 &                 name='aux2',                         & 
 &                 routine='define_overlap_phitilde' )

  endsubroutine define_overlap_phitilde

!
!
!

  subroutine define_invsqover_phitilde( index_manifold )

    integer, intent(in)       :: index_manifold

!
!   Internal variables
!
    integer mu       ! Counter for loop on atomic orbitals
    integer nu       ! Counter for loop on atomic orbitals
    integer lambda   ! Counter for loop on atomic orbitals

!
!   Variables required for the diagonalization of overlaptilde
!
    integer  :: ierror        ! Code for error message from cdiag
    integer  :: npsi_tilde    ! Variable to dimension the coefficient vector

    complex(dp), pointer, save ::   aux(:,:)
    complex(dp), pointer, save ::   unity(:,:)
    complex(dp), pointer, save ::   invsqrtd(:,:)
    real(dp), dimension(:), pointer :: epsilonstilde ! Eigenvalues of 
                                                     !   overlaptilde
    real(dp), pointer :: psi_tilde(:) => null()

!!   For debugging
!    complex(dp), pointer, save ::   normalization(:,:)
!    complex(dp), pointer, save ::   overlaptilde_save(:,:)
!!   End debugging

    integer number_of_orbitals_to_project
    integer number_of_bands_in_manifold

    external cdiag

    number_of_orbitals_to_project =                                      &
 &        manifold_bands_lowdin(index_manifold)%number_of_bands
    number_of_bands_in_manifold =                                        &
 &        manifold_bands_lowdin(index_manifold)%number_of_bands


    nullify( epsilonstilde )
    call re_alloc( epsilonstilde,                                      &
  &                1, number_of_orbitals_to_project,                   &
  &                name='epsilonstilde',                               &
  &                routine='define_invsqover_phitilde' )

    nullify( unity )
    call re_alloc( unity,                                              &
 &                 1, number_of_orbitals_to_project,                   &
 &                 1, number_of_orbitals_to_project,                   &
 &                 name='unity', routine='define_invsqover_phitilde',  &
 &                 shrink=.false., copy=.false.)
    unity = cmplx(0.0_dp,0.0_dp,kind=dp)

    do mu = 1, number_of_orbitals_to_project
      unity(mu,mu) = cmplx(1.0_dp,0.0_dp,kind=dp)
    enddo

    nullify( aux )
    call re_alloc( aux,                                                       &
 &                 1, number_of_orbitals_to_project,                          &
 &                 1, number_of_orbitals_to_project,                          &
 &                 name='aux', routine='define_invsqover_phitilde',           &
 &                 shrink=.false., copy=.false.)
    aux = cmplx(0.0_dp, 0.0_dp, kind=dp)


    nullify( invsqrtd )
    call re_alloc( invsqrtd,                                                  &
 &                 1, number_of_orbitals_to_project,                          &
 &                 1, number_of_orbitals_to_project,                          &
 &                 name='invsqrtd', routine='define_invsqover_phitilde',      &
 &                 shrink=.false., copy=.false.)
    invsqrtd = cmplx(0.0_dp,0.0_dp,kind=dp)

    npsi_tilde = 2*number_of_orbitals_to_project*number_of_orbitals_to_project
    call re_alloc( psi_tilde,  1, npsi_tilde,  name='psi_tilde',              &
 &                 routine='define_invsqover_phitilde' )

!!   For debugging
!    nullify( normalization )
!    call re_alloc( normalization,                                             &
! &                 1, number_of_orbitals_to_project,                          &
! &                 1, number_of_orbitals_to_project,                          &
! &                 name='normalization', routine='define_invsqover_phitilde', &
! &                 shrink=.false., copy=.false.)
!    normalization = cmplx(0.0_dp, 0.0_dp, kind=dp)
!
!    nullify( overlaptilde_save )
!    call re_alloc( overlaptilde_save,                                         &
! &                 1, number_of_orbitals_to_project,                          &
! &                 1, number_of_orbitals_to_project,                          &
! &                 name='overlaptilde_save',                                  &
! &                 routine='define_invsqover_phitilde',                       &
! &                 shrink=.false., copy=.false.)
!    overlaptilde_save = overlaptilde
!!   End debugging

!
!   Diagonalize overlaptilde.
!   We need it to compute the inverse of the root square
!
!   In the input, overlaptilde is the matrix to be diagonalized
!   After calling cdiag, overlaptilde contains the unitary transformation
!   matrix that transforms the overlap matrix into its diagonal form
!   L. F. Mattheiss, Phys. Rev. B 2, 3918 (1970)


    call cdiag( overlaptilde, unity,                  &
 &              number_of_orbitals_to_project,        &
 &              number_of_orbitals_to_project,        &
 &              number_of_orbitals_to_project,        &
 &              epsilonstilde,                        &
 &              psi_tilde,                            &
 &              number_of_orbitals_to_project,        &
 &              1, ierror, BlockSize )

    do mu = 1, number_of_orbitals_to_project
      invsqrtd(mu,mu) =      &
 &      cmplx( epsilonstilde(mu)**(-1.0_dp/2.0_dp), 0.0_dp, kind=dp )
    enddo

!!   For debugging
!    do mu = 1, no_u
!      write(6,'(i5,3f12.5)')mu, epsilonstilde(mu), invsqrtd(mu,mu)
!    enddo 
!!   End debugging

    do mu = 1, number_of_orbitals_to_project
      do nu = 1, number_of_orbitals_to_project
        do lambda = 1, number_of_orbitals_to_project
          aux(mu,nu) = aux(mu,nu) +                    &
 &          overlaptilde(mu,lambda) * invsqrtd(lambda,nu)
        enddo 
      enddo
    enddo

    invsqrtover = cmplx(0.0_dp, 0.0_dp, kind=dp)
    do mu = 1, number_of_orbitals_to_project
      do nu = 1, number_of_orbitals_to_project
        do lambda = 1, number_of_orbitals_to_project
          invsqrtover(mu,nu) = invsqrtover(mu,nu) +           &
 &          aux(mu,lambda) * conjg(overlaptilde(nu,lambda))
        enddo 
      enddo
    enddo

!!   For debugging
!!   Check that the computation of the inverse of the root square has been
!!   properly done.
!!
!    aux = cmplx(0.0_dp, 0.0_dp, kind=dp)
!    do mu = 1, number_of_orbitals_to_project
!      do nu = 1, number_of_orbitals_to_project
!        do lambda = 1, number_of_orbitals_to_project
!          aux(mu,nu) = aux(mu,nu) +                    &
! &          invsqrtover(mu,lambda) * overlaptilde_save(lambda,nu)
!        enddo 
!      enddo
!    enddo
!
!    do mu = 1, number_of_orbitals_to_project
!      do nu = 1, number_of_orbitals_to_project
!        do lambda = 1, number_of_orbitals_to_project
!          normalization(mu,nu) = normalization(mu,nu) +     &
! &          aux(mu,lambda) * invsqrtover(lambda,nu) 
!        enddo 
!      enddo
!    enddo
!
!    write(6,'(a)') ' Normalization: '
!    do mu = 1, number_of_orbitals_to_project
!      do nu = 1, number_of_orbitals_to_project
!        write(6,'(2i5,2f12.5)') mu, nu, normalization(mu,nu)
!      enddo
!    enddo
!!   End debugging


    call de_alloc( unity,                              &
 &                 name='unity',                       & 
 &                 routine='define_invsqover_phitilde' )

    call de_alloc( epsilonstilde,                      &
 &                 name='epsilonstilde',               & 
 &                 routine='define_invsqover_phitilde' )

    call de_alloc( invsqrtd,                           &
 &                 name='invsqrtd',                    & 
 &                 routine='define_invsqover_phitilde' )

    call de_alloc( aux,                                &
 &                 name='aux',                         & 
 &                 routine='define_invsqover_phitilde' )

    call de_alloc( psi_tilde,                          &
 &                 name='psi_tilde',                   & 
 &                 routine='define_invsqover_phitilde' )

!!   For debugging
!    call de_alloc( normalization,                      &
! &                 name='normalization',               & 
! &                 routine='define_invsqover_phitilde' )
!
!    call de_alloc( overlaptilde_save,                  &
! &                 name='overlaptilde_save',           & 
! &                 routine='define_invsqover_phitilde' )
!!   End debugging



  endsubroutine define_invsqover_phitilde

!
!
!  
  subroutine define_coeffshatphi( index_manifold, ik )

    integer, intent(in)       :: index_manifold
    integer, intent(in)       :: ik

!
!   Internal variables
!
    integer mu       ! Counter for loop on atomic orbitals
    integer nu       ! Counter for loop on atomic orbitals
    integer lambda   ! Counter for loop on atomic orbitals
    integer nu_index ! Counter for loop on atomic orbitals
    integer jband    ! Counter for loop on bands


    complex(dp), pointer, save ::   aux(:,:)

    integer number_of_orbitals_to_project
    integer number_of_bands_in_manifold_local

    number_of_orbitals_to_project =                                  &
 &        manifold_bands_lowdin(index_manifold)%number_of_bands
    number_of_bands_in_manifold_local =                              &
 &        manifold_bands_lowdin(index_manifold)%nincbands_loc_lowdin

    nullify( aux )
    call re_alloc( aux,                                              &
 &                 1, number_of_bands_in_manifold_local,             &
 &                 1, number_of_orbitals_to_project,                 &
 &                 name='aux', routine='define_coeffshatphi',        &
 &                 shrink=.false., copy=.false. )
    aux = cmplx(0.0_dp,0.0_dp,kind=dp)

    do nu  = 1, number_of_orbitals_to_project
      nu_index = manifold_bands_lowdin(index_manifold)%orbital_indices(nu)
      do jband = 1, number_of_bands_in_manifold_local 
        do lambda = 1, no_u
          aux(jband,nu) = aux(jband,nu) +                    &
 &          conjg(coeffs(lambda,jband,ik)) * overlap_sq(lambda,nu_index)
        enddo 
      enddo
    enddo

    coeffshatphi = cmplx(0.0_dp,0.0_dp,kind=dp)
    do mu  = 1, number_of_orbitals_to_project
      do jband = 1, number_of_bands_in_manifold_local
        do nu = 1, number_of_orbitals_to_project
          coeffshatphi(jband,mu) = coeffshatphi(jband,mu) +       &
 &          aux(jband,nu) * invsqrtover(nu,mu)
        enddo 
      enddo
    enddo

!!   For debugging
!    do jband = 1, number_of_bands_in_manifold_local
!      do mu  = 1, number_of_orbitals_to_project
!        write(6,'(a,2i5,2f12.5)')' jband, mu, coeffshat = ',     &
! &        jband, mu, coeffshatphi(jband,mu)
!      enddo
!    enddo
!!   End debugging

    call de_alloc( aux,                                   &
 &                 name='aux',                            & 
 &                 routine='define_invsqover_phitilde' )


  endsubroutine define_coeffshatphi

  subroutine compute_tight_binding_param( ispin, ik, index_manifold, & 
 &                                        kvector, epsilon )

    integer,     intent(in)  :: ispin
    integer,     intent(in)  :: ik
    integer,     intent(in)  :: index_manifold
    real(dp),    intent(in)  :: kvector(3)
    real(dp),    intent(in)  :: epsilon(no_u)  ! Eigenvalues of the Hamiltonian

    real(dp), pointer, save  :: tight_binding_param_k(:,:,:) 
    complex(dp), pointer     :: aux(:) 
    complex(dp), pointer     :: lpsi(:)


    integer  :: BNode         !
    integer  :: BTest         !
    integer  :: iuo
    integer  :: io
    integer  :: juo
    integer  :: ie 
    integer  :: iie 
    integer  :: ind
    integer  :: jbandini 
    real(dp) :: ckxij
    real(dp) :: skxij
    real(dp) :: kxij
    real(dp) :: qp1
    real(dp) :: qp2
    real(dp) :: eqp1
    real(dp) :: eqp2
    real(dp) :: ee

    integer  :: number_of_orbitals_to_project
    integer  :: number_of_bands_in_manifold_local


    number_of_orbitals_to_project =                                  &
 &        manifold_bands_lowdin(index_manifold)%number_of_bands
    number_of_bands_in_manifold_local =                              &
 &        manifold_bands_lowdin(index_manifold)%nincbands_loc_lowdin


    nullify( aux )
    call re_alloc(aux, 1, 5*number_of_orbitals_to_project,       &
 &                name='aux', routine='diagonalizeHk', &
 &                shrink=.false., copy=.false.)
    aux = cmplx(0.0_dp, 0.0_dp, kind=dp)

!    nullify( lpsi )
!    call re_alloc(lpsi, 1, number_of_orbitals_to_project,                 &
! &                name='lpsi', routine='diagonalizeHk',                   &
! &                shrink=.false., copy=.false.)
!    lpsi = cmplx(0.0_dp, 0.0_dp, kind=dp)

    nullify( tight_binding_param_k )
    call re_alloc( tight_binding_param_k,                                 &
 &                 1, 2,                                                  &
 &                 1, number_of_orbitals_to_project,                      &
 &                 1, number_of_orbitals_to_project,                      &
 &                 name='tight_binding_param_k',                          &
 &                 routine='compute_tight_binding_param',                 &
 &                 shrink=.false., copy=.false. )
    tight_binding_param_k = 0.0_dp

    call define_overlap_phitilde( index_manifold, ik )

    call define_invsqover_phitilde( index_manifold )

    call define_coeffshatphi( index_manifold, ik )


    jbandini = manifold_bands_lowdin(index_manifold)%initial_band
 
!   Add contribution to tight-binding-parameters of unit-cell orbitals
!   WARNING: Dk and Ek may be EQUIVALENCE'd to Haux and Saux
!$OMP parallel do default(shared), collapse(2)
!$OMP&private(iuo,juo)
    do iuo = 1, number_of_orbitals_to_project
      do juo = 1, number_of_orbitals_to_project
        tight_binding_param_k(1,juo,iuo) = 0.0_dp
        tight_binding_param_k(2,juo,iuo) = 0.0_dp
      enddo
    enddo
!$OMP end parallel do

! Global operation to form new density matrix
    BNode = 0
    iie = 0
    do ie = 1, number_of_bands_in_manifold_local
      if ( Node == BNode ) iie = iie + 1

      if ( Node == BNode ) then
         lpsi => coeffshatphi(iie,:)
      else
         lpsi => aux(:)
      endif
#ifdef MPI
      call MPI_Bcast( lpsi(1),number_of_orbitals_to_project,          &
     &                MPI_double_complex,                             &
     &                BNode,MPI_Comm_World,MPIerror )
#endif

      ee = epsilon(jbandini+ie-1)

!$OMP parallel do default(shared),
!$OMP&private(iuo,ind,juo,qp1,qp2,eqp1,eqp2)
      do iuo = 1, number_of_orbitals_to_project

        qp1 = ee * realpart(lpsi(iuo))
        qp2 = ee * imagpart(lpsi(iuo))

        do juo = 1, number_of_orbitals_to_project
          eqp1 = qp1 * realpart(lpsi(juo)) + qp2 * imagpart(lpsi(juo))
          eqp2 = qp1 * imagpart(lpsi(juo)) - qp2 * realpart(lpsi(juo))

          tight_binding_param_k(1,juo,iuo) = tight_binding_param_k(1,juo,iuo) + eqp1
          tight_binding_param_k(2,juo,iuo) = tight_binding_param_k(2,juo,iuo) + eqp2
!!     For debugging
!        write(6,'(a,3i5,1f12.5,2i5,6f12.5)')                              &
! &         ' k, i_man, band_local, band_global, eigen, iuo, juo = ',      &
! &           index_manifold, ie, jbandini+ie-1, ee*13.6058_dp, juo, iuo,  &
! &           realpart(lpsi(iuo)), imagpart(lpsi(iuo)), realpart(lpsi(juo)), imagpart(lpsi(juo)), &
! &           tight_binding_param_k(1,juo,iuo), tight_binding_param_k(2,juo,iuo)
!!     End debugging
        enddo
      enddo
!$OMP end parallel do
      BTest = ie/BlockSize
      if (BTest*BlockSize.eq.ie) then
        BNode = BNode + 1
        if (BNode .gt. Nodes-1) BNode = 0
      endif
    enddo

!!   For debugging
!    do iuo = 1, number_of_orbitals_to_project
!      do juo = 1, number_of_orbitals_to_project
!        write(6,'(a,2i5,2f15.7)')'iuo, juo, tbk1, tbk2 = ',       &
! &        iuo, juo,                                               &
! &        tight_binding_param_k(1,juo,iuo), tight_binding_param_k(2,juo,iuo)
!      enddo 
!    enddo
!!   End debugging

!$OMP parallel do default(shared),
!$OMP&private(iuo,ind,juo,kxij,ckxij,skxij)
    do io = 1, number_of_orbitals_to_project
      iuo = manifold_bands_lowdin(index_manifold)%orbital_indices(io)
      
      do ind = listhptr(iuo) + 1, listhptr(iuo) + numh(iuo)
        juo = indxuo(listh(ind))
        if( .not. manifold_bands_lowdin(index_manifold)%orbexcluded(juo) ) then
           juo = manifold_bands_lowdin(index_manifold)%orb_in_manifold(juo)
           kxij = kvector(1) * xijo(1,ind) +  &
 &                kvector(2) * xijo(2,ind) +  &
 &                kvector(3) * xijo(3,ind)
           kxij = -1.0_dp * kxij
           ckxij = cos(kxij)
           skxij = sin(kxij)
           tight_binding_param(index_manifold,ind,ispin) =      &
 &              tight_binding_param(index_manifold,ind,ispin) + &
 &              tight_binding_param_k(1,juo,io)*ckxij -        &
 &              tight_binding_param_k(2,juo,io)*skxij
!!        For debugging      
!           write(6,'(a,5i7,8f12.5)')   &
! &          'i_man,o_man,o_global,ind,o_neig,xijo,tb,tbk1,cos,tbk2,sin',  &
! &           index_manifold, io, iuo, ind, juo, xijo(:,ind),              &
! &           tight_binding_param(index_manifold,ind,ispin),               &
! &           tight_binding_param_k(1,juo,io), ckxij,                      &
! &           tight_binding_param_k(2,juo,io), skxij
!!        End debugging      
        endif
      enddo
    enddo
!$OMP end parallel do

    call de_alloc( aux,                                   &
 &                 name='aux',                            & 
 &                 routine='compute_tight_binding_param' )

!    call de_alloc( lpsi,                                   &
! &                 name='lpsi',                            & 
! &                 routine='compute_tight_binding_param' )

    call de_alloc( tight_binding_param_k,                  &
 &                 name='tight_binding_param_k',           & 
 &                 routine='compute_tight_binding_param' )


  end subroutine compute_tight_binding_param

  subroutine write_tight_binding_param( ispin, numkpoints )

    integer,     intent(in)  :: ispin
    integer,     intent(in)  :: numkpoints

    character(len=30) :: sname 
    character(len=41) :: filenameparam ! Name of the file where the 
                                       !   tight-binding parameters  
                                       !   will be written
    integer           :: fileunitparam ! Logical unit of the file

    integer           :: io
    integer           :: iuo
    integer           :: j
    integer           :: juo
    integer           :: ind
    integer           :: jo
    integer           :: index_manifold
    integer           :: number_of_orbitals_to_project
    integer      :: eof

    external     :: io_assign ! Assign a logical unit
    external     :: io_close  ! Close a logical unit

    if (Node.eq.0) then
      sname = fdf_string('SystemLabel','siesta')
    endif

    do index_manifold = 1, n_lowdin_manifolds
      number_of_orbitals_to_project =                              &
 &        manifold_bands_lowdin(index_manifold)%number_of_bands

      write(filenameparam,"(a,'.',i1.1,'.tb.param')") &
 &      trim(sname), index_manifold
!     For debugging
!      write(6,*)filenameparam
!     End debugging

      call io_assign(fileunitparam)
    
      open( unit=fileunitparam, err=199, file=filenameparam,       &
 &          status='replace', form='formatted', iostat=eof )

      do io = 1, number_of_orbitals_to_project
        iuo = manifold_bands_lowdin(index_manifold)%orbital_indices(io)
        do ind = listhptr(iuo) + 1, listhptr(iuo) + numh(iuo)
          juo = indxuo(listh(ind))
          if( .not. manifold_bands_lowdin(index_manifold)%orbexcluded(juo) ) then
            write(fileunitparam,'(4i7,4f15.7)') io, iuo,                      & 
 &           manifold_bands_lowdin(index_manifold)%orb_in_manifold(juo), juo, &
 &           xijo(:,ind),                                                     &
 &           tight_binding_param(index_manifold,ind,ispin)/numkpoints/2.0_dp
          endif
        enddo
      enddo
      call io_close(fileunitparam)

    enddo 
    
    return

199 call die('write_tight_binding_param: Error creating output parameter files')
     
  end subroutine write_tight_binding_param

  subroutine compute_position( ispin, index_manifold )

    use m_switch_local_projection, only: seedname 
                                         ! Seed for the name of the file
                                         !   where the matrix elements of the
                                         !   position operator will be written.
    use m_switch_local_projection, only: nncount
    use m_switch_local_projection, only: numkpoints
    use m_switch_local_projection, only: bvectorsfrac
    use m_switch_local_projection, only: Mmnkb
    use m_switch_local_projection, only: Amnmat
    use m_switch_local_projection, only: Amnmat_man
    use siesta_options, only: n_lowdin_manifolds
    use lowdin_types, only: kmeshlowdin
    use w90_constants,  only : cmplx_i
    use w90_constants,  only : cmplx_0
    use w90_constants,  only : cmplx_1
    use w90_constants,  only : twopi
    use w90_constants,  only : eps7
    use w90_parameters, only : wb
    use w90_parameters, only : bk
    use w90_parameters, only : lenconfac
    use w90_parameters, only : u_matrix
    use w90_parameters, only : m_matrix
    use w90_parameters, only : num_bands
    use w90_parameters, only : num_wann
    use w90_parameters, only : num_kpts
    use w90_parameters, only : nntot
    use w90_parameters, only : nnlist
    use w90_parameters, only : kpt_latt
    use w90_parameters, only : real_metric
    use w90_overlap,    only : overlap_project

    integer :: ispin
    integer :: index_manifold
    integer :: idir
    integer :: nn
    integer :: n
    integer :: n1
    integer :: n2
    integer :: n3
    integer :: i1
    integer :: i
    integer :: j 
    integer :: i2
    integer :: i3
    integer :: icnt
    integer :: ndiff(3)
    real(kind=dp) :: dist(125), dist_min, tot
    integer :: m
    integer :: nkp
    integer :: ik
    integer :: iw
    integer :: ind
    integer :: iband
    integer :: iproj
    integer :: iorb
    integer :: ix
    integer :: number_of_bands_in_manifold_local
    integer :: m_reset 
    integer :: n_reset 

    integer              :: nrpts      !! Number of Wigner-Seitz grid points
    integer              :: irpts      !! Counter for the Wigner-Seitz grid poin
    integer, allocatable :: irvec(:,:) !! The irpt-th Wigner-Seitz grid point 
                                       !!   has components irvec(1:3,irpt) 
                                       !!   in the basis of the lattice vectors

    complex(kind=dp) :: pos_r(3)
    complex(kind=dp) :: fac
    real(kind=dp)    :: delta
    real(kind=dp)    :: rdotk

    complex(kind=dp), allocatable :: csheet(:,:,:)
    real(kind=dp),    allocatable :: sheet (:,:,:)
    real(kind=dp),    allocatable :: rave(:,:)
    real(kind=dp),    allocatable :: ln_tmp(:,:,:)

    character(len=len_trim(seedname)+6) :: posfilename
                                           ! Name of the file where the 
                                           !   matrix elements of the position
                                           !   operator will be written
    integer      :: posunit   
                                           ! Logical unit assigned to the file 
                                           !   where the position matrices   
                                           !   will be written
    integer      :: eof

    external     :: io_assign              ! Assign a logical unit
    external     :: io_close               ! Close a logical unit


    number_of_bands_in_manifold_local =                                  &
 &        manifold_bands_lowdin(index_manifold)%nincbands_loc_lowdin

    num_bands = number_of_bands_in_manifold_local 
    num_wann  = number_of_bands_in_manifold_local 
    num_kpts  = numkpoints
    nntot     = nncount

    nnlist    = nnlist_lowdin

    allocate( csheet (number_of_bands_in_manifold_local, nncount, numkpoints) )
    allocate( sheet  (number_of_bands_in_manifold_local, nncount, numkpoints) )
    allocate( rave   (3, number_of_bands_in_manifold_local) )
    allocate( ln_tmp (number_of_bands_in_manifold_local, nncount, numkpoints) )

    csheet = cmplx_1 
    sheet  = 0.0_dp
    rave   = 0.0_dp

    if ( allocated(u_matrix) ) deallocate(u_matrix)
    allocate( u_matrix(number_of_bands_in_manifold_local,               &
 &                     number_of_bands_in_manifold_local,               &
 &                     numkpoints) )
    u_matrix = cmplx(0.0_dp,0.0_dp,kind=dp)
    if ( allocated(m_matrix) ) deallocate(m_matrix)
    allocate( m_matrix(number_of_bands_in_manifold_local,               &
 &                     number_of_bands_in_manifold_local,               &
 &                     nncount,                                         &
 &                     numkpoints) )
    m_matrix = cmplx(0.0_dp,0.0_dp,kind=dp)

    call amn( ispin )

    Amnmat_man(index_manifold,:,:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
    do ik = 1, numkpoints
      do n = 1, number_of_bands_in_manifold_local
        do m = 1, number_of_bands_in_manifold_local
          Amnmat_man(index_manifold,m,n,ik) = Amnmat(m,n,ik) 
!         write(6,*) index_manifold, m, n, ik, Amnmat_man(index_manifold,m,n,ik)
        enddo 
      enddo 
    enddo 

    call writeeig( ispin )
     
    u_matrix = Amnmat

    if( index_manifold .eq. 3) then
      u_matrix = cmplx(0.0_dp,0.0_dp,kind=dp)
      do ik = 1, numkpoints
        do n = 1, manifold_bands_lowdin(1)%nincbands_loc_lowdin
          do m = 1, manifold_bands_lowdin(1)%nincbands_loc_lowdin
            u_matrix(m,n,ik) = Amnmat_man(1,m,n,ik)
          enddo
        enddo
      enddo 

      do ik = 1, numkpoints
        n_reset = 0
        do n = manifold_bands_lowdin(1)%nincbands_loc_lowdin + 1,  &
 &             manifold_bands_lowdin(1)%nincbands_loc_lowdin +     & 
 &             manifold_bands_lowdin(2)%nincbands_loc_lowdin 
          n_reset = n_reset + 1
          m_reset = 0
          do m = manifold_bands_lowdin(1)%nincbands_loc_lowdin + 1,  &
 &               manifold_bands_lowdin(1)%nincbands_loc_lowdin +     & 
 &               manifold_bands_lowdin(2)%nincbands_loc_lowdin 
            m_reset = m_reset + 1
            u_matrix(m,n,ik) = Amnmat_man(2,m_reset,n_reset,ik)
          enddo
        enddo
      enddo 

!!     For debugging
!      do ik = 1, numkpoints
!        do n = 1, number_of_bands_in_manifold_local
!          do m = 1, number_of_bands_in_manifold_local
!           write(6,fmt="(3i5,1x,f24.16,2x,f24.16)")            &
! &              m, n, ik,                                      &
! &              real(u_matrix(m,n,ik)),aimag(u_matrix(m,n,ik))
!          enddo
!        enddo
!      enddo
!!     End debugging

    endif

!   Compute the matrix elements of the plane wave,
!   for all the wave vectors that connect a given k-point to its nearest
!   neighbours
    call compute_pw_matrix( nncount, bvectorsfrac )

    call mmn( ispin )

    do nkp = 1, numkpoints
       do nn = 1, nncount
          do n = 1, number_of_bands_in_manifold_local
            do m = 1, number_of_bands_in_manifold_local
              m_matrix(m,n,nn,nkp) = Mmnkb(m,n,nkp,nn)
            enddo 
          enddo 
       enddo 
    enddo 

    call overlap_project()

    do nkp = 1, numkpoints
       do nn = 1, nncount
          do n = 1, number_of_bands_in_manifold_local
             ! Note that this ln_tmp is defined differently wrt the one in wann_domega
             ln_tmp(n,nn,nkp)=( aimag(log(csheet(n,nn,nkp) &
                     * m_matrix(n,n,nn,nkp))) - sheet(n,nn,nkp) )
          end do
      end do
    end do

    rave  = 0.0_dp
    do iw = 1, number_of_bands_in_manifold_local
       do ind = 1, 3
          do nkp = 1, numkpoints
             do nn = 1, nncount
                rave(ind,iw) = rave(ind,iw) + wb(nn) * bk(ind,nn,nkp) &
                      *ln_tmp(iw,nn,nkp)
             enddo
          enddo
       enddo
    enddo
    rave = -rave/real(numkpoints,dp)

    do iw=1, number_of_bands_in_manifold_local
       write(6,1000) iw,(rave(ind,iw)*lenconfac,ind=1,3)
    end do

1000 format(2x,'WF centre ', &
         &       i5,2x,'(',f10.6,',',f10.6,',',f10.6,' )')

    irpts = 0
    do n1 = -kmeshlowdin(1), kmeshlowdin(1)
      do n2 = -kmeshlowdin(2), kmeshlowdin(2)
        do n3 = -kmeshlowdin(3), kmeshlowdin(3)
          ! Loop over the 125 points R. R=0 corresponds to
          ! i1=i2=i3=0, or icnt=63
          icnt = 0
          do i1 = -2, 2
            do i2 = -2, 2
              do i3 = -2, 2
                icnt = icnt + 1
                ! Calculate distance squared |r-R|^2
                ndiff(1) = n1 - i1 * kmeshlowdin(1)
                ndiff(2) = n2 - i2 * kmeshlowdin(2)
                ndiff(3) = n3 - i3 * kmeshlowdin(3)
                dist(icnt) = 0.0_dp
                do i = 1, 3
                  do j = 1, 3
                    dist(icnt) = dist(icnt) +                &
 &                    real(ndiff(i),dp) * real_metric(i,j) * real(ndiff(j),dp)
                  enddo ! j
                enddo   ! i
              enddo     ! i3
            enddo       ! i2
          enddo         ! i1
          dist_min = minval(dist)
          if (abs(dist(63) - dist_min ) .lt. eps7 ) then
            irpts = irpts + 1
          endif 
        enddo ! n3
      enddo   ! n2
    enddo     ! n1
    nrpts = irpts

    allocate(irvec(3,nrpts))
    irvec=0

    irpts = 0
    do n1 = -kmeshlowdin(1), kmeshlowdin(1)
      do n2 = -kmeshlowdin(2), kmeshlowdin(2)
        do n3 = -kmeshlowdin(3), kmeshlowdin(3)
          ! Loop over the 125 points R. R=0 corresponds to
          ! i1=i2=i3=0, or icnt=63
          icnt = 0
          do i1 = -2, 2
            do i2 = -2, 2
              do i3 = -2, 2
                icnt = icnt + 1
                ! Calculate distance squared |r-R|^2
                ndiff(1) = n1 - i1 * kmeshlowdin(1)
                ndiff(2) = n2 - i2 * kmeshlowdin(2)
                ndiff(3) = n3 - i3 * kmeshlowdin(3)
                dist(icnt) = 0.0_dp
                do i = 1, 3
                  do j = 1, 3
                    dist(icnt) = dist(icnt) +                &
 &                    real(ndiff(i),dp) * real_metric(i,j) * real(ndiff(j),dp)
                  enddo ! j
                enddo   ! i
              enddo     ! i3
            enddo       ! i2
          enddo         ! i1
          dist_min = minval(dist)
          if (abs(dist(63) - dist_min ) .lt. eps7 ) then
            irpts = irpts + 1
            irvec(1,irpts) = n1
            irvec(2,irpts) = n2
            irvec(3,irpts) = n3
          endif 
        enddo ! n3
      enddo   ! n2
    enddo     ! n1

!!   For debugging
!    do irpts = 1, nrpts
!      write(6,'(a,4i5)') ' irpts, irvec = ', irpts, irvec(:,irpts)
!    enddo 
!!   End debugging

!   Compute and write the matrix elements of the position operator
!   Assign a name to the file where the position matrices will be written
    posfilename = trim( seedname )// "_r.dat"

!   Open the output file where the position matrices will be written
    open( unit=posunit, err=1992, file=posfilename, status="replace", &
 &        iostat=eof )

    do irpts = 1, nrpts
      write(posunit,'(/,3i5)') irvec(:,irpts)
      do i = 1, number_of_bands_in_manifold_local
        do j = 1, number_of_bands_in_manifold_local
          delta=0._dp
          if (i==j) delta=1._dp
          pos_r(:)=0._dp
          do ik = 1, numkpoints
            rdotk = twopi*dot_product(kpt_latt(:,ik),real(irvec(:,irpts),dp))
            fac   = exp(-cmplx_i*rdotk)/real(numkpoints,dp)
            do idir=1,3
              do nn=1,nntot
                if(i==j) then
                  ! For irpts==rpt_origin, this reduces to
                  ! Eq.(31) of Marzari and Vanderbilt PRB 56,
                  ! 12847 (1997). Otherwise, is is Eq.(44)
                  ! Wang, Yates, Souza and Vanderbilt PRB 74,
                  ! 195118 (2006), modified according to
                  ! Eqs.(27,29) of Marzari and Vanderbilt
                  pos_r(idir) = pos_r(idir) - &
 &                   wb(nn)*bk(idir,nn,ik)*aimag(log(m_matrix(i,i,nn,ik)))*fac
                else
                  ! Eq.(44) Wang, Yates, Souza and Vanderbilt 
                  ! PRB 74, 195118 (2006)
                  pos_r(idir) = pos_r(idir) + &
 &                 cmplx_i*wb(nn)*bk(idir,nn,ik)*(m_matrix(j,i,nn,ik)-delta)*fac
                endif
              enddo ! nn
            enddo   ! idir
          enddo     ! ik
          write(posunit,'(2i5,3x,6(e15.8,1x))') j, i, pos_r(:)
        enddo       ! j
      enddo         ! i
    enddo           ! irpts

!   Close the output file where the position matrices will be written
    call io_close(posunit)


    deallocate( csheet )
    deallocate( sheet  )
    deallocate( rave   )
    deallocate( ln_tmp )

    return
1992 call die("compute_position: Error writing to _r.dat file")

  end subroutine compute_position

  subroutine define_trial_orbitals( i_manifold )
    use trialorbitalclass,  only: trialorbital  ! Derived type to define the
                                                !    localized trial
                                                !    orbitals
    use atmfuncs,    only: lofio      ! Returns angular momentum number
    use atmfuncs,    only: mofio      ! Returns magnetic quantum number
    use atmfuncs,    only: mofio      ! Returns magnetic quantum number
    use atmfuncs,    only: rcut       ! Returns orbital cutoff radius
    use atmfuncs,    only: orb_gindex ! Returns the global index of a 
                                      !    basis orbital
    use atomlist,    only: iaorb      ! Atomic index of each orbital
    use atomlist,    only: iphorb     ! Orbital index of each orbital 
                                      !    in its atom
    use siesta_geom, only: xa         ! Atomic coordinates
    use siesta_geom, only: isa        ! Species index of each atom

    integer,  intent(in)  :: i_manifold
    integer  :: number_projections
    integer  :: iproj
    integer  :: iorb
    integer  :: ia
    integer  :: iao
    integer  :: is
    integer  :: l
    integer  :: r
    integer  :: m
    real(dp) :: rc
    real(dp) :: zaxis(3)
    real(dp) :: xaxis(3)
    real(dp) :: yaxis(3)
    real(dp) :: zovera

    xaxis(1) = 1.0_dp
    xaxis(2) = 0.0_dp
    xaxis(3) = 0.0_dp

    yaxis(1) = 0.0_dp
    yaxis(2) = 1.0_dp
    yaxis(3) = 0.0_dp

    zaxis(1) = 0.0_dp
    zaxis(2) = 0.0_dp
    zaxis(3) = 1.0_dp

    zovera = 1.0_dp/0.529177_dp

    r = 1

    number_projections = manifold_bands_lowdin(i_manifold)%number_of_bands
    if (allocated(manifold_bands_lowdin(i_manifold)%proj_lowdin)) &
 &     deallocate( manifold_bands_lowdin(i_manifold)%proj_lowdin )
    allocate(manifold_bands_lowdin(i_manifold)%proj_lowdin(number_projections))
    do iproj = 1, number_projections
      iorb = manifold_bands_lowdin(i_manifold)%orbital_indices(iproj)
      ia = iaorb(iorb)        
      is = isa(ia)
      iao = iphorb( iorb )           ! Orbital index within atom
      l  = lofio( is, iao )          ! Orbital's angular mumentum number
      m  = mofio( is, iao )          ! (Real) orbital's magnetic quantum number
      rc = rcut(  is, iao )          ! Orbital's cutoff radius
      manifold_bands_lowdin(i_manifold)%proj_lowdin(iproj)%center = xa(:,ia)
      manifold_bands_lowdin(i_manifold)%proj_lowdin(iproj)%zaxis  = zaxis
      manifold_bands_lowdin(i_manifold)%proj_lowdin(iproj)%xaxis  = xaxis
      manifold_bands_lowdin(i_manifold)%proj_lowdin(iproj)%yaxis  = yaxis
      manifold_bands_lowdin(i_manifold)%proj_lowdin(iproj)%zovera = zovera
      manifold_bands_lowdin(i_manifold)%proj_lowdin(iproj)%r      = r
      manifold_bands_lowdin(i_manifold)%proj_lowdin(iproj)%l      = l
      manifold_bands_lowdin(i_manifold)%proj_lowdin(iproj)%mr     = m
      manifold_bands_lowdin(i_manifold)%proj_lowdin(iproj)%rcut   = rc
      manifold_bands_lowdin(i_manifold)%proj_lowdin(iproj)%lmax   = l
      manifold_bands_lowdin(i_manifold)%proj_lowdin(iproj)%from_basis_orbital =&
 &                                                         .true.
      manifold_bands_lowdin(i_manifold)%proj_lowdin(iproj)%iorb   = iorb
      manifold_bands_lowdin(i_manifold)%proj_lowdin(iproj)%iorb_gindex = &
 &                                     orb_gindex(is,iao) 
    enddo

  end subroutine define_trial_orbitals

endmodule m_lowdin
