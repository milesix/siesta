! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
subroutine diagonalizeHk( ispin )
!
! Solves H(k)c(k) = e(k)S(k)c(k), with the SCF hamiltonian and
! overlap matrices in k points of the Wannier90 Monkhorst-Pack.
!
! The output are the coefficients of the wavefunctions at the Wannier90 mesh.
! Stores the coefficients in coeffs(orb,band,kpoint)
! It also computes the eigenvalues
! WARNING: these arrays are relabeled in the band indices because
! user specified bands are excluded.
!
!
! Notes about the use of psi and coeffs:
! the array psi is allocated with a size equal to 2 * no_u * no_l
! 2 is for the complex nature of the coefficients (we have to store the
! real and the imaginary parts).
! no_u comes from the fact that for every band we have as many coefficients
! as basis orbitals in the unit cell.
! no_l is due to the fact that each node knows the coefficients of no_l bands.
! When running in serial, no_l = no_u.
! When running in parallel, no_l < no_u, and sum_{nodes} no_l = no_u.
! In diagonalizeHk.F90, psi is a one dimensional array.
! 
! In diagpol, as in cdiag, this array has three dimensions,
! (2, nuotot=no_u, nuo=no_l).
! The first index refers to the real or imaginary part.
! The second index refers to the atomic orbital.
! The third index is the band index.
! 
! Since only a few bands might be used for wannierization, 
! and those bands might not be the 
! numincbands (number of included bands)lowest bands,
! a first reordering of the coefficients is done in the subroutine reordpsi.
! After calling it, the bands whose number ranges from 1 to numincbands(ispin)
! correspond to those bands that will be included for wannierization.
! Only the bands considered for wannierization are included in coeffs.
!
! An example: bulk BaTiO3 (five atoms per unit cell).
! With the DZP basis set plus semicore in valence, there are no_u = 72 bands per
! k-point.
! Those bands are distributed among the different Nodes (for instance,
! if we run with 4 nodes, each node knows the coefficients for no_l = 18 bands:
! Node 0: bands from  1 to 18
! Node 1: bands from 19 to 36
! Node 2: bands from 37 to 54
! Node 3: bands from 55 to 72
! Since bulk BaTiO3 has 40 electrons per unit cell and is non spin-polarized
! material, there are nocc = 20 occupied bands: 
! 18 stored on Node 0, 
! 2 stored on Node 1. 
!
! Now, might be that, for wannierization, we are only interested in the 
! top of the valence bands that come mostly from the O2p orbitals:
! the nine bands (3 p orbitals x 3 O atoms in the unit cell) that 
! range between 12 and 20.
!
! In the matrix of coefficients of the wave functions we want to store
! only the coefficients of those bands and, in the final matrix those
! band should be numbered between 1 and 9. 
! The maximum number of bands considered for wannierization is included
! in the variable nincbands
! 
! Moreover, those bands should be distributed among all the available Nodes,
! so for the matrix of coefficients:
! Node 0: bands 12, 13, and 20 (this Node takes care of three bands)
! Node 1: bands 14 and 15      (this Node takes care of two bands)
! Node 2: bands 16 and 17      (this Node takes care of two bands)
! Node 3: bands 18 and 19      (this Node takes care of two bands)
!
! The number of bands that each Node takes care is defined in the
! variable nincbands_loc
!
! ------------------------------ INPUT -----------------------------------------
! integer ispin        : Spin component
! ------------------------------ OUTPUT ----------------------------------------
! Coefficients of the wavefunctions at the Wannier90 mesh. 
! ------------------------------------------------------------------------------
! Used module variables
  use precision,          only: dp           ! Real douple precision type 
  use densematrix,        only: Haux         ! Hamiltonian matrix in dense form
  use densematrix,        only: Saux         ! Overlap matrix in dense form
  use densematrix,        only: psi          ! Coefficients of the wave function
                                             ! NOTE: while running in parallel,
                                             ! each core knows about the 
                                             ! wave functions of no_l bands
  use m_spin,             only: nspin        ! Number of spin components
  use atomlist,           only: no_l         ! Number of orbitals in local node
                                             ! NOTE: When running in parallel,
                                             !   this is core dependent
                                             !   Sum_{cores} no_l = no_u
  use atomlist,           only: no_u         ! Number of orbitals in unit cell
                                             ! NOTE: When running in parallel,
                                             !   this is core independent
  use atomlist,           only: no_s         ! Number of orbitals in supercell
                                             ! NOTE: When running in parallel,
                                             !   this is core independent
  use atomlist,           only: indxuo       ! Index of equivalent orbital in 
                                             !   the unit cell
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
  use sparse_matrices,    only: H            ! Hamiltonian matrix in sparse form
  use sparse_matrices,    only: S            ! Overlap matrix in sparse form
  use sparse_matrices,    only: xijo         ! Vectors between orbital centers
  use sparse_matrices,    only: tight_binding_param ! Parameters of the 
                                             !   tight-binding Hamiltonian
  use units,              only: eV           ! Conversion factor from Ry to eV
  use m_siesta2wannier90, only: numkpoints   ! Total number of k-points
                                             !   for which the overlap of the
                                             !   periodic part of the wavefunct
                                             !   with a neighbour k-point will
                                             !   be computed
  use m_siesta2wannier90, only: kpointsfrac  ! List of k points relative to the 
                                             !   reciprocal lattice vectors.
                                             !   First  index: component
                                             !   Second index: k-point index 
                                             !      in the list
  use m_siesta2wannier90, only: numbands     ! Number of bands for wannierizatio
                                             !   before excluding bands 
  use m_siesta2wannier90, only: numincbands  ! Number of bands for wannierizatio
                                             !   after excluding bands 
  use m_siesta2wannier90, only: nincbands_loc! Number of bands for wannierizatio
                                             !   after excluding bands 
                                             !   in the local node
  use m_siesta2wannier90, only: blocksizeincbands ! Maximum number of bands
                                             !   considered for wannierization
                                             !   per node
  use m_siesta2wannier90, only: isexcluded   ! Masks excluded bands
  use m_siesta2wannier90, only: coeffs       ! Coefficients of the wavefunctions
                                             !   First  index: orbital
                                             !   Second index: band
                                             !   Third  index: k-point
  use m_siesta2wannier90, only: eo           ! Eigenvalues of the Hamiltonian 
                                             !   at the numkpoints introduced in
                                             !   kpointsfrac 
                                             !   First  index: band index
                                             !   Second index: k-point index
!
! Allocation/Deallocation routines
!
  use alloc,              only: re_alloc     ! Reallocation routines
  use alloc,              only: de_alloc     ! Deallocation routines
  use parallelsubs,       only: LocalToGlobalOrb

#ifdef MPI
  use parallelsubs,       only : set_blocksizedefault
! 
! Subroutine to order the indices of the different bands after 
! excluding some of them for wannierization 
! 
  use m_orderbands,       only: order_index
#endif

! For debugging
  use parallel,           only: Node, Nodes, IOnode, BlockSize
! End debugging


  implicit none

  integer,intent(in)    :: ispin

! Internal variables
  integer  :: ik            ! Counter for loops on k-points
  integer  :: iband         ! Counter for loops on bands
  integer  :: iuo           ! Counter for loops on atomic orbitals
  integer  :: juo           ! Counter for loops on atomic orbitals
  integer  :: j             ! Counter for loops on atomic orbitals
  integer  :: jo            ! Counter for loops on atomic orbitals
  integer  :: ind           ! Counter for loops on atomic orbitals
  integer  :: io            ! Counter for loops on atomic orbitals
  integer  :: n             ! Counter for loops on Nodes
  integer  :: nincbands     ! Total number of included bands for wannierization 
  real(dp) :: kvector(3)    ! k-point vector where the Hamiltonian
                            !    will be diagonalized
  integer  :: nhs           ! Variable to dimension the Hamiltonian and Overlap
  integer  :: npsi          ! Variable to dimension the coefficient vector
  integer  :: nsave         ! Variable to dimension the coefficient vector
  integer  :: ierror        ! Code for error message from cdiag
  integer  :: BNode         ! 
  integer  :: BTest         ! 
  integer  :: ie            ! 
  integer  :: iie            ! 

  real(dp), dimension(:), pointer :: epsilon ! Eigenvalues of the Hamiltonian
  real(dp), dimension(:), pointer :: epsilonstilde 
                                             ! Eigenvalues of the Overlap matrix
! jjunquer
  integer  :: iio   
  real(dp) :: ckxij
  real(dp) :: skxij
  real(dp) :: kxij
  real(dp) :: qp1
  real(dp) :: qp2
  real(dp) :: eqp1
  real(dp) :: eqp2
  real(dp) :: ee
  complex(dp), pointer, save ::   coeffsnew(:,:)
  complex(dp), pointer, save ::   overlap_sq(:,:)
  complex(dp), pointer, save ::   normalization(:,:)
  complex(dp), pointer, save ::   phitilde(:,:)
  complex(dp), pointer, save ::   overlaptilde(:,:)
  complex(dp), pointer, save ::   overlapaux(:,:)
  complex(dp), pointer, save ::   invsqrtover(:,:)
  complex(dp), pointer, save ::   invsqrtd(:,:)
  complex(dp), pointer, save ::   coeffshatphi(:,:)
  real(dp), pointer, save ::   tight_binding_param_k(:,:,:)
  complex(dp), pointer, save ::   aux(:)
  complex(dp), pointer, save ::   lpsi(:)
  logical, allocatable :: done_juo(:)
! end jjunquer

#ifdef MPI
  integer, external :: numroc
#endif
  external cdiag

  call timer('diagonalizeHk',1)

! jjunquer
  allocate(done_juo(no_u))

  nullify( lpsi )
  call re_alloc(lpsi, 1, no_u,                         &
 &              name='lpsi', routine='diagonalizeHk',  &
 &              shrink=.false., copy=.false.)   
  lpsi = cmplx(0.0_dp, 0.0_dp, kind=dp)

  nullify( tight_binding_param )
  call re_alloc(tight_binding_param, 1, maxnh, 1, nspin,              &
 &              name='tight_binding_param', routine='diagonalizeHk',  &
 &              shrink=.false., copy=.false.)   
  tight_binding_param = 0.0_dp

  nullify( coeffsnew )
  call re_alloc(coeffsnew, 1, no_u, 1, no_u,        &
 &              name='coeffsnew', routine='diagonalizeHk',  &
 &              shrink=.false., copy=.false.)   
  coeffsnew = cmplx(0.0_dp, 0.0_dp, kind=dp)

  nullify( overlap_sq )
  call re_alloc(overlap_sq, 1, no_u, 1, no_u,       &
 &              name='overlap_sq', routine='diagonalizeHk', &
 &              shrink=.false., copy=.false.)
  overlap_sq = cmplx(0.0_dp, 0.0_dp, kind=dp)

  nullify( normalization )
  call re_alloc(normalization, 1, no_u, 1, no_u,       &
 &              name='normalization', routine='diagonalizeHk', &
 &              shrink=.false., copy=.false.)
  normalization = cmplx(0.0_dp, 0.0_dp, kind=dp)

  nullify( phitilde )
  call re_alloc(phitilde, 1, no_u, 1, no_u,       &
 &              name='phitilde', routine='diagonalizeHk', &
 &              shrink=.false., copy=.false.)
  phitilde = cmplx(0.0_dp, 0.0_dp, kind=dp)

  nullify( overlaptilde )
  call re_alloc(overlaptilde, 1, no_u, 1, no_u,       &
 &              name='overlaptilde', routine='diagonalizeHk', &
 &              shrink=.false., copy=.false.)
  overlaptilde = cmplx(0.0_dp, 0.0_dp, kind=dp)

  nullify( overlapaux )
  call re_alloc(overlapaux, 1, no_u, 1, no_u,       &
 &              name='overlapaux', routine='diagonalizeHk', &
 &              shrink=.false., copy=.false.)
  overlapaux = cmplx(0.0_dp, 0.0_dp, kind=dp)

  nullify( aux )
  call re_alloc(aux, 1, 5*no_u,       &
 &              name='aux', routine='diagonalizeHk', &
 &              shrink=.false., copy=.false.)
  aux = cmplx(0.0_dp, 0.0_dp, kind=dp)

  nullify( invsqrtover )
  call re_alloc(invsqrtover, 1, no_u, 1, no_u,       &
 &              name='invsqrtover', routine='diagonalizeHk', &
 &              shrink=.false., copy=.false.)
  invsqrtover = cmplx(0.0_dp, 0.0_dp, kind=dp)

  nullify( coeffshatphi )
  call re_alloc(coeffshatphi, 1, no_u, 1, no_u,       &
 &              name='coeffshatphi', routine='diagonalizeHk', &
 &              shrink=.false., copy=.false.)
  coeffshatphi = cmplx(0.0_dp, 0.0_dp, kind=dp)

  nullify( invsqrtd )
  call re_alloc(invsqrtd, 1, no_u, 1, no_u,       &
 &              name='invsqrtd', routine='diagonalizeHk', &
 &              shrink=.false., copy=.false.)
  invsqrtd = cmplx(0.0_dp, 0.0_dp, kind=dp)

  nullify( tight_binding_param_k )
  call re_alloc(tight_binding_param_k, 1, 2, 1, no_u, 1, no_l,         &
 &              name='tight_binding_param_k', routine='diagonalizeHk', &
 &              shrink=.false., copy=.false.)
  tight_binding_param_k = 0.0_dp

! end jjunquer

! Initialize the number of occupied bands
  nincbands = numincbands( ispin )

! Allocate memory related with the dense matrices that will be used in
! the diagonalization routines.
! (Hamiltonian, Overlap and Coefficient vectors)
! These matrices are defined in the module dense matrix
  nhs  = 2 * no_u * no_l
  npsi = 2 * no_u * no_l

  nullify( Saux, Haux, psi )
  call re_alloc( Haux,     1, nhs,   name='Haux',    routine='diagonalizeHk' )
  call re_alloc( Saux,     1, nhs,   name='Saux',    routine='diagonalizeHk' )
  call re_alloc( psi,      1, npsi,  name='psi',     routine='diagonalizeHk' )

! Allocate memory related with the eigenvalues of the Hamiltonian (epsilon)
! and with a local variable where the coefficients of the eigenvector at the
! k-point will be stored
  nullify( epsilon )
  call re_alloc( epsilon, 1, no_u, name='epsilon', routine='diagonalizeHk' )

  nullify( epsilonstilde )
  call re_alloc( epsilonstilde, 1, no_u, name='epsilonstilde', routine='diagonalizeHk' )

! Allocate memory related with the coefficients of the wavefunctions
#ifdef MPI
! Find the number of included bands for Wannierization that will be stored 
! per node. Use a block-cyclic distribution of nincbands over Nodes.
!
     call set_blocksizedefault(Nodes,nincbands,blocksizeincbands)

!     write(6,'(a,3i5)')' diagonalizeHk: Node, Blocksize = ', &
! &                                      Node, BlockSizeincbands

     nincbands_loc = numroc(nincbands,blocksizeincbands,node,0,nodes)
#else
     nincbands_loc = nincbands
#endif
  nullify( coeffs )
  call re_alloc( coeffs,             &
 &               1, no_u,            &
 &               1, nincbands_loc,   &
 &               1, numkpoints,      &
 &               'coeffs',           &
 &               'diagonalizeHk' )


! Allocate memory related with the eigenvalues of the Hamiltonian
  nullify( eo )
  call re_alloc( eo,                 &
 &               1, nincbands,       &
 &               1, numkpoints,      &
 &               'eo',               &
 &               'diagonalizeHk' )

! Initialise psi
  do io = 1, npsi
    psi(io)     = 0.0_dp
  enddo


#ifdef MPI
! Set up the arrays that control the indices of the bands to be 
! considered after excluding some of them for wannierization
! This is done once and for all the k-points
  call order_index( no_l, no_u, nincbands )
#endif

!
! Solve for eigenvectors of H(k) for the k's given in the .nnkp
!
! Loop on k-points
kpoints:                                                             &  
  do ik = 1, numkpoints

!   Compute the wave vector in bohr^-1 for every vector in the list
!   (done in the subroutine getkvector).
!   Remember that kpointsfrac are read from the .nnkp file in reduced units, 
!   so we have to multiply then by the reciprocal lattice vector.
    call getkvector( kpointsfrac(:,ik), kvector )

!   Diagonalize the Hamiltonian for the k-point.
!   Here, we obtain $\psi_{n} (\vec{k})$, where n runs between 1 and no_l
    call diagpol( ispin, nspin, no_l, no_s, no_u,                             &
 &                maxnh, numh, listhptr, listh, H, S, xijo, indxuo, kvector,  &
 &                epsilon, psi, 2, Haux, Saux )

! jjunquer
    write(6,*)' ispin, ik    = ', ispin, ik
!   Coeffsnew is a complex square matrix of dimension:
!   (number_of_atomic_orbitals_in_the_unit_cell)^2.
!   It has two indices:
!   The first one refers to the atomic orbital 
!   The second one refers to the band index
!   
    coeffsnew = cmplx(0.0_dp,0.0_dp,kind=dp)
    ind = 0
    do iband = 1, no_l
      do io = 1, no_u
        coeffsnew(io,iband) =                                  &
 &         cmplx( psi(ind+1), psi(ind+2), kind=dp )
        ind = ind + 2
!        write(6,*)' iband, io, eo, psi = ',                    &
! &                  iband, io, eo(iband,ik),                   &
! &                  coeffsnew(io, iband)
      enddo
    enddo

! 
!   In overlap_sq we store the overlap matrix between two Bloch orbitals
!   that enters in the diagonalization routines
!   It is computed inside diagpol
!   This is a complex square matrix of dimension
!   (number_of_atomic_orbitals_in_the_unit_cell)^2.
!
    overlap_sq     = cmplx(0.0_dp,0.0_dp,kind=dp)
    ind = 0
    do iuo = 1, no_u
      do juo = 1, no_u
        overlap_sq(juo,iuo) =                                  &
 &         cmplx(Saux(ind+1),Saux(ind+2),kind=dp)
        ind = ind + 2
      enddo
    enddo
! end jjunquer


!!   NOTE OF CAUTION: beware when comparing the coefficients of the
!!   wave function obtained in differente machines, specially for 
!!   degenerate states. 
!!   There is a phase that is arbitrary and might change from one machine
!!   to the other. Also, any linear combination of eigenvectors with 
!!   the same eigenvalue is also a solution of the Hamiltonian, nincbands
!!   and the coefficients of the linear combination might be different.


!   Store the eigenvalues, while skipping the excluded bands
    eo(1:numincbands(ispin),ik) = pack(epsilon/eV,.not.isexcluded)

!
!   Keep only the coefficients of the included bands for wannierization,
!   while skipping the excluded eigenvectors.
!   The coefficients of the included eigenvectors will be stored in coeffs,
!   bands whose band index ranges from 1 to nincbands correspond
!   to the bands included for wannierization.
!   
    call reordpsi( coeffs(1:no_u,1:nincbands_loc,ik), psi, no_l, &
                   no_u, numbands(ispin), nincbands_loc )

! jjunquer
!
!   For debugging: check the normalization of the eigenfunctions
!   that come out of the diagonalization 
    normalization = cmplx(0.0_dp,0.0_dp,kind=dp)
    normalization = matmul( overlap_sq, coeffsnew )

    normalization = matmul( transpose(conjg(coeffsnew)),normalization)

!    write(6,'(a)') ' Normalization: '
!    do iuo = 1, no_u
!      do juo = 1, no_u
!        write(6,'(2i5,2f12.5)') iuo, juo, normalization(iuo,juo)
!      enddo
!    enddo
!   End debugging: check the normalization of the eigenfunctions

!
!   Define phitilde in the basis of Blochs made of numerical atomic orbitals
!   It is a complex square matrix of dimension
!   (number_of_atomic_orbitals_in_the_unit_cell)^2.
!   And this is independent of number of bands considered
!

    phitilde = cmplx(0.0_dp,0.0_dp,kind=dp)
    phitilde = matmul( transpose(conjg(coeffsnew)), overlap_sq )
    phitilde = matmul( coeffsnew, phitilde )

!    do iuo = 1, no_u
!      do juo = 1, no_u
!        write(6,'(2i5,2f12.5)') iuo, juo, phitilde(iuo,juo)
!      enddo
!    enddo

!
!   Define the overlap between phitildes
!   It is a complex square matrix of dimension
!   (number_of_atomic_orbitals_in_the_unit_cell)^2.
!   And this is independent of number of bands considered.
!

    overlaptilde = cmplx(0.0_dp,0.0_dp,kind=dp)
    overlaptilde = matmul( overlap_sq, coeffsnew )
    overlaptilde = matmul( overlaptilde, transpose(conjg(coeffsnew)))
    overlaptilde = matmul( overlaptilde, overlap_sq )

    overlapaux = cmplx(0.0_dp,0.0_dp,kind=dp)
    do iuo = 1, no_u
      overlapaux(iuo,iuo) = cmplx(1.0_dp,0.0_dp,kind=dp)
!      do juo = 1, no_u
!        write(6,'(2i5,6f12.5)') iuo, juo, overlaptilde(iuo,juo), overlap_sq(iuo,juo), overlapaux(iuo,juo)
!      enddo
    enddo

!
!   Diagonalize overlaptilde.
!   We need it to compute the inverse of the root square
!

    call cdiag( overlaptilde, overlapaux, no_u, no_l, no_u, epsilonstilde, &
 &              psi, no_u, 1, ierror, BlockSize )


    overlapaux = matmul( conjg(transpose(overlaptilde)), overlap_sq )
    overlapaux = matmul( overlapaux, overlaptilde )

!    do iuo = 1, no_u
!      write(6,*)epsilonstilde(iuo)
!      do juo = 1, no_u
!        write(6,'(2i5,6f12.5)') iuo, juo, overlaptilde(iuo,juo), overlap_sq(iuo,juo), overlapaux(iuo,juo)
!      enddo
!    enddo

    invsqrtd = cmplx( 0.0_dp, 0.0_dp, kind=dp )
    do iuo = 1, no_u
      invsqrtd(iuo, iuo) = cmplx( epsilonstilde(iuo)**(-1.0_dp/2.0_dp), 0.0_dp, kind=dp )
    enddo

    invsqrtover = matmul( overlaptilde, invsqrtd )
    invsqrtover = matmul( invsqrtover, conjg(transpose(overlaptilde)) )

    normalization = matmul( invsqrtover, overlap_sq ) 
    normalization = matmul( normalization, invsqrtover ) 

!    do iuo = 1, no_u
!      write(6,*)epsilonstilde(iuo)
!      do juo = 1, no_u
!        write(6,'(2i5,6f12.5)') iuo, juo, normalization(iuo,juo)
!      enddo
!    enddo

!   
!   coeffshatphi is a complex square matrix of dimension
!   (number_of_atomic_orbitals_in_the_unit_cell)^2.
!   the first index refers to a band
!   the second index refer to a basis function
!
!   For debugging:
!   If we want the tight-binding Hamiltonian to be the same as the 
!   Hamiltonian in real space between SIESTA atomic orbitals
!
    invsqrtover = cmplx( 0.0_dp, 0.0_dp, kind=dp )
    do iuo = 1, no_u
      invsqrtover(iuo, iuo) = cmplx( 1.0_dp, 0.0_dp, kind=dp )
    enddo

    coeffshatphi = matmul( transpose(conjg(coeffsnew)), overlap_sq )
    coeffshatphi = matmul( coeffshatphi, invsqrtover )

    call timer( 'tight-binding-compute', 1 )

!   Add contribution to tight-binding-parameters of unit-cell orbitals
!   WARNING: Dk and Ek may be EQUIVALENCE'd to Haux and Saux
!$OMP parallel do default(shared), collapse(2)
!$OMP&private(iuo,juo)
    do iuo = 1, no_l
      do juo = 1, no_u
        tight_binding_param_k(1,juo,iuo) = 0.0_dp
        tight_binding_param_k(2,juo,iuo) = 0.0_dp
      enddo
    enddo
!$OMP end parallel do

! Global operation to form new density matrix
    BNode = 0
    iie = 0
    do ie = 1, no_u  ! Loop over bands
      if ( Node == BNode ) iie = iie + 1

      if ( Node == BNode ) then
         lpsi => coeffshatphi(iie,:)
      else
         lpsi => aux(:)
      endif
#ifdef MPI
      call MPI_Bcast(lpsi(1),no_u,MPI_double_complex,
     &               BNode,MPI_Comm_World,MPIerror)
#endif

      ee = epsilon(ie)

!$OMP parallel do default(shared),
!$OMP&private(iuo,iio,done_juo,ind,juo,qp1,qp2,eqp1,eqp2)
      do iuo = 1, no_l
        call LocalToGlobalOrb(iuo,Node,Nodes,iio)

        qp1 = ee * realpart(lpsi(iio))
        qp2 = ee * imagpart(lpsi(iio))

       ! Process only the elements really needed
       ! Dk and Ek are re-used storage, so
       ! no memory is really wasted  [AG, Aug 2009]
        done_juo(1:no_u) = .false.
        do ind = listhptr(iuo) + 1, listhptr(iuo) + numh(iuo)
          juo = indxuo(listh(ind))
          if (done_juo(juo)) cycle
          eqp1 = qp1 * realpart(lpsi(juo)) + qp2 * imagpart(lpsi(juo))
          eqp2 = qp1 * imagpart(lpsi(juo)) - qp2 * realpart(lpsi(juo))

          tight_binding_param_k(1,juo,iuo) = tight_binding_param_k(1,juo,iuo) + eqp1 
          tight_binding_param_k(2,juo,iuo) = tight_binding_param_k(2,juo,iuo) + eqp2 
          done_juo(juo) = .true.
        enddo
      enddo
!$OMP end parallel do
      BTest = ie/BlockSize
      if (BTest*BlockSize.eq.ie) then
        BNode = BNode + 1
        if (BNode .gt. Nodes-1) BNode = 0
      endif
    enddo

!$OMP parallel do default(shared),
!$OMP&private(iuo,ind,juo,kxij,ckxij,skxij)
    do iuo = 1,no_l
      do ind = listhptr(iuo) + 1, listhptr(iuo) + numh(iuo)
        juo = indxuo(listh(ind))
        kxij = kvector(1) * xijo(1,ind) +  &
 &             kvector(2) * xijo(2,ind) +  &
 &             kvector(3) * xijo(3,ind)
        kxij = -1.0_dp * kxij
        ckxij = cos(kxij)
        skxij = sin(kxij)
        tight_binding_param(ind,ispin) = tight_binding_param(ind,ispin) + &
 &           tight_binding_param_k(1,juo,iuo)*ckxij -                     &
 &           tight_binding_param_k(2,juo,iuo)*skxij
      enddo
    enddo
!$OMP end parallel do

    call timer( 'tight-binding-compute', 2 )

! end jjunquer


  enddo kpoints

  do iuo = 1, no_u
    do j = 1, numh(iuo)
      ind = listhptr(iuo) + j
      jo  = listh(ind)
      write(6,'(2i5,4f12.5)') iuo, jo,  &
 &       H(ind,:), S(ind), tight_binding_param(ind,ispin)/numkpoints
    enddo 
  enddo 


  call de_alloc( Haux,    name='Haux',         routine='diagonalizeHk' )
  call de_alloc( Saux,    name='Saux',         routine='diagonalizeHk' )
  call de_alloc( psi,     name='psi',          routine='diagonalizeHk' )
  call de_alloc( epsilon, name='epsilon',      routine='diagonalizeHk' )
  call de_alloc( coeffsnew,     name='coeffsnew',      routine='diagonalizeHk' )
  call de_alloc( overlap_sq,    name='overlap_sq',     routine='diagonalizeHk' )
  call de_alloc( normalization, name='normalization',  routine='diagonalizeHk' )
  call de_alloc( phitilde,      name='phitilde',       routine='diagonalizeHk' )
  call de_alloc( epsilonstilde, name='epsilonstilde',  routine='diagonalizeHk' )
  call de_alloc( overlaptilde,  name='overlaptilde',   routine='diagonalizeHk' )
  call de_alloc( overlapaux,    name='overlapaux',     routine='diagonalizeHk' )
  call de_alloc( coeffshatphi,  name='coeffshatphi',   routine='diagonalizeHk' )
  call de_alloc( invsqrtover,   name='invsqrtover',    routine='diagonalizeHk' )
  call de_alloc( invsqrtd,      name='invsqrtd',       routine='diagonalizeHk' )
  call de_alloc( tight_binding_param_k, name='tight_binding_param_k', routine='diagonalizeHk' )
  call de_alloc( lpsi,          name='lpsi',           routine='diagonalizeHk' )
  call de_alloc( tight_binding_param, name='tight_binding_param', routine='diagonalizeHk' )
  call de_alloc( aux,           name='aux',            routine='diagonalizeHk' )
  deallocate(done_juo)

  call timer('diagonalizeHk',2)

end subroutine diagonalizeHk
