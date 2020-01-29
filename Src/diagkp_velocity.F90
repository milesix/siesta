! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

#ifdef MPI
subroutine diagkp_velocity( spin, no_l, no_u, no_s, nnz, &
    ncol, ptr, col, H, S, getD, &
    fixspin, qtot, qs, temp, e1, e2, xij, &
    nk, kpoint, wk, eo, qo, DM, EDM, ef, efs, &
    Entropy, Haux, Saux, psi, aux, &
    occtol, iscf, neigwanted)
  
  ! *********************************************************************
  ! Subroutine to calculate the eigenvalues and eigenvectors, density
  ! and energy-density matrices, and occupation weights of each 
  ! eigenvector, for given Hamiltonian and Overlap matrices (including
  ! spin polarization). K-sampling version.
  ! This routine also calculates the velocities for the bands
  ! and shifts the eigenvalues according to the 
  ! Created from diagkp, written by Nick Papior, 2018
  ! Uses parallelisation over K points instead of parallelisation 
  ! within them.
  ! **************************** INPUT **********************************
  ! type(t_spin) spin           : Spin type
  ! integer no_l                : Number of basis orbitals in unit cell
  !                               local to this processor
  ! integer no_u                : Number of basis orbitals in unit cell
  ! integer no_s                : Number of basis orbitals in supercell
  ! integer nnz                 : Maximum number of orbitals interacting
  ! integer ncol(no_l)          : Number of nonzero elements of each row 
  !                               of hamiltonian/density matrix locally
  ! integer ptr(no_l)           : Pointer to each row (-1) of the
  !                               hamiltonian matrix locally
  ! integer col(nnz)            : Nonzero hamiltonian-matrix element  
  !                               column indexes for each matrix row
  ! real*8  H(nnz,spin%H)       : Hamiltonian in sparse form
  ! real*8  S(nnz)              : Overlap in sparse form
  ! logical getD                : Find occupations and density matrices?
  ! logical fixspin             : Spin-dependent fermi-levels
  ! real*8  qtot                : Number of electrons in unit cell
  ! real*8  qs                  : Number of electrons in unit cell per spin
  ! real*8  temp                : Electronic temperature 
  ! real*8  e1, e2              : Energy range for density-matrix states
  !                               (to find local density of states)
  !                               Not used if e1 > e2
  ! real*8  xij(3,nnz)          : Vectors between orbital centers (sparse)
  !                               (not used if only gamma point)
  ! integer nk                  : Number of k points
  ! real*8  kpoint(3,nk)        : k point vectors
  ! real*8  wk(nk)              : k point weights (must sum one)
  ! real*8  occtol              : Occupancy threshold for DM build
  ! integer neigwanted          : maximum number of eigenvalues wanted
  ! *************************** OUTPUT **********************************
  ! real*8 eo(maxo,spin%spinor,nk)   : Eigenvalues
  ! ******************** OUTPUT (only if getD=.true.) *******************
  ! real*8 qo(maxo,spin%spinor,nk)   : Occupations of eigenstates
  ! real*8 Dnew(maxnd,spin%DM)    : Output Density Matrix
  ! real*8 Enew(maxnd,spin%EDM)    : Output Energy-Density Matrix
  ! real*8 ef                   : Fermi energy
  ! real*8 Entropy              : Electronic entropy
  ! *************************** AUXILIARY *******************************
  ! real*8 Haux(2,no_u,no_l) : Auxiliary space for the hamiltonian matrix
  ! real*8 Saux(2,no_u,no_l) : Auxiliary space for the overlap matrix
  ! real*8 psi(2,no_u,no_l)  : Auxiliary space for the eigenvectors
  ! real*8 aux(2,no_u*5)     : Extra auxiliary space
  ! *************************** UNITS ***********************************
  ! xij and kpoint must be in reciprocal coordinates of each other.
  ! temp and H must be in the same energy units.
  ! eo, Enew and ef returned in the units of H.
  ! *************************** PARALLEL ********************************
  ! The auxiliary arrays are now no longer symmetric and so the order
  ! of referencing has been changed in several places to reflect this.
  ! Note : It is assumed in a couple of places that the sparsity of
  ! H/S and Dscf are the same (which is the case in the code at present)
  ! and therefore the only one set of pointer arrays are globalised.
  ! *********************************************************************
  !
  !  Modules
  !
  use precision
  use sys
  use units, only: eV, Ang
  use parallel,     only : Node, Nodes
  use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb
  use mpi_siesta
  use m_fermid,     only : fermid, fermispin, stepf
  use alloc,        only : re_alloc, de_alloc
  use t_spin, only: tSpin
  use m_velocity_shift, only: velocity_shift, calc_velocity_current, velocity_current
  use siesta_geom, only: ucell, cell_periodic
  use intrinsic_missing, only: MODP, VNORM

  implicit none

  integer :: MPIerror

  type(tSpin), intent(in) :: spin

  ! Pass system size information
  ! Local, unit-cell, supercell
  integer, intent(in) :: no_l, no_u, no_s

  ! K-point information
  integer, intent(in) :: nk
  real(dp), intent(in) :: kpoint(3,nk), wk(nk)

  ! Now pass sparse patterns
  integer, intent(in) :: nnz
  integer, intent(in) :: ncol(no_l), ptr(no_l), col(nnz)

  ! Matrices
  real(dp), intent(in) :: H(nnz,spin%H), S(nnz), xij(3,nnz)
  real(dp), intent(inout) :: DM(nnz,spin%DM), EDM(nnz,spin%EDM)

  ! Sought charges per spin
  real(dp), intent(in) :: Qs(spin%spinor), Qtot
  ! Fermi-levels
  real(dp), intent(inout) :: Ef, Efs(spin%spinor)
  ! Energy range of DM
  real(dp), intent(in) :: E1, E2
  
  ! Eigenvalues and occupations (charges)
  real(dp), intent(inout) :: qo(no_u,spin%spinor,nk), eo(no_u,spin%spinor,nk)
  real(dp), intent(in) :: Occtol, Temp
  ! Calculated quantities
  real(dp), intent(inout) :: Entropy

  ! Control whether DM/EDM is calculated and for fixed spin
  logical, intent(in) :: getD, fixspin

  ! Current SCF step
  integer, intent(in) :: iSCF

  ! Number of calculated eigenstates
  integer, intent(in) :: neigwanted

  ! Auxiliary arrays
  real(dp), intent(inout), target :: Haux(2,no_u,no_u), Saux(2,no_u,no_u)
  real(dp), intent(inout), target :: psi(2,no_u,no_u), aux(2,no_u*5)

  ! Internal variables
  integer :: BNode, ie, ierror, ik, is
  integer :: io, iio, jo, ind, j, neigneeded
  real(dp) :: ee, qe, kxij, ckxij, skxij, t
  real(dp) :: pipj1, pipj2

  ! Current calculation
  real(dp) :: I
  real(dp), parameter :: Coulomb = 1.6021766208e-19_dp
  real(dp), parameter :: hbar_Rys = 4.8377647940592375e-17_dp

  ! Arrays for figuring out the degenerate states
  real(dp), parameter :: deg_EPS = 7.3498067e-06_dp ! 1e-4 eV
  integer :: ndeg
  real(dp), pointer :: vdeg(:,:) => null(), Identity(:,:) => null()

  ! Globalized matrices
  integer :: g_nnz
  integer, pointer :: g_ncol(:), g_ptr(:), g_col(:)
  real(dp), pointer :: g_H(:,:), g_S(:), g_xij(:,:)
  real(dp), pointer :: g_DM(:,:), g_EDM(:,:)
  real(dp), pointer :: v(:,:)

  ! Not allocated, only pointers
  real(dp), pointer :: Dk(:,:,:), Ek(:,:,:)
  real(dp), pointer :: dH(:,:,:), dS(:,:,:)
  real(dp), pointer :: paux(:,:), eig_aux(:)

  real(dp) :: vcross(3)
  real(dp), external :: volcel

#ifdef DEBUG
  call write_debug( '    PRE diagkp_velocity' )
#endif

  ! Globalise sparsity pattern
  call MPI_AllReduce(nnz, g_nnz, 1, MPI_Integer, MPI_Sum, &
      MPI_Comm_World, MPIerror)

  call diagkp_velocity_2d1d(no_u, aux(1,no_u*2+1), eig_aux)

  ! Nullify arrays
  nullify(g_ncol, g_ptr, g_col, g_H, g_S, g_xij, g_DM, g_EDM)
  nullify(v)
  
  ! Allocate local memory for global list arrays
  call re_alloc( g_ncol, 1, no_u, name='g_ncol', routine= 'diagkp_velocity' )
  call re_alloc( g_ptr, 1, no_u, name='g_ptr', routine= 'diagkp_velocity' )
  call re_alloc( g_col, 1, g_nnz, name='g_col', routine= 'diagkp_velocity' )
  call re_alloc( g_H, 1, g_nnz, 1, spin%spinor, name='g_H', routine= 'diagkp_velocity' )
  call re_alloc( g_S, 1, g_nnz, name='g_S', routine= 'diagkp_velocity' )
  call re_alloc( g_xij, 1, 3, 1, g_nnz, name='g_xij', routine= 'diagkp_velocity' )
  call re_alloc( v, 1, 3, 1, no_u, name='v', routine= 'diagkp_velocity' )

  ! Create pointers for d* arrays
  ! Note that these arrays may then not be inter-used
  dH => Haux
  dS => Saux
  Dk => Haux
  Ek => Saux

  ! Globalize arrays
  g_ptr(1) = 0
  do io = 1, no_u
    
    call WhichNodeOrb(io,Nodes,BNode)
    
    if ( Node == BNode ) then
      call GlobalToLocalOrb(io,Node,Nodes,iio)
      g_ncol(io) = ncol(iio)
      do jo = 1, ncol(iio)
        g_col(g_ptr(io)+jo) = col(ptr(iio)+jo)
        g_H(g_ptr(io)+jo,:) = H(ptr(iio)+jo,:)
        g_S(g_ptr(io)+jo) = S(ptr(iio)+jo)
        g_xij(:,g_ptr(io)+jo) = xij(:,ptr(iio)+jo)
      end do
    end if
    
    call MPI_Bcast(g_ncol(io),1,MPI_Integer, BNode, &
        MPI_Comm_World,MPIerror)
    call MPI_Bcast(g_col(g_ptr(io)+1),g_ncol(io),MPI_Integer, &
        BNode,MPI_Comm_World,MPIerror)
    do is = 1, spin%spinor
      call MPI_Bcast(g_H(g_ptr(io)+1,is),g_ncol(io),MPI_Double_Precision, &
          BNode,MPI_Comm_World,MPIerror)
    end do
    call MPI_Bcast(g_S(g_ptr(io)+1),g_ncol(io),MPI_Double_Precision, &
        BNode,MPI_Comm_World,MPIerror)
    call MPI_Bcast(g_xij(1,g_ptr(io)+1),3*g_ncol(io),MPI_Double_Precision, &
        BNode,MPI_Comm_World,MPIerror)

    ! Update list-pointer
    if ( io < no_u ) g_ptr(io+1) = g_ptr(io) + g_ncol(io)
    
  end do

  ! Perform diagonalization loop
  do ik = 1 + Node, nk, Nodes
    
    do is = 1, spin%spinor
      
      call setup_k(kpoint(:,ik))

      ! Since we want to calculate the velocities as well we do need the eigenstates
      call cdiag(Haux,Saux,no_u,no_u,no_u,eo(1,is,ik),psi, &
          neigwanted,iscf,ierror, -1)

      ! Check error flag and take appropriate action
      if ( ierror > 0 ) then
        call die('Terminating due to failed diagonalisation')
      else if ( ierror < 0 ) then
        call setup_k(kpoint(:,ik))
        call cdiag(Haux,Saux,no_u,no_u,no_u,eo(1,is,ik),psi, &
            neigwanted,iscf,ierror, -1)
      end if

      ! Figure out the degenerate states
      ndeg = 1
      do ie = 2, neigwanted
        
        ! Eigenvalues are sorted, so this will work fine
        if ( eo(ie,is,ik) - eo(ie-1,is,ik) < deg_EPS ) then

          ndeg = ndeg + 1
          
        else if ( ndeg > 1 ) then

          ! Decouple ndeg from ie - 1 and ndeg back
          call degenerate_decouple(neigwanted, eo(1,is,ik), ndeg, ie - 1)

          ! To debug, one may calculate the velocities
          ! before and after the degenerate decoupling
          ndeg = 1
          
        end if
        
      end do

      ! Decouple the last elements
      if ( ndeg > 1 ) &
          call degenerate_decouple(neigwanted, eo(1,is,ik), ndeg, neigwanted)

      ! After having decoupled the degenerate states we can calculate
      ! the velocities
      call calculate_velocity(neigwanted, eo(:,is,ik))

      ! Shift eigenvalues according to the velocity projection onto the field
      call velocity_shift(+1, neigwanted, eo(:,is,ik), v)

    end do
    
  end do

  ! Globalise eigenvalues
  if ( neigwanted /= no_u ) then
    do ik = 1,nk
      do is = 1,spin%spinor
        BNode = mod(ik-1, Nodes)
        call MPI_Bcast(eo(1,is,ik), neigwanted, MPI_Double_Precision, &
            BNode,MPI_Comm_World,MPIerror)
      end do
    end do
  else
    io = no_u * 2
    do ik = 1,nk
      BNode = mod(ik-1, Nodes)
      call MPI_Bcast(eo(1,1,ik), io, MPI_Double_Precision, &
          BNode,MPI_Comm_World,MPIerror)
    end do
  end if

  ! Check if we are done ................................................
  if ( .not. getD ) goto 999

  ! Find new Fermi energy and occupation weights ........................
  if ( fixspin ) then
    call fermispin( spin%spinor, spin%spinor, nk, wk, no_u, &
        neigwanted, eo, Temp, qs, qo, efs, Entropy )
  else
    call fermid( spin%spinor, spin%spinor, nk, wk, no_u, &
        neigwanted, eo, Temp, qtot, qo, ef, Entropy )
    ! This just means we can use efs independently on the used method
    efs(:) = ef
  end if

  ! Initialize current
  I = 0._dp

  ! Allocate globalized DM and EDM
  call re_alloc( g_DM, 1, g_nnz, 1, spin%DM, name='g_DM', routine= 'diagkp_velocity' )
  call re_alloc( g_EDM, 1, g_nnz, 1, spin%EDM, name='g_EDM', routine= 'diagkp_velocity' )

!$OMP parallel default(shared), private(t,ik,is,io)

  ! Find weights for local density of states ............................
  if ( e1 < e2 ) then
    
    t = max( temp, 1.d-6 )
!$OMP do collapse(3)
    do ik = 1,nk
      do is = 1,spin%spinor
        do io = 1,neigwanted
          qo(io,is,ik) = wk(ik) * &
              ( stepf((eo(io,is,ik)-e2)/t) - &
              stepf((eo(io,is,ik)-e1)/t)) * 2.0d0 / spin%spinor
        end do
      end do
    end do
!$OMP end do nowait

  end if

  ! Initialize to 0
!$OMP workshare
  g_DM(:,:) = 0._dp
  g_EDM(:,:) = 0._dp
!$OMP end workshare nowait

!$OMP end parallel

  do ik = 1 + Node, nk, Nodes
    do is = 1, spin%spinor

      ! Find maximum eigenvector that is required for this k point and spin
      ! Note that since eo has averaged out the degeneracy eigenvalues this below
      ! block will also group *all* degenerate eigenvalues!
      neigneeded = 1
      do ie = neigwanted, 1, -1
        if ( abs(qo(ie,is,ik)) > occtol ) then
          neigneeded = ie
          exit
        end if
      end do

      ! Find eigenvectors
      call setup_k(kpoint(:,ik))
      call cdiag(Haux,Saux,no_u,no_u,no_u,eig_aux,psi,neigneeded,iscf,ierror, -1)

      ! Check error flag and take appropriate action
      if ( ierror > 0 ) then
        call die('Terminating due to failed diagonalisation')
      else if ( ierror < 0 ) then
        call setup_k(kpoint(:,ik))
        call cdiag(Haux,Saux,no_u,no_u,no_u,eig_aux,psi,neigneeded,iscf,ierror, -1)
      end if

      ! Before expanding eigenvectors we need to decouple the degenerate
      ! states so that we don't populate the degeneracies wrongly
      ! This assumes that the decoupling is stable.
      ndeg = 1
      do ie = 2, neigneeded
        
        ! Eigenvalues are sorted, so this will work fine
        if ( eig_aux(ie) - eig_aux(ie-1) < deg_EPS ) then

          ndeg = ndeg + 1
          
        else if ( ndeg > 1 ) then

          ! Decouple ndeg from ie - 1 and ndeg back
          call degenerate_decouple(neigneeded, eig_aux, ndeg, ie - 1)

          ! To debug, one may calculate the velocities
          ! before and after the degenerate decoupling
          ndeg = 1
          
        end if
        
      end do

      ! Decouple the last elements
      if ( ndeg > 1 ) &
          call degenerate_decouple(neigneeded, eig_aux, ndeg, neigneeded)

      ! If we want to calculate the current and print it, we can do so here
      if ( calc_velocity_current ) then
        call calculate_velocity(neigneeded, eig_aux)
        call velocity_current(neigneeded, eig_aux, v, Efs(is), wk(ik), Temp, I)
      end if

      ! Expand the eigenvectors to the density matrix
      
!$OMP parallel default(shared), &
!$OMP&private(ie,qe,ee,paux,io,jo,pipj1,pipj2,j,ind), &
!$OMP&private(kxij,ckxij,skxij)

      ! Add contribution to density matrices of unit-cell orbitals
      
!$OMP workshare
      Dk = 0._dp
      Ek = 0._dp
!$OMP end workshare

      ! Global operation to form new density matrix
      do ie = 1, neigneeded
        
        qe = qo(ie,is,ik)
        ee = qe * eig_aux(ie)
        
        ! Point to the wavefunction
        paux => psi(:,:,ie)
        
!$OMP do collapse(2)
        do io = 1, no_u
          do jo = 1, no_u
            pipj1 = paux(1,io) * paux(1,jo) + paux(2,io) * paux(2,jo)
            pipj2 = paux(1,io) * paux(2,jo) - paux(2,io) * paux(1,jo)
            Dk(1,jo,io) = Dk(1,jo,io) + qe * pipj1
            Dk(2,jo,io) = Dk(2,jo,io) + qe * pipj2
            Ek(1,jo,io) = Ek(1,jo,io) + ee * pipj1
            Ek(2,jo,io) = Ek(2,jo,io) + ee * pipj2
          end do
        end do
!$OMP end do

      end do
      
!$OMP do
      do io = 1, no_u
        do j = 1, g_ncol(io)
          ind = g_ptr(io) + j
          jo = modp(g_col(ind), no_u)
          kxij = kpoint(1,ik) * g_xij(1,ind) + kpoint(2,ik) * g_xij(2,ind) + kpoint(3,ik) * g_xij(3,ind)
          ckxij = cos(kxij)
          skxij = sin(kxij)
          
          g_DM(ind,is) = g_DM(ind,is) + Dk(1,jo,io)*ckxij - Dk(2,jo,io)*skxij
          g_EDM(ind,is) = g_EDM(ind,is) + Ek(1,jo,io)*ckxij - Ek(2,jo,io)*skxij
          
        end do
      end do
!$OMP end do
    
!$OMP end parallel

    end do

  end do

  ! Now do bcast of the respective arrays

  do io = 1, no_u
    
    call WhichNodeOrb(io,Nodes,BNode)
    call GlobalToLocalOrb(io,BNode,Nodes,iio)

    if ( BNode == Node ) then
      
      do is = 1, spin%spinor
        call MPI_Reduce(g_DM(g_ptr(io)+1,is), &
            DM(ptr(iio)+1,is),g_ncol(io),MPI_double_precision, &
            MPI_sum,BNode,MPI_Comm_World,MPIerror)
        call MPI_Reduce(g_EDM(g_ptr(io)+1,is), &
            EDM(ptr(iio)+1,is),g_ncol(io),MPI_double_precision, &
            MPI_sum,BNode,MPI_Comm_World,MPIerror)
      end do

    else

      ! This is because this is *NOT* the receiving node
      ! and hence the recv_buff is not important
      do is = 1, spin%spinor
        call MPI_Reduce(g_DM(g_ptr(io)+1,is), &
            aux(1,1),g_ncol(io),MPI_double_precision, &
            MPI_sum,BNode,MPI_Comm_World,MPIerror)
        call MPI_Reduce(g_EDM(g_ptr(io)+1,is), &
            aux(1,1),g_ncol(io),MPI_double_precision, &
            MPI_sum,BNode,MPI_Comm_World,MPIerror)
      end do
      
    end if
  end do

  call de_alloc( g_DM, name='g_DM', routine= 'diagkp_velocity' )
  call de_alloc( g_EDM, name='g_EDM', routine= 'diagkp_velocity' )

  if ( calc_velocity_current ) then
    call MPI_Reduce(I, ee, 1, MPI_Double_Precision, &
        MPI_sum,0,MPI_Comm_World,MPIerror)
    ! Current I is in [e Bohr Ry]
    I = ee / hbar_Rys
    ! Now I is in [e Bohr/s]
    ! Now figure out whether we are in a 1D, 2D or 3D system
    ! and convert to micro-amps
    I = I * Coulomb * 1.e6_dp * 2._dp / spin%spinor

    if ( Node == 0 ) then
      select case ( count( cell_periodic(:) ) )
      case ( 3 )
        ! We are dealing with the volume
        I = I / volcel(ucell) * Ang ** 2
        write(*,'(a,e14.6,tr1,a)') 'Bulk current: ', I, ' uA/Ang^2'
      case ( 2 )
        if ( .not. cell_periodic(1) ) then
          call cross(ucell(:, 2), ucell(:, 3), vcross)
        else if ( .not. cell_periodic(2) ) then
          call cross(ucell(:, 1), ucell(:, 3), vcross)
        else if ( .not. cell_periodic(3) ) then
          call cross(ucell(:, 1), ucell(:, 2), vcross)
        end if
        I = I / VNORM(vcross) * Ang
        write(*,'(a,e14.6,tr1,a)') 'Bulk current: ', I, ' uA/Ang'
      case ( 1 )
        if ( cell_periodic(1) ) then
          I = I / VNORM(ucell(:,1))
        else if ( cell_periodic(2) ) then
          I = I / VNORM(ucell(:,2))
        else if ( cell_periodic(3) ) then
          I = I / VNORM(ucell(:,3))
        end if
        write(*,'(a,e14.6,tr1,a)') 'Bulk current: ', I, ' uA'
      end select
    end if
  end if

  ! Exit point 
999 continue

  ! Clean-up the arrays required for degenerate decoupling
  call de_alloc(vdeg, name='vdeg', routine= 'diagkp_velocity' )
  call de_alloc(Identity, name='identity', routine= 'diagkp_velocity' )

  ! Clean-up the velocities
  call de_alloc( v, name='v', routine= 'diagkp_velocity' )

  ! Clean up memory
  call de_alloc( g_ncol, name='g_ncol', routine= 'diagkp_velocity' )
  call de_alloc( g_ptr, name='g_ptr', routine= 'diagkp_velocity' )
  call de_alloc( g_col, name='g_col', routine= 'diagkp_velocity' )
  call de_alloc( g_H, name='g_H', routine= 'diagkp_velocity' )
  call de_alloc( g_S, name='g_S', routine= 'diagkp_velocity' )
  call de_alloc( g_xij, name='g_xij', routine= 'diagkp_velocity' )

#ifdef DEBUG
  call write_debug( '    POS diagkp_velocity' )
#endif

contains
  
  subroutine setup_k(k)
    real(dp), intent(in) :: k(3)

!$OMP parallel default(shared), private(io,jo,j,ind,kxij,ckxij,skxij)

!$OMP do
    do io = 1,no_u
      Saux(:,:,io) = 0._dp
      Haux(:,:,io) = 0._dp
      do j = 1, g_ncol(io)
        ind = g_ptr(io) + j
        jo = modp(g_col(ind), no_u)
        kxij = k(1) * g_xij(1,ind) + k(2) * g_xij(2,ind) + k(3) * g_xij(3,ind)
        
        ckxij = cos(kxij)
        skxij = sin(kxij)
        
        ! Note : sign of complex part changed to match change in order of io/jo
        Saux(1,jo,io) = Saux(1,jo,io) + g_S(ind) * ckxij
        Saux(2,jo,io) = Saux(2,jo,io) - g_S(ind) * skxij
        Haux(1,jo,io) = Haux(1,jo,io) + g_H(ind,is) * ckxij
        Haux(2,jo,io) = Haux(2,jo,io) - g_H(ind,is) * skxij
        
      end do
    end do
!$OMP end do nowait

!$OMP end parallel

  end subroutine setup_k
  
  subroutine setup_dk(ix, k)
    integer, intent(in) :: ix
    real(dp), intent(in) :: k(3)

!$OMP parallel default(shared), private(io,jo,j,ind,kxij,ckxij,skxij)

!$OMP do
    do io = 1, no_u
      dS(:,:,io) = 0._dp
      dH(:,:,io) = 0._dp
      do j = 1, g_ncol(io)

        ind = g_ptr(io) + j
        jo = modp(g_col(ind), no_u)
        kxij = k(1) * g_xij(1,ind) + k(2) * g_xij(2,ind) + k(3) * g_xij(3,ind)
        
        ckxij = cos(kxij) * g_xij(ix,ind)
        skxij = sin(kxij) * g_xij(ix,ind)
        
        ! Note : sign of complex part changed to match change in order of io/jo
        ! Since the phases are
        !    exp(-i k xij)
        ! we get:
        !  - i xij exp(-i k xij)
        dS(1,jo,io) = dS(1,jo,io) - g_S(ind) * skxij
        dS(2,jo,io) = dS(2,jo,io) - g_S(ind) * ckxij
        dH(1,jo,io) = dH(1,jo,io) - g_H(ind,is) * skxij
        dH(2,jo,io) = dH(2,jo,io) - g_H(ind,is) * ckxij
        
      end do
    end do
!$OMP end do nowait

!$OMP end parallel

  end subroutine setup_dk

  ! This subroutine decouples consecutive degenerate states
  !
  ! NOTE this routine may not use the global variable ie, is or ik
  subroutine degenerate_decouple(neig, eig, ndeg, io_end)
    integer, intent(in) :: neig
    real(dp), intent(inout) :: eig(neig)
    integer, intent(in) :: ndeg, io_end
    integer :: io_start, ix
    real(dp) :: e_avg
    
    call re_alloc(vdeg, 1, 2, 1, ndeg**2, name='vdeg', routine= 'diagkp_velocity', &
        copy=.false., shrink=.false.)
    call re_alloc(Identity, 1, 2, 1, ndeg**2, name='identity', routine= 'diagkp_velocity', &
        copy=.false., shrink=.false.)

    ! Figure out the start and end of the orbitals
    io_start = io_end - ndeg + 1

    ! Calculate the average eigenvalue
    e_avg = 0._dp
    do io = io_start, io_end
      e_avg = e_avg + eig(io)
    end do
    
    ! Average the eigenvalues
    e_avg = e_avg / ndeg

    ! Copy them back to the energy so that they have a unified energy
    ! A decoupling forces them to be aligned since we "mix" the eigenvalues
    do io = io_start, io_end
      eig(io) = e_avg
    end do

    do ix = 1, 3
      
      call setup_dk(ix, kpoint(:,ik))

      ! Calculate < psi_i | dH - e dS | psi_j >
      Identity(:,1:ndeg**2) = 0._dp
      do io = 1, ndeg
        jo = io_start + io - 1

        ! dH | psi_i >
        call zgemv('N', no_u, no_u, dcmplx(1._dp, 0._dp), &
            dH, no_u, psi(1,1,jo), 1, dcmplx(0._dp, 0._dp), aux, 1)
        ! dH - e dS | psi_i >
        call zgemv('N', no_u, no_u, dcmplx(-e_avg, 0._dp), &
            dS, no_u, psi(1,1,jo), 1, dcmplx(1._dp, 0._dp), aux, 1)

        ! Calculte < psi_: | dH - e dS | psi_i >
        call zgemv('C', no_u, ndeg, dcmplx(1._dp, 0._dp), &
            psi(1,1,io_start), no_u, aux, 1, dcmplx(0._dp, 0._dp), vdeg(1,(io-1)*ndeg+1), 1)

        ! Initialize identity matrix
        Identity(1,(io-1)*ndeg+io) = 1._dp

      end do

      ! Now we have the < psi_j | dH - e dS | psi_i > matrix, we need to diagonalize to find
      ! the linear combinations of the new states
      call cdiag(vdeg, Identity, ndeg, ndeg, ndeg, aux, Haux, ndeg, 1, ierror, -1)
      if ( ierror > 0 ) then
        call die('Terminating due to failed diagonalisation in degenerate sub-space')
      end if
      
      ! Re-create the new states
      ! Now Haux contains the linear combination of the states that
      ! should decouple the degenerate subspace
      call zgemm('N', 'N', no_u, ndeg, ndeg, dcmplx(1._dp, 0._dp), &
          psi(1,1,io_start), no_u, Haux, ndeg, dcmplx(0._dp, 0._dp), Saux, no_u)

      ! Copy back the new states
      call zcopy(ndeg*no_u, Saux, 1, psi(1,1,io_start), 1)
      
    end do

  end subroutine degenerate_decouple

  ! This routine calculates the velocities of all the states
  ! It uses the variable ie, so don't use this outside
  subroutine calculate_velocity(neig, eig)
    integer, intent(in) :: neig
    real(dp), intent(in) :: eig(neig)
    integer :: ix, ie
    
    do ix = 1, 3
      
      ! Calculate velocities using the analytic d H(k) / d k_alpha where
      ! alpha is the Cartesian direction
      call setup_dk(ix, kpoint(:,ik))
      
      do ie = 1, neig

        ! dH | psi >
        call zgemv('N', no_u, no_u, dcmplx(1._dp, 0._dp), &
            dH, no_u, psi(1,1,ie), 1, dcmplx(0._dp, 0._dp), aux, 1)
        ! dH - e dS | psi >
        call zgemv('N', no_u, no_u, dcmplx(-eig(ie), 0._dp), &
            dS, no_u, psi(1,1,ie), 1, dcmplx(1._dp, 0._dp), aux, 1)

        ! Now calculate the velocity along ix
        ! < psi | H - e dS | psi >
        v(ix,ie) = 0._dp
        do io = 1, no_u
          v(ix,ie) = v(ix,ie) + psi(1,io,ie) * aux(1,io) + psi(2,io,ie) * aux(2,io)
        end do

      end do
      
    end do

!!$    do ie = 1, neigwanted
!!$      print '(i5,3(tr1,f10.4))', ik, v(:,ie) / Ang / hbar_Rys * 1.e-12_dp
!!$    end do

  end subroutine calculate_velocity

  subroutine diagkp_velocity_2d1d(n,from,to)
    integer, intent(in) :: n
    real(dp), intent(in), target :: from(2*n)
    real(dp), pointer :: to(:)
    to => from(:)
  end subroutine diagkp_velocity_2d1d

end subroutine diagkp_velocity


#else

subroutine diagkp_velocity_dummy()
end subroutine diagkp_velocity_dummy

#endif
