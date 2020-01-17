! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_lib_omm

use MatrixSwitch

use fdf,            only : fdf_integer
use parallel,       only : BlockSize, Node, Nodes, ionode
use precision,      only : dp
use sys,            only : die
#ifdef MPI
use mpi_siesta
use parallelsubs,   only : set_BlockSizeDefault
#endif

implicit none

!**** PUBLIC ************************************!

public :: omm_min

!************************************************!

contains

!================================================!
! use the orbital minimization method (OMM) from libOMM to   !
! solve the eigenvalue problem (double precision !
! routine, for Gamma point-only calculations)    !
! not parallel yet
!================================================!
subroutine omm_min(CalcE,PreviousCallDiagon,iscf,istp,nbasis,nspin,h_dim,nhmax,numh,listhptr,listh,d_sparse,eta,qs,h_sparse,&
    s_sparse,t_sparse)
  

  implicit none

  !**** INPUT ***********************************!

  logical, intent(in) :: CalcE              ! Calculate the energy-density matrix from the existing coeffs.?
  logical, intent(in) :: PreviousCallDiagon ! Previous SCF iteration solved by diagonalization?

  integer, intent(in) :: iscf               ! SCF iteration num.
  integer, intent(in) :: istp               ! MD iteration num.
  integer, intent(in) :: nbasis             ! dimension of numh and listhptr
  integer, intent(in) :: nspin              ! num. of spins
  integer, intent(in) :: h_dim              ! num. of AOs (global)
  integer, intent(in) :: nhmax              ! first dimension of listh and sparse matrices
  integer, intent(in) :: numh(1:nbasis)     ! num. of nonzero elements of each row of sparse matrices
  integer, intent(in) :: listhptr(1:nbasis) ! pointer to start of row in listh
  integer, intent(in) :: listh(1:nhmax)     ! list of nonzero elements of each row of sparse matrices

  real(dp), intent(in) :: qs(1:2)                             ! num. of electrons per spin
  real(dp), intent(in) :: eta(1:2)                            ! chemical potential for Kim functional
  real(dp), intent(in), optional :: h_sparse(1:nhmax,1:nspin) ! hamiltonian matrix (sparse)
  real(dp), intent(in), optional :: t_sparse(1:nhmax)         ! kinetic energy matrix (sparse)
  real(dp), intent(in), optional :: s_sparse(1:nhmax)         ! overlap matrix (sparse)

  !**** OUTPUT **********************************!

  real(dp), intent(out) :: d_sparse(1:nhmax,1:nspin) ! (energy-)density matrix (sparse)

  !**** LOCAL ***********************************!

  character(5) :: m_storage
  character(3) :: m_operation

  logical :: init_C

#ifdef MPI
  integer :: MPIerror, mpi_err, BlockSize_c_default, BlockSize_c
  integer :: nbasis_i, nhmax_i
  integer, allocatable, save :: h_dim_l2g(:) ! local-to-global index transform for AOs
  integer, allocatable :: numh_i(:), listhptr_i(:), listh_i(:), h_dim_l2g_i(:)
  real(dp), allocatable :: h_sparse_i(:), s_sparse_i(:), C_old_i(:)
#endif

  integer :: i, j, io, jo, ind, k, l, N_occ
  real(dp) :: he, se, e_min, co, de
  type(matrix) :: H, S, D_min, T, C_min  
  real(dp), allocatable, save :: C_old(:,:)
  logical, save :: first_call = .true.

  !**********************************************!

  if (nspin == 1) then
    N_occ = nint(0.5_dp*qs(1))
  else
    N_occ = nint(qs(1))
  end if
  
#ifdef MPI
! calculate the ScaLAPACK blocking factor for distributing the WF coeffs. matrix
  call set_blocksizedefault(Nodes,N_occ,BlockSize_c_default)
  BlockSize_c=fdf_integer('OMM.BlockSize',BlockSize_c_default)

  if (first_call) call ms_scalapack_setup(mpi_comm_world,1,'c',BlockSize_c)
 
  m_storage ='pddbc'
  m_operation ='lap'
#else
  m_storage = 'sdden'
  m_operation = 'ref'
#endif
  
  if(first_call) then
    allocate(C_old(1:N_occ, 1:nbasis))
    C_old(:,:) = 0.0_dp
    if (ionode) then
      print'(a)','OMM with libOMM'
    end if

#ifdef MPI
    allocate(h_dim_l2g(1:nbasis))
    h_dim_l2g(:) = 0
    j = 0
    k = 0
    l = 0
    do i = 1, h_dim
      k = k+1
      if (j == Node) then
        l = l+1
        h_dim_l2g(l) = i
      end if
      if (k == BlockSize) then
        k = 0
        j = j+1
        if (j == Nodes) j = 0
      end if
    end do
#endif
  end if

  do i = 1, nbasis
!    print'(a,i5,a,i5,a,i5)','Node = ', Node, ' i = ', i, ' l2g = ', h_dim_l2g(i) 
  end do

  if (.not. H%is_initialized) call m_allocate(H,h_dim,h_dim,m_storage)
  if (.not. S%is_initialized) call m_allocate(S,h_dim,h_dim,m_storage)
  if (.not. D_min%is_initialized) call m_allocate(D_min,h_dim,h_dim,m_storage)
  if (.not. C_min%is_initialized) call m_allocate(C_min,N_occ,h_dim,m_storage)

  init_C = .true.  

  call m_set(H,'a',0.0_dp,0.0_dp,m_operation)
  call m_set(S,'a',0.0_dp,0.0_dp,m_operation)

#ifdef MPI
  do i = 0, Nodes - 1
    if(Node == i) then
      nbasis_i = nbasis
      nhmax_i = nhmax
    end if
    call MPI_Bcast(nbasis_i,1,MPI_Integer,i,MPI_Comm_World,MPIerror)
    call MPI_Bcast(nhmax_i,1,MPI_Integer,i,MPI_Comm_World,MPIerror)
    allocate(numh_i(1:nbasis_i))
    allocate(listhptr_i(1:nbasis_i))
    allocate(h_dim_l2g_i(1:nbasis_i))
    allocate(listh_i(1:nhmax_i))
    allocate(h_sparse_i(1:nhmax_i))
    allocate(s_sparse_i(1:nhmax_i))
    if(Node == i) then
      numh_i(:) = numh(:)
      h_dim_l2g_i(:) = h_dim_l2g(:)
      listhptr_i(:) = listhptr(:)
      listh_i(:) = listh(:)
      h_sparse_i(:) = h_sparse(:,1)
      s_sparse_i(:) = s_sparse(:)
    end if
    call MPI_Bcast(numh_i,nbasis_i,MPI_Integer,i,MPI_Comm_World,MPIerror)
    call MPI_Bcast(h_dim_l2g_i,nbasis_i,MPI_Integer,i,MPI_Comm_World,MPIerror)
    call MPI_Bcast(listhptr_i,nbasis_i,MPI_Integer,i,MPI_Comm_World,MPIerror)
    call MPI_Bcast(listh_i,nhmax_i,MPI_Integer,i,MPI_Comm_World,MPIerror)
    call MPI_Bcast(s_sparse_i,nhmax_i,MPI_Double_Precision,i,MPI_Comm_World,MPIerror)  
    call MPI_Bcast(h_sparse_i,nhmax_i,MPI_Double_Precision,i,MPI_Comm_World,MPIerror)
    do io = 1, nbasis_i
      do j = 1, numh_i(io)
        ind = listhptr_i(io)+j
        jo = listh_i(ind)
        he = h_sparse_i(ind)
        se = s_sparse_i(ind)
        call m_set_element(H,jo,h_dim_l2g_i(io),he,0.0_dp, m_operation)
        call m_set_element(S,jo,h_dim_l2g_i(io),se,0.0_dp, m_operation)
      end do
    end do
    deallocate(numh_i)
    deallocate(h_dim_l2g_i)
    deallocate(listh_i)
    deallocate(listhptr_i)
    deallocate(h_sparse_i)
    deallocate(s_sparse_i)
  end do
#else
  do io = 1, nbasis
    do j = 1, numh(io)
      ind = listhptr(io)+j
      jo = listh(ind)
      he = h_sparse(ind, 1)
      se = s_sparse(ind)
      call m_set_element(H,jo,io,he,0.0_dp, m_operation)
      call m_set_element(S,jo,io,se,0.0_dp, m_operation)
    end do
  end do
#endif
  
  call m_set(D_min,m_operation,0.0_dp,0.0_dp)  
  call m_set(C_min,m_operation,0.0_dp,0.0_dp)

#ifdef MPI
  do i = 0, Nodes - 1
    if(Node == i) then
      nbasis_i = nbasis
    end if
    call MPI_Bcast(nbasis_i,1,MPI_Integer,i,MPI_Comm_World,MPIerror)
    allocate(h_dim_l2g_i(1:nbasis_i))
    allocate(C_old_i(1:nbasis_i))
    if(Node == i) h_dim_l2g_i(:) = h_dim_l2g(:)
    call MPI_Bcast(h_dim_l2g_i,nbasis_i,MPI_Integer,i,MPI_Comm_World,MPIerror)
    do io = 1, N_occ
      if(Node == i) C_old_i(:) = C_old(io, :)
      call MPI_Bcast(C_old_i(:),nbasis_i,MPI_Double_Precision,i,MPI_Comm_World,MPIerror)
      do jo = 1, nbasis_i
         call m_set_element(C_min,io,h_dim_l2g_i(jo),C_old_i(jo),0.0_dp)
      end do
    end do
    deallocate(h_dim_l2g_i)
    deallocate(C_old_i)
  end do
#else
  do io = 1, N_occ
    do jo = 1, nbasis
      call m_set_element(C_min,io,jo,C_old(io,jo),0.0_dp)
    end do
  end do
#endif

  if(.not. calcE) then
    if(first_call) init_C = .false.
    if(.not. init_C) print'(a,i5,a)','Node = ', Node, ' init_C = false '
    call omm(h_dim,N_occ,H,S,.true.,e_min,D_min,.false.,0.0_dp,&
      C_min,init_C,T,0.0_dp,0,nspin,1,-1.0_dp,.true.,.true.,&
      m_storage,m_operation)
#ifdef MPI
    do i = 0, Nodes - 1
      if(Node == i) then
        nbasis_i = nbasis
      end if
      call MPI_Bcast(nbasis_i,1,MPI_Integer,i,MPI_Comm_World,MPIerror)
      allocate(h_dim_l2g_i(1:nbasis_i))
      if(Node == i) h_dim_l2g_i(:) = h_dim_l2g(:)
      call MPI_Bcast(h_dim_l2g_i,nbasis_i,MPI_Integer,i,MPI_Comm_World,MPIerror)
      do io = 1, N_occ
        do jo = 1, nbasis_i
           call m_get_element(C_min,io,h_dim_l2g_i(jo),co)
           if(Node == i) C_old(io, jo) = co
        end do
      end do
    deallocate(h_dim_l2g_i)
  end do
#else
  do io = 1, N_occ
    do jo = 1, nbasis
      call m_get_element(C_min,io,jo,C_old(io,jo))
    end do
  end do
#endif
  else
    call omm(h_dim,N_occ,H,S,.true.,e_min,D_min,.false.,0.0_dp,&
      C_min,init_C,T,0.0_dp,0,nspin,1,-1.0_dp,.true.,.false.,&
      m_storage,m_operation)
    call omm(h_dim,N_occ,H,S,.false.,e_min,D_min,.true.,0.0_dp,&
      C_min,init_C,T,0.0_dp,0,nspin,1,-1.0_dp,.true.,.true.,&
      m_storage,m_operation)
  end if	  
  if (ionode) then
    print'(a, f13.7)','e_min = ', e_min
  end if 

#ifdef MPI
  do i = 0, Nodes - 1
    if(Node == i) then
      nbasis_i = nbasis
      nhmax_i = nhmax
    end if
    call MPI_Bcast(nbasis_i,1,MPI_Integer,i,MPI_Comm_World,MPIerror)
    call MPI_Bcast(nhmax_i,1,MPI_Integer,i,MPI_Comm_World,MPIerror)
    allocate(numh_i(1:nbasis_i))
    allocate(listhptr_i(1:nbasis_i))
    allocate(h_dim_l2g_i(1:nbasis_i))
    allocate(listh_i(1:nhmax_i))
    if(Node == i) then
      numh_i(:) = numh(:)
      h_dim_l2g_i(:) = h_dim_l2g(:)
      listhptr_i(:) = listhptr(:)
      listh_i(:) = listh(:)
    end if
    call MPI_Bcast(numh_i,nbasis_i,MPI_Integer,i,MPI_Comm_World,MPIerror)
    call MPI_Bcast(h_dim_l2g_i,nbasis_i,MPI_Integer,i,MPI_Comm_World,MPIerror)
    call MPI_Bcast(listhptr_i,nbasis_i,MPI_Integer,i,MPI_Comm_World,MPIerror)
    call MPI_Bcast(listh_i,nhmax_i,MPI_Integer,i,MPI_Comm_World,MPIerror)
    do io = 1, nbasis_i
      do j = 1, numh_i(io)
        ind = listhptr_i(io) + j
        jo = listh_i(ind)
        call m_get_element(D_min,jo,h_dim_l2g_i(io),de)
        if(i == Node) then
          d_sparse(ind, 1) = 2.0*de 
          if(nspin == 2) d_sparse(ind, 2) = d_sparse(ind, 1)
        end if
      end do
    end do
    deallocate(numh_i)
    deallocate(h_dim_l2g_i)
    deallocate(listh_i)
    deallocate(listhptr_i)
  end do
#else
  do io = 1, nbasis
    do j = 1, numh(io)
      ind = listhptr(io) + j
      jo = listh(ind)
      call m_get_element(D_min,jo,io,de)
      d_sparse(ind, 1) = 2.0*de
      if(nspin == 2) d_sparse(ind, 2) = d_sparse(ind, 1)
    end do
  end do
#endif

  if(first_call) first_call = .false.
   
  call m_deallocate(D_min)
  call m_deallocate(S)
  call m_deallocate(H)
  call m_deallocate(C_min)

end subroutine omm_min

end module m_lib_omm
