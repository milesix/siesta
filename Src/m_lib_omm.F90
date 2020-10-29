! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_lib_omm

use MatrixSwitch
use omm_rand
!use dbcsr_api
!use dbcsr_csr_conversions, only : csr_create_new 


use atomlist,       only : qa, lasto
use fdf,            only : fdf_integer, fdf_boolean, fdf_get
use files,          only : slabel
use parallel,       only : BlockSize, Node, Nodes, ionode, ProcessorY
use precision,      only : dp
use sys,            only : die
#ifdef MPI
use mpi_siesta
#endif

implicit none

!**** PUBLIC ************************************!

public :: omm_min, omm_min_block

!************************************************!

contains

!================================================!
! use the orbital minimization method (OMM) from libOMM to   !
! solve the eigenvalue problem (double precision !
! routine, for Gamma point-only calculations)    !
! not parallel yet
!================================================!
subroutine omm_min(CalcE,PreviousCallDiagon,iscf,istp,nbasis,nspin,h_dim,nhmax,numh,listhptr,listh,d_sparse,eta0,qs,h_sparse,&
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
  real(dp), intent(in) :: eta0(1:2)                            ! chemical potential for Kim functional
  real(dp), intent(in) :: h_sparse(1:nhmax,1:nspin) ! hamiltonian matrix (sparse)
  real(dp), intent(in), optional :: t_sparse(1:nhmax)         ! kinetic energy matrix (sparse)
  real(dp), intent(in) :: s_sparse(1:nhmax)         ! overlap matrix (sparse)

  !**** OUTPUT **********************************!

  real(dp), intent(out) :: d_sparse(1:nhmax,1:nspin) ! (energy-)density matrix (sparse)

  !**** LOCAL ***********************************!

  character(5) :: m_storage
  character(3) :: m_operation

  logical :: init_C, dealloc

#ifdef MPI
  integer :: MPIerror
  integer, save :: ictxt_1D, ictxt_2D, ictxt_c, desc_1D(9), desc_2D(9), desc_c(9), info
  integer :: k, l 
  logical, save :: Use2D
  integer, external :: numroc
#endif
  
  integer :: i, j, io, jo, ind, N_occ, bl, h_dim_loc(2), flavour
  integer, dimension(2) :: dims
  real(dp) :: e_min
  type(matrix), save :: H, S, D_min, T, C_min, C_old
  logical, save :: first_call = .true.
  logical :: found, ReadCoeffs, file_exist
  logical :: new_S, precon
  logical, save :: long_out, WriteCoeffs, io_coeff, C_extrapol
  integer, save :: istp_prev, precon_st, precon_st1, N_occ_loc
  integer, save :: blk_c, blk_h, BlockSize_c, nmax
  real(dp), save :: cg_tol, g_tol, eta, tau
  character(len=100) :: WF_COEFFS_filename

  real(dp), allocatable, save :: h_dense_1D(:,:), s_dense_1D(:,:), d_dense_1D(:,:)
  real(dp), allocatable, save :: h_dense_2D(:,:), s_dense_2D(:,:), d_dense_2D(:,:)
  real(dp), allocatable, save :: t_dense_2D(:,:), t_dense_1D(:,:), c_dense(:,:)
  !**********************************************!

  if (nspin == 1) then
    N_occ = nint(0.5_dp*qs(1))
  else
    N_occ = nint(qs(1))
  end if
  
#ifdef MPI
  m_storage ='pddbc'
  m_operation ='lap'
  if(first_call) Use2D = fdf_boolean('OMM.Use2D',.true.)
#else
  m_storage = 'sdden'
  m_operation = 'ref'
  if(ionode) then
    write(6,*) 'omm_min: Serial version is not available at the moment. Please compile with MPI'
  endif
  call die()
#endif
 
  call timer('libomm',1)
  if(first_call) then
    print'(a, i5, a, i8, a, i8, a, i8)','Node =', Node,' hdim =', h_dim, ' nbasis =', nbasis, ' N_occ =', N_occ
    io_coeff = .false.
    long_out = fdf_boolean('OMM.LongOutput',.true.)
    cg_tol = fdf_get('OMM.RelTol',1.0d-9)
    g_tol = fdf_get('OMM.GTol',1.0d-5)
    precon_st1 = fdf_integer('OMM.PreconFirstStep',-1)
    precon_st = fdf_integer('OMM.Precon',-1)
    tau = fdf_get('OMM.TPreconScale',10.0_dp,'Ry')
    eta = fdf_get('OMM.Eta',0.0_dp,'eV')
    WriteCoeffs=fdf_boolean('OMM.WriteCoeffs',.false.)
    ReadCoeffs=fdf_boolean('OMM.ReadCoeffs',.false.)
    if(WriteCoeffs .or. ReadCoeffs) io_coeff = .true.
    C_extrapol = fdf_boolean('OMM.Extrapolate',.false.)
    nmax = fdf_integer('OMM.MaxIter',100000000)
  end if     

  new_S = .false.
  if(first_call) istp_prev = 0
  if(first_call .or. (istp .ne. istp_prev)) then
    new_S = .true.
    istp_prev = istp
  end if
  init_C = .true.
  if(first_call) init_C = .false.

  if(istp==2 .and. new_S .and. WriteCoeffs) then
    call timer('WriteCoeffs',1)
    WF_COEFFS_filename=trim(slabel)//'.WF_LIBOMM_ST1'
    call m_write(C_min,WF_COEFFS_filename)
    if(ionode) print'(a)', 'File for C is written'
    call timer('WriteCoeffs',2)
  end if

  if(new_S .and. (istp>1) .and. C_extrapol) then
    call timer('c_extrapol',1)
    if(.not. C_old%is_initialized) call m_allocate(C_old,N_occ,h_dim,&
       label=m_storage)

    if(C_extrapol .and. (istp>2)) then  
      call m_add(C_old,'n',C_min,-1.0_dp,2.0_dp,m_operation)
    end if

    call m_copy(C_old,C_min)
    call timer('c_extrapol',2)
  end if

  precon = .false.
  if(present(t_sparse)) then
    if((istp .le. 1) .and. (precon_st1 .ge. iscf)) precon = .true.
    if(precon_st .ge. iscf) precon = .true.
  end if
  
  if(new_S) then
    if(ionode) print'(a)','OMM with libOMM'
    if(ionode) print'(a, i8, a, i8)',' h_dim = ', h_dim, '     N_occ = ', N_occ
  end if

  if(first_call) then
    call timer('m_allocate',1)
    call blacs_get(-1,0,ictxt_1D)
    call blacs_gridinit(ictxt_1D,'C',1,Nodes) 
    if(Use2D) then
      call blacs_get(ictxt_1D,0,ictxt_2D)
      call blacs_gridinit(ictxt_2D,'C',ProcessorY,Nodes/ProcessorY)      
      call ms_scalapack_setup(MPI_Comm_world,ProcessorY,'c',BlockSize,icontxt=ictxt_2D)
      call blacs_gridinfo(ictxt_2D,i,j,k,l)
      h_dim_loc(1) = numroc(h_dim,BlockSize,k,0,processorY)
      h_dim_loc(2) = numroc(h_dim,BlockSize,l,0,Nodes/processorY)
      print'(a, i5, a, i5, a,i5)','Node =', Node,' h_dim1 = ', h_dim_loc(1),' h_dim2 = ', h_dim_loc(2)
      call descinit(desc_2D,h_dim,h_dim,BlockSize,BlockSize,0,0,ictxt_2D,h_dim_loc(1),info)
      allocate(h_dense_2D(1:h_dim_loc(1),1:h_dim_loc(2)))
      allocate(d_dense_2D(1:h_dim_loc(1),1:h_dim_loc(2)))
      allocate(s_dense_2D(1:h_dim_loc(1),1:h_dim_loc(2)))
      if(precon) allocate(t_dense_2D(1:h_dim_loc(1),1:h_dim_loc(2)))
      d_dense_2D(:,:) = 0.0_dp
    else 
      call ms_scalapack_setup(MPI_Comm_world,1,'c',BlockSize,icontxt=ictxt_1D)
    endif
    call descinit(desc_1D,h_dim,h_dim,BlockSize,BlockSize,0,0,ictxt_1D,h_dim,info)
  
    if(info .ne. 0) print'(a)','qerror in descinit'
     
    allocate(h_dense_1D(1:h_dim,1:nbasis))
    allocate(s_dense_1D(1:h_dim,1:nbasis))
    allocate(d_dense_1D(1:h_dim,1:nbasis))  
    if(precon) allocate(t_dense_1D(1:h_dim,1:nbasis)) 
    d_dense_1D(:,:) = 0.0_dp
    if (.not. C_min%is_initialized) call m_allocate(C_min,N_occ,h_dim,m_storage)
    call m_set(C_min,m_operation,0.0_dp,0.0_dp)
    call timer('m_allocate',2)
  end if

   call timer('m_copy',1)
   h_dense_1D(:,:) = 0.0_dp
   if(new_S) s_dense_1D(:,:) = 0.0_dp
   if(precon) t_dense_1D(:,:) = 0.0_dp
   do io = 1, nbasis
     do j = 1, numh(io)
       ind = listhptr(io) + j
       jo = listh(ind)
       h_dense_1D(jo,io) = h_dense_1D(jo,io) + h_sparse(ind,1)
       if(new_S) s_dense_1D(jo,io) = s_dense_1D(jo,io) + s_sparse(ind)
       if(precon) t_dense_1D(jo,io) = t_dense_1D(jo,io) + t_sparse(ind)
     end do
   end do
   if(Use2D) then
      if(new_S) call pdgemr2d(h_dim,h_dim,s_dense_1D,1,1,desc_1D,s_dense_2D,1,1,desc_2D,ictxt_2D)
      if(precon) call pdgemr2d(h_dim,h_dim,t_dense_1D,1,1,desc_1D,t_dense_2D,1,1,desc_2D,ictxt_2D)
      call pdgemr2d(h_dim,h_dim,h_dense_1D,1,1,desc_1D,h_dense_2D,1,1,desc_2D,ictxt_2D)
      d_dense_2D(:,:) = 0.0_dp
   end if
   call timer('m_copy',2)

   if(first_call) then
     if(ReadCoeffs) then
       call timer('ReadCoeffs',1)
       WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_LIBOMM'
       call m_read(C_min,WF_COEFFS_filename,file_exist)
       if(file_exist .and. ionode) print'(a)','File for C is read'
       if(file_exist) init_C = .true.
       call timer('ReadCoeffs',2)
     end if
     call timer('m_register',1)
     if(Use2D) then
       call m_register_pdbc(H,h_dense_2D,desc_2D)
       call m_register_pdbc(S,s_dense_2D,desc_2D)
       call m_register_pdbc(D_min,d_dense_2D,desc_2D)
       if(precon) call m_register_pdbc(T,t_dense_2D,desc_2D)
     else
       call m_register_pdbc(H,h_dense_1D,desc_1D)
       call m_register_pdbc(S,s_dense_1D,desc_1D)
       call m_register_pdbc(D_min,d_dense_1D,desc_1D)
       if(precon) call m_register_pdbc(T,t_dense_1D,desc_1D)
     end if
     call timer('m_register',2)
  end if
  flavour = 0
  if(precon) flavour = 3
  dealloc = .false.

  if(.not. calcE) then
    call timer('omm_density',1)
    call omm(h_dim,N_occ,H,S,new_S,e_min,D_min,.false.,eta,&
      C_min,init_C,T,tau,flavour,nspin,1,cg_tol,g_tol,long_out,dealloc,&
      m_storage,m_operation,nmax)
    call timer('omm_density',2)
    if(ionode) print'(a, f13.7)','e_min = ', e_min
  else
    call timer('omm_energy',1)
    call omm(h_dim,N_occ,H,S,new_S,e_min,D_min,.true.,eta,&
      C_min,init_C,T,tau,flavour,nspin,1,cg_tol,g_tol,long_out,dealloc,&
      m_storage,m_operation,nmax)
    call timer('omm_energy',2)
  end if

  call timer('d_copy',1)
  if(Use2D) then
    call pdgemr2d(h_dim,h_dim,d_dense_2D,1,1,desc_2D,d_dense_1D,1,1,desc_1D,ictxt_2D)
  end if

  do io = 1, nbasis 
    do j = 1, numh(io)
      ind = listhptr(io) + j
      jo = listh(ind) 
      d_sparse(ind,1) = d_dense_1D(jo,io)
      if(nspin == 2) then 
        d_sparse(ind,2) = d_sparse(ind,1)
      else
        d_sparse(ind,1) = 2.0_dp*d_sparse(ind,1)
      end if
    end do
  end do
  call timer('d_copy',2)

  if(.not. CalcE) then
    if(WriteCoeffs) then
      call timer('WriteCoeffs',1)
      WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_LIBOMM'
      call m_write(C_min,WF_COEFFS_filename)
      if(ionode) print'(a)','File for C is written '
      call timer('WriteCoeffs',2)
    end if
  end if
 
  if(first_call) first_call = .false.
  call timer('libomm',2)
end subroutine omm_min

subroutine omm_min_block(CalcE,PreviousCallDiagon,iscf,istp,nbasis,nspin,h_dim,nhmax,numh,listhptr,listh,d_sparse,&
    eta0,qs,h_sparse,s_sparse,t_sparse)
  
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
  real(dp), intent(in) :: eta0(1:2)                            ! chemical potential for Kim functional
  real(dp), intent(in) :: h_sparse(1:nhmax,1:nspin) ! hamiltonian matrix (sparse)
  real(dp), intent(in), optional :: t_sparse(1:nhmax)         ! kinetic energy matrix (sparse)
  real(dp), intent(in) :: s_sparse(1:nhmax)         ! overlap matrix (sparse)

  !**** OUTPUT **********************************!

  real(dp), intent(out) :: d_sparse(1:nhmax,1:nspin) ! (energy-)density matrix (sparse)

  !**** LOCAL ***********************************!

  character(5) :: m_storage
  character(3) :: m_operation

  logical :: init_C, dealloc
  integer :: i, j, ind, io, jo, N_occ, MPIerror, flavour, is
  logical, save :: first_call=.true.
  type(matrix), save :: H, C_min, S, D_min, T, C_old, C_old2
  real(dp) :: e_min, tau, qout(2), qtmp(2)
  real(dp) :: block_data(1,1), c_occ
  integer, save :: istp_prev, BlockSize_c, nmax
  logical :: new_S, ReadCoeffs, file_exist
  logical, save :: long_out, Use2D, WriteCoeffs, C_extrapol
  real(dp), save :: cg_tol, g_tol, eta
  character(len=100) :: WF_COEFFS_filename

#ifdef MPI
  m_storage ='pdcsr'
  m_operation ='lap'
#else
  if(ionode) then
    write(6,*) 'omm_min: The use of DBCSR matrices requires compilation with MPI'
  endif
  call die()  
#endif

  call timer('blomm',1)

  if (nspin == 1) then
    N_occ = nint(0.5_dp*qs(1))
  else
    N_occ = nint(qs(1))
  end if

 
  if(first_call) then
    BlockSize_c = fdf_integer('OMM.BlockSizeC',BlockSize)
    if(ionode) print'(a, i8, a, i8, a, i5, a, i5)','hdim =', h_dim, '    N_occ =', N_occ, &
      '    BlockSize =', BlockSize, '    BlockSize_c =', BlockSize_c
    Use2D = fdf_boolean('OMM.Use2D',.true.)
    C_extrapol = fdf_boolean('OMM.Extrapolate',.false.)
    long_out = fdf_boolean('OMM.LongOutput',.true.)
    cg_tol = fdf_get('OMM.RelTol',1.0d-9)
    g_tol = fdf_get('OMM.GTol',1.0d-3)
    WriteCoeffs=fdf_boolean('OMM.WriteCoeffs',.false.)
    ReadCoeffs=fdf_boolean('OMM.ReadCoeffs',.false.)
    eta = fdf_get('OMM.Eta',0.0_dp,'eV')
    nmax = fdf_integer('OMM.MaxIter',100000000)
  end if

  new_S = .false.
  if(first_call) istp_prev = 0
  if(first_call .or. (istp .ne. istp_prev)) then
    new_S = .true. 
    istp_prev=istp
  end if

  if(new_S .and. istp>1) then
    call timer('c_extrapol',1)
    if(.not. C_old%is_initialized) call m_allocate(C_old,N_occ,h_dim,&
       BlockSize_c,BlockSize,label=m_storage,use2D=Use2D)

    if(C_extrapol .and. istp>2) then
      if(.not. C_old2%is_initialized) call m_allocate(C_old2,N_occ,&
        h_dim,BlockSize_c,BlockSize,label=m_storage,use2D=Use2D)
      call m_copy(C_old2,C_old)
    end if
    call m_copy(C_old,C_min)
    if(C_min%is_initialized) call m_deallocate(C_min)
    call init_c_matrix(C_min, N_occ, h_dim, BlockSize_c, BlockSize, Use2D, .false.)

    if(C_extrapol .and. (istp>2)) then 
      call m_add(C_old,'n',C_old2,2.0_dp,-1.0_dp)
      call m_copy(C_min,C_old2,keep_sparsity=.true.)
    else
      call m_copy(C_min,C_old,keep_sparsity=.true.)
    end if  
    call timer('c_extrapol',2)
  end if

  if(first_call) then
    call timer('m_allocate',1)
    if(ionode) print'(a)','sparse OMM with libOMM'
 
    call ms_dbcsr_setup(MPI_Comm_World) 
    
    if (.not. H%is_initialized) call m_allocate(H,h_dim,h_dim,BlockSize,BlockSize,&
       label=m_storage,use2D=Use2D)
    if (.not. S%is_initialized) call m_allocate(S,h_dim,h_dim,BlockSize,BlockSize,&
        label=m_storage,use2D=Use2D)
    if (.not. D_min%is_initialized) &
      call m_allocate(D_min,h_dim,h_dim,BlockSize,BlockSize,label=m_storage,use2D=Use2D)
    call timer('m_allocate',2)
    call init_c_matrix(C_min, N_occ, h_dim, BlockSize_c, BlockSize, Use2D, .true.)
    if(ReadCoeffs) then
      call timer('ReadCoeffs',1)
      WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_BLOMM'
      call m_read(C_min,WF_COEFFS_filename,file_exist=file_exist, keep_sparsity=.true.)
      if(file_exist) then
        if(ionode) print'(a)','File for C is read'
      end if
      call timer('ReadCoeffs',2)
    end if
    call m_dbcsr_occupation(C_min, c_occ)
    if(ionode) print'(a,f10.8)','C occupation = ', c_occ
    if(WriteCoeffs) then
      call timer('WriteCoeffs',1)
      WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_BLOMM'
      call m_write(C_min,WF_COEFFS_filename)
      if(ionode) print'(a)', 'File for C is written'
      call timer('WriteCoeffs',2)
    end if
  end if

  call timer('m_register',1)
  if(new_S) then    
    call m_register_pdcsr(H, nbasis,BlockSize,listhptr,listh,numh,h_sparse(:,1))
    call m_register_pdcsr(S, nbasis,BlockSize,listhptr,listh,numh,s_sparse)
    d_sparse(:,1) = 0.0_dp
    call m_register_pdcsr(D_min, nbasis,BlockSize,listhptr,listh,numh,d_sparse(:,1))    
   if(ionode) print'(a)','CSR matrix created'
    call m_convert_csrdbcsr(S,threshold=1.0e-14_dp, bl_size=BlockSize)
  end if
  call m_convert_csrdbcsr(H,threshold=1.0e-14_dp, bl_size=BlockSize)
  call timer('m_register',2)

  first_call  = .false.
  init_C = .true.
  dealloc = .false.
  
  flavour = 0
  tau = 0.0_dp

  if(.not. calcE) then
    call timer('omm_density',1)
    call omm(h_dim,N_occ,H,S,new_S,e_min,D_min,.false.,eta,&
      C_min,init_C,T,tau,flavour,nspin,1,cg_tol,g_tol,long_out,dealloc,&
      m_storage,m_operation,nmax)
    call timer('omm_density',2)
    if(ionode) print'(a, f13.7)','e_min = ', e_min
  else
    call timer('m_register',1)
    call m_register_pdcsr(D_min,d_sparse(:,1))
    call timer('m_register',2)
    call timer('omm_energy',1)
    call omm(h_dim,N_occ,H,S,new_S,e_min,D_min,.true.,eta,&
      C_min,init_C,T,tau,flavour,nspin,1,cg_tol,g_tol,long_out,dealloc,&
      m_storage,m_operation,nmax)
    call timer('omm_energy',2)
  end if

  call m_dbcsr_occupation(C_min, c_occ)
  if(ionode) print'(a,f10.8)','C occupation = ', c_occ

  call timer('d_copy',1)
  call m_convert_dbcsrcsr(D_min, bl_size=BlockSize)

  qout(1:2) = 0.0_dp
  do io = 1, nbasis
    do j = 1, numh(io)
      ind = listhptr(io) + j
      jo = listh(ind)
      if(nspin==2) then
        d_sparse(ind,2) = d_sparse(ind,1)
      else
        d_sparse(ind,1) = 2.0_dp * d_sparse(ind,1)
      end if
      do is = 1, nspin
        qout(is) = qout(is) + d_sparse(ind,is) * s_sparse(ind)
      end do
    end do
  end do
  call timer('d_copy',2)

  if(.not. calcE) then
    call timer('d_charge',1)
    qtmp(1:2) = qout(1:2)
    call MPI_AllReduce(qtmp, qout, nspin, MPI_double_precision, MPI_sum, MPI_Comm_World, MPIerror)
    do is = 1, nspin
      if(abs(qout(is)) .gt. 1.d-10) qout(is) = qs(is)/qout(is)
    end do

    do io = 1, nbasis
      do j = 1, numh(io)
        ind = listhptr(io) + j
        jo = listh(ind)
        do is = 1, nspin
          d_sparse(ind,is) = qout(is)*d_sparse(ind,is)
        end do
      end do
    end do
    call timer('d_charge',2)

    if(WriteCoeffs) then
      call timer('WriteCoeffs',1)
      WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_BLOMM'
      call m_write(C_min,WF_COEFFS_filename)
      call timer('WriteCoeffs',2)
      if(ionode) print'(a)', 'File for C is written'
    end if

  end if
  call timer('blomm',2)

  contains
   
  subroutine init_c_matrix(C_min, N_occ, h_dim, BlockSize_c, BlockSize_h, Use2D, set_rand)

    use neighbour,      only : jan, mneighb
    use siesta_geom,    only : na_u, xa, ucell
    use siesta_options, only : rcoor
    use alloc,          only : re_alloc
  
    type(matrix), intent(inout) :: C_min    
    integer, intent(in) :: N_occ, h_dim
    integer, intent(in) :: BlockSize_c, BlockSize_h
    logical, intent(in) :: Use2D, set_rand

    integer :: i, j, k, io, jo, nna, MPIerror    
    integer :: seed, iorb, norb, ja, neib0, index, indexi
    integer :: iwf, iwf1, iwf2, ia
    integer :: nblks, iblks
    integer :: nrows, nze_c
    real(dp), allocatable :: c_loc(:)
    integer, dimension(:), pointer:: id_col_p
    integer, allocatable :: id_col(:), id_row(:), nze_row(:)
    integer, allocatable :: neib(:), iatom(:)
    real(dp) :: rn(2), coef, rmax, rr(3), rrmod, cgval
    logical :: found
    character(5) :: m_storage

    m_storage ='pdcsr'
    call timer('init_c_matrix', 1)
    rcoor = fdf_get('OMM.RcLWF',9.5_dp,'Bohr')
    rmax = 0.0_dp
    do i = -1,1
      do j = -1,1        
        do k = -1,1
          rr(1) = i*ucell(1,1) + j*ucell(1,2) + k*ucell(1,3)
          rr(2) = i*ucell(2,1) + j*ucell(2,2) + k*ucell(2,3)
          rr(3) = i*ucell(3,1) + j*ucell(3,2) + k*ucell(3,3)
          rrmod = sqrt( rr(1)**2 + rr(2)**2 + rr(3)**2 )
          if (rrmod .gt. rmax) rmax = rrmod
        enddo
      enddo
    enddo
    if(ionode) print'(a,f13.7)',    'rmax = ', rmax

    call mneighb(ucell,rcoor,na_u,xa,0,0,nna)
    if(ionode) print'(a,f13.7)',    'rcoor = ', rcoor

   
    coef = 1.0d-2/sqrt(real(h_dim,dp))
    
    allocate(iatom(1:N_occ))
    iwf = 0
    do ia = 1, na_u
      call get_number_of_lwfs_on_atom(ia, indexi)
      do i=1, indexi
        iwf = iwf + 1
        iatom(iwf) = ia
      end do
    end do

    if(set_rand) seed = omm_rand_seed()

    io = 0
    jo = 0
    nze_c = 0

    nblks = N_occ/(Nodes * BlockSize_c)
    iwf1 = nblks * BlockSize_c * Nodes + BlockSize_c * Node + 1
    if(iwf1 .gt. N_occ) then
      nrows = nblks * BlockSize_c
    else    
      iwf2 = nblks * BlockSize_c * Nodes + BlockSize_c * (Node + 1)
      if(iwf2 .gt. N_occ) iwf2 = N_occ  
      nrows = nblks * BlockSize_c + iwf2 - iwf1 + 1
      nblks = nblks  + 1
    end if

    allocate(id_row(1:nrows))
    allocate(nze_row(1:nrows))
    nze_row(:) = 0
    id_row(:) = 0    
    allocate(id_col_p(1:1))
    do iblks = 0, nblks
      iwf1 = iblks * BlockSize_c * Nodes + BlockSize_c * Node + 1
      iwf2 = iblks * BlockSize_c * Nodes + BlockSize_c * (Node + 1)
      if(iwf2 .gt. N_occ) iwf2 = N_occ
      do iwf = iwf1, iwf2
        if((iwf==iwf1) .or. (iatom(iwf) .ne. iatom(iwf-1))) then
          if(allocated(neib)) deallocate(neib)
          if(2.0 * rcoor .lt. rmax) then
            call mneighb(ucell,rcoor,na_u,xa,iatom(iwf),0,nna) !
          else
            nna = na_u
            do j = 1, na_u
              jan(j) = j
            end do
          end if
          allocate(neib(1:nna))
          neib(:) = 0
          index = 0
          do i = 1, nna
            found = .false.
            do j = 1, index
              if(neib(j) == jan(i)) found = .true.
            end do
            if(.not. found) then
              index = index + 1
              neib(index) = jan(j)
            end if
          end do
          nna = index
          do i = 2, nna
            neib0 = neib(i)
            j = i - 1
            do while (neib0 .lt. neib(j))
              neib(j + 1) = neib(j)
              neib(j) = neib0
              j = j - 1
            end do
          end do
        end if

        io = io + 1
        id_row(io) = jo
        do j = 1, nna
          ja = neib(j)
          norb = lasto(ja) - lasto(ja-1) 
          nze_c = nze_c + norb
          nze_row(io) = nze_row(io) + norb
          if(nze_c .gt. size(id_col_p)) call re_alloc(id_col_p, 1, nze_c, 'id_col_p', 'omm_min_block')
          do iorb = 1, norb
            jo = jo + 1
            id_col_p(jo) = lasto(ja-1)+iorb
          end do
        end do      
      end do
    end do 
    if(allocated(neib)) deallocate(neib) 
    allocate(c_loc(1:nze_c))
    allocate(id_col(1:nze_c))
    iwf=1
    do i=1, nze_c
      if(set_rand) then
        do k = 1, 2
          call omm_bsd_lcg(seed, rn(k))
        end do
        cgval = sign(0.5_dp * rn(1), rn(2) - 0.5_dp)
      else
        cgval=1.0_dp
      end if
      c_loc(i) = cgval * coef
      id_col(i) = id_col_p(i)
    end do
    if(allocated(iatom)) deallocate(iatom)
    nullify(id_col_p)

    if (.not. C_min%is_initialized) call m_allocate(C_min,N_occ,h_dim,&
       BlockSize_c,BlockSize_h,label=m_storage,use2D=Use2D)

    call m_register_pdcsr(C_min, nrows, BlockSize_c,&
       id_row,id_col,nze_row,c_loc)

    if(ionode) print'(a)','C_min CSR matrix created'

    call m_convert_csrdbcsr(C_min, threshold=1.0d-14, bl_size=BlockSize_h)
     
    if(ionode) print'(a)','C_min CSR matrix converted to DBCSR'

    call timer('init_c_matrix',2)

  end subroutine init_c_matrix

end subroutine omm_min_block


subroutine get_number_of_lwfs_on_atom(ia, indexi)
  implicit none
  
  integer, intent(in)  :: ia
  integer, intent(out) :: indexi
  
  real(dp) :: tiny
  integer  :: nelectr
  logical, save :: secondodd = .false.
 
  tiny = 1.d-10
  nelectr = qa(ia) + tiny 
 
  if (abs(nelectr - qa(ia) + tiny) .gt. 1e-3) then
    if(ionode) then
      write(6,*) 'omm_min: Wrong atomic charge for atom ',ia
      write(6,*) '      qa = ',qa(ia),' must be an integer'
    endif
    call die()
  endif
  if ( (nelectr / 2) * 2 .ne. nelectr) then
    if (secondodd) then
      indexi = ( ( nelectr - 1 ) / 2 )
      secondodd = .false.
    else
      indexi = ( ( nelectr + 1 ) / 2 )
      secondodd = .true.
    endif
  else
    indexi = ( ( nelectr ) / 2 )
  endif

end subroutine get_number_of_lwfs_on_atom

end module m_lib_omm
