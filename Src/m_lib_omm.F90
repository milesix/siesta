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

public :: omm_min_block

!************************************************!

contains

!===============================================================================!
! Use the orbital minimization method (OMM) to solve the eigenvalue problem
! (double precision routine, Gamma point-only calculations).
! Linear scaling is achieved when sparse matrices are used (Ordejon-Mauri
! and Kim method). Cubic-scaling OMM is performed when working with dense
! matrices. The libOMM library is used for conjugate gradient (CG) minimization
! of the energy functional. The MatrixSwitch library serves as an interface
! low-level routines for algebraic operations, which are performed using
! the DBCSR library for sparse matrices and ScaLAPACK for dense ones.
!===============================================================================!

subroutine omm_min_block(CalcE,PreviousCallDiagon,iscf,istp,nbasis,nspin,h_dim,nhmax,numh,listhptr,listh,d_sparse,&
    eta0,qs,h_sparse,s_sparse,t_sparse)

  use siesta_options, only : rcoor
  
  implicit none

  !**** INPUT ***********************************!

  logical, intent(in) :: CalcE              ! Calculate the energy-density matrix from the existing coeffs.?
  logical, intent(in) :: PreviousCallDiagon ! Previous SCF iteration solved by diagonalization?

  integer, intent(in) :: iscf               ! SCF iteration number
  integer, intent(in) :: istp               ! MD iteration number
  integer, intent(in) :: nbasis             ! Dimension of numh and listhptr arrays
                                            ! (the number of local rows in sparse Hamiltonian and overlap matrices)
  integer, intent(in) :: nspin              ! Number of spin components
  integer, intent(in) :: h_dim              ! Global number of atomic orbitals
  integer, intent(in) :: nhmax              ! First dimension of listh and sparse matrices
  integer, intent(in) :: numh(1:nbasis)     ! Number of nonzero elements of each local row of sparse matrices
  integer, intent(in) :: listhptr(1:nbasis) ! Pointers to start of local rows in listh
  integer, intent(in) :: listh(1:nhmax)     ! List of nonzero elements of each local row for sparse matrices

  real(dp), intent(in) :: qs(1:2)                   ! Number of electrons per spin
  real(dp), intent(in) :: eta0(1:2)                 ! Chemical potential for Kim functional
  real(dp), intent(in) :: h_sparse(1:nhmax,1:nspin) ! Hamiltonian matrix (sparse)
  real(dp), intent(in), optional :: t_sparse(1:nhmax)  ! Kinetic energy matrix (sparse)
  real(dp), intent(in) :: s_sparse(1:nhmax)         ! Overlap matrix (sparse)

  !**** OUTPUT **********************************!

  real(dp), intent(out) :: d_sparse(1:nhmax,1:nspin) ! (Energy-)density matrix (sparse)

  !**** LOCAL ***********************************!

  type(matrix) :: csr_mat              ! csr matrix h_dim x h_dim in MatrixSwitch (MS) format
                                       ! with nbasis local rows and block size BlockSize
  type(matrix), allocatable, save :: H(:)     ! Hamiltonian matrix in MS format
  type(matrix), allocatable, save :: C_min(:) ! Coefficient matrix (localized wavefunctions expanded
                                       ! in basis functions) in MS format
  type(matrix), save :: S              ! Overlap matrix in MS format
  type(matrix), allocatable, save :: D_min(:) ! (Energy-)density matrix in MS format
  type(matrix), save :: T              ! Kinetic energy matrix in MS format
  type(matrix), allocatable, save :: C_old(:) ! Coefficient matrix from the previous MD step in MS format
  type(matrix), allocatable, save :: C_old2(:)! Coefficient matrix from the last but one MD step
  type(matrix), save :: brd_mat        ! Matrix h_dim x h_dim in MS format distributed on a 1D MPI grid
                                       ! with nbasis local rows and block size BlockSize

  integer :: i, j, io, jo, MPIerror
  integer :: ind                       ! Index of sparse matrices
  integer :: flavour                   ! Flavour of OMM calculation with dense matrices: 0 for basic,
                                       ! 1 for Cholesky factorization, 3 for preconditioning (see libOMM)
  integer :: is                        ! Spin component
  integer :: nze                       ! Estimated number of nonzero elements in the coefficient matrix
                                       ! C_min when distributed on the 1D MPI grid
  integer, save :: istp_prev           ! MD step at the previous call
  integer, save :: BlockSize_c         ! Block size for wavefunctions (rows of coefficient matrix C_min)
  integer, save :: N_occ               ! Number of occupied states
  integer, save :: wf_dim              ! Number of wavefunctions considered
  integer, save :: precon_st           ! Number of SCF steps at which preconditioning is applied
                                       ! (all for negative values) for dense matrices
  integer, save :: precon_st1          ! Number of SCF steps at which preconditioning is applied
                                       ! at the first MD step for dense matrices
  integer, allocatable, save :: ind_ordered(:) ! indices of sparse matrices ordered in such a way that
                                               ! column indices are in the growing order for each row
  real(dp) :: e_min                    ! Band structure energy at each SCF step
  real(dp) :: qout(2), qtmp(2)         ! Charges for each spin component
  real(dp) :: c_occ                    ! Matrix occupation (fraction of nonempty elements)
  real(dp) :: rcoor0                   ! Cutoff radius for wavefunctions (intermediate variable)
  real(dp) :: rcoor_init               ! Initial cutoff radius for wavefunctions used for initialization
  real(dp), save :: cg_tol             ! Tolerance for CG minimization. The CG minimization is stopped
                                       ! when 2(E_n - E_(n-1))/(E_n + E_(n-1)) becomes less than cg_tol,
                                       ! where E_n is the energy at CG iteration n
  real(dp), save :: eta                ! Chemical potential for Kim method
  real(dp), save :: tau                ! Kinetic energy scale for preconditioning

  logical, save :: sparse              ! Use sparse matrices?
  logical, save :: use_kim             ! Use Kim functional? If false, Ordejon-Mauri functional is used
  logical :: new_S                     ! Is the overlap matrix new?
  logical :: ReadCoeffs                ! Read wavefunctions (C_min matrix) from the file?
  logical :: ReadUseLib                ! Use DBCSR library for reading wavefunctions?
  logical, save :: WriteCoeffs         ! Write wavefunctions (C_min matrix) to the file at each SCF step?
  logical, save :: WriteUseLib         ! Use DBCSR library for writing the restart for wavefunctions?
  logical :: file_exist(1:nspin)       ! Has the file with wavefunctions been found?
  logical :: precon                    ! Apply preconditioning at this SCF step (for dense matrices only)?
  logical, save :: use_cholesky        ! Apply Cholesky factorization (for dense matrices only)?
  logical :: init_C                    ! Are wavefunctions (C_min matrix) initialized?
  logical :: my_order                  ! Are the indices of the csr matrix ordered in such a way that
                                       ! column indices for each row are in the growing order?
  logical, save :: long_out            ! Print detailed information on CG iterations to libOMM.log?
  logical, save :: Use2D               ! Use a 2D MPI grid for distribution of matrix blocks?
  logical, save :: C_extrapol          ! Extrapolate coefficient matrix C_min linearly based on the results
                                       ! of the two previous MD steps?
  logical :: dealloc                   ! Deallocate matrices at each SCF step?
  logical, save :: first_call=.true.   ! Is this the first call to this subroutine?

  character(len=100) :: WF_COEFFS_filename ! Path to the restart file for wavefunctions (C_min matrix)
  character(len=2) :: sfile            ! Spin-related suffix
  character(5) :: m_storage            ! MS format used for matrices
  character(3) :: m_operation          ! MS implementation of the operations performed

#ifndef MPI
  if(ionode) then
    write(6,*) 'omm_min: OMM requires compilation with MPI'
  endif
  call die()  
#endif

  call timer('blomm',1)
  m_operation = 'lap'
 
  if(first_call) then
    if (nspin == 1) then
      N_occ = nint(0.5_dp*qs(1)) ! Number of occupied states
    else
      N_occ = nint(qs(1))
    end if
    if(nspin > 2) then
      if(ionode) write(6,*) 'omm_min: not implemented for spinors '
      call die()
    endif
    ! Setting input parameters
    BlockSize_c = fdf_integer('OMM.BlockSizeC', 0)        ! Block size for wavefunctions
                                                          ! (rows of coefficient matrix C_min)
    sparse = fdf_boolean('OMM.UseSparse', .true.)         ! Use sparse matrices?
    use_kim = fdf_boolean('OMM.UseKimFunctional', .true.) ! Use Kim (or Ordejon-Mauri) functional?
    if(use_kim) then
      if(ionode) print'(a)','Using the Kim functional'
    else
      if(ionode) print'(a)','Using the Ordejon-Mauri functional'
    end if
    Use2D = fdf_boolean('OMM.Use2D', .true.)              ! Distribute matrices on a 2D MPI grid?
    C_extrapol = fdf_boolean('OMM.Extrapolate', .false.)  ! Extrapolate the coefficient C_min based on the
                                                          ! results of two previous MD steps?
    long_out = fdf_boolean('OMM.LongOutput', .true.)      ! Print detailed information on CG iterations to libOMM.log?
    cg_tol = fdf_get('OMM.RelTol', 1.0d-9)                ! Tolerance for CG iterations
    WriteCoeffs=fdf_boolean('OMM.WriteCoeffs', .false.)   ! Write wavefunctions (C_min) to the file?
    ReadCoeffs=fdf_boolean('OMM.ReadCoeffs', .false.)     ! Read wavefunctions (C_min) from the file?
    WriteUseLib=fdf_boolean('OMM.WriteUseLib', .false.)   ! Use DBCSR library for writing wavefunctions?
    ReadUseLib=fdf_boolean('OMM.ReadUseLib', .false.)     ! Use DBCSR library for reading wavefunctions?
    eta = fdf_get('OMM.Eta', 0.0_dp, 'Ry')                ! Chemical potential for Kim method
    if(ionode) print'(a, i15, a, i15, a, i8)','hdim =', h_dim, '    N_occ =', N_occ, &
      '    BlockSize =', BlockSize
    wf_dim = N_occ                                        ! The number of wavefunctions considered is set
                                                          ! to the number of occupied states by default
    precon_st = fdf_integer('OMM.Precon', -1)             ! Number of SCF steps at which preconditioning is applied
                                                          ! (all for negative values) for dense matrices
    precon_st1 = fdf_integer('OMM.PreconFirstStep', precon_st)  ! Number of SCF steps at which preconditioning
                                                          ! at the first MD step for dense matrices
    tau = fdf_get('OMM.TPreconScale', 10.0_dp, 'Ry')      ! Kinetic energy scale for preconditioning
    use_cholesky = .false.                                ! Apply Cholesky factorization (for dense matrices only)?
    if((.not. use_kim) .and. (.not. sparse)) use_cholesky=fdf_boolean('OMM.UseCholesky', .false.)
    rcoor = fdf_get('OMM.RcLWF', 9.5_dp)                  ! Cutoff radius for wavefunctions
    rcoor_init = fdf_get('OMM.RcLWFInit', 0.0_dp)         ! Initial cutoff radius for wavefunctions
  end if

  if(sparse) then
    m_storage ='pdcsr'   ! pdcsr MS format is used for sparse matrices
  else
    m_storage ='pddbc'   ! pddbc MS format is used for dense matrices
  end if

  dealloc = .true.
  new_S = .false.
  if(first_call) istp_prev = 0
  if(first_call .or. (istp .ne. istp_prev)) then
    new_S = .true. 
    istp_prev=istp
  end if

  if(new_S .and. istp>1) then ! If this is a new (but not the first) MD step
    call timer('c_extrapol', 1)
    ! Creating C_old matrix for wavefunctions at the previous MD step
    ! Creating C_old2 matrix for wavefunctions two MD steps before if extrapolation is used
    if(C_extrapol .and. istp > 2) then
      do is = 1, nspin
        call m_copy(C_old2(is), C_old(is)) ! Creating C_old copy
      end do
    end if
    do is = 1, nspin
      if(C_old(is)%is_initialized) call m_deallocate(C_old(is)) ! Deallocating old C_old
      call m_copy(C_old(is), C_min(is)) ! Creating C_min copy
    end do
    call timer('c_extrapol', 2)

    if(sparse) then
      ! Updating the sparsity of the coefficient matrix according to the current system geometry
      ! when sparse matrices are used. Nonempty elements are set to zero
      call init_c_matrix(C_min(1), wf_dim, h_dim, BlockSize_c, BlockSize, &
        use_kim, .false., m_storage)
      if(nspin > 1) then
         do is = 2, nspin
           if(C_min(is)%is_initialized) call m_deallocate(C_min(is))
           call m_copy(C_min(is), C_min(1)) ! Setting the sparsity of the coefficient matrix for
                                            ! the second spin component the same as of the first one
         end do
      end if
    else
      ! Just deallocate the existing coefficient matrix if dense matrices are used
      if(nspin > 1) then
         do is = 2, nspin
           if(C_min(is)%is_initialized) call m_deallocate(C_min(is))
         end do
      end if
    end if

    call timer('c_extrapol', 1)
    if(C_extrapol .and. (istp > 2)) then
      do is = 1, nspin
        ! Extrapolation: C_min = C_old + (C_old - C_old2)
        ! Here C_old2 = 2C_old - C_old2
        call m_add(C_old(is), 'n', C_old2(is), 2.0_dp, -1.0_dp)
        ! Now copying C_old2 to C_min maintaining the sparsity of the latter. That is the elements that
        ! should be now empty are empty. The new nonempty elements are initialized to zero.
        call m_copy(C_min(is), C_old2(is))
        call m_deallocate(C_old2(is)) ! Deallocating C_old2
      end do
    else
      if(.not. first_call) then
        do is = 1, nspin
          ! Copying C_old to C_min maintaining the sparsity of the latter.
          call m_copy(C_min(is), C_old(is))
          call m_deallocate(C_old(is)) ! Deallocating C_old
        end do
      end if
    end if  
    if(sparse) then
      do is = 1, nspin
        call m_occupation(C_min(is), c_occ)  ! Calculating the occupation of C_min matrix
        if(ionode) print'(a, i1, a, f10.8)','C occupation (', is,') = ', c_occ
      end do
    end if
    call timer('c_extrapol', 2)
  end if

  precon = .false. ! Is preconditioning applied at this SCF step?
  if((.not. use_cholesky) .and. present(t_sparse) .and. (.not. sparse)) then
    if(istp .eq. 1) then
      if((precon_st1 .lt. 0) .or. (precon_st1 .ge. iscf)) precon = .true.
    else
      if((precon_st .lt. 0) .or. (precon_st .ge. iscf)) precon = .true.
    end if
  end if


  if(.not. allocated(H)) allocate(H(1:nspin))
  if(.not. allocated(D_min)) allocate(D_min(1:nspin))

  if(first_call) then ! If this is the first call to this subroutine
    if(sparse) then
      if(ionode) print'(a)','sparse OMM with libOMM'
    else
      if(ionode) print'(a)','dense OMM with libOMM'
    end if

    if(sparse) then
      call ms_dbcsr_setup(MPI_Comm_World, BlockSize, Use2D) ! Setting up the DBCSR library for sparse matrices
    else
      if(Use2D) then
        call ms_scalapack_setup(MPI_Comm_world, ProcessorY, 'c', BlockSize) ! Setting up the ScaLAPACK library
                                                                            ! for dense matrices and 2D MPI grid
      else
        call ms_scalapack_setup(MPI_Comm_world, 1, 'c', BlockSize) ! Setting up the ScaLAPACK library
                                                                   ! for dense matrices in the case of 1D MPI grid
      end if
    end if

    if(.not. allocated(C_min)) allocate(C_min(1:nspin))
    if(.not. allocated(C_old)) allocate(C_old(1:nspin))
    if(C_extrapol) then
      if(.not. allocated(C_old2)) allocate(C_old2(1:nspin))
    end if

    if(sparse .and. (rcoor_init .gt. 0.0_dp)) then
      ! Initializing the first spin component for the coefficient matrix.
      ! It is later copied to initialize the second component.

      ! If the initial cutoff radius is defined for sparse matrices
      rcoor0 = rcoor ! Temporarily saving the cutoff radius
      ! Initializing C_old using the initial cutoff radius
      rcoor = rcoor_init
      call init_c_matrix(C_old(1), wf_dim, h_dim, BlockSize_c, BlockSize, use_kim, &
        .true., m_storage)
      ! Computing the sparsity of C_min using the cutoff radius used for futher calculations.
      ! Setting nonempty elements to zero.
      rcoor = rcoor0
      call init_c_matrix(C_min(1), wf_dim, h_dim, BlockSize_c, BlockSize, use_kim, &
        .false., m_storage)
      ! Copying C_old to C_min keeping the sparsity of the latter.
      call m_copy(C_min(1), C_old(1))
      call m_deallocate(C_old(1))
    else
      ! If the initial cutoff radius is not defined, the usual value is used for initialization.
      ! For dense matrices, the cutoff radius (initial or usual one) is also used for initialization
      ! but the sparsity is not maintained in CG iterations.
      if((.not. sparse) .and. (rcoor_init .gt. 0.0_dp)) rcoor=rcoor_init
      call init_c_matrix(C_min(1), wf_dim, h_dim, BlockSize_c, BlockSize, use_kim, &
        .true., m_storage)
    end if

    nze = (nhmax * wf_dim)/h_dim ! Estimated number of nonzero elements in local rows of a coefficient matrix
                                 ! distributed on a 1D MPI grid

    do is = 1, nspin
      file_exist(is) = .false.
    end do

    ! Reading wavefunctions
    if(ReadCoeffs) then
      call timer('ReadCoeffs', 1)
      do is = 1, nspin
        write(sfile,'(a, i1)') '.', is
        WF_COEFFS_filename = trim(slabel)//'.WF_COEFFS_BLOMM'//trim(sfile)
        if(nspin==1) then
          WF_COEFFS_filename = trim(slabel)//'.WF_COEFFS_BLOMM' ! The wavefunctions are read from the file *.WF_COEFFS_BLOMM
        else
          if(is > 1) then
            if(C_min(is)%is_initialized) call m_deallocate(C_min(is))
            call m_copy(C_min(is), C_min(1)) ! Setting the sparsity pattern of the second component of the coefficient
                                             ! matrix the same as of the first one
          end if
        end if
        call m_read(C_min(is), WF_COEFFS_filename, file_exist=file_exist(is),&
          keep_sparsity=.true., use_dbcsrlib=ReadUseLib,nze=nze)
        if(file_exist(is)) then
          if(ionode) print'(a,i1)','File for C is read ', is
        end if
      end do
      call timer('ReadCoeffs', 2)
    end if

    ! Printing the occupation of the coefficient matrix C_min
    if(sparse) then
      call m_occupation(C_min(1), c_occ)
      if(ionode) print'(a,f10.8)','C occupation (1)= ', c_occ
    end if
    if(ionode) print'(a,i15, a, i8)','WF_dim = ', wf_dim, '   BlockSizeC = ', BlockSize_c

    ! Writing the initial guess for the wavefunctions
    if(WriteCoeffs) then
      if((.not. ReadCoeffs) .or. (.not. file_exist(1))) then
        call timer('WriteCoeffs', 1)
        if(nspin == 2) then
          WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_BLOMM.1'
        else
          WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_BLOMM' ! The wavefunctions are read from the file *.WF_COEFFS_BLOMM
        end if
        call m_write(C_min(1),WF_COEFFS_filename,use_dbcsrlib=WriteUseLib,nze=nze)
        if(ionode) print'(a)', 'File for C(1) is written'
        call timer('WriteCoeffs', 2)
      end if
    end if
  end if

  call timer('m_register', 1)
  if(.not. calcE) then
    ! If computing the density matrix, not the energy density
    if(new_S) then
      ! If sparsity has changed, the indices of the csr matrices should be reordered again
      ! (in such a way that the column indices are in the growing order)
      my_order = .false.
      if(sparse) then
        if(allocated(ind_ordered)) then
          if(size(ind_ordered) .lt. nhmax) then
            deallocate(ind_ordered)
            allocate(ind_ordered(nhmax))
          end if
        else
          allocate(ind_ordered(nhmax))
        end if
        ind_ordered(:) = 0
        my_order = .true.
      end if
      ! Passing pointers to the overlap matrix in the csr format to MS and ordering the indices
      ! of the csr matrix using the array ind_ordered
      call m_register_csr(csr_mat, h_dim, h_dim, nbasis, listhptr, listh, numh, s_sparse, &
        ind_ordered=ind_ordered, order=my_order)
      ! Allocating overlap matrix in the MS format
      if (.not. S%is_initialized) call m_allocate(S, h_dim, h_dim, label=m_storage)
      ! Converting the overlap matrix from csr to MS format using m_sp as an intermediate matrix
      ! distributed on a 1D MPI grid
      call m_copy(S, csr_mat, m_sp=brd_mat) ! Converting and creating the intermediate matrix
      ! Removing the pointers to the overlap matrix in the csr format from MS
      call m_deallocate(csr_mat)
      if(precon) then ! If preconditioning is used
        ! Allocating the kinetic energy matrix in the MS format
        if (.not. T%is_initialized) call m_allocate(T, h_dim, h_dim, label=m_storage)
        ! Passing the pointers to the arrays of the kinetic energy matrix in the csr format to MS
        call m_register_csr(csr_mat, h_dim, h_dim, nbasis, listhptr, listh, numh, t_sparse)
        ! Converting the format of the kinetic energy matrix from csr to MS
        call m_copy(T, csr_mat)
        ! Removing the pointers to the kinetic energy matrix in the csr format from MS
        call m_deallocate(csr_mat)
      end if
    end if
    do is = 1, nspin
      if (.not. H(is)%is_initialized) then
        ! Allocating the Hamiltonian matrix
        if(Use2D .or. (.not. sparse)) then
          call m_allocate(H(is), h_dim, h_dim, label=m_storage)
        else
          call m_copy(H(is), S) ! This is needed to reuse the sparsity pattern if
                                ! the overlap and Hamiltonian matrices are distributed on
                                ! a 1D MPI grid
        end if
      end if
      ! Passing pointers to the Hamiltonian matrix in the csr format to MS and using the
      ! previously prepared array of ordered indices ind_ordered
      call m_register_csr(csr_mat, h_dim, h_dim, nbasis, listhptr, listh, numh, h_sparse(:,is), &
        ind_ordered=ind_ordered, order=.false.)
      ! Converting the Hamiltonian matrix from csr to MS format reusing m_sp as an intermediate matrix
      call m_copy(H(is), csr_mat, m_sp=brd_mat)
      ! Removing the pointers to the Hamiltonian matrix in the csr format from MS
      call m_deallocate(csr_mat)
    end do
  end if
  ! Allocating density matrix in the MS format
  do is = 1, nspin
    if (.not. D_min(is)%is_initialized) call m_allocate(D_min(is), h_dim, h_dim, label=m_storage)
  end do
  call timer('m_register',2)

  init_C = .true.
  
  ! Setting up OMM flavours
  flavour = 0
  if(precon) flavour = 3
  if(use_cholesky) flavour = 1

  do is = 1, nspin
    if(.not. calcE) then
      call timer('omm_density', 1)
    else
      call timer('omm_energy', 1)
    end if
    if(first_call .and. (is > 1)) then
      if((.not. ReadCoeffs) .or. (.not. file_exist(is))) then
        if(C_min(is)%is_initialized) call m_deallocate(C_min(is))
        if(ionode) print'(a)','Copy C_min(1)'
        call m_copy(C_min(is), C_min(1))
      end if
    end if

    ! Calling the libOMM library for CG minimization of the energy functional and calculation
    ! of the density matrix
    call omm(h_dim, wf_dim, n_occ, H(is), S, new_S, e_min, D_min(is), calcE, eta,&
      C_min(is), init_C, T, tau, flavour, nspin, is, cg_tol, long_out, dealloc,&
      m_storage, m_operation)
    if(.not. calcE) then
       call timer('omm_density', 2)
    else
       call timer('omm_energy', 2)
    end if
    if((.not. calcE) .and. ionode) print'(a, f20.7)','e_min = ', e_min
    if(.not. calcE) then
      call timer('d_copy', 1)
    else
      call timer('e_copy', 1)
    end if
    ! Passing pointers to the density matrix in the csr format to MS and using the
    ! previously prepared array of ordered indices ind_ordered
    call m_register_csr(csr_mat, h_dim, h_dim, nbasis, listhptr, listh, numh, d_sparse(:,is), &
      ind_ordered=ind_ordered, order=.false.)
    ! Converting the format of the density matrix from MS format to csr
    call m_copy(csr_mat, D_min(is))
    ! Removing the pointers to the density matrix in the csr format from MS
    call m_deallocate(csr_mat)
    if(.not. calcE) then
      call timer('d_copy', 2)
    else
      call timer('e_copy', 2)
    end if
  end do
  ! Deallocating matrices
  do is = 1, nspin
    if(D_min(is)%is_initialized) call m_deallocate(D_min(is))
    if(H(is)%is_initialized) call m_deallocate(H(is))
  end do
  deallocate(D_min)
  deallocate(H)
  if(T%is_initialized) call m_deallocate(T)
  if(calcE) then
    if(S%is_initialized) call m_deallocate(S)
    ! Deallocating the intermediate 1D-distributed matrix used for format conversion of
    ! the overlap and Hamiltonian matrices
    if(brd_mat%is_initialized) call m_deallocate(brd_mat)
  end if

  first_call = .false.

  ! Taking into account spin components
  if(.not. calcE) then
    call timer('d_copy', 1)
  else
    call timer('e_copy', 1)
  end if
  do io = 1, nbasis
    do j = 1, numh(io)
      ind = listhptr(io) + j
      jo = listh(ind)
      if(nspin == 1) then
        d_sparse(ind, 1) = 2.0_dp * d_sparse(ind, 1)
      end if
    end do
  end do
  if(.not. calcE) then
    call timer('d_copy', 2)
  else
    call timer('e_copy', 2)
  end if

  if(.not. calcE) then
    call timer('d_charge', 1)
    ! Computing the charge
    qout(1:2) = 0.0_dp
    do io = 1, nbasis
      do j = 1, numh(io)
        ind = listhptr(io) + j
        jo = listh(ind)
        do is = 1, nspin
          qout(is) = qout(is) + d_sparse(ind, is) * s_sparse(ind)
        end do
      end do
    end do
    qtmp(1:2) = qout(1:2)
    call MPI_AllReduce(qtmp, qout, nspin, MPI_double_precision, MPI_sum, MPI_Comm_World, MPIerror)

    ! Correcting the charge
    do is = 1, nspin
      qout(is) = qs(is)/qout(is)
    end do

    do io = 1, nbasis
      do j = 1, numh(io)
        ind = listhptr(io) + j
        jo = listh(ind)
        do is = 1, nspin
          d_sparse(ind,is) = qout(is) * d_sparse(ind, is)
        end do
      end do
    end do
    call timer('d_charge', 2)
  end if

  ! Writing wavefunctions (coefficient matrix C_min) to the restart file
  if(WriteCoeffs) then
    call timer('WriteCoeffs', 1)
    nze = (nhmax * wf_dim)/h_dim ! Estimated number of nonzero elements in the local rows
                                 ! of the coefficient matrix distributed on a 1D MPI grid
    do is = 1, nspin
      write(sfile,'(a, i1)') '.', is
      WF_COEFFS_filename = trim(slabel)//'.WF_COEFFS_BLOMM'//trim(sfile)
      if(nspin == 1) then
        WF_COEFFS_filename = trim(slabel)//'.WF_COEFFS_BLOMM' ! Writing to the file *.WF_COEFFS_BLOMM
      end if
      call m_write(C_min(is),WF_COEFFS_filename,use_dbcsrlib=WriteUseLib,nze=nze)
    end do
    call timer('WriteCoeffs', 2)
    if(ionode) print'(a)', 'Files for C are written'
  end if

  call timer('blomm', 2)

  contains

!===============================================================================!
! Analysis of sparsity of the coefficient matrix and initialization by random
! numbers.
!===============================================================================!
  subroutine init_c_matrix(C_min, wf_dim, h_dim, BlockSize_c, BlockSize_h, &
    use_kim, set_rand, m_storage)

    use neighbour,      only : jan, mneighb
    use siesta_geom,    only : na_u, xa, ucell
    use alloc,          only : re_alloc
  
    type(matrix), intent(inout) :: C_min   ! Coefficient matrix initialized
    integer, intent(in) :: h_dim           ! Size of the basis of localized atomic orbitals
    integer, intent(inout) :: BlockSize_c  ! Block size for wavefunctions (rows of C_min)
    integer, intent(in) :: BlockSize_h     ! Blocks size for atomic orbitals (columns of C_min)
    logical, intent(in) :: use_kim         ! Is the Kim method used?
    logical, intent(in) :: set_rand        ! Are nonempty matrix elements filled in with random numbers?
    integer, intent(out) :: wf_dim         ! Number of wavefunctions considered (global number of rows in C_min)
    character(5), intent(in) :: m_storage  ! MS format used for the coefficient matrix

    type(matrix) :: C_csr                  ! Coefficient matrix in the csr format

    integer :: i, j, k, io, jo
    integer :: nna                         ! Number of neighbours found through neighbour search with mneighb
    integer :: seed                        ! Seed for generator of random numbers
    integer :: iorb                        ! Atomic orbital index for the atom considered
    integer :: norb                        ! Number of orbitals for the atom considered
    integer :: ja                          ! Neighbour index
    integer ::  index                      ! Neighbour index in the revised list of neighbours excluding
                                           ! repeating atoms
    integer :: indexi                      ! Number of wavefunctions centered on the atom considered
    integer :: nelectr                     ! Atomic charge
    integer :: iwf, iwf1, iwf2             ! Wavefunction indices
    integer :: ia                          ! Atom considered
    integer :: nblks                       ! Number of blocks of rows in the coefficient matrix in the csr format
    integer :: iblks                       ! Row index of the block in the coefficient matrix
    integer :: nrows                       ! Number of local rows in the coefficient matrix
    integer :: nze_c                       ! Number of nonempty elements in the local rows of the coefficient matrix
    integer :: dn                          ! Step at which the size of the array of column indices of nonempty
                                           ! elements of the coefficient matrix in the csr format, id_col_tmp,
                                           ! is increased if needed
    integer :: ncol                        ! Number of defined elements in the array id_col_tmp
    integer :: nmax                        ! Size of the array of column indices of nonempty
                                           ! elements of the coefficient matrix in the csr format, id_col_tmp

    integer, allocatable :: id_col_tmp(:)  ! Array of column indices of nonempty elements of the coefficient
                                           ! matrix in the csr format
    integer, allocatable :: id_col(:)      ! Final array of column indices of nonempty elements of the coefficient
                                           ! matrix in the csr format
    integer, allocatable :: id_row(:)      ! Array with indices of start of rows for the coefficient matrix in the
                                           ! csr format
    integer, allocatable :: nze_row(:)     ! Number of nonempty elements in each row of the coefficient matrix in
                                           ! the csr format
    integer, allocatable :: neib(:)        ! Corrected list of neighbour of the atom considered excluding repeating
                                           ! neighbours
    integer, allocatable :: iatom(:)       ! Array of atoms to which wavefunctions belong to

    real(dp) :: rn(2)                      ! Random numbers
    real(dp) :: coef                       ! Coefficient for initialization of matrix values
    real(dp) :: rmax                       ! Maximal dimension of the simulation cell
    real(dp) :: rr(3), rrmod               ! Variables used to compute simulation cell dimensions
    real(dp) :: cgval                      ! Random value used for coefficient matrix element

    real(dp), allocatable :: c_loc(:)      ! Local nonempty coefficient matrix elements
    real(dp), allocatable :: fact(:)       ! Occupation of states

    logical :: found                       ! Has the neighbour been already included in the neighbour list
    logical :: include_all                 ! Are all atoms included in the neighbour list?
    logical :: set_neib                    ! Is it needed to make a new neighbour list?

    call timer('init_c_matrix', 1)
    ! Calculating the largest dimension of the simulation cell
    rmax = 0.0_dp
    do i = -1,1
      do j = -1,1        
        do k = -1,1
          rr(1) = i*ucell(1, 1) + j*ucell(1, 2) + k*ucell(1, 3)
          rr(2) = i*ucell(2, 1) + j*ucell(2, 2) + k*ucell(2, 3)
          rr(3) = i*ucell(3, 1) + j*ucell(3, 2) + k*ucell(3, 3)
          rrmod = sqrt( rr(1)**2 + rr(2)**2 + rr(3)**2 )
          if (rrmod .gt. rmax) rmax = rrmod
        enddo
      enddo
    enddo
    if(ionode) print'(a, f13.7)',    'rmax = ', rmax
    if(ionode) print'(a, f13.7)',    'rcoor = ', rcoor

    ! If the simulation cell is larger than the cutoff radius, all atoms are neighbours
    ! and all matrix elements of the coefficient matrix are treated as nonempty
    include_all = .true.
    if(2.0 * rcoor .lt. rmax) include_all = .false.
    ! If the simulation cell is smaller than the cutoff radius, the subroutine for
    ! neighbour search is initialized
    if(.not. include_all) call mneighb(ucell, rcoor, na_u, xa, 0, 0, nna)
   
    ! Coefficient for normalization of the coefficient matrix in order to avoid instabilities
    coef = 1.0d-2/sqrt(real(h_dim, dp))

    ! Assigning wavefunctions to atoms depending on their atomic charge.
    ! Computing the total number of wavefunctions (the number of rows of the coefficient matrix)
    wf_dim = 0
    do ia = 1, na_u
      ! Getting the number of wavefunctions, indexi, on atom ia
      call get_number_of_lwfs_on_atom(ia, use_kim, nelectr, indexi)
      wf_dim = wf_dim + indexi ! Total number of wavefunctions
    end do

    ! Setting the block size for wavefunctions (rows of the coefficient matrix) if not defined
    if(BlockSize_c .le. 0) then
      BlockSize_c = wf_dim * BlockSize_h  ! By default, it equals to the block size for basis functions
      BlockSize_c = BlockSize_c/h_dim     ! (columns of the coefficient matrix) multiplied by the ratio
      if(BlockSize_c==0) BlockSize_c=1    ! of the total numbers of wavefunctions and basis functions
    end if

    ! Setting the occupations and assigning atom to each wavefunction
    allocate(iatom(1:wf_dim))
    allocate(fact(1:wf_dim))
    iwf = 0
    do ia = 1, na_u
      ! Getting the number of wavefunctions, indexi, on atom ia depending on its charge
      call get_number_of_lwfs_on_atom(ia, use_kim, nelectr, indexi)
      do i = 1, indexi  ! For each wavefunction
        iwf = iwf + 1
        iatom(iwf) = ia  ! Assigning the atom to which the wavefunction belongs to
        fact(iwf) = 1.0  ! Occupation
      end do
      if(use_kim) then  !! Revise for nspin=2
        if(2*(nelectr/2) .eq. nelectr) then
           fact(iwf) = sqrt(1.d-6)  ! If the Kim method is used, there is an additional wavefunction
         else                       ! on each atom and this state is not occupied or partially occupied
           fact(iwf) = sqrt(0.5)
         end if 
      end if
    end do

    if(set_rand) seed = omm_rand_seed() ! Setting the seed for random number generator

    io = 0
    jo = 0
    nze_c = 0

    ! Computing the number of local rows for the coefficient matrix in the csr format
    ! and the number of blocks for rows
    nblks = wf_dim/(Nodes * BlockSize_c) ! Number of blocks for rows
    iwf1 = nblks * BlockSize_c * Nodes + BlockSize_c * Node + 1
    if(iwf1 .gt. wf_dim) then
      nrows = nblks * BlockSize_c ! Number of local rows
    else    
      iwf2 = nblks * BlockSize_c * Nodes + BlockSize_c * (Node + 1)
      if(iwf2 .gt. wf_dim) iwf2 = wf_dim  
      nrows = nblks * BlockSize_c + iwf2 - iwf1 + 1
      nblks = nblks  + 1
    end if

    ! Step at which the size of the array of column indices of nonempty elements of
    ! the coefficient matrix in the csr format, id_col_tmp, is increased if needed
    dn = 10 * (h_dim/wf_dim) * nrows
    if(dn .lt. 100) dn = 100
    ncol = 0
    nmax = 3 * dn ! The size of the array of column indices
    if(ionode) print'(a,i8)', 'dn  = ', dn

    allocate(id_row(1:nrows)) ! Row indices of nonempty elements
    allocate(nze_row(1:nrows)) ! Number of nonempty elements in local rows
    if(.not. include_all) then ! If not all atoms in the system are neighbours
      nze_row(:) = 0
      id_row(:) = 0
      allocate(id_col_tmp(1:nmax))
      do iblks = 0, nblks  ! Checking all local rows of blocks
        iwf1 = iblks * BlockSize_c * Nodes + BlockSize_c * Node + 1
        iwf2 = iblks * BlockSize_c * Nodes + BlockSize_c * (Node + 1)
        if(iwf2 .gt. wf_dim) iwf2 = wf_dim
        do iwf = iwf1, iwf2 ! Checking all rows within the block of rows
          set_neib = .false.
          if(iwf == iwf1) then
            set_neib = .true.
          else
            if(iatom(iwf) .ne. iatom(iwf-1)) set_neib = .true.
          end if
          ! If the wavefunction that corresponds to this row belows to a new atom, it is
          ! needed to find its neighbours
          if(set_neib) then
            if(allocated(neib)) deallocate(neib) ! Deallocating the old list of neighbours
            call mneighb(ucell,rcoor,na_u,xa,iatom(iwf),0,nna) ! Looking for neighbours
            allocate(neib(1:nna)) ! Allocating the new list of neighbours
            neib(:) = 0
            index = 0
            ! Setting the list of neighbours excluding the repeating atoms
            do i = 1, nna
              found = .false.
              do j = 1, index
                if(neib(j) == jan(i)) found = .true.
              end do
              if(.not. found) then
                index = index + 1
                neib(index) = jan(i)
              end if
            end do
            nna = index ! Corrected number of neighbours
          end if

          io = io + 1
          id_row(io) = jo
          ! Matrix elements that correspond to atomic orbitals centered on neighbours are set nonempty
          do j = 1, nna
            ja = neib(j) ! One of the neighbours
            norb = lasto(ja) - lasto(ja-1) ! The number of atomic orbitals corresponding to this atom
            nze_c = nze_c + norb  ! The total number of nonempty matrix elements
            nze_row(io) = nze_row(io) + norb ! The number of nonempty matrix elements in the row considered
            ! Adjusting the size of the array with column indices of nonempty matrix elements if necessary
            if(nze_c .gt. nmax) then
              if(ncol .gt. 0) then
                allocate(id_col(1:ncol))
                do i = 1, ncol
                  id_col(i) = id_col_tmp(i)
                end do
              end if
              deallocate(id_col_tmp)
              nmax = nze_c + dn ! New size of the array of column indices
              allocate(id_col_tmp(1:nmax))
              if(ncol .gt. 0) then
                do i = 1, ncol
                  id_col_tmp(i) = id_col(i)
                end do
                deallocate(id_col)
              end if
            end if
            ncol = ncol + norb ! The number of defined column indices
            do iorb = 1, norb
              jo = jo + 1
              id_col_tmp(jo) = lasto(ja-1)+iorb ! Setting up column inidices
            end do
          end do
        end do
      end do 
      if(allocated(neib)) deallocate(neib) ! Deallocating the list of neighbours
    else ! If all atoms are neighbours, then all matrix elements are nonempty
      do io = 1, nrows
        id_row(io) = (io-1) * h_dim  ! Row inidices of nonempty matrix elements
        nze_row(io) = h_dim ! Number of nonempty matrix elements in each local row
      end do
      nze_c = h_dim * nrows ! Total number of nonempty matrix elements
    end if      

    allocate(c_loc(1:nze_c))  ! Local nonempty elements of the coefficient matrix
    allocate(id_col(1:nze_c)) ! Final array of column indices of local nonempty
                              ! elements of correct size
    jo = 0
    io = 0
    do iblks = 0, nblks ! For each block of rows
      iwf1 = iblks * BlockSize_c * Nodes + BlockSize_c * Node + 1
      iwf2 = iblks * BlockSize_c * Nodes + BlockSize_c * (Node + 1)
      if(iwf2 .gt. wf_dim) iwf2 = wf_dim
      do iwf = iwf1, iwf2 ! For each row within the block
        io = io + 1
        do i = 1, nze_row(io)
          jo = jo + 1
          if(set_rand) then ! Random initialization of the nonempty matrix elements
            do k = 1, 2
              call omm_bsd_lcg(seed, rn(k))
            end do
            cgval = sign(0.5_dp * rn(1), rn(2) - 0.5_dp) ! The value from -0.5 to 0.5
          else
            cgval = 0.0_dp ! If only sparsity is checked, the matrix elements are initialized by zeros
          end if
          ! The matrix element value with account of occupation and normalization coefficients
          c_loc(jo) = cgval * coef * fact(iwf)
          if(include_all) then
            id_col(jo) = i ! Column inidices if all atoms in the system are neighbours
          else
            id_col(jo) = id_col_tmp(jo) ! Column inidices if not all atoms are neighbours
          end if
        end do
      end do 
    end do
    if(allocated(iatom)) deallocate(iatom)
    if(allocated(fact)) deallocate(fact)
    if(.not. include_all) then
      if(allocated(id_col_tmp)) deallocate(id_col_tmp)
    end if

    ! Passing pointers to the arrays of the coefficient matrix in the csr format to MS
    ! Ordering the indices in such a way that the column indices are in the growing order
    ! for each row
    call m_register_csr(C_csr, wf_dim, h_dim, nrows, id_row, id_col, &
      nze_row, c_loc, order=.true., blk_size=BlockSize_c)

    ! Allocating the coefficient matrix in the MS format
    if (.not. C_min%is_initialized) call m_allocate(C_min, wf_dim, h_dim,&
      label=m_storage, blocksize1=BlockSize_c, blocksize2=BlockSize)

    if(ionode) print'(a)','C_min CSR matrix created'
    ! Converting the coefficient matrix from the csr format to MS
    call m_copy(C_min, C_csr)
    if(ionode) print'(a)','C_min CSR matrix converted to DBCSR'
    ! Removing the pointers to the arrays of the csr matrix
    call m_deallocate(C_csr)

    ! Deallocating arrays of the csr matrix
    deallocate(id_col)
    deallocate(id_row)
    deallocate(nze_row)
    deallocate(c_loc)
    call timer('init_c_matrix', 2)
  end subroutine init_c_matrix

end subroutine omm_min_block

!===============================================================================!
! Calculation of the number of wavefunctions corresponding to the atom
! considered. For Ordejon-Mauri method, this number equals half of the
! atomic charge (qa/2). If the atomic charge is odd, that for one atom it is
! (qa+1)/2, for the next one (qa-1)/2, then again (qa+1)/2 and (qa-1)/2, etc.
! For Kim method, the number of wavefunctions is (qa+2)/2. Correspondingly,
! the highest state is either unoccupied or half-occupied.
!===============================================================================!
subroutine get_number_of_lwfs_on_atom(ia, use_kim, nelectr, indexi)
  implicit none
  
  integer, intent(in)  :: ia
  logical, intent(in)  :: use_kim
  integer, intent(out) :: nelectr
  integer, intent(out) :: indexi
  
  real(dp) :: tiny
  logical, save :: secondodd = .false.
 
  tiny = 1.d-10
  nelectr = qa(ia) + tiny 
 
  if (abs(nelectr - qa(ia) + tiny) .gt. 1e-3) then
    if(ionode) then
      write(6,*) 'omm_min: Wrong atomic charge for atom ',ia
      write(6,*) '      qa = ',qa(ia),' must be an integer'
    endif
    call die()
  end if
  if(use_kim) then
    indexi = ( (nelectr + 2) / 2 )
  else
    if ( (nelectr / 2) * 2 .ne. nelectr) then
      if (secondodd) then
        indexi = ( ( nelectr - 1 ) / 2 )
        secondodd = .false.
      else
        indexi = ( ( nelectr + 1 ) / 2 )
        secondodd = .true.
      end if
    else
      indexi = ( ( nelectr ) / 2 )
    end if
  end if

end subroutine get_number_of_lwfs_on_atom

end module m_lib_omm
