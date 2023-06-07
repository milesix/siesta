! ---
! Copyright (C) 1996-2023      The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!==================================================================
!
program fat_atomic
!
!
! Computes the projection of each eigenvector (nk) on a given Ilm subspace,
! where I is an atom index and lm are the angular momentum quantum numbers.
! 
! Two ways of computing these projections are offered:
! a) The traditional way, modelled after the computation of the pDOS  
!    with mprop, using
!          g_mu = Real(SUM_{ c_mu  S(k)_mu_nu c_nu})
!
!     The 'real part' is inherited from pDOS calculations in which a sum  
!     over k-points is done).  This sometimes results in negative, although
!     typically very small, projections in some cases. To get the Ilm projection,
!     the contributions of all the mu's in the Ilm set are added.
!
! b) An orthogonal projection, using implicitly the Lowdin  
!    orthogonal basis associated the the basis set, and taking advantage
!    of the "optimal resemblance" feature of the Lowdin orthonormalization
!    to compute the projection over a given orbital mu as the square modulus of the
!    coefficient of wfs_nk in the Lowdin basis' mu. To get the Ilm projection,
!    the contributions of all the mu's in the Ilm set are added.
!
! The program computes both projections, using '&n' and '&o' markers that
! allow selection of the 'new' (Lowdin) and 'old' methods with grep.
!

  use main_vars
  use io_hs, only: read_hs_file

  use iso_c_binding, only: c_loc, c_f_pointer

  implicit none

  logical :: gamma_wfsx

  complex(dp), allocatable  :: S_sqroot(:,:)
  complex(dp), allocatable, target  :: S_inv_sqroot(:,:)
  complex(dp), allocatable  :: S_k(:,:)
  complex(dp), pointer  :: psi_coeffs(:) => null()
  complex(dp), pointer  :: coeff(:,:) => null()
  complex(dp), allocatable, target  :: wf_gamma(:)   ! For Gamma-point calculations

  real(dp) :: projs(9), old_projs(9)
 
  real(dp) :: sum_projs
  complex(dp) :: proj_new, proj_old, norm, phase
  integer :: ja


  integer  :: nwfmx, nwfmin
  integer  :: min_band = 1
  integer  :: max_band =  huge(1)
  logical  :: min_band_set = .false.
  logical  :: max_band_set = .false.
  logical  :: band_interval_set = .false.
  real(dp) :: min_eigval
  real(dp) :: max_eigval
  real(dp) :: min_eigval_in_file
  real(dp) :: min_eigval_in_band_set
  real(dp) :: max_eigval_in_band_set
  real(dp) :: minimum_spec_eigval = -huge(1.0_dp)
  real(dp) :: maximum_spec_eigval = huge(1.0_dp)

  integer  :: ib, nbands, nspin_blocks
  logical  :: non_coll

  integer :: proj_u
  logical :: write_lowdin_basis = .false.
  character(len=2) :: mark

  !
  !     Process options
  !
  n_opts = 0
  do
     call getopts('dhlR:b:B:s',opt_name,opt_arg,n_opts,iostat)
     if (iostat /= 0) exit
     select case(opt_name)
     case ('d')
        debug = .true.
     case ('l')
        energies_only = .true.
     case ('R')
        ref_line_given = .true.
        ref_line = opt_arg
     case ('b')
        read(opt_arg,*) min_band
        min_band_set = .true.
     case ('B')
        read(opt_arg,*) max_band
        max_band_set = .true.
     case ('p')
        write_lowdin_basis = .true.
     case ('h')
        call manual()
     case ('?',':')
        write(0,*) "Invalid option: ", opt_arg(1:1)
        write(0,*) "Usage: fat [ options ] SystemLabel "
        write(0,*) "Use -h option for manual"
        STOP
     end select
  enddo

  nargs = command_argument_count()
  nlabels = nargs - n_opts + 1
  if (nlabels /= 1)  then
     write(0,*) "Usage: fat [ options ] DescriptionFile "
     write(0,*) "Use -h option for manual"
     STOP
  endif

  call get_command_argument(n_opts,value=sflnm,status=iostat)
  if (iostat /= 0) then
     STOP "Cannot get SystemLabel string prefix"
  endif

  band_interval_set = (min_band_set .or. max_band_set)

  !==================================================
  ! Read WFSX file with wave-function information

  write(6,"(a)") "Reading wf file: " // trim(sflnm) // ".WFSX"

  open(wfs_u,file=trim(sflnm)//'.WFSX',status='old',form='unformatted')
  read(wfs_u) nkp, gamma_wfsx
  allocate (wk(nkp), pk(3,nkp))

  read(wfs_u) nsp
  non_coll = (nsp >= 4)
  read(wfs_u) nao
  read(wfs_u)        !! Symbols, etc
  if (debug) print *, "WFSX read: nkp, nsp, nnao: ", nkp, nsp, nao

  nwfmx = -huge(1)
  nwfmin = huge(1)
  min_eigval = huge(1.0_dp)
  max_eigval = -huge(1.0_dp)
  min_eigval_in_band_set = huge(1.0_dp)
  max_eigval_in_band_set = -huge(1.0_dp)

  if (non_coll) then
     !----
     write(0,*) 'This program does not work yet with non-collinear spin'
     error stop 1
     !----
     nspin_blocks = 1
  else
     nspin_blocks = nsp
  endif
  
  do ik=1,nkp
     do is=1,nspin_blocks

        read(wfs_u) idummy, pk(1:3,ik), wk(ik)
        if (idummy /= ik) stop "ik index mismatch in WFS file"
        read(wfs_u) is0
        read(wfs_u) number_of_wfns
        nwfmx = max(nwfmx,number_of_wfns)
        nwfmin = min(nwfmin,number_of_wfns)

        do iw=1,number_of_wfns
           read(wfs_u) iw0
           read(wfs_u) eigval
           min_eigval = min(min_eigval,eigval)
           max_eigval = max(max_eigval,eigval)
           ! 
           !
           if ((iw>=min_band).and.(iw<=max_band)) then
              min_eigval_in_band_set = min(min_eigval_in_band_set,eigval)
              max_eigval_in_band_set = max(max_eigval_in_band_set,eigval)
           endif
           read(wfs_u)
        enddo
     enddo
  enddo

  print "(a,2i5)", " Minimum/Maximum number of wfs per k-point: ", nwfmin, nwfmx
  print "(a,2f12.4)", "Min_eigval, max_eigval on WFS file: ",  &
       min_eigval, max_eigval
  print "(a,2f12.4)", "Min_eigval, max_eigval in band set : ",  &
          min_eigval_in_band_set, max_eigval_in_band_set

  if (.not. band_interval_set) then

     min_band = 1      ! Already set by default
     max_band = nwfmin

  else

     if (min_band_set .and. (min_band < 1)) then
        print "(a)", " ** Min_band implicitly reset to 1..."
        min_band = 1
     endif
     if (min_band_set .and. (min_band > nwfmin)) then
        print "(a,2i5)", " ** Min_band is too large for some k-points: (min_band, nwfmin):", min_band, nwfmin
        STOP
     endif
     if (max_band_set .and. (max_band > nwfmin)) then
        print "(a,2i5)", " ** Max_band is too large for some k-points: (max_band, nwfmin):", max_band, nwfmin
        print "(a)", " ** Max_band will be effectively reset to its maximum allowed value"
        max_band = nwfmin
     endif
     if (max_band_set .and. (max_band < min_band)) then
        print "(a,2i5)", " ** Max_band is less than min_band: (max_band, min_band):", max_band, min_band
        STOP
     endif

     min_eigval = min_eigval_in_band_set
     max_eigval = max_eigval_in_band_set
  endif
  print "(a,3i4)", "Band set used: (min, max):",  min_band, max_band

  ! Read HSX file
  ! We do not really neeed H here; only S, but there is other relevant
  ! information in the file.
  ! Will pick up atoms, zval, and thus the nominal number of electrons,
  ! but the total charge is read as qtot.

  call read_hs_file(trim(sflnm)//".HSX")

  if (energies_only) STOP

  !====================

  ! * Orbital list

  allocate(za(no_u), zc(no_u), zn(no_u), zl(no_u), zx(no_u), zz(no_u))
  nao = 0
  do ia=1,na_u
     it = isa(ia)
     io = 0
     do 
        io = io + 1
        if (io > no(it)) exit
        lorb = lquant(it,io)
        do ko = 1, 2*lorb + 1
           nao = nao + 1
           za(nao)=ia
           zc(nao)=it
           zn(nao)=nquant(it,io)
           zl(nao)=lorb
           zx(nao)=ko
           zz(nao)=zeta(it,io)
        enddo
        io = io + 2*lorb
     enddo
  enddo
  if (nao /= no_u) STOP "nao /= no_u"

  ! ==================================

  write(6,"('Writing files: ',a,'.stt ...')") trim(sflnm)
  open(stt_u,file=trim(sflnm)//'.info')
  write(stt_u,"(/'UNIT CELL ATOMS:')")
  write(stt_u,"(3x,i4,2x,i3,2x,a20)") (i, isa(i), label(isa(i)), i=1,na_u)
  write(stt_u,"(/'BASIS SET:')")
  write(stt_u,"(5x,a20,3(3x,a1))") 'spec', 'n', 'l', 'z'
  do it=1,nspecies
     write(stt_u,"(5x,a20)") trim(label(it))
     io = 0
     do 
        io = io + 1
        if (io > no(it)) exit
        write(stt_u,"(3(2x,i2))") nquant(it,io), lquant(it,io), zeta(it,io)
        io = io + 2*lquant(it,io)
     enddo
  enddo

  if ( nsp == 8 ) then
    write(stt_u,"(/'SPIN (spin-orbit): ',i2)") nspin_blocks
  else if ( nsp == 4 ) then
    write(stt_u,"(/'SPIN (non-coll): ',i2)") nspin_blocks
  else
    write(stt_u,"(/'SPIN: ',i2)") nspin_blocks
  endif
  
  write(stt_u,"(/'AO LIST:')")
  taux=repeat(' ',len(taux))
  do io=1,no_u
     taux(1:30)=repeat(' ',30)
     ik=1
     if (zl(io).eq.0) then
        taux(ik:ik)='s'
        ik=0
     elseif (zl(io).eq.1) then
        if (zx(io).eq.1) taux(ik:ik+1)='py'
        if (zx(io).eq.2) taux(ik:ik+1)='pz'
        if (zx(io).eq.3) taux(ik:ik+1)='px'
        ik=0
     elseif (zl(io).eq.2) then
        if (zx(io).eq.1) taux(ik:ik+2)='dxy'
        if (zx(io).eq.2) taux(ik:ik+2)='dyz'
        if (zx(io).eq.3) taux(ik:ik+2)='dz2'
        if (zx(io).eq.4) taux(ik:ik+2)='dxz'
        if (zx(io).eq.5) taux(ik:ik+5)='dx2-y2'
        ik=0
     elseif (zl(io).eq.3) then
        taux(ik:ik)='f'
     elseif (zl(io).eq.4) then
        taux(ik:ik)='g'
     elseif (zl(io).eq.5) then
        taux(ik:ik)='h'
     endif
     write(stt_u,"(3x,i5,2x,i3,2x,a20)",advance='no')  &
                          io, za(io), trim(label(zc(io)))
     if (ik.eq.0) then
        write(stt_u,"(3x,i2,a)") zn(io), trim(taux)
     else
        write(stt_u,"(3x,i2,a,i2.2)") zn(io), trim(taux), zx(io)
     endif
  enddo

     write(stt_u,"(/'KPOINTS:',i7)") nkp
     do ik=1,nkp
        write(stt_u,"(3x,3f9.6)") pk(:,ik)
     enddo

  close(stt_u)

  !==================================


  allocate(S_k(no_u,no_u))
  allocate(S_sqroot(no_u,no_u), S_inv_sqroot(no_u,no_u))

  ! * Fatband weights

  ! nspin has been read in iohs
  nbands = max_band - min_band + 1

     ! The first dimension is the number of real numbers per orbital
     ! 1 for real wfs, 2 for complex, and four for the two spinor components

     if (non_coll) then
        allocate(wf_single(4,1:no_u))
        allocate(wf(4,1:no_u))
     else
        if (gamma_wfsx) then
           allocate(wf_single(1,1:no_u))
           ! Use complex array to fit the complex-case machinery
           !!  allocate(wf(1,1:no_u))     
           allocate(wf_gamma(1:no_u))
        else
           allocate(wf_single(2,1:no_u))
           allocate(wf(2,1:no_u))
        endif
     endif

     ! Second pass for computation of projections and Writing in a
     ! format appropriate for Pyprocar and others

        open(newunit=proj_u,file=trim(sflnm)// '.projs')
        write(proj_u,"(a,2i5)") "# " // trim(sflnm) // " min_band, max_band: ", min_band, max_band
        write(proj_u,"(6i6)")   nkp, nbands, nspin_blocks, na_u, 9

        rewind(wfs_u)

        read(wfs_u) 
        read(wfs_u) 
        read(wfs_u) 
        read(wfs_u) 
     
        do ik=1,nkp

           write(proj_u,"(//,a,i4,3(1x,f10.5),2x,a,1x,f12.8)") 'K-point: ', ik, pk(1:3,ik), ' weight: ', wk(ik)

           !-------------------------------------
           ! Process S_k and friends here

           compute_Sk: block
             complex(dp), parameter :: ii = (0.0_dp, 1.0_dp)
             complex(dp) :: phase
             integer jo, ind

             do io = 1,no_u
                S_k(:,io) = 0.0_dp
                do j = 1,numh(io)
                   ind = listhptr(io) + j
                   phase = dot_product(pk(:,ik),xij(:,ind))
                   jo = listh(ind)
                   jo = MODP(jo,no_u)  ! To allow auxiliary supercells
                   ! Maybe transpose?
                   S_k(jo,io) = S_k(jo,io) + Sover(ind)*exp(-ii*phase)  
                enddo
             enddo
             if (debug) then
                write(6,*) "Real(S(k)):"
                do i = 1, min(8,no_u)
                   do j = 1, min(8,no_u)
                      write(6,fmt="(1x,f10.4)",advance="no") real(S_k(i,j),kind=dp)
                   enddo
                   write(6,*)
                enddo
             endif
           end block compute_Sk

           call get_sk_funcs(S_k, S_sqroot, S_inv_sqroot)

              ! Use S_inv_sqroot to get the orthogonal basis
              ! The coefficients of the ith orthogonal vector will be
              ! those in the ith COLUMN of S_inv_sqroot

           if (write_lowdin_basis) then
              coeff => S_inv_sqroot
              do i = 1, no_u
                 print "(/,a,i0)", "coeffs of orthog orbital number: ", i
                 do ja = 1, no_u
                       write(*,"(1x,a1,f6.3,f7.3,a1)",advance="no") '(', coeff(ja,i), ')'
                 enddo
                 write(*,*)
                 !! norm = inner_prod(coeff(:,i),coeff(:,i),S_k)
                 !! print *, "NORM:", abs(norm)
              enddo
           endif
              !-------------------------------------

        
           do is=1,nspin_blocks
              if (nspin_blocks > 1) then
                 write(proj_u,"(/a,i1)") 'Spin: ', is
              endif
              ib = 0 ! for counting wfs in set...
              read(wfs_u) 
              read(wfs_u) 
              read(wfs_u)  number_of_wfns
              do iw=1,number_of_wfns
                 read(wfs_u) 
                 read(wfs_u) eigval

                 ! Use only the specified band set
                 if ( (iw<min_band) .or. (iw>max_band)) then
                    read(wfs_u)   ! Still need to read this
                    CYCLE
                 endif

                 ib = ib + 1
                 write(proj_u,"(/,a,i4,f12.6,a,i3,a)") 'Wf: ', iw, eigval, " (in set: ", ib, ")"

                 read(wfs_u) (wf_single(:,io), io=1,no_u)
                 
                 ! Use a double precision form in what follows
                    if (gamma_wfsx) then
                       wf_gamma(:) = cmplx(wf_single(1,:),0.0_dp,kind=dp)
                       psi_coeffs => wf_gamma
                    else
                       wf(:,:) = real(wf_single(:,:), kind=dp)
                       call c_f_pointer(c_loc(wf), psi_coeffs, [no_u]) ! Note: only for collinear spin...
                    endif

                    ! For spinors:
                    !call c_f_pointer(c_loc(wf), psi_coeffs, [2*no_u])
                    ! or, for more clarity
                    !call c_f_pointer(c_loc(wf), psi_spinor, [2,no_u])
                    ! Then dispatch to the appropriate contraction below
                    
                    write(proj_u,"(/,3x,9(1x,a7))") "s", "py", "pz", "px", &
                      "dxy", "dyz", "dz2", "dxz", "dx2-z2"

                    ! Aggregate the different orbital projections into lm bins
                    nao = 0
                    do ia = 1, na_u
                       projs(:) = 0.0_dp
                       old_projs(:) = 0.0_dp
                       it = isa(ia)
                       io = 0
                       do
                          io = io + 1
                          if (io > no(it)) exit
                          lorb = lquant(it,io)
                          do ko = 1, 2*lorb + 1
                             nao = nao + 1

                             proj_new = abs(sqroot_contraction(nao))**2
                             proj_old = real(old_contraction(nao),kind=dp)
                             select case (lorb)
                             case (0)
                                projs(1) = projs(1) + proj_new
                                old_projs(1) = old_projs(1) + proj_old
                             case (1)
                                projs(1+ko) = projs(1+ko) + proj_new
                                old_projs(1+ko) = old_projs(1+ko) + proj_old
                             case (2) 
                                projs(1+3+ko) = projs(1+3+ko) + proj_new
                                old_projs(1+3+ko) = old_projs(1+3+ko) + proj_old
                             end select
                          enddo
                          io = io + 2*lorb
                       enddo
                       write(proj_u,"(i3,9f8.4,2x,f8.4,1x,a2)") ia, (projs(i),i=1,9), sum(projs), "&n"
                       write(proj_u,"(i3,9f8.4,2x,f8.4,1x,a2)") ia, (old_projs(i),i=1,9), sum(old_projs), "&o"
                    enddo



              enddo   ! iwf
           enddo      ! is

        enddo         ! ik
        
 
 CONTAINS

      subroutine manual()

      write(6,"('* FAT_ATOMIC (lm projections)')")
      write(6,"('  Alberto Garcia, ICMAB-CSIC, 2023 ')")
      write(6,*)
      write(6,"('    FAT_ATOMIC calculates eigenvector projections ')")
      write(6,"('    using output files obtained with SIESTA. ')")
      write(6,"('    The Lowdin basis is used implicitly for projections. ')")
      write(6,"('    The output is in Pyprocar style.')")
      write(6,"('  ')")
      write(6,*) "Usage: fat_atomic [ options ] SystemLabel"
      write(6,*) "Options:"
      write(6,*) "           -h:  print manual                    "
      write(6,*) "           -d:  debug                    "
      write(6,*) "           -l:  print summary of energy information         "
      write(6,*) "           -p:  print Lowdin basis information              "
      write(6,*) "    "
      write(6,*) "   Selection of eigenstates to be used: "
      write(6,*) "    "
      write(6,*) "   -b Min_band  :  set minimum band index to be used               "
      write(6,*) "   -B Max_band  :  set maximum band index to be used               "
      write(6,*) "    "
      write(6,*)
      stop

      end subroutine manual

      ! Some of these functions use host association for
      ! relevant variables
      function old_contraction(i) result (res)
        ! This is what Siesta did implicitly for DOS and fat
        integer, intent(in) :: i
        complex(dp) :: res
        integer  j

        res = 0.0_dp
        do j = 1, no_u
           res = res + S_k(i,j) *  psi_coeffs(j)
        enddo
        res =  conjg(psi_coeffs(i)) * res

      end function old_contraction

      function sqroot_contraction(i) result (res)
        integer, intent(in) :: i
        complex(dp) :: res
        integer  j

        res = 0.0_dp
        do j = 1, no_u
           res = res + psi_coeffs(j) * S_sqroot(i,j)
        end do

      end function sqroot_contraction

      function inner_prod(a,b,S) result (res)
        complex(dp), intent(in) :: a(:)
        complex(dp), intent(in) :: b(:)
        complex(dp), intent(in) :: S(:,:)
        complex(dp)             :: res

        complex(dp), allocatable :: wrk(:)
        
        integer i, j, n

        n = size(a)
        allocate(wrk(n))
        wrk = matmul(S,b)
        res = dot_product(a,wrk)

        deallocate(wrk)
        
      end function inner_prod

      ! Compute functions of S_k

      subroutine get_sk_funcs(S_k ,S_sqroot, S_inv_sqroot)
        complex(dp), intent(in)  :: S_k(:,:)
        complex(dp), intent(out) :: S_sqroot(:,:), S_inv_sqroot(:,:)


        integer :: no_u

        complex(dp), allocatable :: S_base(:,:), B(:,:)

        real(dp), allocatable :: w(:)              ! Eigenvalues
        integer :: lwork                  ! Size of the workspace
        integer :: lrwork                 ! Size of the workspace
        integer :: liwork                 ! Size of the workspace

        complex(dp), allocatable :: work(:)
        real(dp), allocatable :: rwork(:)
        integer, allocatable :: iwork(:)

        integer :: n, iwork_query(1)
        complex(dp) :: work_query(1)         ! Workspace query
        real(dp) :: rwork_query(1)         ! Workspace query

        integer :: i, j, jo, info, io, ind
        complex(dp), parameter :: ii = (0.0_dp, 1.0_dp)

        no_u = size(S_k,dim=1)
        
        allocate(S_base(no_u,no_u))
        allocate(w(no_u))
       
        S_base = S_k

        ! Workspace query
        lwork = -1
        lrwork = -1
        liwork = -1

        call zheevd('V', 'U', no_u, S_base, no_u, w, work_query, lwork, rwork_query, lrwork, iwork_query, liwork, info)

        lwork = int(real(work_query(1)))
        lrwork = int(rwork_query(1))
        liwork = iwork_query(1)

        allocate(work(lwork))
        allocate(rwork(lrwork))
        allocate(iwork(liwork))

        call zheevd('V', 'U', no_u, S_base, no_u, w, work, lwork, rwork, lrwork, iwork, liwork, info)

        ! Check for successful execution
        if (info == 0) then
           if (debug) then
              print *, "Range of Eigenvalues of S_k:", w(1), w(no_u)
           endif
        else
           print *, "Error: zheevd failed with info =", info
        endif

        ! We have that S = Z D Z^H
        ! Compute B=ZD
        allocate(B(no_u,no_u))

        ! Compute square root ----------------------
        do j = 1, no_u
           b(1:no_u,j) = S_base(1:no_u,j) * sqrt(w(j))
        enddo
        ! Now, compute B*Z^T
        call zgemm('N','C', no_u,no_u,no_u, (1.0_dp,0.0_dp), B,no_u, S_base,no_u, (0.0_dp,0.0_dp),S_sqroot,no_u)

        B = matmul(S_sqroot,S_sqroot)

        if (debug) then
           write(6,*) "Check of square root: Product (only valid when func=sqrt():"
           if ( any(  abs(B-S_k) > 1.0e-5) ) then
              print *, "S_sqroot is not a proper square root of S_k"
              error stop 1
           endif
        endif
        
        ! Compute inverse square root ----------------------
        do j = 1, no_u
           b(1:no_u,j) = S_base(1:no_u,j) / sqrt(w(j))
        enddo
        ! Now, compute B*Z^T
        call zgemm('N','C', no_u,no_u,no_u, (1.0_dp,0.0_dp), B,no_u, S_base,no_u, (0.0_dp,0.0_dp),S_inv_sqroot,no_u)

        B = matmul(S_inv_sqroot,S_inv_sqroot)
        S_base = matmul(B,S_k)  ! Note reuse

        if (debug) then
           write(6,*) "Check of S_inv_sqroot. It should be identity matrix"
           do i = 1, min(8,no_u)
              do j = 1, min(8,no_u)
                 write(6,fmt="(1x,f10.4)",advance="no")  real(s_base(i,j))
              enddo
              write(6,*)
           enddo
        endif

        deallocate(work,iwork,w,S_base,B)

      end subroutine get_sk_funcs


! A MOD function which behaves differently on the edge.                                                               
! It will NEVER return 0, but instead return the divisor.                                                             
! This makes it useful in FORTRAN do-loops which are not zero-based                                                   
! but 1-based.                                                                                                        
! Thus we have the following scheme:                                                                                  
!  MODP(x',x) == x         for x' == x                                                                                
!  MODP(x',x) == MOD(x',x) for x' /= x                                                                                
  elemental function MODP(a,p)
    integer, intent(in) :: a,p
    integer :: MODP
    MODP = MOD(a-1,p) + 1
  end function MODP


end program fat_atomic


