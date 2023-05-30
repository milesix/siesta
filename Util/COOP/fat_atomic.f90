! ---
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!==================================================================
!
program fatband_atomic
!
!
! Computes the generalized projection of each eigenvector on a given
! orbital set.
! The syntax for orbital sets is the same as in mprop, so this
! program can re-use as description file the same DOS-type file
! as mprop.
!
!

  use main_vars
  use orbital_set, only: get_orbital_set
  use io_hs, only: read_hs_file
  use read_curves, only: read_curve_information, mask_to_arrays

  use iso_c_binding, only: c_loc, c_f_pointer

  implicit none

  logical :: gamma_wfsx, got_qcos
  integer :: ii1, ii2, ind, ind_red, no1, no2, n_int, nnz
  real(dp) :: factor

  real(dp), parameter  :: tol_overlap = 1.0e-10_dp

  logical, allocatable   :: mask2(:)
  logical, allocatable   :: orb_mask_atomic(:,:)
  integer, allocatable   :: num_red(:), ptr(:), list_io2(:), list_ind(:)
  real(dp), allocatable  :: eig(:,:), fat(:,:,:,:)

  complex(dp), allocatable  :: S_sqroot(:,:)
  complex(dp), allocatable  :: S_inv_sqroot(:,:)
  complex(dp), allocatable  :: S_k(:,:)
  complex(dp), pointer  :: psi_coeffs(:) => null()
  complex(dp), allocatable  :: coeff(:,:)

  real(dp) :: sum_projs
  complex(dp) :: proj_new, norm, phase
  integer :: n_orbs_atom, ja


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

  integer :: proj_u, nc_p
  logical :: atom_lm_scheme = .false.
  logical :: all_interactions = .true.
  character(len=2) :: mark

  !
  !     Process options
  !
  n_opts = 0
  do
     call getopts('dhlaR:b:B:s',opt_name,opt_arg,n_opts,iostat)
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
     case ('a','n')
        atom_lm_scheme = .true.
     case ('s')   ! "same sets I and II"
        all_interactions = .false.
     case ('h')
        call manual()
     case ('?',':')
        write(0,*) "Invalid option: ", opt_arg(1:1)
        write(0,*) "Usage: fat [ options ] DescriptionFile "
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

  call get_command_argument(n_opts,value=mflnm,status=iostat)
  if (iostat /= 0) then
     STOP "Cannot get .mpr file root"
  endif

  band_interval_set = (min_band_set .or. max_band_set)

  !==================================================

  ierr=0

  ! Read type of job

  open(mpr_u,file=trim(mflnm) // ".mpr", status='old')
  read(mpr_u,*) sflnm
  read(mpr_u,*) what
  if (trim(what) /= "DOS") STOP "Fatbands needs DOS-style jobfile"
  
  !==================================================
  ! Read WFSX file

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

  ! Process orbital sets

  ! Here we should implement new code to generate "orbital sets" for
  ! each lm pair and for each atom involved 

  if (atom_lm_scheme) then

     ! There are 9 "curves" per atom
     ncb = 9*na_u
     allocate(orb_mask_atomic(no_u,ncb))
     allocate (koc(ncb,2,no_u))
     orb_mask_atomic(:,:) = .false.
     nc_p = 0
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
              print *, "ia, io, nao, lorb, ko: ", ia, io, nao, lorb, ko
              select case (lorb)
              case (0)
                 orb_mask_atomic(nao,nc_p+1) = .true.
                 write(tit(nc_p+1),"(a,i3)") "s on atom", ia
                 print *, "s on atom", ia, " put in set", nc_p+1
              case (1)
                 orb_mask_atomic(nao,nc_p+1+ko) = .true.
                 write(tit(nc_p+1+ko),"(a,i1,a,i3)") "p ", ko," on atom", ia
                 print *, "p ", ko, "on atom", ia, " put in set", nc_p+1+ko
              case (2) 
                 orb_mask_atomic(nao,nc_p+4+ko) = .true.
                 write(tit(nc_p+4+ko),"(a,i1,a,i3)") "d ", ko," on atom", ia
                 print *, "d ", ko, "on atom", ia, " put in set", nc_p+1+ko
              end select
           enddo
           io = io + 2*lorb
        enddo
        nc_p = nc_p + 9
     enddo

     call mask_to_arrays(ncb,orb_mask_atomic(:,:),noc(:,1),koc(:,1,:))
     if (all_interactions) then
        ! Consider all interactions
        orb_mask_atomic(:,:) = .true.
     else
        ! Leave the mask for second set equal to the first
     endif
     call mask_to_arrays(ncb,orb_mask_atomic(:,:),noc(:,2),koc(:,2,:))
        
  else
!
     allocate(orb_mask(no_u,2,ncbmx))
     allocate (koc(ncbmx,2,no_u))
     ! Give values to flags for reuse of the reading routine
     dos = .true.
     coop = .false.
     call read_curve_information(.true.,.false.,  &
                                    mpr_u,no_u,ncbmx,ncb,tit,orb_mask,dtc)

     if (all_interactions) then
        orb_mask(:,2,1:ncb) = .true.       ! All orbitals considered
     else
        orb_mask(:,2,1:ncb) = orb_mask(:,1,1:ncb)  ! Just the same set as I
     endif

     call mask_to_arrays(ncb,orb_mask(:,1,:),noc(:,1),koc(:,1,:))
     call mask_to_arrays(ncb,orb_mask(:,2,:),noc(:,2),koc(:,2,:))
  endif
!!
  write(6,"('Writing files: ',a,'.stt ...')") trim(mflnm)
  open(stt_u,file=trim(mflnm)//'.info')
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

     write(stt_u,"(/'FATBAND ORBITAL SETS:')")
     do ic=1,ncb
        write(stt_u,"(3x,a)") trim(tit(ic))
        write(stt_u,"(3x,'AO set I: ',/,15x,12i5)") (koc(ic,1,j),j=1,noc(ic,1))
        write(stt_u,"(3x,'Number of set II orbs:',i8)") noc(ic,2)
     enddo

  close(stt_u)

  !==================================

  if (ref_line_given) then
     allocate(ref_mask(no_u))
     print *, "Orbital set spec: ", trim(ref_line)
     call get_orbital_set(ref_line,ref_mask)
     do io=1, no_u
        if (ref_mask(io)) write(6,fmt="(i5)",advance="no") io
     enddo
     deallocate(ref_mask)
     write(6,*)
     STOP "bye from ref_line processing"
  endif


  !=====================


  allocate(S_k(no_u,no_u))
  allocate(S_sqroot(no_u,no_u), S_inv_sqroot(no_u,no_u))

  ! Assume binary compound with atoms with the same number of orbitals...
  n_orbs_atom = no_u/2

  allocate(coeff(1:no_u,1:no_u))


  ! * Fatband weights

  ! nspin has been read in iohs
  nbands = max_band - min_band + 1
  allocate(eig(nbands,nspin_blocks), fat(nbands,nspin_blocks,nkp,ncb))

     ! The first dimension is the number of real numbers per orbital
     ! 1 for real wfs, 2 for complex, and four for the two spinor components

     if (non_coll) then
        allocate(wf_single(4,1:no_u))
        allocate(wf(4,1:no_u))
     else
        if (gamma_wfsx) then
           allocate(wf_single(1,1:no_u))
           allocate(wf(1,1:no_u))
        else
           allocate(wf_single(2,1:no_u))
           allocate(wf(2,1:no_u))
        endif
     endif
     allocate (mask2(1:no_u))

     do ic=1,ncb

        no1 = noc(ic,1)
        no2 = noc(ic,2)

        mask2(1:no_u) = .false.
        do i2=1,no2
           io2=koc(ic,2,i2)              ! AO Set II
           mask2(io2) = .true.
        enddo

        ! Create reduced pattern
        ! First pass for checking dimensions

        allocate (num_red(no1))
        do i1=1,no1
           num_red(i1) = 0
           io1=koc(ic,1,i1)              ! AO Set I
           do ii1 = 1,numh(io1)
              ind = listhptr(io1)+ii1
              ii2 = indxuo(listh(ind))      ! Equiv orb in unit cell

              if ( .not. mask2(ii2)) Cycle    ! Is not one of Set II

              ! (what to do with semicore states??)

              num_red(i1) = num_red(i1) + 1
           enddo
        enddo
        allocate (ptr(no1))
        ptr(1)=0
        do i1=2,no1
           ptr(i1)=ptr(i1-1)+num_red(i1-1)
        enddo
        nnz = sum(num_red(1:no1))  

        write(*,"(a,3x,a,2x,a,i6,1x,i12)") 'Fatband coeffs set: ', trim(tit(ic)),  &
                                      'Base orbitals and interactions: ', &
                                       no1, nnz

        allocate (list_io2(nnz))
        allocate (list_ind(nnz))

        n_int = 0
        do i1=1,no1
           io1=koc(ic,1,i1)              ! AO Set I
           do ii1 = 1,numh(io1)
              ind = listhptr(io1)+ii1
              ii2 = indxuo(listh(ind))

              if ( .not. mask2(ii2)) Cycle

              ! (what to do with semicore states??)

              n_int = n_int + 1
              list_io2(n_int) = ii2
              list_ind(n_int) = ind
           enddo
        enddo
        if (n_int .ne. nnz) then
           print *, "n_int, nnz:", n_int, nnz
           STOP "mismatch"
        endif


     !Stream over file, without using too much memory

        rewind(wfs_u)

        read(wfs_u) 
        read(wfs_u) 
        read(wfs_u) 
        read(wfs_u) 

        if (.not. atom_lm_scheme) then
           open(fat_u,file=trim(mflnm)// "." // trim(tit(ic)) // '.EIGFAT')
           write(fat_u,"(a,2i5)") "# " // trim(sflnm) // " min_band, max_band: ", min_band, max_band
           write(fat_u,"(3i6)")   nbands, nspin_blocks, nkp
        endif

        do ik=1,nkp

           if (.not. atom_lm_scheme) then
              write(fat_u,"(i4,3(1x,f10.5))")  ik, pk(1:3,ik)
           endif
           
           do is=1,nspin_blocks

              ib = 0
              fat(:,is,ik,ic) = 0.0_dp

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
                 eig(ib,is) = eigval    ! This will be done for every curve... harmless

                 read(wfs_u) (wf_single(:,io), io=1,no_u)
                 ! Use a double precision form in what follows
                 wf(:,:) = real(wf_single(:,:), kind=dp)


                 do i1 = 1, no1
                    io1=koc(ic,1,i1)              ! AO Set I

                      do i2 = 1,num_red(i1)
                         ind_red = ptr(i1)+i2
                         io2 = list_io2(ind_red)
                         ind = list_ind(ind_red)
                             
                                ! (qcos, qsin) = C_1*conjg(C_2)
                                !AG: Corrected:  (qcos, qsin) = conjg(C_1)*(C_2)
                                ! We might want to avoid recomputing this
                                if (non_coll) then
                                   ! Take as weight the "complete spinor" product
                                   qcos= wf(1,io1)*wf(1,io2) + &
                                        wf(2,io1)*wf(2,io2) + &
                                        wf(3,io1)*wf(3,io2) + &
                                        wf(4,io1)*wf(4,io2)
                                   qsin= wf(1,io1)*wf(2,io2) - &
                                        wf(2,io1)*wf(1,io2) + &
                                        wf(3,io1)*wf(4,io2) - &
                                        wf(4,io1)*wf(3,io2) 
                                else
                                   if (gamma_wfsx) then
                                      qcos = wf(1,io1)*wf(1,io2) 
                                      qsin = 0.0_dp
                                   else
                                      qcos= (wf(1,io1)*wf(1,io2) + &
                                           wf(2,io1)*wf(2,io2))
                                      qsin= (wf(1,io1)*wf(2,io2) - &
                                           wf(2,io1)*wf(1,io2))
                                   endif
                                endif
                             ! k*R_12    (r_2-r_1)
                             alfa=dot_product(pk(1:3,ik),xij(1:3,ind))

                             ! Crb = Real(C_1*conjg(C_2)*exp(-i*alfa)) * S_12
                             !AG: This one better --  or Real(conjg(C_1)*C_2)*exp(+i*alfa)) * S_12
                             ! Common factor computed here
                             factor =  (qcos*cos(alfa)-qsin*sin(alfa))

                             fat(ib,is,ik,ic) = fat(ib,is,ik,ic) + Sover(ind)*factor

                        enddo   ! i2
                    enddo  ! i1

                 enddo   ! iwf

                 if (.not. atom_lm_scheme) then

                    write(fat_u,"(4(4x,f10.4,f9.5))")   &
                         (eig(ib,is),fat(ib,is,ik,ic),ib=1,nbands)
                 endif
                 
              enddo      ! is

           enddo         ! ik
           
           deallocate (num_red)
           deallocate (ptr)
           deallocate (list_io2)
           deallocate (list_ind)

           close(fat_u)

        enddo    ! ic

       if (atom_lm_scheme) then
           ! Write in a format appropriate for Pyprocar and others
           
           ! Use the structure of the file to get the wavefunctions per k-point, etc

        open(newunit=proj_u,file=trim(mflnm)// '.projs')
        write(proj_u,"(a,2i5)") "# " // trim(sflnm) // " min_band, max_band: ", min_band, max_band
        write(proj_u,"(3i6)")   nbands, nspin_blocks, nkp

        rewind(wfs_u)

        read(wfs_u) 
        read(wfs_u) 
        read(wfs_u) 
        read(wfs_u) 
     
        do ik=1,nkp

           write(proj_u,"(//,a,i4,3(1x,f10.5))") 'K-point: ', ik, pk(1:3,ik)

           !-------------------------------------
           ! Process S_k and friends here

           compute_Sk: block
             complex(dp), parameter :: ii = (0.0_dp, 1.0_dp)
             complex(dp) :: phase
             integer jo

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
             write(6,*) "Real(S(k)):"
             do i = 1, min(8,no_u)
                do j = 1, min(8,no_u)
                   write(6,fmt="(1x,f10.4)",advance="no") real(S_k(i,j),kind=dp)
                enddo
                write(6,*)
             enddo

           end block compute_Sk

           call get_sk_funcs(S_k, S_sqroot, S_inv_sqroot)

              ! Use S_inv_sqroot to get the orthogonal basis
              ! The coefficients of the ith orthogonal vector will be
              ! those in the ith COLUMN of S_inv_sqroot

              do i = 1, no_u
                 coeff(1:no_u,i) = S_inv_sqroot(1:no_u,i)
                 print "(/,a,i0)", "coeffs of orthog orbital number: ", i
                 do ia = 1, 2
                    do j = 1, n_orbs_atom
                       ja = (ia-1)*n_orbs_atom + j
                       write(*,"(1x,a1,2f8.4,a1)",advance="no") '(', coeff(ja,i), ')'
                    enddo
                    write(*,*)
                 enddo
                 norm = inner_prod(coeff(:,i),coeff(:,i),S_k)
                 print *, "NORM:", abs(norm)
              enddo
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
                    wf(:,:) = real(wf_single(:,:), kind=dp)
                    call c_f_pointer(c_loc(wf), psi_coeffs, [no_u])

                    write(proj_u,"(/,3x)",advance="no")
                    do j=1, n_orbs_atom
                       write(proj_u,"(1x,i7)",advance="no") j
                    enddo
                    write(proj_u,*)
                    
                    ! Use inner product of wf with projecting orbital,
                    ! which is the matching one in the Lowdin orthogonal basis
                    do ia = 1, 2
                       write(proj_u,"(i3)",advance="no") ia
                       sum_projs = 0.0_dp
                       do i = 1, n_orbs_atom
                          j = (ia-1)* n_orbs_atom + i
                          proj_new = inner_prod(coeff(:,j),psi_coeffs(:),S_k)
                          write(proj_u,"(f8.4)",advance="no") abs(proj_new)**2
                          sum_projs = sum_projs + abs(proj_new)**2
                       enddo
                       write(proj_u,"(2x,f8.4)") sum_projs
                    enddo

                    ! Alternative method with contraction with sqrt(S)
                    do ia = 1, 2
                       write(proj_u,"(i3)",advance="no") ia
                       sum_projs = 0.0_dp
                       do i = 1, n_orbs_atom
                          j = (ia-1)* n_orbs_atom + i
                          proj_new = sqroot_contraction(j)
                          write(proj_u,"(f8.4)",advance="no") abs(proj_new)**2
                          sum_projs = sum_projs + abs(proj_new)**2
                       enddo
                       write(proj_u,"(2x,f8.4)") sum_projs
                    enddo

                    write(proj_u,"(/,3x,9(1x,a7))") "s", "py", "pz", "px", &
                      "dxy", "dyz", "dz2", "dxz", "dx2-z2"
                    nc_p = 0
                    do ia = 1, na_u

                       ! There are sometimes very small negative numbers
                       ! We mark them with '**'. Clients can choose what to do with them
                       if (any(fat(ib,is,ik,nc_p+1:nc_p+9) < -0.1)) then
                          mark="!!"
                       else  if (any(fat(ib,is,ik,nc_p+1:nc_p+9) < -0.00005)) then
                          mark="**"
                       else
                          mark ="  "
                       endif

                       write(proj_u,"(i3,9f8.4,2x,f8.4,1x,a2)") ia, &
                            (fat(ib,is,ik,ic),ic=nc_p+1,nc_p+9), sum(fat(ib,is,ik,nc_p+1:nc_p+9)), mark

                       nc_p = nc_p + 9
                    end do

              enddo   ! iwf
           enddo      ! is

        enddo         ! ik
        
     endif

 
 CONTAINS

      subroutine manual()

      write(6,"('* FAT(BANDS) PROGRAM')")
      write(6,"('  Alberto Garcia, ICMAB-CSIC, 2012 ')")
      write(6,*)
      write(6,"('    FAT calculates eigenvector projections ')")
      write(6,"('    using output files obtained with SIESTA. The atomic orbital (AO)')")
      write(6,"('    sets are defined in an input file (MLabel.mpr).')")
      write(6,"('  ')")
      write(6,*) "Usage: fat [ options ] MPROP_FILE_BASENAME"
      write(6,*) "Options:"
      write(6,*) "           -h:  print manual                    "
      write(6,*) "           -d:  debug                    "
      write(6,*) "           -l:  print summary of energy information         "
      write(6,*) "           -a:  use atom-lm scheme and generate projections file "
      write(6,*) "           -s:  ignore non-local interactions (for testing only) "
      write(6,*) "    "
      write(6,*) "   Selection of eigenstates to be used: "
      write(6,*) "    "
      write(6,*) "   -b Min_band  :  set minimum band index to be used               "
      write(6,*) "   -B Max_band  :  set maximum band index to be used               "
      write(6,*) "    "
      write(6,*)
      write(6,"('* .mpr FILE STRUCTURE')")
      write(6,"('         SLabel                   # Name of the siesta output files')")
      write(6,"('         DOS                      # DOS option of mprop is mandatory')")
      write(6,"('    /-  As many blocks as projections wanted ]')")
      write(6,"('    |    projection_name         # DOS projection name')")
      write(6,"('    \-   Subset of AO (*)        # Subset of orbitals included')")
      write(6,"('     (*) See below how to define subsets of AO')")
      write(6,"('     A final line with leading chars  ----  can signal the end of the input')")
      write(6,*)
      write(6,"('* INPUT FILES')")
      write(6,"('    [output files from SIESTA >=  2.4.1]')")
      write(6,"('    SLabel.WFSX and SLabel.HSX (new format)')")
      write(6,*)
      write(6,"('* OUTPUT FORMAT')")
      write(6,*) 
      write(6,*) " MLabel.CurveName.EIGFAT    :  File with eigenvalue and projection info "
      write(6,"('    [An information file with extension .info will always be generated]')")
      write(6,*)
      write(6,"('* PROJECTION AND CURVES NAMES')")
      write(6,"('    Alphanumerical string up to 30 char. with no spaces')")
      write(6,"('* SUBSET OF AO USING ORDER NUMBERS')")
      write(6,"('    List of integer numbers preceeded by a + symbol')")
      write(6,"('    Each number refers to one AO in the final list of AO of SIESTA')")
      write(6,"('    Example: + 23 65 78')")
      write(6,"('* SUBSET OF AO USING ATOM_SHELL NOTATION')")
      write(6,"('    List of atoms and shell groups of AO')")
      write(6,"('    General notation: ATOM_SHELL')")
      write(6,"('     > ATOM:  Atomic symbol refers to all the atoms of that type')")
      write(6,"('              Integer number refers to the N-th atom in unit cell')")
      write(6,"('     > SHELL: Integer1+Letter+Integer2')")
      write(6,"('               > Integer1 refers to the n quantum number')")
      write(6,"('               > Letter   refers to the l quantum number (s,p,d,f,g,h)')")
      write(6,"('               > Integer2 refers to a single AO into the n-l shell')")
      write(6,"('                   Alternatively, alphanumerical strings can be used')")
      write(6,"('                     p-shells   1  y    d-shells   1  xy   4  xz')")
      write(6,"('                                2  z               2  yz   5  x2-y2')")
      write(6,"('                                3  x               3  z2')")
      write(6,"('    Particular cases:')")
      write(6,"('     > Just ATOM is indicated: all the AO of the atom will be included')")
      write(6,"('     > No value for Integer2:  all the AO of the shell will be included')")
      write(6,"('    Example: Ca_3p Al 4_4d3 5 O_2py')")
      stop

      end subroutine manual

      ! Some of these functions use host association for
      ! relevant variables
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
!!$        do i = 1, n
!!$           wrk(i) = dot_product(a,S(:,i))
!!$        enddo
!!$        res  = dot_product(conjg(wrk),b)
        deallocate(wrk)
        
      end function inner_prod

      function sqroot(x) result(res)
        real(dp), intent(in) :: x
        real(dp)             :: res
        res = sqrt(x)
      end function sqroot

      function inv_sqroot(x) result(res)
        real(dp), intent(in) :: x
        real(dp)             :: res
        res = 1.0_dp/sqrt(x)
      end function inv_sqroot

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
           print *, "Range of Eigenvalues of S_k:", w(1), w(no_u)
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

        write(6,*) "Check of square root: Product (only valid when func=sqrt():"
        if ( any(  abs(B-S_k) > 1.0e-5) ) then
           print *, "S_sqroot is not a proper square root of S_k"
           error stop 1
        endif

        ! Compute inverse square root ----------------------
        do j = 1, no_u
           b(1:no_u,j) = S_base(1:no_u,j) / sqrt(w(j))
        enddo
        ! Now, compute B*Z^T
        call zgemm('N','C', no_u,no_u,no_u, (1.0_dp,0.0_dp), B,no_u, S_base,no_u, (0.0_dp,0.0_dp),S_inv_sqroot,no_u)

        B = matmul(S_inv_sqroot,S_inv_sqroot)
        S_base = matmul(B,S_k)  ! Note reuse

        write(6,*) "Check of S_inv_sqroot. It should be identity matrix"
        do i = 1, min(8,no_u)
           do j = 1, min(8,no_u)
              write(6,fmt="(1x,f10.4)",advance="no")  real(s_base(i,j))
           enddo
           write(6,*)
        enddo

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


end program fatband_atomic


