! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
!******************************************************************************
! MODULE m_fft
! 3-D fast complex Fourier transform
! Written by J.D.Gale (July 1999)
!******************************************************************************
!
!   PUBLIC procedures available from this module:
! fft   : 3-D complex FFT. Old version of J.D.Gale
!
!   PUBLIC parameters, types, and variables available from this module:
! none
!
!******************************************************************************
!
!   USED module procedures:
! use sys,          only: die           ! Termination routine
! use alloc,        only: de_alloc      ! De-allocation routine
! use alloc,        only: re_alloc      ! Re-allocation routine
! use fft1d,        only: gpfa          ! 1-D FFT routine
! use fft1d,        only: setgpfa       ! Sets gpfa routine
! use m_timer,      only: timer_start   ! Start counting CPU time
! use m_timer,      only: timer_stop    ! Stop counting CPU time
!
!   USED MPI procedures
! use mpi_siesta
!
!   USED module parameters:
! use precision,    only: dp            ! Real double precision type
! use precision,    only: gp=>grid_p    ! Real precision type of mesh arrays
! use parallel,     only: Node, Nodes, ProcessorY
! use mesh,         only: nsm
!
!   EXTERNAL procedures used:
! gpfa    : 1D complex FFT
! setgpfa : Initializes gpfa
! timer   : CPU time counter
!
!******************************************************************************

      MODULE m_fft

      use precision,    only : dp, grid_p
      use parallel,     only : Node, Nodes, processorY
      use moreMeshSubs, only : UNIFORM, getMeshBox
      use sys,          only : die
      use alloc,        only : re_alloc, de_alloc
      use mesh,         only : nsm
      use gpfa_fft,     only : gpfa        ! 1-D FFT routine
      use gpfa_fft,     only : setgpfa=>setgpfa_check     ! Sets gpfa routine
      use m_timer,      only : timer_start ! Start counting CPU time
      use m_timer,      only : timer_stop  ! Stop counting CPU time
#ifdef MPI
      use mpi_siesta
#endif
      use fftw3

      implicit none

      PUBLIC :: fft
      PUBLIC :: fft_init
      PUBLIC :: getFFTLims

      PRIVATE ! Nothing is declared public beyond this point
      logical :: frsttime = .TRUE.
#ifdef _FFTW_
      type(C_PTR), target :: plan_fw(3)
      type(C_PTR), target :: plan_bw(3)
#else
      integer :: maxtrigs = 256
      real(dp), pointer :: trigs(:,:)
#endif
      integer :: OldMesh(3) = (/ 0, 0, 0 /)

      integer :: gmesh(3)
      integer :: ProcessorZ
      integer :: PY
      integer :: PZ
      integer :: COMM_Y
      integer :: COMM_Z
      logical :: k_ori_distr = .false.
      integer, pointer :: xzdisp(:) => null()
      integer, pointer :: ydispl(:) => null()
      integer, pointer :: zdispl(:) => null()
      integer, pointer :: xdispl(:) => null()
      integer, pointer :: yzdisp(:) => null()

      !type FFT_ARR
      !  real(grid_p), pointer :: f1(:)
      !  real(grid_p), pointer :: f2(:)
      !end type FFT_ARR

      CONTAINS

!******************************************************************************

      subroutine fft_init( ori_distr, nMesh, size1, size2 )
      implicit none
      logical :: ori_distr                      ! Keep original distribution
      integer :: nMesh(3), size1, size2
      integer, pointer :: box(:,:,:)
      integer :: i, di, mo, s1, s2, s3
      integer :: ierr, maxmaxtrigs, ntrigs
      logical :: lredimension
#ifdef _FFTW_
      integer :: rank, n(1), howmany, inembed(1), istride, idist
      complex(C_DOUBLE_COMPLEX) :: in(1)
#endif
      gmesh(1:3) = nMesh(1:3)

      if (frsttime) then
#ifndef _FFTW_
        nullify(trigs)
        call re_alloc( trigs, 1, maxtrigs, 1, 3, 'trigs', 'fft' )
#endif
#ifdef MPI
        ProcessorZ = NODES/ProcessorY
        PY = mod(Node,ProcessorY)+1
        PZ = Node/ProcessorY+1
        call MPI_Comm_Split( MPI_Comm_World, PZ, PY, COMM_Y, ierr )
        call MPI_Comm_Split( MPI_Comm_World, PY, PZ, COMM_Z, ierr )
        call re_alloc( xdispl, 1, ProcessorY+1, 'xdispl', 'fft' )
        call re_alloc( ydispl, 1, ProcessorY+1, 'ydispl', 'fft' )
        call re_alloc( zdispl, 1, ProcessorZ+1, 'zdispl', 'fft' )
        call re_alloc( xzdisp, 1, ProcessorZ+1, 'xzdisp', 'fft' )
        yzdisp => xzdisp
#endif
        frsttime = .false.
      endif
#ifdef MPI
#else
        size1 = gMesh(1)*gMesh(2)*gMesh(3)
        size2 = 0
#endif

      if ( any(OldMesh/=gmesh) .or. ori_distr/=k_ori_distr ) then
C
C       Compute the limits of the differents pencils in each possible direction
C
        call getMeshBox( UNIFORM, box )
        di = gmesh(1)/ProcessorY
        mo = mod(gmesh(1),ProcessorY)
        xdispl(1) = 1
        do i= 1, ProcessorY
          xdispl(i+1) = xdispl(i) + di + merge(1,0,i<=mo)
        enddo
        ydispl(1) = 1
        do i= 1, ProcessorY
          ydispl(i+1) = box(2,2,i)*NSM+1
        enddo
        zdispl(1) = 1
        do i= 1, ProcessorZ
          zdispl(i+1) = box(2,3,(i-1)*ProcessorY+1)*NSM+1
        enddo
        if (ori_distr) then
          di = gmesh(1)/ProcessorZ
          mo = mod(gmesh(1),ProcessorZ)
          xzdisp(1) = 1
          do i= 1, ProcessorZ
            xzdisp(i+1) = xzdisp(i) + di + merge(1,0,i<=mo)
          enddo
        else
          di = gmesh(2)/ProcessorZ
          mo = mod(gmesh(2),ProcessorZ)
          yzdisp(1) = 1
          do i= 1, ProcessorZ
            yzdisp(i+1) = yzdisp(i) + di + merge(1,0,i<=mo)
          enddo
        endif
#ifdef _FFTW_
        if (any(OldMesh/=0)) then
          call fftw_destroy_plan( plan_fw(1) )
          call fftw_destroy_plan( plan_fw(2) )
          call fftw_destroy_plan( plan_fw(3) )
          call fftw_destroy_plan( plan_bw(1) )
          call fftw_destroy_plan( plan_bw(2) )
          call fftw_destroy_plan( plan_bw(3) )
        endif
        rank = 1
        n = (/ gmesh(1) /)
        howmany = (ydispl(py+1)-ydispl(py))*(zdispl(pz+1)-zdispl(pz))
        inembed = (/ gmesh(1) /)
        istride = 1
        idist   = gmesh(1)
        plan_fw(1) = fftw_plan_many_dft( rank, n, howmany,
     &                                    in, inembed, istride, idist,
     &                                    in, inembed, istride, idist,
     &                                    FFTW_FORWARD, FFTW_ESTIMATE )
        plan_bw(1) = fftw_plan_many_dft( rank, n, howmany,
     &                                    in, inembed, istride, idist,
     &                                    in, inembed, istride, idist,
     &                                    FFTW_BACKWARD, FFTW_ESTIMATE )
        n = (/ gmesh(2) /)
        howmany = (xdispl(py+1)-xdispl(py))
        inembed = (/ gmesh(2) /)
        istride = (xdispl(py+1)-xdispl(py))
        idist   = 1
        plan_fw(2) = fftw_plan_many_dft( rank, n, howmany,
     &                                    in, inembed, istride, idist,
     &                                    in, inembed, istride, idist,
     &                                    FFTW_FORWARD, FFTW_ESTIMATE )
        plan_bw(2) = fftw_plan_many_dft( rank, n, howmany,
     &                                    in, inembed, istride, idist,
     &                                    in, inembed, istride, idist,
     &                                    FFTW_BACKWARD, FFTW_ESTIMATE )
 
        n = (/ gmesh(3) /)
        if (ori_distr) then
          howmany = (xzdisp(pz+1)-xzdisp(pz))*(ydispl(py+1)-ydispl(py))
        else
          howmany = (xdispl(py+1)-xdispl(py))*(yzdisp(pz+1)-yzdisp(pz))
        endif
        istride = howmany
        idist   = 1
        inembed = (/ gmesh(3) /)
        plan_fw(3) = fftw_plan_many_dft( rank, n, howmany,
     &                                    in, inembed, istride, idist,
     &                                    in, inembed, istride, idist,
     &                                    FFTW_FORWARD, FFTW_ESTIMATE )
        plan_bw(3) = fftw_plan_many_dft( rank, n, howmany,
     &                                    in, inembed, istride, idist,
     &                                    in, inembed, istride, idist,
     &                                    FFTW_BACKWARD, FFTW_ESTIMATE )
        OldMesh = gmesh
#else
C
C       Initialise the tables for the FFT if the mesh has changed
C
        do
          lredimension = .false.
          maxmaxtrigs = maxtrigs
          do i=1,3
            if (OldMesh(i).ne.gmesh(i)) then
              call setgpfa( trigs(:,i), maxtrigs, ntrigs, gmesh(i) )
              if (ntrigs.gt.maxmaxtrigs) then
                lredimension = .true.
                maxmaxtrigs = ntrigs
              else
                OldMesh(i) = gmesh(i)
              endif
            endif
          enddo
          if (.not. lredimension) exit
C
C         Resize FFT array for trig factors and set OldMesh to 0 to force recalculation
C
          maxtrigs = maxmaxtrigs
          call re_alloc( trigs, 1, maxtrigs, 1, 3, 'trigs', 'fft' )
          OldMesh(1:3) = 0
        enddo
#endif
      endif
#ifdef MPI
        s1 = gMesh(1)*(ydispl(py+1)-ydispl(py))*
     &                (zdispl(pz+1)-zdispl(pz))
        s2 = gMesh(2)*(xdispl(py+1)-xdispl(py))*
     &                (zdispl(pz+1)-zdispl(pz))
      if (ori_distr) then
        s3 = gMesh(3)*(ydispl(py+1)-ydispl(py))*
     &                (xzdisp(pz+1)-xzdisp(pz))
        size1 = s1
        size2 = max(s2,s3)
      else
        s3 = gMesh(3)*(xdispl(py+1)-xdispl(py))*
     &       (yzdisp(pz+1)-yzdisp(pz))
        size1 = max(s1,s3)
        size2 = s2
      endif
#else
        size1 = gMesh(1)*gMesh(2)*gMesh(3)
        size2 = 0
#endif
      k_ori_distr = ori_distr
      end subroutine fft_init

      subroutine getFFTLims( lbox )
      implicit none
      integer :: lbox(2,3)
#ifdef MPI
      if (k_ori_distr) then
        lbox(1,1) = 1
        lbox(2,1) = gmesh(1)
        lbox(1,2) = ydispl(py) 
        lbox(2,2) = ydispl(py+1)-1
        lbox(1,3) = zdispl(pz)
        lbox(2,3) = zdispl(pz+1)-1
      else
        lbox(1,1) = xdispl(py)
        lbox(2,1) = xdispl(py+1)-1
        lbox(1,2) = yzdisp(pz)
        lbox(2,2) = yzdisp(pz+1)-1
        lbox(1,3) = 1
        lbox(2,3) = gmesh(3)
      endif
#else
      lbox(1,1:3) = 1
      lbox(2,1:3) = gMehs(1:3)
#endif
      end subroutine getFFTLims

      subroutine fft( f1, f2, isn )
      implicit none
      real(grid_p), target :: f1(*), f2(*)
      integer :: isn
      integer :: ierr, n, n1, n2, n3, n2a, n3a, n1b, n3b, n1c, n2c
      integer :: i, j, k, ng, ioffset
      real(dp) :: scale
#     ifdef _FFTW_
      complex(grid_p), pointer :: in(:)
      type(C_PTR), pointer :: plan(:)
#     endif
      call timer( 'fft', 1 )

      n1 = gmesh(1)
      n2 = gmesh(2)
      n3 = gmesh(3)
#ifdef MPI
#     ifdef _FFTW_
      if (isn<0) then
        plan => plan_fw
      else
        plan => plan_bw
      endif
#     endif
      if (k_ori_distr) then
        n2a = ydispl(py+1)-ydispl(py)
        n3a = zdispl(pz+1)-zdispl(pz)
        n1b = xdispl(py+1)-xdispl(py)
        n3b = n3a
        n1c = xzdisp(pz+1)-xzdisp(pz)
        n2c = n2a
C
C       FFT in X direction
C
        n = n1*n2a*n3a
        if (n>0) then  
#ifdef _FFTW_
          call c_f_pointer( c_loc(f1), in, shape=[n] )
          call fftw_execute_dft( plan(1), in, in )
#else
          call gpfa( f1(1:2*n), f1(2:2*n), trigs(:,1), 2, 2*n1, n1,
     &               n2a*n3a, -isn )
#endif
        endif
C
C       Transpose data from X to Y direction
C
        call redistributePencil( 1, n1, n2a, n3a, f1,
     &                           2, n1b, n2, n3b, f2,
     &                           py, ProcessorY, xdispl, ydispl,
     &                           COMM_Y )
C
C       FFT in Y direction
C
        n = n1b*n2*n3b
        if (n>0) then
          IOffSet = 1
          do i= 1, n3b
#           ifdef _FFTW_
            call c_f_pointer( c_loc(f2(IOffSet:IOffSet+n1b*n2-1)),
     &                  in, shape=[n1b*n2] )
            call fftw_execute_dft( plan(2), in, in )
#           else
            call gpfa(f2(IOffSet:2*n),f2(IOffSet+1:2*n), trigs(:,2),
     .            2*n1b, 2, n2, n1b, -isn )
#           endif
            IOffSet = Ioffset + 2*n1b*n2
          enddo
        endif
C
C       Transpose data from Y to X and then to Z direction
C
        call redistributePencil( 2, n1b, n2, n3b, f2,
     &                           1, n1, n2a, n3a, f1,
     &                           py, ProcessorY, ydispl, xdispl,
     &                           COMM_Y )
        call redistributePencil( 1, n1, n2a, n3a, f1,
     &                           3, n1c, n2c, n3, f2,
     &                           pz, ProcessorZ, xzdisp, zdispl,
     &                           COMM_Z )
C
C       FFT in Z direction
C
        n = n3*n1c*n2c
        if (n>0) then
#         ifdef _FFTW_
          call c_f_pointer( c_loc(f2), in, shape=[n] )
          call fftw_execute_dft( plan(3), in, in )
#         else
          call gpfa( f2(1:2*n), f2(2:2*n), trigs(:,3), 2*n1c*n2c, 2,
     &    n3, n1c*n2c,-isn )
#         endif
        endif
C
C       Transpose data back from Z to X
C
        call redistributePencil( 3, n1c, n2c, n3, f2,
     &                           1, n1, n2a, n3a, f1,
     &                           pz, ProcessorZ, zdispl, xzdisp,
     &                           COMM_Z )
C
C       Scale values
C
        if (isn.gt.0) then
          ng = n1*n2*n3
          scale=1.0_dp/dble(ng)
          do i=1,2*n
            f1(i)=f1(i)*scale
          enddo
        endif
      else ! .not. k_ori_distr
        n2a = ydispl(py+1)-ydispl(py)
        n3a = zdispl(pz+1)-zdispl(pz)
        n1b = xdispl(py+1)-xdispl(py)
        n3b = n3a
        n1c = n1b
        n2c = yzdisp(pz+1)-yzdisp(pz)
        if (isn==-1) then
C
C         FFT in X direction
C
          n = n1*n2a*n3a
          if (n>0) then
#ifdef _FFTW_
            call c_f_pointer( c_loc(f1), in, shape=[n] )
            call fftw_execute_dft( plan(1), in, in )
#else
            call gpfa( f1(1:2*n), f1(2:2*n), trigs(:,1), 2, 2*n1, n1,
     &                  n2a*n3a, -isn )
#endif
           endif
C
C         Transpose data from X to Y direction
C
          call redistributePencil( 1, n1, n2a, n3a, f1,
     &                             2, n1b, n2, n3b, f2,
     &                             py, ProcessorY, xdispl, ydispl,
     &                             COMM_Y )
C
C         FFT in Y direction
C
          n = n1b*n2*n3b
          if (n>0) then
            IOffSet = 1
            do i= 1, n3b
#ifdef _FFTW_
              call c_f_pointer( c_loc(f2(IOffSet:IOffSet+n1b*n2-1)),
     &                  in, shape=[n1b*n2] )
              call fftw_execute_dft( plan(2), in, in )
#else
              call gpfa(f2(IOffSet:2*n),f2(IOffSet+1:2*n), trigs(:,2),
     .              2*n1b, 2, n2, n1b, -isn )
#endif
              IOffSet = Ioffset + 2*n1b*n2
            enddo
          endif
C
C         Transpose data from Y to Z direction
C
          call redistributePencil( 2, n1b, n2, n3b, f2,
     &                             3, n1c, n2c, n3, f1,
     &                             pz, ProcessorZ, yzdisp, zdispl,
     &                             COMM_Z )
C
C         FFT in Z direction
C
          n = n3*n1c*n2c
          if (n>0) then
#ifdef 	    _FFTW_
            call c_f_pointer( c_loc(f1), in, shape=[n] )
            call fftw_execute_dft( plan(3), in, in )
#else
            call gpfa( f1(1:2*n), f1(2:2*n), trigs(:,3), 2*n1c*n2c, 2,
     &      n3, n1c*n2c,-isn )
#endif
          endif
        else ! isn==1
C
C         FFT in Z direction
C
          n = n3*n1c*n2c
          if (n>0) then
#ifdef      _FFTW_
            call c_f_pointer( c_loc(f1), in, shape=[n] )
            call fftw_execute_dft( plan(3), in, in )
#else
            call gpfa( f1(1:2*n), f1(2:2*n), trigs(:,3), 2*n1c*n2c, 2,
     &      n3, n1c*n2c,-isn )
#endif
          endif
C
C         Transpose data from Z to Y direction
C
          call redistributePencil( 3, n1c, n2c, n3, f1,
     &                             2, n1b, n2, n3b, f2,
     &                             pz, ProcessorZ, zdispl, yzdisp,
     &                             COMM_Z )
C
C         FFT in Y direction
C
          n = n1b*n2*n3b
          if (n>0) then
            IOffSet = 1
            do i= 1, n3b
#ifdef        _FFTW_
              call c_f_pointer( c_loc(f2(IOffSet:IOffSet+n1b*n2-1)),
     &                  in, shape=[n1b*n2] )
              call fftw_execute_dft( plan(2), in, in )
#else
              call gpfa(f2(IOffSet:2*n),f2(IOffSet+1:n), trigs(:,2),
     .              2*n1b, 2, n2, n1b, -isn )
#endif
              IOffSet = Ioffset + 2*n1b*n2
            enddo
          endif
C
C         Transpose data from Y to X direction
C
          call redistributePencil( 2, n1b, n2, n3b, f2,
     &                             1, n1, n2a, n3a, f1,
     &                             py, ProcessorY, ydispl, xdispl,
     &                             COMM_Y )
C
C         FFT in X direction
C
          n = n1*n2a*n3a
          if (n>0) then
#ifdef _FFTW_
            call c_f_pointer( c_loc(f1), in, shape=[n] )
            call fftw_execute_dft( plan(1), in, in )
#else
            call gpfa( f1(1:2*n), f1(2:2*n), trigs(:,1), 2, 2*n1, n1,
     &                  n2a*n3a, -isn )
#endif
C
C           SCALE the result
C
            ng = n1*n2*n3
            scale=1.0_dp/dble(ng)
            do i=1,2*n
              f1(i)=f1(i)*scale
            enddo
          endif
        endif
      endif
#else
        n = n1*n2*n3
C
C       FFT in X direction
C
#ifdef _FFTW_
        call c_f_pointer( c_loc(f1), in, shape=[n] )
        call fftw_execute_dft( plan(1), in, in )
#else
        call gpfa( f1(1:2*n), f1(2:2*n), trigs(:,1), 2, 2*n1, n1,
     &             n2*n3, -isn )
#endif
C
C  FFT in Y direction
C
      IOffSet = 1
      do i=0,n3-1
#       ifdef  _FFTW_
        call c_f_pointer( c_loc(f1(IOffSet:IOffSet+n1*n2-1)),
     &              in, shape=[n1*n2] )
        call fftw_execute_dft( plan(2), in, in )
#       else
        call gpfa(f1(IOffSet:2*n),f1(IOffSet:2*n),trigs(:,2),
     .             2*n1,2,n2,n1,-isn)
#       endif
        IOffSet = Ioffset + 2*n1*n2
      enddo
C
C  FFT in Z direction
C
#     ifdef  _FFTW_
      call c_f_pointer( c_loc(f1), in, shape=[n] )
      call fftw_execute_dft( plan(3), in, in )
#     else
      call gpfa(f1(1:2*n),f1(2:2*n),trigs(:,3),2*n1*n2,2,n3,n1*n2,-isn)
#     endif
C
C  Scale values
C
      if (isn.gt.0) then
        scale=1.0_dp/dble(n)
        do i=1,2*n
          f(i)=f(i)*scale
        enddo
      endif
#endif /* MPI */
      call timer( 'fft', 2 )
      end subroutine fft

#ifdef MPI
      subroutine redistributePencil( src, s1, s2, s3, fsrc,
     &                               dst, d1, d2, d3, fdst,
     &                               p, np, srcd, dstd, comm )
      implicit none
      integer :: src, dst, s1, s2, s3, d1, d2, d3
      integer :: p, np, srcd(np+1), dstd(np+1), comm
      real*8 :: fsrc(2,s1,s2,s3), fdst(2,d1,d2,d3)
      integer :: i, j, k, l, ii, jj, kk, pi, pl, pr, nsend, nrecv
      integer :: ierr, r_recv, r_send, status(MPI_Status_Size)
      integer :: bsize_s, bsize_r
      integer :: sbox(2,3), dbox(2,3)
      real(grid_p), pointer  :: sbuf(:,:), rbuf(:,:)
      sbox(1,:) = 1
      sbox(2,:) = (/ s1, s2, s3 /)
      dbox(1,:) = 1
      dbox(2,:) = (/ d1, d2, d3 /)
C     Compute the buffer sizes
      bsize_s = (s1*s2*s3/sbox(2,src))*(srcd(2)-srcd(1))
      bsize_r = (d1*d2*d3/dbox(2,dst))*(dstd(2)-dstd(1))
      call re_alloc( sbuf, 1, 2, 1, bsize_s, 'sbuf', 'fft' )
      call re_alloc( rbuf, 1, 2, 1, bsize_r, 'rbuf', 'fft' )
C
C     Handle transfer of terms which are purely local
C
      sbox(1,src) = srcd(p)
      sbox(2,src) = srcd(p+1)-1
      dbox(1,dst) = dstd(p)
      dbox(2,dst) = dstd(p+1)-1
      ii = dbox(1,3)
      do i= sbox(1,3), sbox(2,3)
        jj = dbox(1,2)
        do j= sbox(1,2), sbox(2,2)
          kk = dbox(1,1)
          do k= sbox(1,1), sbox(2,1)
            fdst(:,kk,jj,ii) = fsrc(:,k,j,i)
            kk = kk+1
          enddo
          jj = jj+1
        enddo
        ii =ii+1
      enddo
C
C  Loop over all Node-Node vectors exchanging local data
C
      do pi=1, NP-1
        ! Let's receive from left node and send to the right node
        pl = mod(p-1+NP-pi,NP)+1
        pr = mod(p-1+pi,NP)+1
        ! Destination boxes
        dbox(1,dst) = dstd(pl)
        dbox(2,dst) = dstd(pl+1)-1
        sbox(1,src) = srcd(pr)
        sbox(2,src) = srcd(pr+1)-1
        ! Number of elements to send/receive
        nsend = 2*PRODUCT( sbox(2,:)-sbox(1,:)+1 )
        nrecv = 2*PRODUCT( dbox(2,:)-dbox(1,:)+1 )
C
C       Collect data to send
C
        l = 1
        do i= sbox(1,3), sbox(2,3)
          do j= sbox(1,2), sbox(2,2)
            do k= sbox(1,1), sbox(2,1)
              sbuf(:,l) = fsrc(:,k,j,i)
              l = l+1
            enddo
          enddo
        enddo
C
C       Exchange data - send to right and receive from left
C
        call MPI_IRecv( rbuf(1,1), nrecv,
     .    MPI_grid_real, pl-1, 1, COMM, r_recv, ierr )
        call MPI_ISend( sbuf(1,1), nsend,
     .    MPI_grid_real, pr-1, 1, COMM, r_send, ierr )
C
C       Wait for receive to complete
C
        call MPI_Wait( r_recv, Status, ierr )
C
C       Place received data into correct array
C
        l = 1
        do i= dbox(1,3), dbox(2,3)
          do j= dbox(1,2), dbox(2,2)
            do k= dbox(1,1), dbox(2,1)
              fdst(:,k,j,i) = rbuf(:,l)
              l = l+1
            enddo
          enddo
        enddo
C
C  Wait for send to complete
C
        call MPI_Wait( r_send, Status, ierr )
      enddo

      call de_alloc( sbuf, 'sbuf', 'fft' ) 
      call de_alloc( rbuf, 'rbuf', 'fft' )
      end subroutine redistributePencil
#endif /* MPI */

!******************************************************************************

      END MODULE m_fft
