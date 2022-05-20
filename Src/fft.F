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

      implicit none

      PUBLIC:: 
     .  fft      ! 3D complex FFT

      PRIVATE ! Nothing is declared public beyond this point
      logical :: frsttime = .TRUE.
      integer :: maxtrigs = 256
      real(dp), pointer :: trigs(:,:)
      integer :: OldMesh(3) = (/ 0, 0, 0 /)

      integer :: ProcessorZ
      integer :: PY
      integer :: PZ
      integer :: COMM_Y
      integer :: COMM_Z
      integer, pointer :: xydisp(:)
      integer, pointer :: xzdisp(:)
      integer, pointer :: ydispl(:)
      integer, pointer :: zdispl(:)

      CONTAINS

!******************************************************************************

      subroutine fft( f, nMesh, isn )
C
C  Parallel FFT based on FFT routine of Templeton
C
C  On input :
C
C  real*8 f()       : contains data to be Fourier transformed
C  integer nMesh(3) : contains global dimensions of grid
C  integer isn      : indicates the direction of the transform
C
C  On exit :
C
C  real*8 f()      : contains data Fourier transformed
C
C  Julian Gale, July 1999
C
      implicit none
C
C  Passed arguments
C
      real(grid_p) :: f(*)
      integer :: nMesh(3), isn
C
C  Local variables
C 
      integer :: n1, n2, n3, ng, n1l, n2l, n3l, n, n1l_y, n1l_z
      integer :: i, di, mo, ntrigs, bsize, nft, ioffset, j, ii
      integer :: maxmaxtrigs, ierr
      integer, pointer :: box(:,:,:)
      logical :: LocalPoints, lredimension
      real(grid_p), pointer :: ft(:)
      real(dp) :: scale

      call timer( 'fft', 1 )
C
C     Global dimensions
C
      n1 = nMesh(1)
      n2 = nMesh(2)
      n3 = nMesh(3)
      ng = n1*n2*n3
C
C     Local dimensions
C
      call getMeshBox( UNIFORM, box )
      !n1l = (box(2,1,Node+1)-box(1,1,Node+1)+1)*NSM
      n2l = (box(2,2,Node+1)-box(1,2,Node+1)+1)*NSM
      n3l = (box(2,3,Node+1)-box(1,3,Node+1)+1)*NSM
      n = n1*n2l*n3l
C
C  Set logical as to whether there are any locally stored points
C
      LocalPoints = (n.gt.0) 
      if (frsttime) then
        nullify(trigs)
        call re_alloc( trigs, 1, maxtrigs, 1, 3, 'trigs', 'fft' )
        frsttime = .FALSE.
#ifdef MPI
        ProcessorZ = NODES/ProcessorY
        PY = mod(Node,ProcessorY)+1
        PZ = Node/ProcessorY+1
        call MPI_Comm_Split( MPI_Comm_World, PZ, PY, COMM_Y, ierr )
        call MPI_Comm_Split( MPI_Comm_World, PY, PZ, COMM_Z, ierr )
        nullify(xydisp,xzdisp,ydispl,zdispl)
        call re_alloc( xydisp, 1, ProcessorY+1, 'xydisp', 'fft' )
        call re_alloc( xzdisp, 1, ProcessorZ+1, 'xzdisp', 'fft' )
        call re_alloc( ydispl, 1, ProcessorY+1, 'ydispl', 'fft' )
        call re_alloc( zdispl, 1, ProcessorZ+1, 'zdispl', 'fft' )
        
        di = n1/ProcessorY
        mo = mod(n1,ProcessorY)
        xydisp(1) = 1
        do i= 1, ProcessorY
          xydisp(i+1) = xydisp(i) + di + merge(1,0,i<=mo)
        enddo
        di = n1/ProcessorZ
        mo = mod(n1,ProcessorZ)
        xzdisp(1) = 1
        do i= 1, ProcessorZ
          xzdisp(i+1) = xzdisp(i) + di + merge(1,0,i<=mo)
        enddo
        ydispl(1) = 1
        do i= 1, ProcessorY
          ydispl(i+1) = box(2,2,i)*NSM+1
        enddo
        zdispl(1) = 1
        do i= 1, ProcessorZ
          zdispl(i+1) = box(2,3,(i-1)*ProcessorY+1)*NSM+1
        enddo
#endif
      endif
C
C  Initialise the tables for the FFT if the mesh has changed
C
      do
        lredimension = .false.
        maxmaxtrigs = maxtrigs
        do i=1,3
          if (OldMesh(i).ne.nMesh(i)) then
            call setgpfa( trigs(:,i), maxtrigs, ntrigs, nMesh(i) )
            if (ntrigs.gt.maxmaxtrigs) then
              lredimension = .true.
              maxmaxtrigs = ntrigs
            else
              OldMesh(i) = nMesh(i)
            endif
          endif
        enddo
        if (.not. lredimension) exit
C
C  Resize FFT array for trig factors and set OldMesh to 0 to force recalculation
C
        maxtrigs = maxmaxtrigs
        call re_alloc( trigs, 1, maxtrigs, 1, 3, 'trigs', 'fft' )
        OldMesh(1:3) = 0
      enddo
C
C  FFT in X direction
C
      if (LocalPoints) then
        call gpfa(f(1:2*n),f(2:2*n),trigs(:,1),2,2*n1,n1,n2l*n3l,-isn)
      endif

#ifdef MPI
C***********************
C  2-D Processor Grid  *
C***********************
      n1l_y = xydisp(py+1)-xydisp(py)
      n1l_z = xzdisp(pz+1)-xzdisp(pz)
C
C  Allocate local memory
C
      nullify( ft )
      bsize = 2*max(n1l_y*n2*n3l,n1l_z*n2l*n3)
      call re_alloc( ft, 1, bsize, 'ft', 'fft' )
      n1l = n1l_y

C
C  Redistribute data to be distributed by X and Z
C
      call redistribXZ( f, n1, n2l, n3l, ft, n1l, n2, 1,
     &    nsm, Node, Nodes )
C
C  FFT in Y direction
C
      if (LocalPoints) then
        nft = 2*n1l*n2*n3l
        do i= 0, n3l-1
          IOffSet=2*n1l*n2*i
          call gpfa(ft(IOffSet+1:nft),ft(IOffSet+2:nft),trigs(:,2),
     .            2*n1l, 2, n2, n1l, -isn )
        enddo
      endif
C
C  Redistribute data back to original form
C
      call redistribXZ( f, n1, n2l, n3l, ft, n1l, n2,
     &    -1, nsm, Node, Nodes )
C
C  Redistribute data to be distributed by X and Y
C
      n1l = n1l_z
      call redistribXY( f, n1, n2l, n3l, ft, n1l, n3,
     &    1, nsm, Node, Nodes )
C
C  FFT in Z direction
C
      if (LocalPoints) then
        nft = 2*n1l*n2l*n3
        call gpfa(ft(1:nft),ft(2:nft),trigs(:,3),
     $            2*n1l*n2l,2,n3,n1l*n2l,-isn)
      endif
C
C  Redistribute data to be distributed by Z again
C
      call redistribXY(f,n1,n2l,n3l,ft,n1l,n3,-1,nsm,Node,Nodes)
C
C  Free local memory
C
      call de_alloc( ft, 'ft', 'fft' )
#else
C
C  FFT in Y direction
C
      do i=0,n3-1
         IOffSet=2*n1*n2*i
         call gpfa(f(IOffSet+1:nf),f(IOffSet+2:nf),trigs(:,2),
     .             2*n1,2,n2,n1,-isn)
      enddo
C
C  FFT in Z direction
C
      call gpfa(f(1:nf),f(2:nf),trigs(:,3),2*n1*n2,2,n3,n1*n2,-isn)
#endif /* MPI */
C
C  Scale values
C
      if (LocalPoints.and.isn.gt.0) then
        scale=1.0_dp/dble(ng)
        do i=1,2*n
          f(i)=f(i)*scale
        enddo
      endif

      call timer( 'fft', 2 )
      return

      end subroutine fft





!******************************************************************************

#ifdef MPI
      subroutine redistribXZ(f,n1,n2l,n3l,ft,n1l,n2,idir,nsm,
     .  Node,Nodes)
C
C  This routine redistributes the data over the Nodes as needed
C  for the FFT routines between the arrays f and ft
C
C  Array f is distributed in the Y/Z direction while ft is distributed
C  in the X/Z direction
C
C  idir = direction of redistribution : > 0 => f->ft
C                                       < 0 => ft->f
C
C  Use the processor grid to divide communicator according to Z
C
C  Julian Gale, July 1999
C

      implicit none

C
C  Passed arguments
C
      integer
     .  n1, n2, n2l, n1l, n3l, Node, Nodes, idir, nsm
      real(grid_p)
     .  f(2,n1,n2l,n3l), ft(2,n1l,n2,n3l)
C
C  Local variables
C
      integer :: i, j, k, j1, k1, pl, pr, l, p
      integer :: maxn1l, maxn2l, n1t, n2t
      integer :: ierr, r_send, r_recv, Status(MPI_Status_Size)
      real(grid_p), pointer  :: fbuf(:,:), tbuf(:,:)
C
C  Allocate local memory
C
      maxn1l = xydisp(2)-xydisp(1)
      maxn2l = ydispl(2)-ydispl(1)
      nullify(fbuf,tbuf)
      call re_alloc( fbuf, 1, 2, 1, maxn1l*n2l*n3l, 'fbuf', 'fft' ) 
      call re_alloc( tbuf, 1, 2, 1, n1l*maxn2l*n3l, 'tbuf', 'fft' )

      if (idir.gt.0) then
C***********
C  F -> FT *
C***********
C
C  Handle transfer of terms which are purely local
C
        do i= 1, n3l
          j1 = ydispl(py)
          do j= 1, n2l
            k1 = 1
            do k=xydisp(py), xydisp(py+1)-1
              ft(:,k1,j1,i) = f(:,k,j,i)
              k1 = k1 + 1
            enddo
            j1 = j1 + 1
          enddo
        enddo
C  Loop over all Node-Node vectors exchanging local data
C
        do p=1, ProcessorY-1
          ! Let's receive from left node and send to the right node
          pl = mod(py-1+ProcessorY-p,ProcessorY)+1
          pr = mod(py-1+p,ProcessorY)+1
C
C  Collect data to send
C
          l = 1
          do i= 1, n3l
            do j= 1, n2l
              do k= xydisp(pr), xydisp(pr+1)-1
                fbuf(:,l) = f(:,k,j,i)
                l = l+1
              enddo
            enddo
          enddo
C
C  Exchange data - send to right and receive from left
C
          n2t = ydispl(pl+1)-ydispl(pl)
          n1t = xydisp(pr+1)-xydisp(pr)
          call MPI_IRecv( tbuf(1,1), 2*n1l*n2t*n3l,
     .      MPI_grid_real, pl-1, 1, COMM_Y, r_recv, ierr )
          call MPI_ISend( fbuf(1,1), 2*n1t*n2l*n3l,
     .      MPI_grid_real, pr-1, 1, COMM_Y, r_send, ierr )
C
C  Wait for receive to complete
C
          call MPI_Wait( r_recv, Status, ierr )
C
C  Place received data into correct array
C
          l = 1
          do i= 1, n3l
            do j= ydispl(pl), ydispl(pl+1)-1
              do k= 1, n1l
                ft(:,k,j,i) = tbuf(:,l)
                l = l+1
              enddo
            enddo
          enddo
C  Wait for send to complete
C
          call MPI_Wait( r_send, Status, ierr )
        enddo
      else
C***********
C  FT -> F *
C***********
C
C  Handle transfer of terms which are purely local
        do i= 1, n3l
          j1 = ydispl(py)
          do j= 1, n2l
            k1 = 1
            do k=xydisp(py), xydisp(py+1)-1
              f(:,k,j,i) = ft(:,k1,j1,i)
              k1 = k1 + 1
            enddo
            j1 = j1 + 1
          enddo
        enddo
C  Loop over all Node-Node vectors exchanging local data
C
        do p=1, ProcessorY-1
          ! Let's receive from left node and send to the right node
          pl = mod(py-1+ProcessorY-p,ProcessorY)+1
          pr = mod(py-1+p,ProcessorY)+1
C
C  Collect data to send
C
          l = 1
          do i= 1, n3l
            do j= ydispl(pr), ydispl(pr+1)-1
              do k= 1, n1l
                tbuf(:,l) = ft(:,k,j,i)
                l = l+1
              enddo
            enddo
          enddo

C  Exchange data - send to right and receive from left
C
          n2t = ydispl(pr+1)-ydispl(pr)
          n1t = xydisp(pl+1)-xydisp(pl)
          call MPI_IRecv( fbuf(1,1),  2*n1t*n2l*n3l,
     .      MPI_grid_real, pl-1, 1, COMM_Y, r_recv, ierr )
          call MPI_ISend( tbuf(1,1), 2*n1l*n2t*n3l,
     .      MPI_grid_real, pr-1, 1, COMM_Y, r_send, ierr )
C
C  Wait for receive to complete
C
          call MPI_Wait( r_recv, Status, ierr )
C
C  Place received data into correct array
          l = 1
          do i= 1, n3l
            do j= 1, n2l
              do k= xydisp(pl), xydisp(pl+1)-1
                f(:,k,j,i) = fbuf(:,l)
                l = l+1
              enddo
            enddo
          enddo
C  Wait for send to complete
C
          call MPI_Wait( r_send, Status, ierr )
        enddo
      endif
C
C  Free local memory
C
      call de_alloc( fbuf, 'fbuf', 'fft' )
      call de_alloc( tbuf, 'tbuf', 'fft' )
      end subroutine redistribXZ

!******************************************************************************

      subroutine redistribXY(f,n1,n2l,n3l,ft,n1l,n3,idir,nsm,
     .  Node,Nodes)
C
C  This routine redistributes the data over the Nodes as needed
C  for the FFT routines between the arrays f and ft
C
C  Array f is distributed in the Y/Z direction while ft is distributed
C  in the X/Y direction
C
C  idir = direction of redistribution : > 0 => f->ft
C                                       < 0 => ft->f
C
C  Use the processor grid to divide communicator according to Y
C
C  Currently written in not the most efficient but simple way! 
C  Need to improve communication later.
C
C  Julian Gale, July 1999
C

      implicit none

C
C  Passed arguments
C
      integer
     .  n1, n3, n2l, n1l, n3l, nsm, Node, Nodes, idir
      real(grid_p)
     .  f(2,n1,n2l,n3l), ft(2,n1l,n2l,n3)
C
C  Local variables
C
      integer :: i, j, k, i1, k1, pl, pr, l, p
      integer :: maxn1l, maxn3l, n1t, n3t
      integer :: ierr, r_send, r_recv, Status(MPI_Status_Size)
      real(grid_p), pointer  :: fbuf(:,:), tbuf(:,:)
C
C  Allocate local memory
C
      maxn1l = xzdisp(2)-xzdisp(1)
      maxn3l = zdispl(2)-zdispl(1)
      nullify(fbuf,tbuf)
      call re_alloc( fbuf, 1, 2, 1, maxn1l*n2l*n3l, 'fbuf', 'fft' ) 
      call re_alloc( tbuf, 1, 2, 1, n1l*n2l*maxn3l, 'tbuf', 'fft' )
      if (idir.gt.0) then
C***********
C  F -> FT *
C***********
C
C  Handle transfer of terms which are purely local
        i1 = zdispl(pz)
        do i= 1, n3l
          do j= 1, n2l
            k1 = 1
            do k=xzdisp(pz), xzdisp(pz+1)-1
              ft(:,k1,j,i1) = f(:,k,j,i)
              k1 = k1 + 1
            enddo
          enddo
            i1 = i1 + 1
        enddo
C  Loop over all Node-Node vectors exchanging local data
C
        do p=1, ProcessorZ-1
          ! Let's receive from left node and send to the right node
          pl = mod(pz-1+ProcessorZ-p,ProcessorZ)+1
          pr = mod(pz-1+p,ProcessorZ)+1
C
C  Collect data to send
          l = 1
          do i= 1, n3l
            do j= 1, n2l
              do k= xzdisp(pr), xzdisp(pr+1)-1
                fbuf(:,l) = f(:,k,j,i)
                l = l+1
              enddo
            enddo
          enddo
C  Exchange data - send to right and receive from left
C
          n3t = zdispl(pl+1)-zdispl(pl)
          n1t = xzdisp(pr+1)-xzdisp(pr)
          call MPI_IRecv( tbuf(1,1), 2*n1l*n2l*n3t,
     .      MPI_grid_real, pl-1, 1, COMM_Z, r_recv, ierr )
          call MPI_ISend( fbuf(1,1), 2*n1t*n2l*n3l,
     .      MPI_grid_real, pr-1, 1, COMM_Z, r_send, ierr )
C
C  Wait for receive to complete
C
          call MPI_Wait( r_recv, Status, ierr )
C
C  Place received data into correct array
C
          l = 1
          do i= zdispl(pl), zdispl(pl+1)-1
            do j=1, n2l 
              do k= 1, n1l
                ft(:,k,j,i) = tbuf(:,l)
                l = l+1
              enddo
            enddo
          enddo
C  Wait for send to complete
C
          call MPI_Wait( r_send, Status, ierr )
        enddo
      elseif (idir.lt.0) then
C***********
C  FT -> F *
C***********
C
C  Handle transfer of terms which are purely local
        i1 = zdispl(pz)
        do i= 1, n3l
          do j= 1, n2l
            k1 = 1
            do k=xzdisp(pz), xzdisp(pz+1)-1
              f(:,k,j,i) = ft(:,k1,j,i1)
              k1 = k1 + 1
            enddo
          enddo
          i1 = i1 + 1
        enddo
C  Loop over all Node-Node vectors exchanging local data
C
        do p=1, ProcessorZ-1
          ! Let's receive from left node and send to the right node
          pl = mod(pz-1+ProcessorZ-p,ProcessorZ)+1
          pr = mod(pz-1+p,ProcessorZ)+1
C
C  Collect data to send
C
          l = 1
          do i= zdispl(pr), zdispl(pr+1)-1
            do j= 1, n2l
              do k= 1, n1l
                tbuf(:,l) = ft(:,k,j,i)
                l = l+1
              enddo
            enddo
          enddo
C
C  Exchange data - send to right and receive from left
C
          n3t = zdispl(pr+1)-zdispl(pr)
          n1t = xzdisp(pl+1)-xzdisp(pl)
          call MPI_IRecv( fbuf(1,1),  2*n1t*n2l*n3l,
     .      MPI_grid_real, pl-1, 1, COMM_Z, r_recv, ierr )
          call MPI_ISend( tbuf(1,1), 2*n1l*n2l*n3t,
     .      MPI_grid_real, pr-1, 1, COMM_Z, r_send, ierr )
C
C  Wait for receive to complete
C
          call MPI_Wait( r_recv, Status, ierr )
C
C  Place received data into correct array
C
          l = 1
          do i= 1, n3l
            do j= 1, n2l
              do k= xzdisp(pl), xzdisp(pl+1)-1
                f(:,k,j,i) = fbuf(:,l)
                l = l+1
              enddo
            enddo
          enddo
C  Wait for send to complete
C
          call MPI_Wait( r_send, Status, ierr )
        enddo
      endif
C
C  Free local memory
C
      call de_alloc( fbuf, 'fbuf', 'fft' )
      call de_alloc( tbuf, 'tbuf', 'fft' )
      end subroutine redistribXY
#endif /* MPI */

!******************************************************************************

      END MODULE m_fft
