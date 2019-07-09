 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_dscfComm
C ----------------------------------------------------------------------
C Stores quantities used to transfer data between a matrix that is
C distributed by orbitals and one distributed by mesh
C
C Written by Rogeli Grima (BSC) Dec.2007
C
C ----------------------------------------------------------------------
C integer DCsize        : Number of communications
C integer DCself        : Number of orbitals that the current process should keep
C integer DCtotal       : Total Number of orbitals of the local matrix
C integer DCmax         : Maximum number of orbitals that we should send/receive
C                         from/to another process
C integer DCmaxnd       : Maximum size of matrix that we should send/receive
C                         from/to another process
C integer DCBuffer(:)   : Buffer that contains all the data
C integer DCsrc(DCsize) : Process ID of the source.
C integer DCdst(DCsize) : Process ID of the destiny.
C integer DCsic(DCsize) : Number of orbitals to send/receive
C integer DCinvp(DCtotal) : Permutation vector used to transform the received
C                           data to the local ordering
C ----------------------------------------------------------------------
C
C     Modules
      use precision
      implicit none
      
      integer, public,  save :: DCsize, DCmax, DCmaxnd
      integer, public, pointer, save :: DCsrc(:), DCdst(:),
     $                                  DCsic(:)
      integer, pointer, save :: DCBuffer(:)

      integer, public          :: DCself, DCncom, DCtotal
      integer, public, pointer :: DCsnd(:), DCrcv(:), DCinvp(:),
     &                            DCpid(:), DCreq(:,:), DCsne(:),
     &                            DCrce(:)

      public :: dscfComm, resetdscfComm, CommWaitAll, CommWaitSnd,
     &          CommWaitRcv
      private

      CONTAINS

      subroutine dscfComm( nuotot, nrowsDscfL, NeedDscfL )
C ==================================================================
C Precomputes the needed communications to transform a matrix
C that is distributed by orbitals to one distributed by mesh
C ==================================================================
C subroutine dscfComm( nuotot, nrowsDscfL, NeedDscfL )
C
C INPUT:
C integer nuotot            : Total number of basis orbitals in unit cell
C integer nrowsDscfL        : Local number of rows of matrix DSCF
C integer NeedDscfL(nuotot) : Permutation vector of orbitals from global to local.
C
C OUTPUT:
C The output values are in the module m_dscfComm
C
C BEHAVIOR:
C Count the number of orbitals that we should receive from every process.
C Create a list of communications that every process should do. The process
C 0 receive all the communications and calls scheduleComm to planify
C the communication schedule. The rest of processes receive the schedule
C from process 0.
C
C ==================================================================
      use precision
      use parallel,     only : NODE, NODES
      use parallelsubs, only : GlobalToLocalOrb, WhichNodeOrb
      use alloc
      use scheComm
#ifdef MPI
      use mpi_siesta
      use m_mpi_utils, only: globalize_sum
#endif
      implicit none
C     Input variables
      integer            :: nuotot, nrowsDscfL
      integer            :: NeedDscfL(nuotot)
C     Local variables
      integer            :: ii, ind, io, oNode, ncom
      integer            :: MPIerr
      integer, pointer   :: send(:), recv(:), tinv(:), tsrc(:), xnei(:),
     &                      invp(:)
!------------------------------------------------------------------------- BEGIN
      call timer( "dscfComm", 1 )
      ! allocate memory for temporal buffers
      nullify(send,recv,tinv,tsrc)
      call re_alloc( send, 0, NODES-1, 'send', 'dscfComm' )
      call re_alloc( recv, 0, NODES-1, 'recv', 'dscfComm' )
      call re_alloc( tinv, 1, nrowsDscfL, 'tinv', 'dscfComm' )
      call re_alloc( tsrc, 1, nrowsDscfL, 'tsrc', 'dscfComm' )
      do ii = 0, NODES-1
        send(ii) = 0
        recv(ii) = 0
      enddo
      ii = 0
      do io = 1, nuotot
        if ( NeedDscfL( io ).ne.0 ) then
          call WhichNodeOrb( io, Nodes, oNode )
          recv(oNode) = recv(oNode) + 1
          ii = ii + 1
          tinv(ii) = io
          tsrc(ii) = oNode
        endif
      enddo
      ! 
      nullify(xnei)
      call re_alloc( xnei, 0, NODES, 'xnei', 'dscfComm' )
      xnei(0) = 0
      do ii= 1, NODES
        xnei(ii) = xnei(ii-1) + recv(ii-1)
      enddo
      ! 
      nullify(invp)
      call re_alloc( invp, 1, nrowsDscfL, 'invp', 'dscfComm' )
      do io= 1, nrowsDscfL
        ii    = tinv(io)
        oNode = tsrc(io)
        xnei(oNode) = xnei(oNode) + 1
        invp(xnei(oNode)) = ii
      enddo
      call de_alloc( tinv, 'tinv', 'dscfComm' )
      call de_alloc( tsrc, 'tsrc', 'dscfComm' )
      xnei(0) = 0
      do ii= 1, NODES
        xnei(ii) = xnei(ii-1) + recv(ii-1)
      enddo

#ifdef MPI
      call MPI_Alltoall( recv, 1, MPI_Integer, send, 1, MPI_Integer,
     &                   MPI_COMM_WORLD, MPIerr )
#endif

      ncom = 0
      do ii = 0, NODES-2
        oNode = MOD((ii/2+1)*
     &             MERGE(1,-1,IAND(NODE,ishft(1,ii/2))==0)*
     &             MERGE(1,-1,IAND(ii,1)==0)+NODES+NODE,NODES)
        if (send(oNode)/=0 .or. recv(oNode)/=0) ncom = ncom + 1
      enddo
      DCncom = ncom

      call re_alloc( DCpid, 1, DCncom, 'DCpid', 'dscfComm' )
      call re_alloc( DCsnd, 1, DCncom, 'DCsnd', 'dscfComm' )
      call re_alloc( DCrcv, 1, DCncom, 'DCrcv', 'dscfComm' )
      call re_alloc( DCsne, 1, DCncom, 'DCsne', 'dscfComm' )
      call re_alloc( DCrce, 1, DCncom, 'DCrce', 'dscfComm' )
      call re_alloc( DCreq, 1, 2, 1, DCncom, 'DCreq', 'dscfComm' )
      ncom = 0
      do ii = 0, NODES-2
        oNode = MOD((ii/2+1)*
     &             MERGE(1,-1,IAND(NODE,ishft(1,ii/2))==0)*
     &             MERGE(1,-1,IAND(ii,1)==0)+NODES+NODE,NODES)
        if (send(oNode)/=0 .or. recv(oNode)/=0) then
          ncom = ncom + 1
          DCpid(ncom) = oNode
          DCsnd(ncom) = send(oNode)
          DCrcv(ncom) = recv(oNode)
        endif
      enddo
      DCself = send(Node)
      DCtotal = DCself + SUM(DCsnd) + SUM(DCrcv)

      call de_alloc( send, 'send', 'dscfComm' )
      call de_alloc( recv, 'recv', 'dscfComm' )

      nullify(DCinvp)
      call re_alloc( DCinvp, 1, DCtotal, 'DCinvp', 'dscfComm' )
      DCinvp(:DCself) = invp(xnei(node)+1:xnei(node+1))
      ind = DCself
#ifdef MPI
      do ii= 1, DCncom
        oNode = DCpid(ii)
        if (DCrcv(ii)/=0) then
          DCinvp(ind+1:ind+DCrcv(ii)) =
     &      invp(xnei(oNode)+1:xnei(oNode+1))
          call MPI_ISEND( DCinvp(ind+1:), DCrcv(ii), MPI_INTEGER,
     &      oNode, 0, MPI_COMM_WORLD, DCreq(1,ii), MPIerr )
          ind = ind + DCrcv(ii)
        endif
        if (DCsnd(ii)/=0) then
          call MPI_Irecv( DCinvp(ind+1:),  DCsnd(ii), MPI_INTEGER,
     &      oNode, 0, MPI_COMM_WORLD, DCreq(2,ii), MPIerr )
          ind = ind + DCsnd(ii)
        endif
      enddo

      call CommWaitAll( )
#endif
      call de_alloc( invp, 'invp', 'dscfComm' )
      call de_alloc( xnei, 'xnei', 'dscfComm' )

      call timer( "dscfComm", 2 )
!--------------------------------------------------------------------------- END
      end subroutine dscfComm

      subroutine CommWaitAll( )
#ifdef MPI
      use mpi_siesta, only: MPI_Wait, MPI_Status_Size
#endif
      implicit none
      integer :: ii, MPIerr, Status(MPI_Status_Size)
      ! Wait until all communications are done
#ifdef MPI
      do ii= 1, DCncom
        if (DCrcv(ii)/=0) call MPI_Wait( DCreq(1,ii), status, MPIerr )
        if (DCsnd(ii)/=0) call MPI_Wait( DCreq(2,ii), status, MPIerr )
      enddo
#endif
      end subroutine CommWaitAll

      subroutine CommWaitSnd( )
#ifdef MPI
      use mpi_siesta, only: MPI_Wait, MPI_Status_Size
#endif
      implicit none
      integer :: ii, MPIerr, Status(MPI_Status_Size)
      ! Wait until all send communications are done
#ifdef MPI
      do ii= 1, DCncom
        if (DCsnd(ii)/=0) call MPI_Wait( DCreq(2,ii), status, MPIerr )
      enddo
#endif
      end subroutine CommWaitSnd

      subroutine CommWaitRcv( )
#ifdef MPI
      use mpi_siesta, only: MPI_Wait, MPI_Status_Size
#endif
      implicit none
      integer :: ii, MPIerr, Status(MPI_Status_Size)
      ! Wait until all send communications are done
#ifdef MPI
      do ii= 1, DCncom
        if (DCrcv(ii)/=0) call MPI_Wait( DCreq(1,ii), status, MPIerr )
      enddo
#endif
      end subroutine CommWaitRcv

      subroutine resetdscfComm( )
      use alloc,      only : de_alloc
      implicit none
!------------------------------------------------------------------------- BEGIN
      call de_alloc( DCBuffer, 'DCBuffer', 'dscfComm' )
      call de_alloc( DCinvp, 'DCinvp', 'dscfComm' )
!--------------------------------------------------------------------------- END
      end subroutine resetdscfComm

      end module m_dscfComm
