! ---
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

      MODULE m_iodynmat
      
     use siesta_geom,  only: na_u
     use files,        only: slabel, label_length
     use parallel,     only: Node, IOnode, Nodes
     use units,        only: Ang, eV
     use precision,    only : dp
     use sys,          only : die
     use fdf,          only : fdf_boolean, fdf_string
     use alloc,        only : re_alloc, de_alloc
     use m_mpi_utils,  only : broadcast
     use parallelsubs, only : LocalToGlobalOrb
#ifdef MPI
     use mpi_siesta
     use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb
#endif


      implicit none
      character(len=*), parameter :: floatfmt = '(ES22.14)'

      private

      public :: writedynmat, readdynmat, write_dmlr, read_dmlr


CONTAINS

      subroutine writedynmat (iai,iaf,ialr,dynmat, final_call)
! **************************************************************************
! S.Illera: Linres subroutine which stores the dynamical matrix to a file.
! For intermediate "storages", it writes an extra flag to know the last atom 
! which is well calculated.
! *****************************INPUTS***************************************
! INTEGER IAF		:Index of final atom which will be perturbed
! INTEGER IAI           :Index of first atom which is perturbed
! INTEGER IALR          :Index perturbed atom



! **************************************************************************
! ****************************OUTPUTS***************************************
! *****************************none*****************************************
! **************************************************************************
  
      implicit none

      integer, intent(in) :: iai,iaf,ialr 
      real(dp), intent(in) :: dynmat(na_u,3,na_u,3)
      logical, intent (in) :: final_call

! Internal variables ----------------------------------------------------
      character(len=label_length+10) :: fname
      integer    :: unit1, j, i, ix, jx
      real(dp)   :: conv
! Begin.....
     
! If ialr ==iaf, this is the final calculation of Linres and we want to store the 
! Force constant matrix to be read by vibra.
! In the other cases, we are saving the intermediate dynamat elements just in case
! explosive events....
      if (.not. IONode) return 

      if (final_call ) then
        write(*,'(a)') 'writedynmat: Saving into a file the force &
                      &constant matrix'
        fname = trim(slabel)//'.FC'
      else
        write(*,'(a)')
        write(*,'(a)') 'writedynmat: Saving into a file calculated &
                      &dynamical matrix'
        fname = trim(slabel)//'.LRDYNMAT'
      endif 

      call io_assign(unit1)

      if (final_call) then !Store the force constant matrix
         open(unit1, file=fname, status='unknown' )
         rewind(unit1)
         write(unit1,'(a)') 'Force constants matrix'
         conv=-Ang ** 2 / eV
         do i = iai,iaf
           do ix = 1,3
             do j = 1,na_u
                write(unit1,'(3f15.7)') &
                       conv*(dynmat(j,:,i,ix))
             enddo
             do j = 1,na_u
                write(unit1,'(3f15.7)') &
                     conv*(dynmat(j,:,i,ix))
             enddo
           enddo
         enddo
      else ! Store flag + dynamical matrix
         open(unit1, file=fname, status='unknown' )
         rewind(unit1)
         write(unit1,'(a)') 'ialr='
         write(unit1,*) ialr
         do i = iai, iaf
           do ix = 1,3
             do j = 1, na_u
               do jx = 1,3
                 write(unit1,floatfmt) dynmat(j,jx,i,ix)
               enddo
             enddo
           enddo
         enddo
      endif !final_call

      call io_close(unit1)

      end subroutine writedynmat


      subroutine readdynmat(init,iai,dynmat)
! **************************************************************************
! S.Illera: Linres subroutine which reads the dynamical matrix from file.
! It also reads the first flag of the last atom to update the first atom to move
! *****************************INPUTS***************************************

! **************************************************************************
! ****************************OUTPUTS***************************************
! *****************************none*****************************************
! **************************************************************************
      implicit none


      integer, intent(inout) :: init
      real(dp), intent(inout) :: dynmat(na_u,3,na_u,3)
      integer, intent (in) :: iai
! Internal variables

      character(len=label_length+10) :: fname
      logical               :: found
      integer               :: unit1,jx,ix,j,i

      if (.not. IONode) return

      write(6,'(a)')
      write(6,'(a)')'readdynmat: Reading from file the &
     &dynamical matrix'

      fname = trim(slabel)//'.LRDYNMAT'
      inquire( file=fname, exist=found )

      if (found) then
        call io_assign(unit1)
        open( unit1, file=fname, status='old' )
        read( unit1, *)
        read( unit1, *)  init !read the first flag
        do i=1,na_u
          do ix=1,3
            do j=1,na_u
              do jx=1,3
                read(unit1,floatfmt) dynmat(j,jx,i,ix) 
              enddo
            enddo
          enddo
        enddo
        call io_close( unit1 )

        init=init+1
        write(6,'(a,i7)') 'readdynmat: New Initial perturbed atom !!!',&
                            init
      else
        write(6,'(a,a)') 'readdynmat: File ',fname
        write(6,'(a)') 'readdynmat: File not found...&
     & starting from the begining '
         init=iai
      endif

      end subroutine readdynmat


      subroutine write_dmlr (maxnd, no_l, nspin, &
                        numd, listdptr, listd, dm, ialr)
! Subroutine to write the Linres dDscf to disk

      integer, intent(in) :: maxnd
      integer, intent(in) :: no_l
      integer, intent(in) :: nspin
      integer, intent(in) :: numd(1:no_l)
      integer, intent(in) :: listdptr(1:no_l)
      integer, intent(in) :: listd(maxnd)
      real(dp), intent(in) :: dm(maxnd,nspin,3)
      integer, intent(in) :: ialr

      integer :: no_u, m, ml, im, ndmaxg, unit1, is,ix
#ifdef MPI
      integer   MPIerror, Request, Status(MPI_Status_Size)
      integer   BNode
      real(dp), dimension(:), pointer :: buffer
      integer,  dimension(:), pointer :: ibuffer
#endif
      integer, dimension(:), pointer  :: numdg
      character(len=label_length+5+5) :: fname
      character(len=5) :: atomdisp 
      logical :: fmto
      character(len=*), parameter :: intfmt = '(I11)'

      fmto=.false.

      call timer("WriteLRDM",1)

!      call setup_file_modes()

!     Find total number of basis functions over all Nodes
#ifdef MPI
      call MPI_AllReduce(no_l,no_u,1,MPI_integer,MPI_sum, &
          MPI_Comm_World,MPIerror)
#else
      no_u = no_l
#endif
      write(atomdisp,'(i5)') ialr
      fname = slabel

      if (Node.eq.0) then
         call io_assign(unit1)
         open( unit1, file=trim(fname)//'.LRDM'//adjustl(atomdisp), &
             form='unformatted', status='unknown' )
         rewind(unit1)
         if (fmto) then
            write(unit1, intfmt) no_u, nspin
         else
            write(unit1) no_u, nspin
         endif
      endif

      nullify(numdg)
      call re_alloc( numdg, 1, no_u, 'numdg', 'write_dmlr' )

!     Create globalised numd
      do m = 1,no_u
#ifdef MPI
         call WhichNodeOrb(m,Nodes,BNode)
         if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
            ml = m
#endif
            numdg(m) = numd(ml)
#ifdef MPI
         elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
            call MPI_ISend(numd(ml),1,MPI_integer, &
                0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
         elseif (Node.eq.0) then
            call MPI_IRecv(numdg(m),1,MPI_integer, &
                BNode,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
         endif
         if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
         endif
#endif
      enddo
!     Write out numd array
      if (Node.eq.0) then
         ndmaxg = maxval(numdg(1:no_u))
         if (fmto) then
            write(unit1, intfmt) (numdg(m),m=1,no_u)
         else
            write(unit1) (numdg(m),m=1,no_u)
         endif
#ifdef MPI
         nullify(buffer,ibuffer)
         call re_alloc( buffer,  1, ndmaxg, 'buffer',  'write_dmlr' )
         call re_alloc( ibuffer, 1, ndmaxg, 'ibuffer', 'write_dmlr' )
#endif
      endif
!     Write out listd array
      do m = 1,no_u
#ifdef MPI
         call WhichNodeOrb(m,Nodes,BNode)
         if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
            ml = m
#endif
            if (fmto) then
               write(unit1, intfmt) &
                   (listd(listdptr(ml)+im),im=1,numd(ml))
            else
               write(unit1) (listd(listdptr(ml)+im),im=1,numd(ml))
            endif
#ifdef MPI
         elseif (Node.eq.0) then
            call MPI_IRecv(ibuffer,numdg(m),MPI_integer,BNode,1, &
                MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
         elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
            call MPI_ISend(listd(listdptr(ml)+1),numd(ml),MPI_integer, &
                0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
         endif
         if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
            if (Node.eq.0) then
               if (fmto) then
                  write(unit1, intfmt) (ibuffer(im),im=1,numdg(m))
               else
                  write(unit1) (ibuffer(im),im=1,numdg(m))
               endif
            endif
         endif
#endif
      enddo

#ifdef MPI
      if (Node.eq.0) then
         call de_alloc( ibuffer, 'ibuffer', 'write_dmlr' )
      endif
#endif

!     Write density matrix
      do is=1,nspin
         do m=1,no_u
#ifdef MPI
            call WhichNodeOrb(m,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
               call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
               ml = m
#endif
               if (fmto) then
                  write(unit1, floatfmt) &
                      ((dm(listdptr(ml)+im,is,ix), &
                      im=1,numd(ml)),ix = 1,3)
               else
                  write(unit1) ((dm(listdptr(ml)+im,is,ix), &
                 im=1,numd(ml)),ix = 1,3)
               endif
#ifdef MPI
            elseif (Node.eq.0) then
               call MPI_IRecv(buffer,numdg(m),MPI_double_precision, &
                   BNode,1,MPI_Comm_World,Request,MPIerror)
               call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
               call GlobalToLocalOrb(m,Node,Nodes,ml)
               call MPI_ISend(dm(listdptr(ml)+1,is,ix),numd(ml), &
                   MPI_double_precision,0,1,MPI_Comm_World,Request, &
                   MPIerror)
               call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
               call MPI_Barrier(MPI_Comm_World,MPIerror)
               if (Node.eq.0) then
                  if (fmto) then
                     write(unit1, floatfmt) (buffer(im),im=1,numdg(m))
                  else
                     write(unit1) (buffer(im),im=1,numdg(m))
                  endif
               endif
            endif
#endif
         enddo
      enddo

      if (Node.eq.0) then
#ifdef MPI
         call de_alloc( buffer, 'buffer', 'write_dmlr' )
#endif
         call io_close(unit1)
      endif

      call de_alloc( numdg, 'numdg', 'write_dmlr' )

      call timer("WriteLRDM",2)

      end subroutine write_dmlr


      subroutine read_dmlr ( maxnd, no_l, nspin, numd, &
                          listdptr, listd, dm, fname )

!     Reads density matrix from file
!     Written by P.Ordejon and J.M.Soler. May 1997.
!     Made into a separate routine by Alberto Garcia, August 2007
!     ********* INPUT
!     ***************************************************
!     integer   maxnd    : First dimension of listd and dm
!     integer   no_l   : Number of atomic orbitals
!     integer   nspin    : Number of spins (1 or 2)
!     ********* OUTPUT ************
!     integer numd(no_l)     : Control vector of DM matrix
!     (number of nonzero elements of each row)
!     integer listdptr(no_l) : Control vector of DM matrix
!     (pointer to the start of each row)
!     integer listd(maxnd)     : Control vector of DM matrix
!     (list of nonzero elements of each row)
!     real*8  dm(maxnd,nspin)  : Density matrix
!     logical found : Has DM been found in disk? 
!     *******************************************

      ! These are the "input" sparsity parameters
      ! We might want to read the file info into a temporary
      ! set of arrays and change the structure if needed
      integer, intent(in) :: maxnd
      integer, intent(in) :: no_l
      integer, intent(in) :: nspin
      character(len=label_length+5+5),intent(in) :: fname

      integer, intent(out) :: numd(no_l)
      integer, intent(out) :: listdptr(no_l)
      integer, intent(out) :: listd(maxnd)
      real(dp), intent(out) :: dm(maxnd, nspin,3)


!     Internal variables
      logical   exist3
      integer   im, is, unit1, unit2, m, nb, ndmax, ns
      integer   no_u, ml, ndmaxg, ix
      integer, dimension(:), pointer  :: numdg
      character(len=*), parameter :: intfmt = '(I11)'
      logical :: fmti

#ifdef MPI
      integer   MPIerror, Request, Status(MPI_Status_Size)
      integer   BNode

      real(dp), dimension(:), pointer :: buffer
      integer,  dimension(:), pointer :: ibuffer
#endif

      external          chkdim

      fmti=.false.

      if (Node.eq.0) then
         inquire (file=fname, exist=exist3)
      endif
      call broadcast(exist3)

      if ( .not. exist3) RETURN

!     Find total number of basis functions over all Nodes
#ifdef MPI
      call MPI_AllReduce(no_l,no_u,1,MPI_integer,MPI_sum, &
          MPI_Comm_World,MPIerror)
#else
      no_u = no_l
#endif

      if (Node.eq.0) then
         write(6,'(/,a)') 'readdmlr: Reading Density Matrix from files'
         call io_assign(unit1)
         open( unit1, file=fname, form='unformatted',status='old' )
         rewind(unit1)
         if (fmti) then
            read(unit1, intfmt) nb, ns
         else
            read(unit1) nb, ns
         endif
      endif
!     Communicate the values to all Nodes and adjust to allow for
!     distributed memory before checking the dimensions
#ifdef MPI
      call MPI_Bcast(nb,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(ns,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

      ! This checks for equality of spin and number of rows.
      call chkdim( 'iodm', 'no_u', no_u, nb, 0 )
      call chkdim( 'iodm', 'nspin',  nspin,  ns, 0 )

!     Allocate local buffer array for globalised numd
      nullify(numdg)
      call re_alloc( numdg, 1, no_u, 'numdg', 'read_dm' )
      if (Node.eq.0) then
         if (fmti) then
            read(unit1, intfmt) (numdg(m),m=1,no_u)
         else
            read(unit1) (numdg(m),m=1,no_u)
         endif
      endif
      call broadcast(numdg(1:no_u))

      ndmax = 0
      do m = 1,no_l
         call LocalToGlobalOrb(m,Node,Nodes,ml)
         numd(m) = numdg(ml)
!         numd_tmp(m) = numdg(ml)
         ndmax = ndmax + numdg(ml)
         if (m .eq. 1) then
!            listdptr_tmp(1) = 0
            listdptr(1) = 0
         else
!            listdptr_tmp(m) = listdptr_tmp(m-1) + numd_tmp(m-1)
            listdptr(m) = listdptr(m-1) + numd(m-1)
         endif
      enddo
      ndmaxg = maxval(numdg(1:no_u))

!     Check size of first dimension of dm
      call chkdim( 'iodm', 'maxnd', maxnd, ndmax, 1 )

#ifdef MPI
!     Create buffer arrays for transfering density matrix between nodes
!     and lists
      nullify(buffer,ibuffer)
      call re_alloc( buffer,  1, ndmaxg, 'buffer',  'read_dm' )
      call re_alloc( ibuffer, 1, ndmaxg, 'ibuffer', 'read_dm' )
#endif

      do m = 1,no_u
#ifdef MPI
         call WhichNodeOrb(m,Nodes,BNode)
         if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
            ml = m
#endif
            if (fmti) then
               read(unit1, intfmt) &
                   (listd(listdptr(ml)+im),im=1,numd(ml))
            else
               read(unit1) (listd(listdptr(ml)+im),im=1,numd(ml))
            endif
#ifdef MPI
         elseif (Node.eq.0) then
            if (fmti) then
               read(unit1, intfmt) (ibuffer(im),im=1,numdg(m))
            else
               read(unit1) (ibuffer(im),im=1,numdg(m))
            endif
            call MPI_ISend(ibuffer,numdg(m),MPI_integer, &
                 BNode,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
         elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
            call MPI_IRecv(listd(listdptr(ml)+1),numd(ml), &
                 MPI_integer,0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
         endif
         if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
         endif
#endif
      enddo

#ifdef MPI
      call de_alloc( ibuffer, 'ibuffer', 'read_dm' )
#endif

      do is = 1,nspin
         do m = 1,no_u

#ifdef MPI
            call WhichNodeOrb(m,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
               call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
               ml = m
#endif
               if (fmti) then
                  read(unit1, floatfmt) &
                       ((dm(listdptr(ml)+im,is,ix), &
                       im=1,numd(ml)),ix = 1,3)
               else
                  read(unit1) ((dm(listdptr(ml)+im,is,ix), &
                              im=1,numd(ml)),ix=1,3)
               endif
#ifdef MPI
            elseif (Node.eq.0) then
               if (fmti) then
                  read(unit1, floatfmt) (buffer(im),im=1,numdg(m))
               else
                  read(unit1) (buffer(im),im=1,numdg(m))
               endif
               call MPI_ISend(buffer,numdg(m),MPI_double_precision,&
                    BNode,1,MPI_Comm_World,Request,MPIerror)
               call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
               call GlobalToLocalOrb(m,Node,Nodes,ml)
               call MPI_IRecv(dm(listdptr(ml)+1,is,ix),numd(ml), &
                    MPI_double_precision,0,1,MPI_Comm_World,Request, &
                    MPIerror)
               call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
               call MPI_Barrier(MPI_Comm_World,MPIerror)
            endif
#endif

         enddo
      enddo
#ifdef MPI
      call de_alloc( buffer, 'buffer', 'read_dm' )
#endif
      call de_alloc( numdg, 'numdg', 'read_dm' )

      if (Node.eq.0) then
         call io_close(unit1)
      endif

      end subroutine read_dmlr



      end module m_iodynmat
