! ---
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

      MODULE m_iodynmat
      
     use siesta_geom,       only: na_u
     use files,         only: slabel, label_length
     use parallel,      only: Node, IOnode
     use units,         only: Ang, eV
     use precision,  only : dp


      implicit none
      character(len=*), parameter :: floatfmt = '(ES22.14)'

      private

      public :: writedynmat, readdynmat


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
        write(*,'(a)') 'writedynmat: Saving into a file the force constant matrix'
        fname = trim(slabel)//'.FC'
      else
        write(*,'(a)') 'writedynmat: Saving into a file calculated dynamical matrix'
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
 
      write(6,'(a)')'readdynmat: Reading from file previous dynamical matrix'

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
        write(6,'(a,i7)') 'readdynmat: New Initial perturbed atom !!!', init
      else
        write(6,'(/,a,a,a)') 'readdynmat: File ',fname
        write(6,'(a)') 'File not found... starting from the begining '
         init=iai
      endif

      end subroutine readdynmat

end module m_iodynmat
