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

! Begin.....
     
! If ialr ==iaf, this is the final calculation of Linres and we want to store the 
! Force constant matrix to be read by vibra.
! In the other cases, we are saving the intermediate dynamat elements just in case
! explosive events....
      write(*,'(a)') 'writedynmat: Saving into a file calculated dynamical matrix'

      if (final_call ) then
        fname = trim(slabel)//'.FC'
      else
        fname = trim(slabel)//'.temp.DYNMAT'
      endif 

      if (IONode) then
          if (final_call) then !Store the force constant matrix
            call io_assign(unit1)
            open(unit1, file=fname, status='unknown' )
            rewind(unit1)
            write(unit1,'(a)') 'Force constants matrix'
            do i = iai,iaf
              do ix = 1,3
                do j = 1,na_u
                   write(unit1,'(3f15.7)') &
                         (-Ang**2/eV)*(dynmat(j,1,i,ix)), &
                         (-Ang**2/eV)*(dynmat(j,2,i,ix)), &
                         (-Ang**2/eV)*(dynmat(j,3,i,ix))

                enddo
                do j = 1,na_u
                  write(unit1,'(3f15.7)') &
                         (-Ang**2/eV)*(dynmat(j,1,i,ix)), &
                         (-Ang**2/eV)*(dynmat(j,2,i,ix)), &
                         (-Ang**2/eV)*(dynmat(j,3,i,ix))
                enddo
              enddo
            enddo
            call io_close(unit1) 
          else ! Store flag + dynamical matrix
            call io_assign(unit1)
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
            call io_close(unit1)
          endif !final_call
      endif  !Node

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

      if (IONode) then 
      write(6,'(a)')'readdynmat: Reading from file previous dynamical matrix'

      fname = trim(slabel)//'.temp.DYNMAT'
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
        write(6,'(/,a,a,a)') 'readdynmat: File ',fname,' not found... starting from the begining '
         init=iai
      endif
      endif !IONODE

      end subroutine readdynmat


end module m_iodynmat
