! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_ion_io
!
!     Support for reading PAO and KB information from
!     ASCII .ion files
!
!     Abridged version
      
      use atm_types
      use radial

      implicit none

      public :: read_basis_ascii
      public :: read_ion_ascii

      type(species_info), pointer        :: spp

      integer, parameter :: dp = selected_real_kind(10,100)
      private

      CONTAINS

!=======================================================================

      subroutine read_basis_ascii(nspecies,species_label)
      integer, intent(in) :: nspecies
      character(len=*), dimension(:) :: species_label

      integer is

      allocate(species(nspecies))

      do is = 1, nspecies
        spp                => species(is)
        spp%label          = species_label(is)
        spp%read_from_file = .true.

        nullify(spp%orbnl)
        nullify(spp%pjnl)
        call reset_rad_func( spp%vna )
        call reset_rad_func( spp%chlocal )
        call reset_rad_func( spp%reduced_vlocal )
        call reset_rad_func( spp%core )

        call read_ion_ascii(spp)
      enddo

      end subroutine read_basis_ascii
!
!----------------------
      subroutine read_ion_ascii(spp)
      type(species_info), pointer  :: spp

      character(len=20) filename
      integer i, l, m, lun, nor, nk, ispol

      write(filename,'(a,a)') trim(spp%label), ".ion"
      call io_assign(lun)
      open(lun,file=filename,status='old',form='formatted')
      rewind(lun)

      call read_header(spp,lun)
      read(lun,*) 
         do i=1,spp%n_orbnl
            read(lun,*)
     $           spp%orbnl_l(i), spp%orbnl_n(i), spp%orbnl_z(i),
     $           ispol, spp%orbnl_pop(i)
            spp%orbnl_ispol(i) =  ispol.eq.1
            call radial_read_ascii(spp%orbnl(i),lun,
     $                             yp1=0.0_dp,ypn=huge(1.0_dp))
         enddo
!
!        Update indexes
!
         nor = 0
         do i = 1, spp%n_orbnl
            l = spp%orbnl_l(i)
            do m = -l,l
               nor = nor+1
               spp%orb_n(nor) = spp%orbnl_n(i)
               spp%orb_l(nor) = spp%orbnl_l(i)
               spp%orb_m(nor) = m
               spp%orb_pop(nor) = spp%orbnl_pop(i) / (2*l+1)
               spp%orb_index(nor) = i
            enddo
         enddo
         spp%norbs = nor

         if ( spp%z < 0 ) then
            ! Floating species
            ! Set global number of KB projs to zero. spp%n_pjnl was set in read_header
            spp%nprojs = 0
            ! All the other radfuncs will stay reset, with n=0, so will return zero
            ! when evaluated
            goto 9999
            
         endif
! KBs
         read(lun,*)
         do i=1,spp%n_pjnl
            if (spp%lj_projs) then
               read(lun,*)
     $           spp%pjnl_l(i), spp%pjnl_j(i), spp%pjnl_n(i),
     $           spp%pjnl_ekb(i)
            else
               read(lun,*)
     $           spp%pjnl_l(i), spp%pjnl_n(i), spp%pjnl_ekb(i)
            endif
            call radial_read_ascii(spp%pjnl(i),lun,
     $                             yp1=0.0_dp,ypn=huge(1.0_dp))
         enddo
!
!        Update indexes
!

         nk = 0
         do i = 1, spp%n_pjnl
            l = spp%pjnl_l(i)
            do m = -l,l
               nk = nk+1
               spp%pj_n(nk) = spp%pjnl_n(i)
               spp%pj_l(nk) = spp%pjnl_l(i)
               spp%pj_j(nk) = spp%pjnl_j(i)
               spp%pj_m(nk) = m
               spp%pj_index(nk) = i
            enddo
         enddo
         spp%nprojs = nk
!
!Vna
         read(lun,*)
         call radial_read_ascii(spp%vna,lun,
     $                          yp1=0.0_dp,ypn=huge(1.0_dp))

!
!Chlocal
         read(lun,*)
         call radial_read_ascii(spp%chlocal,lun,
     $                          yp1=0.0_dp,ypn=huge(1.0_dp))
!
!Core
         read(lun,*,end=9999)
         call radial_read_ascii(spp%core,lun,
     $                          yp1=0.0_dp,ypn=huge(1.0_dp))

 9999    continue
         call io_close(lun)

      CONTAINS

      subroutine read_header(p,unit)

      type(species_info), pointer :: p
      integer, intent(in)         :: unit
      
      character(len=78) line
      integer :: iostat

      read(unit,'(a)') line
      if (trim(line) .eq. '<preamble>') then
 1       continue
         read(unit,'(a)') line
         if (trim(line) .ne. '</preamble>') goto 1
      endif
         
      read(unit,'(a2)') p%symbol
      read(unit,'(a20)') p%label
      read(unit,*) p%z
      read(unit,*) p%zval
      read(unit,*) p%mass
      read(unit,*) p%self_energy
      read(unit,*) p%lmax_basis, p%n_orbnl
      read(unit,fmt=*,iostat=iostat) p%lmax_projs, p%n_pjnl, p%lj_projs
      if (iostat /=0 ) then
         backspace(unit)
         read(unit,*) p%lmax_projs, p%n_pjnl
         p%lj_projs = .false.
      endif

      allocate(p%orbnl(p%n_orbnl))
      allocate(p%pjnl(p%n_pjnl))

      end subroutine read_header

      end subroutine read_ion_ascii

      end module m_ion_io












