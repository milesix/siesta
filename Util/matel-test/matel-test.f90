! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
program matel_test

! The PAO and KB information can optionally be read from ASCII files
! (those produced by a standard run of Siesta or Base, but with the extension 
! renamed to '.input' instead of '.dump'), or from NetCDF files (if NetCDF
! is available). Note that there must be files for *all* species.
!
      use precision

      use m_ion_io
      use atm_types, only: species, nspecies
      use m_matel_registry, only: show_pool
      use matel_mod, only: init_matel_main_tables
      use matel_mod, only: init_matel_thermal_transport
      use matel_mod, only: init_matel_orb_xyz_orb
      use matel_mod, only: new_matel
      use m_old_matel, only: old_matel
      
      implicit none

      
      character(len=20), allocatable, dimension(:) :: species_label
      integer :: io, jpx, jpy, jpz, i
      integer :: ko, kpx, kpy, kpz
      integer :: i_r, i_grad
      real(dp) :: val, grad(3)
      real(dp) :: x, y, z

      logical :: xtables = .true.
      
      character(len=1), dimension(3), parameter :: coord_table = [ 'X', 'Y', 'Z' ]

      ! ----------------------
      
      nspecies = 1
      allocate(species_label(1))
      species_label(1) = "Si"
      
      write(6,'(/a)') 'Reading PAOs and KBs from ascii files...'
      call read_basis_ascii(nspecies,species_label)

      write(*,*) species(1)%label
      
      call register_rfs()
      call show_pool()
      call init_matel_main_tables()
      call init_matel_thermal_transport()
      call init_matel_thermal_transport()
      call init_matel_orb_XYZ_orb()
      
      io = species(1)%orb_gindex(1)
      jpy = species(1)%orb_gindex(2)  
      jpz = species(1)%orb_gindex(3)  
      jpx = species(1)%orb_gindex(4)  
      ko = species(1)%pj_gindex(1)  
      kpy = species(1)%pj_gindex(2)  
      kpz = species(1)%pj_gindex(3)  
      kpx = species(1)%pj_gindex(4)  

      call new_MATEL('S', io, io, [ 1.0_dp, 0.0_dp, 0.0_dp ], val, grad)
      print "(a,f10.4,2x,3f10.4)", "<3s|  | 3s(1.0)> : ", val, grad

      do i_grad = 1, 3
         call new_MATEL(SG(i_grad), io, ko, [ 1.0_dp, 0.0_dp, 0.0_dp ], val, grad)
         print "(a,f10.4,2x,3f10.4)", "<s|" // SG(i_grad) // "|ks(1.0)> : ", val, grad
      
         call new_MATEL(SG(i_grad), io, kpx, [ 1.0_dp, 0.0_dp, 0.0_dp ], val, grad)
         print "(a,f10.4,2x,3f10.4)", "<s|" // SG(i_grad) // "|kpx(1.0)> : ", val, grad
      enddo
   
      do i_r = 1, 3
         do i_grad = 1, 3
         
            call new_MATEL(RG(i_r,i_grad), io, ko, [ 0.0_dp, 0.0_dp, 0.0_dp ], val, grad)
            print "(a,f10.4,2x,3f10.4)", "<s|" // RG(i_r,i_grad) // "|ks(0.0)> : ", val, grad
            call new_MATEL(RG(i_r,i_grad), io, ko, [ 1.0_dp, 0.0_dp, 0.0_dp ], val, grad)
            print "(a,f10.4,2x,3f10.4)", "<s|" // RG(i_r,i_grad) // "|ks(1.0)> : ", val, grad

         enddo
      enddo

      print *, "Done thermal transport"
      
      call new_MATEL('S', io, jpx, [ 1.0_dp, 0.0_dp, 0.0_dp ], val, grad)
      print "(a,f10.4,2x,3f10.4)", "<s|  | px(1.0)> : ", val, grad
      call old_MATEL('S', io, jpx, [ 1.0_dp, 0.0_dp, 0.0_dp ], val, grad)
      print "(a,f10.4,2x,3f10.4)", "OLD<s|  | px(1.0)> : ", val, grad

      call new_MATEL('S', io, jpy, [ 0.0_dp, 1.0_dp, 0.0_dp ], val, grad)
      print "(a,f10.4,2x,3f10.4)", "<s|  | py(1.0)> : ", val, grad
      call old_MATEL('S', io, jpy, [ 0.0_dp, 1.0_dp, 0.0_dp ], val, grad)
      print "(a,f10.4,2x,3f10.4)", "OLD<s|  | py(1.0)> : ", val, grad
      
      call new_MATEL('S', io, jpz, [ 0.0_dp, 0.0_dp, 1.0_dp ], val, grad)
      print "(a,f10.4,2x,3f10.4)", "<s|  | pz(1.0)> : ", val, grad
      call old_MATEL('S', io, jpz, [ 0.0_dp, 0.0_dp, 1.0_dp ], val, grad)
      print "(a,f10.4,2x,3f10.4)", "OLD<s|  | pz(1.0)> : ", val, grad

      if ( xtables ) then

      call new_MATEL('X', io, io, [ 1.0_dp, 0.0_dp, 0.0_dp ], val, grad)
      print "(a,f10.4,2x,3f10.4)", "<3s| X | 3s(1.0)> : ", val, grad

      call new_MATEL('Y', io, io, [ 1.0_dp, 0.0_dp, 0.0_dp ], val, grad)
      print "(a,f10.4,2x,3f10.4)", "<3s| Y | 3s(1.0)> : ", val, grad

      call new_MATEL('Z', io, io, [ 1.0_dp, 0.0_dp, 0.0_dp ], val, grad)
      print "(a,f10.4,2x,3f10.4)", "<3s| Z | 3s(1.0)> : ", val, grad


      call new_MATEL('X', io, jpx, [ 1.0_dp, 0.0_dp, 0.0_dp ], val, grad)
      print "(a,f10.4,2x,3f10.4)", "<s| X | px(1.0)> : ", val, grad
      call new_MATEL('Y', io, jpx, [ 1.0_dp, 0.0_dp, 0.0_dp ], val, grad)
      print "(a,f10.4,2x,3f10.4)", "<s| Y | px(1.0)> : ", val, grad
      call new_MATEL('Z', io, jpx, [ 1.0_dp, 0.0_dp, 0.0_dp ], val, grad)
      print "(a,f10.4,2x,3f10.4)", "<s| Z | px(1.0)> : ", val, grad

      do i = 1, 100
         x = (i-1) * 0.1
         call new_MATEL('X', io, io, [ x, 0.0_dp, 0.0_dp ], val, grad)
         print "(a,2f10.4,2x,3f10.4)", "<s|X|s(x)>: ", x, val, grad
      enddo
   endif

    contains
      function SG(i) result (s)
        integer, intent(in) :: i
        character (len=2)   :: s
        s = 'G' // coord_table(i)
      end function SG

      function RG(i,j) result (s)
        integer, intent(in) :: i, j
        character (len=3)   :: s
        s = 'H' // coord_table(i) // coord_table(j)
      end function RG

    end program  matel_test




