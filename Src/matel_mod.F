!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
!
      !> Re-designed module to allow a pre-computation of the needed
      !> matrix elements in parallel, followed by a globalization of the
      !> data among all the MPI processes.  Once the interpolation
      !> tables are setup, further calls to the matrix-element evaluator
      !> are cheap. This has a dramatic effect in some routines (such as
      !> nlefsm) that had to perform the table-building operations under
      !> conditions that did not scale in parallel.

      !> Concept: Rogeli Grima (BSC) and Alberto Garcia (ICMAB)
      !> Initial implementation:  Rogeli Grima (BSC)

      module matel_mod

      use matel_ylm, only: spher_harm_t
      use matel_table, only: matel
      use matel_table, only: MODE_S, MODE_T, MODE_XYZ

      use precision, only : dp
      use alloc,     only : alloc_default, allocDefaults

      use m_radfft,  only : reset_radfft
      use m_matel_registry, ONLY : EVALUATE,
     &                             EVALUATE_X, EVALUATE_Y, EVALUATE_Z

      PRIVATE

      !> Initialize main tables: S, T, TA
      public :: init_matel_main_tables
      !> Initialize tables X, Y, Z, for use in "molecule" optical calculations and
      !> Berry-phase-based polarization calculations
      public :: init_matel_orb_XYZ_orb
      !> Initialize tables S_opt, X_opt, Y_opt, Z_opt for use in solid-state optical calculations
      public :: init_matel_optical_P
      !> Initialize table S_wann for overlaps of Wannier90 trial functions
      public :: init_matel_wannier
#ifdef SIESTA__THERMAL_FLUX
      !> Initialize tables for thermal transport
      public :: init_matel_thermal_transport
#endif
      ! Main routine that evaluates matrix elements, keeping the old interface
      public :: new_matel

      ! These variables allow to call the initialization routines repeatedly 

      logical :: main_tables_initialized = .false.
      logical :: optical_P_initialized = .false.
      logical :: orb_XYZ_orb_initialized = .false.
      logical :: wannier_initialized = .false.
#ifdef SIESTA__THERMAL_FLUX
      logical :: thermal_transport_initialized = .false.
#endif
      !> These relate to the functions recorded in the 'matel registry';
      !> the number of each kind serve as markers for the different
      !> sections needed in the tables.  'init_matel_main_tables' fills
      !> the essential first four values, so that routine must be
      !> called first.  This is enforced by 'check_main_tables'

      !> A prerequisite for the correct setup is that the
      !> atomic-information tables (in the 'species' data structure) has
      !> been correctly filled up In future, this module might be more
      !> integrated with that data structure.

      integer           :: num_orb  ! Number of different orbitals
      integer           :: num_kb  ! Number of different KB projectors
      integer           :: num_dftu ! Number of different LDA+U functions
      integer           :: num_va  ! Number of different Vna functions
      integer           :: num_wannier_projs  ! Number of different Wannier trial orbitals

      !> k-space spherical-harmonic decomposition of registered
      !> functions and composite functions (such as x*Orb, y*KB, etc)

      type(spher_harm_t), target :: ylmk_val_all   ! orbs, kb, dftuprojs, vna
      type(spher_harm_t), target :: ylmk_x_orbs
      type(spher_harm_t), target :: ylmk_y_orbs
      type(spher_harm_t), target :: ylmk_z_orbs
      type(spher_harm_t), target :: ylmk_x_kbs
      type(spher_harm_t), target :: ylmk_y_kbs
      type(spher_harm_t), target :: ylmk_z_kbs
      type(spher_harm_t), target :: ylmk_val_wannier_projs

#ifdef SIESTA__THERMAL_FLUX
      !> For thermal transport
      type(spher_harm_t), target :: ylmk_val_grad_kbs(3)
      type(spher_harm_t), target :: ylmk_xyz_grad_kbs(3,3)
#endif

      !> There are different tables, each appropriate to a given operation and kind of function.
      !> The indexing is a bit cumbersome due to the one-dimensional nature of the matel_registry,
      !> in which all functions are stored in the same section: first orbitals, then KB projectors,
      !> then LDA+U projectors, and finally Vna.
      !> In addition, when using the Wannier interface, "trial orbitals" (called also "projectors"
      !> in the code (numproj of them) need to be dealt with.

      !> The tables are created once and for all for the duration of the program. There is currently
      !> no support for "coexisting qualities" or "on-the-fly" changes (such as could conceivably be
      !> done one day via the Lua interface)

      !> Overlaps of PAOs, KBs, or DFTUprojs with PAOs
      TYPE(MATEL) :: tab_S           ! Unity (overlap).
      !> Laplacian among PAOs
      TYPE(MATEL) :: tab_T           ! -Laplacian.
      !> Laplacian among Vna functions
      TYPE(MATEL) :: tab_TA         ! -Laplacian.

      !> X, Y, Z among orbitals (for dielectric polarization calculations
      !> and "r" optical calculations
      TYPE(MATEL) :: tab_X           ! X projection.
      TYPE(MATEL) :: tab_Y           ! Y projection.
      TYPE(MATEL) :: tab_Z           ! Z projection.

      !> Overlaps of PAOs with KBs for "momentum" optical calculations (different from table "S" because
      !> the KBs are "function 2")
      TYPE(MATEL) :: tab_S_opt        ! Unity (overlap)

      !> X,Y,Z  PAOs with KBs for "momentum" optical calculations
      TYPE(MATEL) :: tab_X_opt        ! X projection
      TYPE(MATEL) :: tab_Y_opt        ! Y projection
      TYPE(MATEL) :: tab_Z_opt        ! Z projection

      !> Overlaps of PAOs with wannier projectors
      TYPE(MATEL) :: tab_S_wann      ! Unity (overlap).

#ifdef SIESTA__THERMAL_FLUX
      !> For thermal transport
      TYPE(MATEL) :: tab_S_Grad_KB(3)           ! Unity (overlap), with gradients of KB
      TYPE(MATEL) :: tab_XYZ_Grad_KB(3,3)       ! XYZ, with gradients of KB
#endif
      ! convenience variable
      TYPE(allocDefaults) :: OLDEFS

      CONTAINS


      !>  Initialize the main MATEL tables S, T, TA

      subroutine init_matel_main_tables( )

      use atm_types, only : nspecies, species

      implicit none
      ! Local Variables
      integer             :: is, top, top_vna

      if (main_tables_initialized) RETURN

      ! Set allocation defaults
      call alloc_default( old=oldefs, copy=.true., shrink=.false. )

      CALL RESET_RADFFT( )

      ! Get the number of orbitals, KB projs, DFTU projs and Neutral-atom potentials

      NUM_ORB = 0
      NUM_KB = 0
      NUM_DFTU = 0
      NUM_VA = 0
      do is= 1, nspecies
        NUM_ORB  = NUM_ORB  + species(is)%norbs
        NUM_KB   = NUM_KB   + species(is)%nprojs
        NUM_DFTU = NUM_DFTU + species(is)%nprojsdftu
        if (species(is)%z > 0) then ! Not a floating orbital
           NUM_VA = NUM_VA + 1
        endif
      enddo

      ! Find the spherical harmonic decomposition and
      ! Fourier transform of all the registry objects
      ! except the wannier trial functions

      top = num_orb + num_kb + num_dftu
      top_vna = top + num_va

      call ylmk_val_all%compute_spha(1, top_vna, EVALUATE, 0)

      ! Overlap table from PAOs, KBs, and DFTU projs to PAOs.
      call tab_S%init(MODE_S, ylmk_val_all, 1, top,
     $                    ylmk_val_all, 1, num_orb)

      ! Laplacian table from PAOs to PAOs.
      call tab_T%init(MODE_T, ylmk_val_all, 1, num_orb,
     $                    ylmk_val_all, 1, num_orb)

      ! Laplacian table from Vnas to Vnas
      call tab_TA%init(MODE_T, ylmk_val_all, top+1, top_vna,
     $                     ylmk_val_all, top+1, top_vna)

      main_tables_initialized = .true.

      ! Restore allocation defaults
      CALL alloc_default( restore=oldefs )

      end subroutine init_matel_main_tables
!-------------------------------------------------------------
      subroutine init_matel_orb_XYZ_orb()

      if (orb_XYZ_orb_initialized) RETURN

      call check_for_main_tables()

      ! Set allocation defaults
      call alloc_default( old=oldefs, copy=.true., shrink=.false. )

      CALL RESET_RADFFT( )

      ! Find the spherical harmonic decomposition and
      ! Fourier transform of {x,y,z}*\phi, where \phi is a PAO

      ! The X, Y, and Z tables are then overlaps of PAOs and xPAOs


         call ylmk_x_orbs%compute_spha(1, num_orb, EVALUATE_X, 1)
         call tab_X%init( MODE_XYZ, ylmk_val_all, 1, num_orb,
     $                         ylmk_x_orbs, 1, num_orb, name="tabX")


         call ylmk_y_orbs%compute_spha(1, num_orb, EVALUATE_Y, 1)
         call tab_Y%init( MODE_XYZ, ylmk_val_all, 1, num_orb,
     $                         ylmk_y_orbs, 1, num_orb, name="tabY")


         call ylmk_z_orbs%compute_spha(1, num_orb, EVALUATE_Z, 1)
         call tab_Z%init( MODE_XYZ, ylmk_val_all, 1, num_orb,
     $                         ylmk_z_orbs, 1, num_orb, name="tabZ")


      orb_XYZ_orb_initialized = .true.

      ! Restore allocation defaults
      CALL alloc_default( restore=oldefs )

      end subroutine init_matel_orb_XYZ_orb

      subroutine check_for_orb_XYZ_orb()
      if (.not. orb_XYZ_orb_initialized) then
         call init_matel_orb_XYZ_orb()
      endif
      end subroutine check_for_orb_XYZ_orb
!-------------------------------------------------------------

      subroutine init_matel_optical_P()

      if (optical_P_initialized) RETURN

      call check_for_main_tables()

      ! Set allocation defaults
      call alloc_default( old=oldefs, copy=.true., shrink=.false. )

      CALL RESET_RADFFT( )

         ! Overlap table from PAOs to KBs.
         call tab_S_opt%init(MODE_S, ylmk_val_all, 1, num_orb,
     $                    ylmk_val_all, num_orb+1, num_orb + num_kb)

         ! Find the spherical harmonic decomposition and
         ! Fourier transform of {x,y,z}*\chi, where \chi is a KB proj.

         ! The X_opt, Y_opt, and Z_opt tables are then overlaps of PAOs and xKBs


         call ylmk_x_kbs%compute_spha(num_orb+1, num_orb+num_kb,
     $                                EVALUATE_X, 1)

         call tab_X_opt%init( MODE_XYZ, ylmk_val_all, 1, num_orb,
     $                         ylmk_x_kbs, 1, num_kb)


         call ylmk_y_kbs%compute_spha(num_orb+1, num_orb+num_kb,
     $                                EVALUATE_Y, 1)

         call tab_Y_opt%init( MODE_XYZ, ylmk_val_all, 1, num_orb,
     $                         ylmk_y_kbs, 1, num_kb)


         call ylmk_z_kbs%compute_spha(num_orb+1, num_orb+num_kb,
     $                                EVALUATE_Z, 1)

         call tab_Z_opt%init( MODE_XYZ, ylmk_val_all, 1, num_orb,
     $                         ylmk_z_kbs, 1, num_kb)


      optical_P_initialized = .true.

      ! Restore allocation defaults
      CALL alloc_default( restore=oldefs )

      end subroutine init_matel_optical_P

!--------------------------------------------------------------------
      subroutine init_matel_wannier( numproj )
!
!     Initialize the table of overlaps between orbitals and
!     Wannier trial functions, once the number of functions is known

      implicit none
      integer, intent(in) :: numproj

      ! Local Variables
      integer :: top_vna

      if (wannier_initialized) RETURN

      call check_for_main_tables()

      ! Set allocation defaults
      call alloc_default( old=oldefs, copy=.true., shrink=.false. )
      CALL RESET_RADFFT( )

      ! This table refers to overlaps between orbitals and
      ! Wannier trial functions, which are registered in the
      ! pool after every other function

      num_wannier_projs = numproj ! set module variable
      top_vna = num_orb + num_kb + num_dftu + num_va

      call ylmk_val_wannier_projs%compute_spha(top_vna+1,
     $                           top_vna+numproj, EVALUATE, 0)

      call tab_S_wann%init(MODE_S, ylmk_val_all, 1, num_orb,
     $                     ylmk_val_wannier_projs, 1, numproj,
     $                     name = "wannier table")

      wannier_initialized = .true.

      ! Restore allocation defaults
      CALL alloc_default( restore=oldefs )

      end subroutine init_matel_wannier

      subroutine check_for_main_tables()
      if (.not. main_tables_initialized) then
         call init_matel_main_tables()
      endif
      end subroutine check_for_main_tables

#ifdef SIESTA__THERMAL_FLUX
!------------------------------------------------------
      subroutine init_matel_thermal_transport()

      use m_matel_registry, only: idx_r, idx_grad
      use m_matel_registry, only: EVALUATE_GRAD
      use m_matel_registry, only: EVALUATE_R_GRAD

      integer :: i_r, i_grad
      character(len=40) :: name

      if (thermal_transport_initialized) RETURN
      call check_for_main_tables()
      call check_for_orb_XYZ_orb()
      call init_matel_optical_P()

      ! Set allocation defaults
      call alloc_default( old=oldefs, copy=.true., shrink=.false. )

      CALL RESET_RADFFT( )

      do i_grad = 1, 3
         ! Set module variable for evaluation
         idx_grad = i_grad

         call ylmk_val_grad_kbs(i_grad)%compute_spha(num_orb+1,
     $                                          num_orb+num_kb,
     $                                          EVALUATE_GRAD, 1)
      ! Overlap table from PAOs to grad_KBs.
         write(name,"(a,i1)") "PAO-GradKB-", i_grad
         call tab_S_Grad_KB(i_grad)%init(MODE_S, ylmk_val_all,
     $        1, num_orb,
     $        ylmk_val_grad_kbs(i_grad), 1, num_kb, name)
      end do

      do i_r = 1, 3
         do i_grad = 1, 3

            ! Set module variables for evaluation
            idx_r = i_r
            idx_grad = i_grad

            call ylmk_xyz_grad_kbs(i_r,i_grad)%compute_spha(num_orb+1,
     $           num_orb+num_kb,
     $           EVALUATE_R_GRAD, 2)
            ! Overlap table from PAOs to grad_KBs.
            write(name,"(a,2i1)") "PAO-R_GradKB-", i_r, i_grad
            call tab_XYZ_Grad_KB(i_r,i_grad)%init(MODE_XYZ,
     $           ylmk_val_all, 1, num_orb,
     $           ylmk_xyz_grad_kbs(i_r,i_grad), 1, num_kb, name)
         end do
      end do

      thermal_transport_initialized = .true.

      ! Restore allocation defaults
      CALL alloc_default( restore=oldefs )

      end subroutine init_matel_thermal_transport

#endif

!===================================================================
      SUBROUTINE new_MATEL( OPERAT, IG1, IG2, R12, S12, DSDR )

! Finds two-center matrix elements between 'atomic orbitals'
! with finite radial and angular momentum cutoffs.

! The interface is the same as "new_matel" in the m_new_matel
! module, but the routine itself only dispatches execution
! to the new table objects (Alberto Garcia, July 2019, after Rogeli Grima)      

! The tables must be initialized previously by calls to the init_matel_XXXX
! routines.

! INPUT

C CHARACTER OPERAT : Operator to be used. The valid options are:
C   'S' => Unity (overlap). Uppercase required for all values.
C   'T' => -Laplacian
C   'U' => 1/|r'-r| (with evaluate returning charge distributions) *NOT IMPLEMENTED
C   'X' => x, returning <phi1(r-R12)|x|phi2(r)> (**origin on second atom**)
C   'Y' => y, returning <phi1(r-R12)|y|phi2(r)>
C   'Z' => z, returning <phi1(r-R12)|z|phi2(r)>
!
#ifdef SIESTA__THERMAL_FLUX
!   For thermal transport
!   'G{X,Y,Z}': <phi1 | (Grad phi2) >
!   'H{X,Y,Z}{X,Y,Z}': <phi1 | {XYZ}(Grad phi2) > (origin on second atom)
#endif
C INTEGER IG1    : Global index of 1st function (must be positive)
C INTEGER IG2    : Global index of 2nd function
C                    Indexes IG1, IG2 are used only to call
C                    routines LCUT, RCUT and EVALUATE (in matel_registry), and
C                    may have other meanings within those routines
C REAL*8  R12(3) : Vector from first to second atom
C ************************* OUTPUT **********************************
C REAL*8 S12      : Matrix element between orbitals.
C REAL*8 DSDR(3)  : Derivative (gradient) of S12 with respect to R12.

C Length units are arbitrary, but must be consistent in MATEL, RCUT
C   and EVALUATE. The laplacian unit is (length unit)**(-2).
C ************************* BEHAVIOUR *******************************
C If |R12| > RCUT(IS1,IO1) + RCUT(IS2,IO2), returns exactly zero.

      !> Operation code
      character(len=*), intent(in) :: OPERAT
      !> Global index of 1st function
      integer, intent(in)      :: IG1
      !> Global index of 2nd function
      integer, intent(in)      :: IG2

      !> Vector from first to second atom
      real(dp), intent(in)    :: R12(3)

      !> Matrix element
      real(dp), intent(out)    :: S12
      !> Derivative (gradient) of S12 with respect to R12.
      real(dp), intent(out)    :: DSDR(3)

      integer :: top

      top = num_orb + num_kb + num_dftu

      select case (operat(1:1))
      case ( 'S' )

         if ( (ig1 <= top) .and. (is_orb(ig2)) ) then
            call tab_S%get_matel(IG1, IG2, R12, S12, DSDR)
         else if (is_orb(ig1) .and. (is_kb(ig2))) then
            call tab_S_opt%get_matel(IG1, IG2, R12, S12, DSDR)
         else if (is_orb(ig1) .and. (is_wannier_proj(ig2))) then
            call tab_S_wann%get_matel(IG1, IG2, R12, S12, DSDR)
         else
            call die("Cannot process 'S' in matel")
         endif

#ifdef SIESTA__THERMAL_FLUX
      case ( 'G' )

         ! Get index from rest of string
         i_grad = coord_spec_to_index(operat(2:2))
         if (is_orb(ig1) .and. (is_kb(ig2))) then
            call tab_S_Grad_KB(i_grad)%get_matel(IG1, IG2,
     $                         R12, S12, DSDR)
         else
            call die("Cannot process 'G' in matel")
         endif

      case ( 'H' )

         ! Get indexes from rest of string
         i_r = coord_spec_to_index(operat(2:2))
         i_grad = coord_spec_to_index(operat(3:3))

         if (is_orb(ig1) .and. (is_kb(ig2))) then
            call tab_XYZ_Grad_KB(i_r,i_grad)%get_matel(
     $           IG1, IG2, R12, S12, DSDR)
         else
            call die("Cannot process 'H' in matel")
         endif
#endif
      case ( 'T' )

         if ( is_orb(ig1) .and. (is_orb(ig2))) then
            call tab_T%get_matel(IG1, IG2, R12, S12, DSDR)
         else if  (is_vna(ig1) .and. (is_vna(ig2))) then
            call tab_TA%get_matel(IG1, IG2, R12, S12, DSDR)
         else
            call die("Cannot process 'T' in matel")
         endif

      case ( 'X' )

         if  (is_orb(ig1) .and. (is_orb(ig2))) then
            call tab_X%get_matel(IG1, IG2, R12, S12, DSDR)
         else if (is_orb(ig1) .and. (is_kb(ig2))) then
            call tab_X_opt%get_matel(IG1, IG2, R12, S12, DSDR)
         else
            call die("Cannot process 'X' in matel")
         endif

      case ( 'Y' )

         if ( is_orb(ig1) .and. (is_orb(ig2))) then
            call tab_Y%get_matel(IG1, IG2, R12, S12, DSDR)
         else if (is_orb(ig1) .and. (is_kb(ig2))) then
            call tab_Y_opt%get_matel(IG1, IG2, R12, S12, DSDR)
         else
            call die("Cannot process 'Y' in matel")
         endif

      case ( 'Z' )

         if  (is_orb(ig1) .and. (is_orb(ig2))) then
            call tab_Z%get_matel(IG1, IG2, R12, S12, DSDR)
         else if (is_orb(ig1) .and. (is_kb(ig2))) then
            call tab_Z_opt%get_matel(IG1, IG2, R12, S12, DSDR)
         else
            call die("Cannot process 'Z' in matel")
         endif

         case ( 'U' )
            call die("Cannot process 'U' in matel")

         case default
            call die("Unrecognized 'operat' in matel")

         end select

      end subroutine new_matel

!----------
!     Convenience functions

      !> Translates {X,Y,Z} to {1, 2, 3}
      function coord_spec_to_index(c) result (i)
      character(len=1) :: c
      integer i

      select case (c)
      case ( 'X' )
         i = 1
      case ( 'Y' )
         i = 2
      case ( 'Z' )
         i = 3
      case default
         call die("wrong coord spec in matel string")
      end select
      end function coord_spec_to_index

      function is_orb(ig) result(res)
      integer, intent(in) :: ig
      logical             :: res

      res = ( ig > 0 .and. ig <= num_orb)
      end function is_orb

      function is_kb(ig) result(res)
      integer, intent(in) :: ig
      logical             :: res

      res = ( ig > num_orb .and. ig <= (num_orb+num_kb))
      end function is_kb

      function is_vna(ig) result(res)
      integer, intent(in) :: ig
      logical             :: res

      integer :: top

      top = num_orb + num_kb + num_dftu
      res = ( ig > top .and. ig <= (top + num_va))
      end function is_vna

      function is_wannier_proj(ig) result(res)
      integer, intent(in) :: ig
      logical             :: res

      integer :: top_vna

      top_vna = num_orb + num_kb + num_dftu + num_va

      res = ( ig > top_vna .and. ig <= (top_vna + num_wannier_projs))
      end function is_wannier_proj

      end module matel_mod
