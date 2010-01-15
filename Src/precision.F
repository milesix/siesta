! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module precision
!
!     Precision handling
!
      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: dp = selected_real_kind(14,100)
!
!     Select precision of certain parts of the program
!     These are set through preprocessor directives. The
!     default behavior is to use single-precision variables
!     for the values of the orbitals on the grid, the Broyden
!     mixing auxiliary arrays, and the O(N) arrays, and 
!     double precision for the grid function arrays.
!
#ifdef BROYDEN_DP
      integer, parameter :: broyden_p = dp
#else
      integer, parameter :: broyden_p = sp
#endif

#ifdef GRID_SP
      integer, parameter :: grid_p = sp
      integer, parameter :: phi_grid_p = sp
#elif defined(GRID_DP)
      integer, parameter :: grid_p = dp
      integer, parameter :: phi_grid_p = dp
#else
      integer, parameter :: grid_p = dp
      integer, parameter :: phi_grid_p = sp
#endif

#ifdef ON_DP
      integer, parameter :: on_p = dp
#else
      integer, parameter :: on_p = sp
#endif

      public

      end module precision

