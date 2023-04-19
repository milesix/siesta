! 
! Copyright (C) 1996-2022	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module units
  
#ifdef SIESTA__UNITS_ORIGINAL
  use units_legacy_m
#else
  use units_codata2018_m
#endif

  implicit none

  public

! Internally, siesta works with length: Bohr.
!                               energy: Rydberg.
!                                 time: femtosecond
!  real(dp), parameter :: Bohr   = 1.0_dp
!  real(dp), parameter :: Rydberg = 1.0_dp
!  real(dp), parameter :: Femtosecond = 1.0_dp
!
!  Ang = Bohr / 0.529177
!   eV = Rydberg / 13.60580
!  Joule = eV / 1.6e-19_dp
!  Meter = Ang / 1.0e-10_dp
!  Pascal = Joule/Meter**3
!   kBar  = Pascal * 1.0e4
!  Ryd^-1 (time) = fs/0.04837769
!   .... and so on.
!
! amu is the conversion from standard atomic mass (Da) to the units
! that SIESTA uses when calculating the kinetic energy of atomic
! nuclei, and for the atomic equations of motion in MD.
!
! The amu is defined 1/12 of the mass of a C-12 isotope (about 1822 times
! the mass of the electron); so in principle:
!
! 1 amu = 1.66053906660e-27 kg
! Thus, amu = 1.66053906660e-27 * Joule * second^2 / meter^2

! For consistency, these should be better defined in terms of the units table.
! Parameter initialization cannot use non-intrinsic functions, though.
! These could be then made into 'protected variables' and initialized
! in a call to a "units initialization function".

  !! real(dp), public, protected :: Ang
  !! ...
  !! subroutine units_initialization()
  !!    Ang = unit_conversion_factor('ang','bohr',...)
  !!    eV  = unit_conversion_factor('eV','ryd',...)
  !!    ...
  !! end subroutine units_initialization
  
end module
