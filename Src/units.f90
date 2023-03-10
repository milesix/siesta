! 
! Copyright (C) 1996-2022	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module units
  
  ! This module defines various unit conversion factors from Siesta's
  ! internal units.

  !=================================================================================
  
  use precision, only : dp
  use units_table, only: inquire_unit

  implicit none

  ! Internally, siesta works with length: Bohr.
  !                               energy: Rydberg.
  !                                 time: femtosecond

!  The easy way to make sense of units conversion:

!  real(dp), parameter :: Bohr   = 1.0_dp
!  real(dp), parameter :: Rydberg = 1.0_dp
!  real(dp), parameter :: Femtosecond = 1.0_dp
!
!  Ang = Bohr / 0.529177
!   eV = Rydberg / 13.60580
!  Joule = eV / 1.6e-19_dp
!  Meter = Ang / 1.0e-10_dp
!  Pascal = Joule/Meter**2
!   kBar  = Pascal * 1.0e4
!  Ryd^-1 (time) = fs/0.04837769
!   .... and so on.

  ! For consistency, these should be better defined in terms of the units table.
  ! Parameter initialization cannot use non-intrinsic functions, though.
  ! These could be then made into 'protected variables' and initialized
  ! in a call to a "units initialization function".

  !   use appropriate module,
  !   perhaps fdf for fdf_convfac?
  
  !!    Ang = unit_conversion_factor('ang','bohr',...)
  !!    eV  = unit_conversion_factor('eV','ryd',...)
  public :: inquire_unit
  public :: units_initialization
  private
  
  real(dp),   public, protected :: Ang   ! = 1._dp / 0.529177_dp
  real(dp),   public, protected :: eV    ! = 1._dp / 13.60580_dp
  real(dp),   public, protected :: kBar  ! = 1._dp / 1.47108e5_dp
  real(dp),   public, protected :: GPa   ! = kBar * 10
  real(dp),   public, protected :: Kelvin ! = eV / 11604.45_dp
  real(dp),   public, protected :: Debye  ! = 0.393430_dp
  real(dp),   public, protected :: amu    ! = 2.133107_dp
  real(dp),   public, protected :: Ryd_time !  = 1._dp/0.04837769_dp

! pi to 50 digits
  real(dp), parameter, public :: pi = 3.14159265358979323846264338327950288419716939937510_dp
  real(dp), parameter, public :: deg = pi / 180.0_dp

CONTAINS

  subroutine units_initialization()
     Ang = unit_conversion_factor('ang','bohr')
     eV  = unit_conversion_factor('eV','ry')
     kBar  = 1000 * unit_conversion_factor('bar','ry/bohr**3')
     GPa   = 10 * kBar
     Kelvin = unit_conversion_factor('kelvin', 'ry')
     Debye = unit_conversion_factor('debye', 'e*bohr')
     amu    = 2.133107_dp
     Ryd_time = 1._dp/0.04837769_dp
  
  end subroutine units_initialization

  function unit_conversion_factor(from,to) result(factor)
    character(len=*), intent(in) :: from, to
    real(dp) :: factor

    integer stat
    character(len=20) :: phys_dim, unit_name
    real(dp) :: val_to, val_from

    call inquire_unit(from, stat, phys_dim, unit_name, val_from)
    if (stat /= 0) stop 'from unit'
    call inquire_unit(to, stat, phys_dim, unit_name, val_to)
    if (stat /= 0) stop 'to unit'

    factor = val_from / val_to

  end function unit_conversion_factor
  
end module units
