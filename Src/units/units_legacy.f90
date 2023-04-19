! 
! Copyright (C) 1996-2022	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module units_legacy_m

  use units_common_m

  implicit none
  private
  
  integer, parameter :: dp = selected_real_kind(14,100)

  ! pi to 50 digits
  real(dp), parameter, public :: pi = 3.14159265358979323846264338327950288419716939937510_dp
  real(dp), parameter, public :: pi2 = pi * 2._dp
  real(dp), parameter, public :: deg = pi / 180.0_dp

  real(dp), parameter, public :: kg     = 1.28460971172331e27_dp
  real(dp), parameter, public :: Joule  = 1._dp/2.1798723611035e-18_dp
  real(dp), parameter, public :: hbar   = 6.62607015e-19_dp * Joule /  pi2 ! in Ry/fs

  real(dp), parameter, public :: Ang    = 1._dp / 0.529177_dp
  real(dp), parameter, public :: eV     = 1._dp / 13.60580_dp
  real(dp), parameter, public :: kBar   = 1._dp / 1.47108e5_dp
  real(dp), parameter, public :: Kelvin = eV / 11604.45_dp
  real(dp), parameter, public :: Debye  = 0.393430_dp
  real(dp), parameter, public :: amu    = 2.133107_dp
  real(dp), parameter, public :: Ryd_time = 1._dp/0.04837769_dp

  real(dp), parameter, public :: GPa = kBar * 10

  public :: inquire_unit

  ! Data generated from table: legacy
  integer, parameter :: nu = 89
  character(8), save :: dimm(nu)
  character(10), save :: name(nu)
  real(dp), save :: unit(nu)
  integer :: iu

  data (dimm(iu), name(iu), unit(iu), iu=1, 3) / &
      'mass    ', 'g         ', 1.d-3, &
      'mass    ', 'kg        ', 1.d0, &
      'mass    ', 'amu       ', 1.66054d-27 /

  data (dimm(iu), name(iu), unit(iu), iu=4, 9) / &
      'length  ', 'm         ', 1.d0, &
      'length  ', 'cm        ', 1.d-2, &
      'length  ', 'nm        ', 1.d-9, &
      'length  ', 'pm        ', 1.d-12, &
      'length  ', 'Ang       ', 1.d-10, &
      'length  ', 'Bohr      ', 0.529177d-10 /

  data (dimm(iu), name(iu), unit(iu), iu=10, 19) / &
      'energy  ', 'J         ', 1.d0, &
      'energy  ', 'kJ        ', 1.d3, &
      'energy  ', 'erg       ', 1.d-7, &
      'energy  ', 'meV       ', 1.60219d-22, &
      'energy  ', 'eV        ', 1.60219d-19, &
      'energy  ', 'mRy       ', 2.17991d-21, &
      'energy  ', 'Ry        ', 2.17991d-18, &
      'energy  ', 'mHa       ', 4.35982d-21, &
      'energy  ', 'mHartree  ', 4.35982d-21, &
      'energy  ', 'Ha        ', 4.35982d-18 /
  data (dimm(iu), name(iu), unit(iu), iu=20, 29) / &
      'energy  ', 'Hartree   ', 4.35982d-18, &
      'energy  ', 'K         ', 1.38066d-23, &
      'energy  ', 'Kelvin    ', 1.38066d-23, &
      'energy  ', 'kcal/mol  ', 6.94780d-21, &
      'energy  ', 'kJ/mol    ', 1.6606d-21, &
      'energy  ', 'Hz        ', 6.6262d-34, &
      'energy  ', 'THz       ', 6.6262d-22, &
      'energy  ', 'cm-1      ', 1.986d-23, &
      'energy  ', 'cm^-1     ', 1.986d-23, &
      'energy  ', 'cm**-1    ', 1.986d-23 /

  data (dimm(iu), name(iu), unit(iu), iu=30, 39) / &
      'time    ', 's         ', 1.d0, &
      'time    ', 'ns        ', 1.d-9, &
      'time    ', 'ps        ', 1.d-12, &
      'time    ', 'fs        ', 1.d-15, &
      'time    ', 'min       ', 60.d0, &
      'time    ', 'mins      ', 60.d0, &
      'time    ', 'hour      ', 3600.d0, &
      'time    ', 'hours     ', 3600.d0, &
      'time    ', 'day       ', 86400.d0, &
      'time    ', 'days      ', 86400.d0 /

  data (dimm(iu), name(iu), unit(iu), iu=40, 43) / &
      'force   ', 'N         ', 1.d0, &
      'force   ', 'eV/Ang    ', 1.60219d-9, &
      'force   ', 'Ry/Bohr   ', 4.11943d-8, &
      'force   ', 'Ha/Bohr   ', 8.23886d-08 /

  data (dimm(iu), name(iu), unit(iu), iu=44, 53) / &
      'pressure', 'Pa        ', 1.d0, &
      'pressure', 'GPa       ', 1.d9, &
      'pressure', 'atm       ', 1.01325d5, &
      'pressure', 'bar       ', 1.d5, &
      'pressure', 'kbar      ', 1.d8, &
      'pressure', 'Mbar      ', 1.d11, &
      'pressure', 'eV/Ang**3 ', 1.60219d11, &
      'pressure', 'eV/Ang^3  ', 1.60219d11, &
      'pressure', 'Ry/Bohr**3', 1.47108d13, &
      'pressure', 'Ry/Bohr^3 ', 1.47108d13 /
  data (dimm(iu), name(iu), unit(iu), iu=54, 55) / &
      'pressure', 'Ha/Bohr^3 ', 2.94216d13, &
      'pressure', 'Ha/Bohr**3', 2.94216d13 /

  data (dimm(iu), name(iu), unit(iu), iu=56, 59) / &
      'surftens', 'N/m       ', 1.d0, &
      'surftens', 'mN/m      ', 1.d3, &
      'surftens', 'dyn/cm    ', 1.d3, &
      'surftens', 'erg/cm**2 ', 1.d3 /

  data (dimm(iu), name(iu), unit(iu), iu=60, 61) / &
      'charge  ', 'c         ', 1.d0, &
      'charge  ', 'e         ', 1.602177d-19 /

  data (dimm(iu), name(iu), unit(iu), iu=62, 66) / &
      'dipole  ', 'c*m       ', 1.d0, &
      'dipole  ', 'D         ', 3.33564d-30, &
      'dipole  ', 'Debye     ', 3.33564d-30, &
      'dipole  ', 'e*Bohr    ', 8.47835d-30, &
      'dipole  ', 'e*Ang     ', 1.602177d-29 /

  data (dimm(iu), name(iu), unit(iu), iu=67, 68) / &
      'mominert', 'kg*m**2   ', 1.d0, &
      'mominert', 'Ry*fs**2  ', 2.17991d-48 /

  data (dimm(iu), name(iu), unit(iu), iu=69, 77) / &
      'efield  ', 'V/m       ', 1.d0, &
      'efield  ', 'V/cm      ', 1.d2, &
      'efield  ', 'V/um      ', 1.d6, &
      'efield  ', 'V/nm      ', 1.d9, &
      'efield  ', 'V/Ang     ', 1.d10, &
      'efield  ', 'V/Bohr    ', 1.8897268d10, &
      'efield  ', 'Ry/Bohr/e ', 2.5711273d11, &
      'efield  ', 'Ha/Bohr/e ', 5.1422546d11, &
      'efield  ', 'Har/Bohr/e', 5.1422546d11 /

  data (dimm(iu), name(iu), unit(iu), iu=78, 79) / &
      'angle   ', 'deg       ', 1.d0, &
      'angle   ', 'rad       ', 5.72957795d1 /

  data (dimm(iu), name(iu), unit(iu), iu=80, 87) / &
      'torque  ', 'meV/deg   ', 1.0d-3, &
      'torque  ', 'meV/rad   ', 1.745533d-5, &
      'torque  ', 'eV/deg    ', 1.0d0, &
      'torque  ', 'eV/rad    ', 1.745533d-2, &
      'torque  ', 'mRy/deg   ', 13.6058d-3, &
      'torque  ', 'mRy/rad   ', 0.237466d-3, &
      'torque  ', 'Ry/deg    ', 13.6058d0, &
      'torque  ', 'Ry/rad    ', 0.237466d0 /

  data (dimm(iu), name(iu), unit(iu), iu=88, 89) / &
      'bfield  ', 'Tesla     ', 1., &
      'bfield  ', 'G         ', 1.d-4 /

contains

  subroutine inquire_unit(unit_str, stat, phys_dim, unit_name, unit_value)
    character(len=*), intent(in)   :: unit_str   ! unit specification
    character(len=*), intent(out)  :: phys_dim   ! physical dimension (e.g. 'mass')
    character(len=*), intent(out)  :: unit_name  ! unit name (e.g. 'g')
    real(dp), intent(out)          :: unit_value ! actual value (e.g. 1.e-3)
    integer, intent(out)           :: stat       ! status code

    call inquire_unit_table(unit_str, stat, phys_dim, unit_name, unit_value, &
        nu, dimm, name, unit)

  end subroutine inquire_unit

end module units_legacy_m
