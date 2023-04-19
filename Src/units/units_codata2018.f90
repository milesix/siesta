! 
! Copyright (C) 1996-2022	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
!> CODATA defined units from the 2018 table
!>
!> These units are unified to follow the CODATA-2018
!> data table.
!> In addition we have some basic definitions such as
!> pi and conversion to degree (from radians).
!> The conversion factors are following this convention:
!>
!>   Ang [Ang] == [Bohr]
!>   [Ang] = [Bohr] / Ang
!>
!> meaning that the units are conversion factors from the named unit
!> to the intrinsic Siesta unit.
module units_codata2018_m

  use units_common_m

  implicit none
  private

  integer, parameter :: dp = selected_real_kind(14,100)

  real(dp), parameter, public :: pi = 3.14159265358979323846264338327950288419716939937510_dp
  real(dp), parameter, public :: pi2 = pi * 2._dp
  real(dp), parameter, public :: deg = pi / 180.0_dp

  ! These values are not present in the *legacy* scheme.
  real(dp), parameter, public :: kg     = 1.28460971172331e27_dp
  real(dp), parameter, public :: Joule  = 1._dp/2.1798723611035e-18_dp
  real(dp), parameter, public :: hbar   = 6.62607015e-19_dp * Joule / pi2 ! in Ry/fs

  real(dp), parameter, public :: Ang    = 1.88972612462577017_dp
  real(dp), parameter, public :: eV     = 7.34986443513115789e-2_dp
  real(dp), parameter, public :: kBar   = 6.79786184348648780e-6_dp
  real(dp), parameter, public :: Kelvin = 6.33362312691136091e-6_dp
  real(dp), parameter, public :: Debye  = 3.93430269519899511e-1_dp
  real(dp), parameter, public :: amu    = 2.13314461165032_dp ! 1.6605390666e-27 *kg
  real(dp), parameter, public :: Ryd_time = 1._dp / hbar

  real(dp), parameter, public :: GPa = kBar * 10

  public :: inquire_unit

  ! Data generated from table: codata-2018
  integer, parameter :: nu = 93
  character(8), save :: dimm(nu)
  character(10), save :: name(nu)
  real(dp), save :: unit(nu)
  integer :: iu

  data (dimm(iu), name(iu), unit(iu), iu=1, 4) / &
      'mass    ', 'g         ', 1.d-3, &
      'mass    ', 'kg        ', 1., &
      'mass    ', 'amu       ', 1.66053906660d-27, &
      'mass    ', 'da        ', 1.99264687991118964e-26_dp /

  data (dimm(iu), name(iu), unit(iu), iu=5, 10) / &
      'length  ', 'm         ', 1., &
      'length  ', 'cm        ', 1.d-2, &
      'length  ', 'nm        ', 1.d-9, &
      'length  ', 'pm        ', 1.d-12, &
      'length  ', 'Ang       ', 1.d-10, &
      'length  ', 'Bohr      ', 0.529177210903d-10 /

  data (dimm(iu), name(iu), unit(iu), iu=11, 20) / &
      'energy  ', 'J         ', 1., &
      'energy  ', 'kJ        ', 1.d3, &
      'energy  ', 'erg       ', 1.d-7, &
      'energy  ', 'meV       ', 1.602176634d-22, &
      'energy  ', 'eV        ', 1.602176634d-19, &
      'energy  ', 'mRy       ', 2.1798723611035d-21, &
      'energy  ', 'Ry        ', 2.1798723611035d-18, &
      'energy  ', 'mHa       ', 4.3597447222071d-21, &
      'energy  ', 'Ha        ', 4.3597447222071d-18, &
      'energy  ', 'Hartree   ', 4.3597447222071d-18 /
  data (dimm(iu), name(iu), unit(iu), iu=21, 29) / &
      'energy  ', 'K         ', 1.380649d-23, &
      'energy  ', 'Kelvin    ', 1.380649d-23, &
      'energy  ', 'kJ/mol    ', 1.6605390671738467d-21, &
      'energy  ', 'kcal/mol  ', 6.94769545705537413e-21_dp, &
      'energy  ', 'Hz        ', 6.62607015d-34, &
      'energy  ', 'THz       ', 6.62607015d-22, &
      'energy  ', 'cm-1      ', 1.986445857d-23, &
      'energy  ', 'cm^-1     ', 1.986445857d-23, &
      'energy  ', 'cm**-1    ', 1.986445857d-23 /

  data (dimm(iu), name(iu), unit(iu), iu=30, 39) / &
      'time    ', 's         ', 1., &
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
      'force   ', 'N         ', 1., &
      'force   ', 'eV/Ang    ', 1.60217663399999979e-09_dp, &
      'force   ', 'Ry/Bohr   ', 4.11936174912694464e-08_dp, &
      'force   ', 'Ha/Bohr   ', 8.23872349825407855e-08_dp /

  data (dimm(iu), name(iu), unit(iu), iu=44, 53) / &
      'pressure', 'Pa        ', 1., &
      'pressure', 'GPa       ', 1.d9, &
      'pressure', 'atm       ', 1.01325d5, &
      'pressure', 'bar       ', 1.d5, &
      'pressure', 'kbar      ', 1.d8, &
      'pressure', 'Mbar      ', 1.d11, &
      'pressure', 'eV/Ang**3 ', 1.60217663399999969e+11_dp, &
      'pressure', 'eV/Ang^3  ', 1.60217663399999969e+11_dp, &
      'pressure', 'Ry/Bohr**3', 1.47105078482607109e+13_dp, &
      'pressure', 'Ry/Bohr^3 ', 1.47105078482607109e+13_dp /
  data (dimm(iu), name(iu), unit(iu), iu=54, 55) / &
      'pressure', 'Ha/Bohr**3', 2.94210156965220977e+13_dp, &
      'pressure', 'Ha/Bohr^3 ', 2.94210156965220977e+13_dp /

  data (dimm(iu), name(iu), unit(iu), iu=56, 59) / &
      'surftens', 'N/m       ', 1.d0, &
      'surftens', 'mN/m      ', 1.d3, &
      'surftens', 'dyn/cm    ', 1.d3, &
      'surftens', 'erg/cm**2 ', 1.d3 /

  data (dimm(iu), name(iu), unit(iu), iu=60, 61) / &
      'charge  ', 'c         ', 1., &
      'charge  ', 'e         ', 1.602176634d-19 /

  data (dimm(iu), name(iu), unit(iu), iu=62, 66) / &
      'dipole  ', 'c*m       ', 1., &
      'dipole  ', 'e*Bohr    ', 8.47835362554076597e-30_dp, &
      'dipole  ', 'e*Ang     ', 1.60217663400000003e-29_dp, &
      'dipole  ', 'D         ', 3.33564095198152075e-30_dp, &
      'dipole  ', 'Debye     ', 3.33564095198152075e-30_dp /

  data (dimm(iu), name(iu), unit(iu), iu=67, 68) / &
      'mominert', 'kg*m**2   ', 1., &
      'mominert', 'Ry*fs**2  ', 2.17987236110350002e-48_dp /

  data (dimm(iu), name(iu), unit(iu), iu=69, 78) / &
      'efield  ', 'V/m       ', 1., &
      'efield  ', 'V/cm      ', 1.d2, &
      'efield  ', 'V/um      ', 1.d6, &
      'efield  ', 'V/nm      ', 1.d9, &
      'efield  ', 'V/Ang     ', 1.d10, &
      'efield  ', 'eV/Ang/e  ', 1.d10, &
      'efield  ', 'V/Bohr    ', 1.88972612462577019e+10_dp, &
      'efield  ', 'Ry/Bohr/e ', 2.57110337381623871e+11_dp, &
      'efield  ', 'Ha/Bohr/e ', 5.14220674763259521e+11_dp, &
      'efield  ', 'Har/Bohr/e', 5.14220674763259521e+11_dp /

  data (dimm(iu), name(iu), unit(iu), iu=79, 80) / &
      'angle   ', 'deg       ', 1., &
      'angle   ', 'rad       ', 5.72957795130823229e+01_dp /

  data (dimm(iu), name(iu), unit(iu), iu=81, 90) / &
      'torque  ', 'N*m       ', 1., &
      'torque  ', 'meV/deg   ', 1.602176634d-22, &
      'torque  ', 'meV/rad   ', 2.79632574618201289e-24_dp, &
      'torque  ', 'eV/deg    ', 1.602176634d-19, &
      'torque  ', 'eV/rad    ', 2.79632574618201262e-21_dp, &
      'torque  ', 'mRy/deg   ', 2.1798723611035d-21, &
      'torque  ', 'mRy/rad   ', 3.80459499744788459e-23_dp, &
      'torque  ', 'Ry/deg    ', 2.1798723611035d-18, &
      'torque  ', 'Ry/rad    ', 3.80459499744788462e-20_dp, &
      'torque  ', 'Ha/deg    ', 4.3597447222071d-18 /
  data (dimm(iu), name(iu), unit(iu), iu=91, 91) / &
      'torque  ', 'Ha/rad    ', 7.60918999489594379e-20_dp /

  data (dimm(iu), name(iu), unit(iu), iu=92, 93) / &
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

end module units_codata2018_m
