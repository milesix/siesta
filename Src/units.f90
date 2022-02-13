! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module units
  ! Define various unit conversion factors from internal units.

  ! internally, siesta works with length: Bohr.
  !                               energy: Rydberg.
  !                                 time: femtosecond

  use precision, only : dp
!     Maybe separate leqi from fdf
  use fdf, only: leqi


  implicit none

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
  
  real(dp), parameter, public :: Ang    = 1._dp / 0.529177_dp
  real(dp), parameter, public :: eV     = 1._dp / 13.60580_dp
  real(dp), parameter, public :: kBar   = 1._dp / 1.47108e5_dp
  real(dp), parameter, public :: GPa    = kBar * 10
  real(dp), parameter, public :: Kelvin = eV / 11604.45_dp
  real(dp), parameter, public :: Debye  = 0.393430_dp
  real(dp), parameter, public :: amu    = 2.133107_dp
  real(dp), parameter, public :: Ryd_time = 1._dp/0.04837769_dp

! pi to 50 digits
  real(dp), parameter, public :: pi = 3.14159265358979323846264338327950288419716939937510_dp
  real(dp), parameter, public :: deg = pi / 180.0_dp

      type, private :: unit
        character(8) :: dimm
        character(10) :: name
        real(dp)      :: value
      end type

      type(unit), dimension(86), private, parameter :: table = [ &
          unit('mass    ', 'g         ', 1.d-3), &
          unit('mass    ', 'kg        ', 1.d0), &
          unit('mass    ', 'amu       ', 1.66054d-27), &
      !
          unit('length  ', 'm         ', 1.d0), &
          unit('length  ', 'cm        ', 1.d-2), &
          unit('length  ', 'nm        ', 1.d-9), &
          unit('length  ', 'pm        ', 1.d-12), &
          unit('length  ', 'ang       ', 1.d-10), &
          unit('length  ', 'bohr      ', 0.529177d-10), &
      !
          unit('energy  ', 'j         ', 1.d0), &
          unit('energy  ', 'kj        ', 1.d3), &
          unit('energy  ', 'erg       ', 1.d-7), &
          unit('energy  ', 'mev       ', 1.60219d-22), &
          unit('energy  ', 'ev        ', 1.60219d-19), &
          unit('energy  ', 'mry       ', 2.17991d-21), &
          unit('energy  ', 'ry        ', 2.17991d-18), &
          unit('energy  ', 'mha       ', 4.35982d-21), &
          unit('energy  ', 'mhartree  ', 4.35982d-21), &
          unit('energy  ', 'ha        ', 4.35982d-18), &
          unit('energy  ', 'hartree   ', 4.35982d-18), &
          unit('energy  ', 'k         ', 1.38066d-23), &
          unit('energy  ', 'kelvin    ', 1.38066d-23), &
          unit('energy  ', 'kcal/mol  ', 6.94780d-21), &
          unit('energy  ', 'kj/mol    ', 1.6606d-21), &
          unit('energy  ', 'hz        ', 6.6262d-34), &
          unit('energy  ', 'thz       ', 6.6262d-22), &
          unit('energy  ', 'cm-1      ', 1.986d-23), &
          unit('energy  ', 'cm^-1     ', 1.986d-23), &
          unit('energy  ', 'cm**-1    ', 1.986d-23), &
      !
          unit('time    ', 's         ', 1.d0), &
          unit('time    ', 'au        ', 2.418884326505e-17), &
          unit('time    ', 'ns        ', 1.d-9), &
          unit('time    ', 'ps        ', 1.d-12), &
          unit('time    ', 'fs        ', 1.d-15), &
          unit('time    ', 'min       ', 60.d0), &
          unit('time    ', 'mins      ', 60.d0), &
          unit('time    ', 'hour      ', 3600.d0), &
          unit('time    ', 'hours     ', 3600.d0), &
          unit('time    ', 'day       ', 86400.d0), &
          unit('time    ', 'days      ', 86400.d0), &
      !
          unit('force   ', 'n         ', 1.d0), &
          unit('force   ', 'ev/ang    ', 1.60219d-9), &
          unit('force   ', 'ry/bohr   ', 4.11943d-8), &
          unit('force   ', 'ha/bohr   ', 8.23886d-8), &
      !
          unit('pressure', 'pa        ', 1.d0), &
          unit('pressure', 'gpa       ', 1.d9), &
          unit('pressure', 'atm       ', 1.01325d5), &
          unit('pressure', 'bar       ', 1.d5), &
          unit('pressure', 'mbar      ', 1.d11), &
          unit('pressure', 'ev/ang**3 ', 1.60219d11), &
          unit('pressure', 'ev/ang^3  ', 1.60219d11), &
          unit('pressure', 'ry/bohr**3', 1.47108d13), &
          unit('pressure', 'ry/bohr^3 ', 1.47108d13), &
          unit('pressure', 'ha/bohr^3 ', 2.94216d13), &
          unit('pressure', 'ha/bohr**3', 2.94216d13), &
      !
          unit('surftens', 'n/m       ', 1.d0), &
          unit('surftens', 'mn/m      ', 1.d3), &
          unit('surftens', 'dyn/cm    ', 1.d3), &
          unit('surftens', 'erg/cm**2 ', 1.d3), &
      !
          unit('charge  ', 'c         ', 1.d0), &
          unit('charge  ', 'e         ', 1.602177d-19), &
      !
          unit('dipole  ', 'c*m       ', 1.d0), &
          unit('dipole  ', 'd         ', 3.33564d-30), &
          unit('dipole  ', 'debye     ', 3.33564d-30), &
          unit('dipole  ', 'e*bohr    ', 8.47835d-30), &
          unit('dipole  ', 'e*ang     ', 1.602177d-29), &
      !
          unit('mominert', 'kg*m**2   ', 1.d0), &
          unit('mominert', 'ry*fs**2  ', 2.17991d-48), &
      !
          unit('efield  ', 'v/m       ', 1.d0), &
          unit('efield  ', 'v/nm      ', 1.d9), &
          unit('efield  ', 'v/ang     ', 1.d10), &
          unit('efield  ', 'v/bohr    ', 1.8897268d10), &
          unit('efield  ', 'ry/bohr/e ', 2.5711273d11), &
          unit('efield  ', 'ha/bohr/e ', 5.1422546d11), &
          unit('efield  ', 'har/bohr/e', 5.1422546d11), &
      !
          unit('angle   ', 'deg       ', 1.d0), &
          unit('angle   ', 'rad       ', 5.72957795d1), &
      !
          unit('torque  ', 'mev/deg   ', 1.0d-3), &
          unit('torque  ', 'mev/rad   ', 1.745533d-5), &
          unit('torque  ', 'ev/deg    ', 1.0d0), &
          unit('torque  ', 'ev/rad    ', 1.745533d-2), &
          unit('torque  ', 'mry/deg   ', 13.6058d-3), &
          unit('torque  ', 'mry/rad   ', 0.237466d-3), &
          unit('torque  ', 'ry/deg    ', 13.6058d0), &
          unit('torque  ', 'ry/rad    ', 0.237466d0), &
      !
      unit('dummy', 'dummy', 0.0) ]

      private
      public :: unit_conversion_factor

    CONTAINS

      function unit_conversion_factor(from,to,stat,msg) result (factor)
        character(len=*), intent(in)   :: from, to
        integer, intent(out)           :: stat
        character(len=96), intent(out) :: msg
        real(dp)                       :: factor
        
        integer :: ifrom, ito, iu

        stat = 0
        factor = 0.0_dp
        
        ifrom = 0
        ito   = 0
        do iu= 1, size(table)
           if (leqi(table(iu)%name, from)) ifrom = iu
           if (leqi(table(iu)%name, to))   ito   = iu
        end do

        if (ifrom .eq. 0) then
           stat = 1
           write(msg,*) 'unknown unit = ', from
           return
        endif
        if (ito .eq. 0) then
           stat = 2
           write(msg,*) 'unknown unit = ', to
           return
        endif
        if (leqi(table(ifrom)%dimm, table(ito)%dimm)) then
           factor = table(ifrom)%value / table(ito)%value
           return
        else
           stat = -1
           write(msg,*) 'unit''s physical dimensions don''t match: ',      &
                     from, ', ', to
           return
        endif
        
      end function unit_conversion_factor

end module units
