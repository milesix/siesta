! 
! Copyright (C) 1996-2022	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module units
  
  ! Comprehensive handling of units (WIP)

  ! This (client-side) module contains the materials-physics "units
  ! table" formerly in the fdf package, and a 'unit_conversion_factor'
  ! function that performs unit conversions in a form tailored to the
  ! domain and the preferred use case.
  !
  ! 'unit_conversion_factor' can be used by fdf (through the use of a
  ! procedure pointer) to provide unit conversion to that package. For
  ! convenience, fdf might re-export the unit converter
  ! ('fdf_convfac'), but this is not strictly necessary, and other
  ! simplified handlers can be written.
  !
  ! fdf should only be concerned with processing two abstract unit
  ! specification strings, and not with the details of the processing,
  ! or the domain.
  !

  ! In addition, this module defines various unit conversion factors from Siesta's
  ! internal units.

  !=================================================================================
  
  use precision, only : dp

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

  !! real(dp), public, protected :: Ang
  !! ...
  !! subroutine units_initialization()
  !!    Ang = unit_conversion_factor('ang','bohr',...)
  !!    eV  = unit_conversion_factor('ang','bohr',...)
  !!    ...
  !! end subroutine units_initialization
  
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

  ! Note re-designed table structure.
  ! Instead of using a 'data' statement, we employ derived-type initializers
  ! and a parameter array.
  !
  ! In the Fortran 95 standard, a maximum of 39 continuation lines are allowed, but most
  ! compilers would allow more.
  ! In Fortran 2003 and Fortran 2008, a maximum of 255 are allowed.

  ! The new format reduces the need for the 'fdf_units.py' script. The extra "consistency checks"
  ! implemented in that script could be part of an initialization sanity check in this module.

      type, private :: unit
        character(8) :: dimm
        character(10) :: name
        real(dp)      :: value
      end type

      type(unit), private, parameter, dimension(85) :: table = [ &
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
          unit('torque  ', 'ry/rad    ', 0.237466d0) ]

      private
      public :: unit_conversion_factor

    CONTAINS

      ! Strings 'from' and 'to' could be anything that this client chooses
      ! For example, to request a specific "quantity", we could use the
      ! form (suggested, and to be implemented below if needed):
      !
      ! B_field = fdf_get("magnetic-field", 0.0_dp, "G@bfield" )
      !
      ! which means "get me G (Gauss) in the bfield group in the table" (not grams)"
      !
      ! The group specification can do away with most ambiguities regarding case,
      ! but if case-sensitivity is allowed, the syntax could be:
      !
      ! B_field = fdf_get("magnetic-field", 0.0_dp, "!G@bfield" )
      !
      
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


!=====================================================================


!
!   Case-insensitive lexical equal-to comparison
!
    FUNCTION leqi(string1, string2)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(len=*) :: string1, string2

!-------------------------------------------------------------- Output Variables
      logical          :: leqi

!--------------------------------------------------------------- Local Variables
      logical          :: completed
      character        :: char1, char2
      integer          :: i, len1, len2, lenc

!------------------------------------------------------------------------- BEGIN
      len1 = len(string1)
      len2 = len(string2)
      lenc = min(len1, len2)

      i = 1
      leqi      = .TRUE.
      completed = .FALSE.
      do while((.not. completed) .and. (i .le. lenc))
        char1 = string1(i:i)
        char2 = string2(i:i)
        call chrcap(char1, 1)
        call chrcap(char2, 1)
        if (char1 .ne. char2) then
          leqi      = .FALSE.
          completed = .TRUE.
        endif

        i = i + 1
      enddo

      if (leqi) then
        if ((len1 .gt. lenc) .and. (string1(lenc+1:len1) .ne. ' '))     &
          leqi = .FALSE.
        if ((len2 .gt. lenc) .and. (string2(lenc+1:len2) .ne. ' '))     &
          leqi = .FALSE.
      endif

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION leqi


!   CHRCAP accepts a STRING of NCHAR characters and replaces
!   any lowercase letters by uppercase ones.
!
    SUBROUTINE chrcap(string, nchar)
      implicit none
!--------------------------------------------------------------- Input Variables
      integer      :: nchar

!-------------------------------------------------------------- Output Variables
      character(*) :: string

!--------------------------------------------------------------- Local Variables
      integer  :: i, itemp, ncopy

!------------------------------------------------------------------------- BEGIN
      if (nchar .le. 0) then
        ncopy = LEN(string)
      else
        ncopy = nchar
      endif

      do i= 1, ncopy
        if (LGE(string(i:i),'a') .and. LLE(string(i:i),'z')) then
          itemp = ICHAR(string(i:i)) + ICHAR('A') - ICHAR('a')
          string(i:i) = CHAR(itemp)
        endif
      enddo

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE chrcap

end module units
