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
  ! table" formerly in the fdf package, and a 'inquire_unit'
  ! function.
  !
  ! 'inquire_unit' can be used by fdf (through the use of a
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
  !!    eV  = unit_conversion_factor('eV','ryd',...)
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

  
  public :: inquire_unit
  private

  integer, parameter :: nu = 80
  integer :: iu
      character(8) :: dimm(nu)
      character(10) :: name(nu)
      real(dp) :: unit(nu)

      data (dimm(iu), name(iu), unit(iu), iu=1, 3) / &
          'mass    ', 'g         ', 1.d-3, &
          'mass    ', 'kg        ', 1.d0, &
          'mass    ', 'amu       ', 1.66054d-27 /

      data (dimm(iu), name(iu), unit(iu), iu=4, 9) / &
          'length  ', 'm         ', 1.d0, &
          'length  ', 'cm        ', 1.d-2, &
          'length  ', 'nm        ', 1.d-9, &
          'length  ', 'pm        ', 1.d-12, &
          'length  ', 'ang       ', 1.d-10, &
          'length  ', 'bohr      ', 0.529177d-10 /

      data (dimm(iu), name(iu), unit(iu), iu=10, 19) / &
          'energy  ', 'j         ', 1.d0, &
          'energy  ', 'kj        ', 1.d3, &
          'energy  ', 'erg       ', 1.d-7, &
          'energy  ', 'mev       ', 1.60219d-22, &
          'energy  ', 'ev        ', 1.60219d-19, &
          'energy  ', 'mry       ', 2.17991d-21, &
          'energy  ', 'ry        ', 2.17991d-18, &
          'energy  ', 'mha       ', 4.35982d-21, &
          'energy  ', 'mhartree  ', 4.35982d-21, &
          'energy  ', 'ha        ', 4.35982d-18 /
      data (dimm(iu), name(iu), unit(iu), iu=20, 29) / &
          'energy  ', 'hartree   ', 4.35982d-18, &
          'energy  ', 'k         ', 1.38066d-23, &
          'energy  ', 'kelvin    ', 1.38066d-23, &
          'energy  ', 'kcal/mol  ', 6.94780d-21, &
          'energy  ', 'kj/mol    ', 1.6606d-21, &
          'energy  ', 'hz        ', 6.6262d-34, &
          'energy  ', 'thz       ', 6.6262d-22, &
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
          'force   ', 'n         ', 1.d0, &
          'force   ', 'ev/ang    ', 1.60219d-9, &
          'force   ', 'ry/bohr   ', 4.11943d-8, &
          'force   ', 'ha/bohr   ', 8.23886d-08 /

      data (dimm(iu), name(iu), unit(iu), iu=44, 52) / &
          'pressure', 'pa        ', 1.d0, &
          'pressure', 'gpa       ', 1.d9, &
          'pressure', 'atm       ', 1.01325d5, &
          'pressure', 'bar       ', 1.d5, &
          'pressure', 'mbar      ', 1.d11, &
          'pressure', 'ev/ang**3 ', 1.60219d11, &
          'pressure', 'ev/ang^3  ', 1.60219d11, &
          'pressure', 'ry/bohr**3', 1.47108d13, &
          'pressure', 'ry/bohr^3 ', 1.47108d13 /

      data (dimm(iu), name(iu), unit(iu), iu=53, 54) / &
          'charge  ', 'c         ', 1.d0, &
          'charge  ', 'e         ', 1.602177d-19 /

      data (dimm(iu), name(iu), unit(iu), iu=55, 59) / &
          'dipole  ', 'c*m       ', 1.d0, &
          'dipole  ', 'd         ', 3.33564d-30, &
          'dipole  ', 'debye     ', 3.33564d-30, &
          'dipole  ', 'e*bohr    ', 8.47835d-30, &
          'dipole  ', 'e*ang     ', 1.602177d-29 /

      data (dimm(iu), name(iu), unit(iu), iu=60, 61) / &
          'mominert', 'kg*m**2   ', 1.d0, &
          'mominert', 'ry*fs**2  ', 2.17991d-48 /

      data (dimm(iu), name(iu), unit(iu), iu=62, 68) / &
          'efield  ', 'v/m       ', 1.d0, &
          'efield  ', 'v/nm      ', 1.d9, &
          'efield  ', 'v/ang     ', 1.d10, &
          'efield  ', 'v/bohr    ', 1.8897268d10, &
          'efield  ', 'ry/bohr/e ', 2.5711273d11, &
          'efield  ', 'ha/bohr/e ', 5.1422546d11, &
          'efield  ', 'har/bohr/e', 5.1422546d11 /

      data (dimm(iu), name(iu), unit(iu), iu=69, 70) / &
          'angle   ', 'deg       ', 1.d0, &
          'angle   ', 'rad       ', 5.72957795d1 /

      data (dimm(iu), name(iu), unit(iu), iu=71, 78) / &
          'torque  ', 'mev/deg   ', 1.0d-3, &
          'torque  ', 'mev/rad   ', 1.745533d-5, &
          'torque  ', 'ev/deg    ', 1.0d0, &
          'torque  ', 'ev/rad    ', 1.745533d-2, &
          'torque  ', 'mry/deg   ', 13.6058d-3, &
          'torque  ', 'mry/rad   ', 0.237466d-3, &
          'torque  ', 'ry/deg    ', 13.6058d0, &
          'torque  ', 'ry/rad    ', 0.237466d0 /

      data (dimm(iu), name(iu), unit(iu), iu=79, 80) / &
          'bfield  ', 'Tesla     ', 1.0d0, &
          'bfield  ', 'G         ', 1.0d-4 /


CONTAINS

  ! Returns information about a unit in the units table
  ! 
  ! Unit specifications might include an optional 'physical dimension'
  ! qualifier (e.g. 'bfield:g')
  ! In this case, 'phys_dim' returns the physical dimension, and the
  ! qualifier is used to match the unit.
  ! This version is case-insensitive (e.g. 'g' and 'G' could stand for 'Gauss').
  ! As the above example indicates, in the absence of a physical dimension qualifier,
  ! 'g' might be ambiguous ('bfield' or 'mass'?). The routine will return 'stat=-1'
  ! in this case.
  ! Units might be ambiguous in a more serious way: 'meV' and 'MeV' could both be
  ! present in the table. In this case, it might be advisable to use a case-sensitive
  ! version of this routine (replacing 'leqi' by 'leqi_strict' below).
  ! If the unit is not found in the table, the routine returns 'stat=-2'.
  
  subroutine inquire_unit(unit_str, stat, phys_dim, unit_name, unit_value)

    character(len=*), intent(in)   :: unit_str   ! unit specification
    character(len=*), intent(out)  :: phys_dim   ! physical dimension (e.g. 'mass')
    character(len=*), intent(out)  :: unit_name  ! unit name (e.g. 'g')
    real(dp), intent(out)          :: unit_value ! actual value (e.g. 1.e-3)
    integer, intent(out)           :: stat       ! status code

    integer           :: idx_colon, iu, idx
    logical           :: phys_dim_specified, match
    
    idx_colon = index(unit_str,":")
    if (idx_colon /= 0) then
       ! spec includes dimension prefix
       phys_dim = unit_str(1:idx_colon-1)
       unit_name = unit_str(idx_colon+1:)
       phys_dim_specified = .true.
    else
       phys_dim = ""
       unit_name = unit_str
       phys_dim_specified = .false.
    endif

    stat = 0
    idx = 0

    do iu= 1, nu
         match = .false.
         if (leqi(name(iu), unit_name)) then
            if (phys_dim_specified) then
               if (leqi(dimm(iu), phys_dim)) then
                  match = .true.
               endif
            else
               match = .true.
            endif
         endif
         if (match) then
            if (idx /= 0) then  ! ambiguous
               stat = 1
               RETURN
            endif
            idx = iu
         endif
      enddo
      
      if (idx == 0) then
         stat = -2    ! not found
      else
         phys_dim = trim(dimm(idx))
         unit_value = unit(idx)
      endif
      
    end subroutine inquire_unit
    
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

!
!   Examples of eq_func's for search function (Case sensitive)
!
    FUNCTION leqi_strict(str1, str2)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(len=*) :: str1, str2

!-------------------------------------------------------------- Output Variables
      logical          :: leqi_strict

!------------------------------------------------------------------------- BEGIN
      leqi_strict = (str1 .eq. str2)
      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION leqi_strict

!
!   CHRCAP accepts a STRING of NCHAR characters and replaces
!   any lowercase letters by uppercase ones.
!
    SUBROUTINE chrcap(string, nchar)
      implicit none
!--------------------------------------------------------------- Input Variables
      integer  :: nchar

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

!
!   CHRLEN accepts a STRING of NCHAR characters and returns LCHAR,
!   the length of the string up to the last NONBLANK, NONNULL.
!     
    SUBROUTINE chrlen(string, nchar, lchar)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*) :: string
      integer  :: nchar

!-------------------------------------------------------------- Output Variables
      integer  :: lchar

!------------------------------------------------------------------------- BEGIN
      lchar = nchar
      if (lchar .le. 0) lchar = LEN(string)

      do while(((string(lchar:lchar) .eq. ' ') .or. (string(lchar:lchar) &
               .eq. CHAR(0))) .and. (lchar .gt. 0))
        lchar = lchar - 1
      enddo

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE chrlen

end module units
