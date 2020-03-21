!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2020, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! Module for accummulating byte counts from array dimensions.

! This module is primarily intended for determining file
! sizes by passing sizes of arrays stored.
!
! Its basic usage may be illustrated by this:
! \code
! character(len=32) :: mem_unit
! type(byte_count_t) :: byte
!
! call byte%add(4, 10, 20, 30) ! add array(10, 20, 30) with each element 4 bytes
!
! call byte%get_string(mem_unit) ! retrieve a string representing the total size (with unit)
!
! \endcode

module byte_count_m

  use precision, only: sp, dp, i4b, i8b

  implicit none

  private

  type :: byte_count_t
    
    !< Bytes accummulated in kB
    real(dp) :: kB = 0._dp

  contains

    procedure, pass :: reset => reset_
    procedure, pass :: get_string => get_string_
#ifdef NCDF_4
    procedure, pass, private :: add_cdf_basic_, add_cdf_complex_
    generic :: add_cdf => add_cdf_basic_, add_cdf_complex_
#endif
    procedure, pass, private :: add_i4_, add_i8_, add_i8i8_
    procedure, pass, private :: add_ai4_, add_ai8_, add_ai8i8_
    generic :: add => add_i4_, add_i8_, add_i8i8_, &
        add_ai4_, add_ai8_, add_ai8i8_
    
  end type byte_count_t

  public :: byte_count_t
  
contains

  subroutine reset_(this)
    class(byte_count_t), intent(inout) :: this
    this%kB = 0._dp
  end subroutine reset_

  subroutine add_i4_(this, bytes, n1, n2, n3, n4, n5)
    class(byte_count_t), intent(inout) :: this
    integer(i4b), intent(in) :: bytes
    integer(i4b), intent(in) :: n1
    integer(i4b), intent(in), optional :: n2, n3, n4, n5

    real(dp) :: kB

    ! Number of elements
    kB = bytes * real(n1, dp) / 1024._dp
    if ( present(n2) ) kB = kB * real(n2, dp)
    if ( present(n3) ) kB = kB * real(n3, dp)
    if ( present(n4) ) kB = kB * real(n4, dp)
    if ( present(n5) ) kB = kB * real(n5, dp)

    this%kB = this%kB + kB

  end subroutine add_i4_

  subroutine add_ai4_(this, bytes, n)
    class(byte_count_t), intent(inout) :: this
    integer(i4b), intent(in) :: bytes
    integer(i4b), intent(in) :: n(:)

    integer :: i
    real(dp) :: kB

    ! Number of elements
    kB = bytes * real(n(1), dp) / 1024._dp
    do i = 2, size(n)
      kB = kB * real(n(i), dp)
    end do
    this%kB = this%kB + kB

  end subroutine add_ai4_

  subroutine add_ai8_(this, bytes, n)
    class(byte_count_t), intent(inout) :: this
    integer(i4b), intent(in) :: bytes
    integer(i8b), intent(in) :: n(:)

    integer :: i
    real(dp) :: kB

    ! Number of elements
    kB = bytes * real(n(1), dp) / 1024._dp
    do i = 2, size(n)
      kB = kB * real(n(i), dp)
    end do
    this%kB = this%kB + kB

  end subroutine add_ai8_

  subroutine add_ai8i8_(this, bytes, n)
    class(byte_count_t), intent(inout) :: this
    integer(i8b), intent(in) :: bytes
    integer(i8b), intent(in) :: n(:)

    integer :: i
    real(dp) :: kB

    ! Number of elements
    kB = bytes * real(n(1), dp) / 1024._dp
    do i = 2, size(n)
      kB = kB * real(n(i), dp)
    end do
    this%kB = this%kB + kB

  end subroutine add_ai8i8_

  subroutine add_i8_(this, bytes, n1, n2, n3, n4, n5)
    class(byte_count_t), intent(inout) :: this
    integer(i4b), intent(in) :: bytes
    integer(i8b), intent(in) :: n1
    integer(i8b), intent(in), optional :: n2, n3, n4, n5

    real(dp) :: kB

    ! Number of elements
    kB = bytes * real(n1, dp) / 1024._dp
    if ( present(n2) ) kB = kB * real(n2, dp)
    if ( present(n3) ) kB = kB * real(n3, dp)
    if ( present(n4) ) kB = kB * real(n4, dp)
    if ( present(n5) ) kB = kB * real(n5, dp)

    this%kB = this%kB + kB

  end subroutine add_i8_
  
  subroutine add_i8i8_(this, bytes, n1, n2, n3, n4, n5)
    class(byte_count_t), intent(inout) :: this
    integer(i8b), intent(in) :: bytes
    integer(i8b), intent(in) :: n1
    integer(i8b), intent(in), optional :: n2, n3, n4, n5

    real(dp) :: kB

    ! Number of elements
    kB = bytes * real(n1, dp) / 1024._dp
    if ( present(n2) ) kB = kB * real(n2, dp)
    if ( present(n3) ) kB = kB * real(n3, dp)
    if ( present(n4) ) kB = kB * real(n4, dp)
    if ( present(n5) ) kB = kB * real(n5, dp)

    this%kB = this%kB + kB

  end subroutine add_i8i8_

#ifdef NCDF_4
  subroutine add_cdf_basic_(this, nf_var, n1, n2, n3, n4, n5)
    use netcdf_ncdf, only: NF90_INT
    use netcdf_ncdf, only: NF90_FLOAT, NF90_DOUBLE

    class(byte_count_t), intent(inout) :: this
    integer, intent(in) :: nf_var
    integer, intent(in) :: n1
    integer, intent(in), optional :: n2, n3, n4, n5

    real(dp) :: kB

    ! Number of elements
    kB = real(n1, dp) / 1024._dp
    if ( present(n2) ) kB = kB * real(n2, dp)
    if ( present(n3) ) kB = kB * real(n3, dp)
    if ( present(n4) ) kB = kB * real(n4, dp)
    if ( present(n5) ) kB = kB * real(n5, dp)

    select case ( nf_var )
    case ( NF90_INT, NF90_FLOAT )
      kB = kB * 4
    case ( NF90_DOUBLE )
      kB = kB * 8
    end select

    this%kB = this%kB + kB

  end subroutine add_cdf_basic_

  subroutine add_cdf_complex_(this, nf_var, n1, n2, n3, n4, n5)
    use netcdf_ncdf, only: NF90_FLOAT_COMPLEX, NF90_DOUBLE_COMPLEX

    class(byte_count_t), intent(inout) :: this
    logical, intent(in) :: nf_var
    integer, intent(in) :: n1
    integer, intent(in), optional :: n2, n3, n4, n5

    real(dp) :: kB

    ! Number of elements
    kB = real(n1, dp) / 1024._dp
    if ( present(n2) ) kB = kB * real(n2, dp)
    if ( present(n3) ) kB = kB * real(n3, dp)
    if ( present(n4) ) kB = kB * real(n4, dp)
    if ( present(n5) ) kB = kB * real(n5, dp)

    if ( nf_var .eqv. NF90_FLOAT_COMPLEX ) then
      kB = kB * 8
    else if ( nf_var .eqv. NF90_DOUBLE_COMPLEX ) then
      kB = kB * 16
    end if

    this%kB = this%kB + kB

  end subroutine add_cdf_complex_
#endif

  subroutine get_string_(this, mem_str, unit)
    use m_char, only: lcase 
    !< Memory class containing the amount of memory
    class(byte_count_t), intent(in) :: this
    !< Output string for the memory with an appropriate unit.
    !!
    !! If `unit` is passed the unit will be in the requested unit.
    !! Otherwise the unit will be automatically determined [kB, MB, GB, TB]
    !! The format will be [es*.2] depending on the length of `mem_str`
    character(len=*), intent(inout) :: mem_str
    character(len=2), intent(in), optional :: unit

    character(len=2) :: lunit
    character(len=24) :: fmt

    real(dp) :: mem
    integer :: len_mem_str
    integer :: fw, fd, fe

    ! Determine length available
    fw = len(mem_str)

    ! Determine default unit
    lunit = 'kB'
    mem = this%kB
    if ( mem > 1024._dp ) then
      lunit = 'MB'
      mem = mem / 1024._dp
      if ( mem > 1024._dp ) then
        lunit = 'GB'
        mem = mem / 1024._dp
        if ( mem > 1024._dp ) then
          lunit = 'TB'
          mem = mem / 1024._dp
        end if
      end if
    end if
    
    if ( present(unit) ) then
      lunit = lcase(unit)
      select case ( lunit )
      case ( 'b' )
        lunit = 'B'
        mem = this%kB * 1024._dp
      case ( 'kb' )
        lunit = 'kB'
        mem = this%kB
      case ( 'mb' )
        lunit = 'MB'
        mem = this%kB / 1024._dp ** 2
      case ( 'gb' )
        lunit = 'GB'
        mem = this%kB / 1024._dp ** 3
      case ( 'tb' )
        lunit = 'TB'
        mem = this%kB / 1024._dp ** 4
      case ( 'pb' )
        lunit = 'PB'
        mem = this%kB / 1024._dp ** 5
      case default
        lunit = 'GB'
        mem = this%kB / 1024._dp ** 3
      end select
    end if

    ! Now determine format
    !  ENw.dEe
    ! w = 3 + [.] + d + [E] + [+-] + e + [ kB]

    if ( mem < 1000._dp ) then

      write(mem_str, '(f9.3,tr1,a)') mem, lunit
      ! Move to the left
      mem_str = adjustl(mem_str)

      return

    end if

    ! Subtract ' *B'
    fw = fw - 3
    
    fe = 1
    if ( abs(exponent(mem)) > 12 ) fe = 2

    ! Maximally 3 digits after '.'
    fd = min(fw - 3 - 1 - 1 - 1 - fe, 3)
    write(fmt,'(''('',a,i0,a,i0,a,i0,'',tr1,a)'')') 'EN',fw,'.',fd,'E',fe
    write(mem_str, trim(fmt)) mem, lunit
    ! Move to the left
    mem_str = adjustl(mem_str)

  end subroutine get_string_

end module byte_count_m
