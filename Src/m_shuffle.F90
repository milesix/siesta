module m_shuffle
!
! Implements the Fisher-Yates shuffling algorithm
! Ref: Wikipedia
! Note that we trust the Fortran intrinsic random number
! generator, which produces numbers in the [0,1) interval,
! to generate the random integers needed.
!
! Alberto Garcia, December 2009

implicit none

public :: shuffle
private


CONTAINS

subroutine shuffle(a)
integer, intent(inout) :: a(:)

integer n, k, itemp

n = size(a)
do while (n > 1)
  k = random_int(n)
  itemp = a(k)
  a(k) = a(n)
  a(n) = itemp
  n = n-1
enddo

end subroutine shuffle
!-----------------------------------------

function random_int(n) result(k)
integer, intent(in) :: n
integer             :: k

real      :: x
intrinsic :: random_number
call random_number(x)
k = int(n*x) + 1
end function random_int
end module m_shuffle
!-----------------------------------------

#ifdef UNIT_TEST
program t

use m_shuffle

implicit none
integer, parameter :: N = 8
integer :: a(N)
integer :: i

do i = 1, N
   a(i) = i
enddo
print "(8i3)", a

do 
   call shuffle(a)
   print "(8i3)", a
   read *
enddo
end program t
#endif /* UNIT_TEST */


