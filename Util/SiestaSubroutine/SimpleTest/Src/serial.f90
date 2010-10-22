program simple

! A very simple driver for Siesta-as-subroutine (or siesta-as-server)
! This version launches a serial version of siesta.

  use fsiesta

  implicit none
  integer,parameter:: dp = kind(1.d0)

  integer,parameter :: na = 3
  integer :: ia
  real(dp):: e, fa(3,na), xa(3,na)

  data xa / 0.0, 0.0, 0.0, &
            0.7, 0.7, 0.0, &
           -0.7, 0.7, 0.0 /

  call siesta_units( 'Ang', 'eV' )

  call siesta_launch( 'h2o' )
  print*, 'siesta launched'

  call siesta_forces( 'h2o', na, xa, energy=e, fa=fa )
  print'(a,/,(3f12.6,3x,3f12.6))', 'xa, fa =', (xa(:,ia),fa(:,ia),ia=1,na)

  xa(1,1) = 0.8
  call siesta_forces( 'h2o', na, xa, energy=e, fa=fa )
  print'(a,/,(3f12.6,3x,3f12.6))', 'xa, fa =', (xa(:,ia),fa(:,ia),ia=1,na)

  call siesta_quit( 'h2o' )

end program simple
