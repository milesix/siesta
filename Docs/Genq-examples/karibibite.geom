
!
real(p), pointer, dimension(:) :: c_As_1_8d
real(p), pointer, dimension(:) :: c_Fe_1_4c
real(p), pointer, dimension(:) :: c_As_2_8d
real(p), pointer, dimension(:) :: c_Fe_2_8d
real(p), pointer, dimension(:) :: c_As_3_8d
real(p), pointer, dimension(:) :: c_O_1_8d 
real(p), pointer, dimension(:) :: c_O_2_8d 
real(p), pointer, dimension(:) :: c_O_3_8d 
real(p), pointer, dimension(:) :: c_O_4_8d 
real(p), pointer, dimension(:) :: c_O_5_8d 
real(p), pointer, dimension(:) :: c_O_6_4c 
real(p), pointer, dimension(:) :: c_O_7_4c 
real(p), pointer, dimension(:) :: c_O_8_8d 
real(p), pointer, dimension(:) :: c_H_1_4c 

! Convenient symbolic names for Wyckoff orbits
!
real(p), pointer, dimension(:,:)  :: As_1_8d, As_2_8d, As_3_8d
real(p), pointer, dimension(:,:)  :: Fe_1_4c, Fe_2_8d
real(p), pointer, dimension(:,:)  :: O_1_8d, O_2_8d, O_3_8d
real(p), pointer, dimension(:,:)  :: O_4_8d, O_5_8d, O_8_8d
real(p), pointer, dimension(:,:)  :: O_6_4c, O_7_4c
real(p), pointer, dimension(:,:)  :: H_1_4c

!  This has to be specified here

natoms = 96
if (size(xa,dim=2) /= natoms) call die("natoms, na mismatch")

c_As_1_8d => q(1:3)
c_Fe_1_4c => q(4:5)
c_As_2_8d => q(6:8)
c_Fe_2_8d => q(9:11)
c_As_3_8d => q(12:14)
c_O_1_8d => q(15:17)
c_O_2_8d => q(18:20)
c_O_3_8d => q(21:23)
c_O_4_8d => q(24:26)
c_O_5_8d => q(27:29)
c_O_6_4c => q(30:31)
c_O_7_4c => q(32:33)
c_O_8_8d => q(34:36)
c_H_1_4c => q(37:38)

if (present(spec_no)) then
   spec_no(1:8)  = 1       ! As
   spec_no(9:12)  = 2      ! Fe
   spec_no(13:20)  = 1
   spec_no(21:28)  = 2
   spec_no(29:36)  = 1
   spec_no(37:92)  = 3     ! O
   spec_no(93:96)  = 4     ! H
endif

! Wyckoff-orbit blocks are convenient, but care should be taken
! to double-check the numbering!!

As_1_8d => xa(:,1:8)
Fe_1_4c => xa(:,9:12)
As_2_8d => xa(:,13:20)
Fe_2_8d => xa(:,21:28)
As_3_8d => xa(:,29:36)
O_1_8d => xa(:,37:44)
O_2_8d => xa(:,45:52)
O_3_8d => xa(:,53:60)
O_4_8d => xa(:,61:68)
O_5_8d => xa(:,69:76)
O_6_4c => xa(:,77:80)
O_7_4c => xa(:,81:84)
O_8_8d => xa(:,85:92)
H_1_4c => xa(:,93:96)
!
!  Now, the coordinates
!  We can use all of Fortran's operators, so this is a more
!  readable format. 

!  For convenience, we could group the blocks by species or
!  Wyckoff-position symbol, and/or use block indices or pointers
!
call generate_8d_orbit(c_As_1_8d,As_1_8d)
call generate_4c_orbit(c_Fe_1_4c,Fe_1_4c,0.25_p)
call generate_8d_orbit(c_As_2_8d,As_2_8d)
call generate_8d_orbit(c_Fe_2_8d,Fe_2_8d)
call generate_8d_orbit(c_As_3_8d,As_3_8d)
call generate_8d_orbit(c_O_1_8d,O_1_8d)
call generate_8d_orbit(c_O_2_8d,O_2_8d)
call generate_8d_orbit(c_O_3_8d,O_3_8d)
call generate_8d_orbit(c_O_4_8d,O_4_8d)
call generate_8d_orbit(c_O_5_8d,O_5_8d)
call generate_4c_orbit(c_O_6_4c,O_6_4c,0.75_p)
call generate_4c_orbit(c_O_7_4c,O_7_4c,0.75_p)
call generate_8d_orbit(c_O_8_8d,O_8_8d)
call generate_4c_orbit(c_H_1_4c,H_1_4c,0.75_p)


CONTAINS

subroutine generate_8d_orbit(xyz,xa)
  real(p), intent(in) :: xyz(3)
  real(p), intent(out) :: xa(3,8)

  real(p) :: x,y,z
  real(p), parameter :: half = 0.5_p
  
x = xyz(1) 
y = xyz(2)
z = xyz(3)

xa(:,1)   = (/  x       ,    y    , z        /)
xa(:,2)   = (/ -x       ,   -y    , -z       /)
xa(:,3)   = (/ -x+half  ,    y    , z        /)
xa(:,4)   = (/ x+half   ,  -y     , -z       /)
xa(:,5)   = (/ x        ,  -y+half, z+half       /)
xa(:,6)   = (/ -x       ,  y+half, -z+half       /)
xa(:,7)   = (/ -x+half  ,  -y+half, z+half       /)
xa(:,8)   = (/ x+half  ,  y+half, -z+half       /)
end subroutine generate_8d_orbit
!
subroutine generate_4c_orbit(yz,xa,x)
  real(p), intent(in) :: yz(2)
  real(p), intent(in) :: x
  real(p), intent(out) :: xa(3,4)

  real(p) :: y, z
  real(p), parameter :: half = 0.5_p
  
y = yz(1)
z = yz(2)

xa(:,1)   = (/  x       ,    y    , z        /)
xa(:,2)   = (/ -x       ,  y+half, -z+half       /)
xa(:,3)   = (/ x + half   ,  -y     , -z       /)
xa(:,4)   = (/ -x + half  ,  -y+half, z+half       /)
end subroutine generate_4c_orbit
!
!-----------------------------------
