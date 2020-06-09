
    ! Coordinates in the P6_3/m space group (No. 176)
    ! There are 4 extra generalized coordinates with respect
    ! to the P6_3/mmc case:
    !
    ! The 12k wpos change into 12i wpos with 3 instead of 2 parameters:
    !     Old:  x1, z1;   x2, z2  ;  x3, z3
    !     New:  x1, y1, z1, with y1 close to 2*x1
    !           ... same for 2 and 3
    !
    ! The 12j wpos split into 2 6h wpos
    !     Old:  x, y
    !     New:  x_a, y_a : close to x, y
    !           x_b, y_b : close to -y, -x
    !
    ! The new coordinates x_b, y_b, and y1, y2, y3 are at the end of
    ! the q(:) array. x_a and y_a substitute x and y at the beginning.
    !
    ! The Nb wpos 2b are now 2a

    !--------------------------------------------------------------

    ! Convenient symbolic names for generalized coordinates
    !
    real(p), pointer        :: x1, y1, z1, x2, y2, z2, x3, y3, z3
    real(p), pointer        :: x_a, y_a, x_b, y_b

    ! work pointers to simplify formulas
    real(p), pointer        :: x, y, z

    ! Convenient symbolic names for Wyckoff orbits
    !
    real(p), pointer, dimension(:,:)  :: Nb_2a, Nb_2c
    real(p), pointer, dimension(:,:)  :: Nb_2d, Nb_6h1, Nb_6h2
    real(p), pointer, dimension(:,:)  :: Se_12i1, Se_12i2, Se_12i3

    ! work pointer
    real(p), pointer, dimension(:,:)  :: Se_12i

    !  This has to be specified here

    natoms = 54
    if (size(xa,dim=2) /= natoms) call die("natoms, na mismatch")

    ! Correspondence of symbols to positions in the q array

    x_a => q(1)    ! Same position, very close values
    y_a => q(2)    ! to the old q(1) and q(2)
    x1 => q(3)
    z1 => q(4)
    x2 => q(5)
    z2 => q(6)
    x3 => q(7)
    z3 => q(8)
    !              Extra coordinates
    x_b => q(9)
    y_b => q(10)
    y1 => q(11)
    y2 => q(12)
    y3 => q(13)


    if (present(spec_no)) then
       spec_no(1:18)  = 1
       spec_no(19:54)  = 2
    endif

    ! Wyckoff-orbit blocks are convenient, but care should be taken
    ! to double-check the numbering!!

    Nb_2a => xa(:,1:2)
    Nb_2c => xa(:,3:4)
    Nb_2d => xa(:,5:6)
!!    Nb_12j => xa(:,7:18)   ! these are split into:
    Nb_6h1 => xa(:,7:12)
    Nb_6h2 => xa(:,13:18)

    Se_12i1 => xa(:,19:30)
    Se_12i2 => xa(:,31:42)
    Se_12i3 => xa(:,43:54)

    !
    !  Now, the coordinates
    !  We can use all of Fortran's operators, so this is a more
    !  readable format. 

    !  For convenience, we could group the blocks by species or
    !  Wyckoff-position symbol, and/or use block indices or pointers
    !
    !  Some atoms might not use any q's
    !  Note change in notation: Nb_2b -> Nb_2a

    Nb_2a(:,1)   = (/  0.0_dp , 0.0_dp, 1.0_dp/4 /)
    Nb_2a(:,2)   = (/  0.0_dp , 0.0_dp, 3.0_dp/4 /)
    !
    Nb_2c(:,1)   = (/  1.0_dp/3 , 2.0_dp/3, 1.0_dp/4 /)
    Nb_2c(:,2)   = (/  2.0_dp/3 , 1.0_dp/3, 3.0_dp/4 /)
    !
    Nb_2d(:,1)   = (/  1.0_dp/3 , 2.0_dp/3, 3.0_dp/4 /)
    Nb_2d(:,2)   = (/  2.0_dp/3 , 1.0_dp/3, 1.0_dp/4 /)
    !
    x => x_a
    y => y_a
    Nb_6h1(:,1)   = (/  x       ,    y    , 1.0_dp/4 /)
    Nb_6h1(:,2)   = (/ -y       ,  x-y    , 1.0_dp/4 /)
    Nb_6h1(:,3)   = (/ -x+y     ,  -x     , 1.0_dp/4 /)
    Nb_6h1(:,4)   = (/ -x       ,  -y     , 3.0_dp/4 /)
    Nb_6h1(:,5)   = (/  y       ,  -x+y   , 3.0_dp/4 /)
    Nb_6h1(:,6)   = (/ x-y      ,   x     , 3.0_dp/4 /)

    x => x_b   ! Note: when re-using old relaxation results,
    y => y_b   ! x_b should be close to -y_old, and y_b
               ! close to -x_old
 
    Nb_6h2(:,1)   = (/  x       ,    y    , 1.0_dp/4 /)
    Nb_6h2(:,2)   = (/ -y       ,  x-y    , 1.0_dp/4 /)
    Nb_6h2(:,3)   = (/ -x+y     ,  -x     , 1.0_dp/4 /)
    Nb_6h2(:,4)   = (/ -x       ,  -y     , 3.0_dp/4 /)
    Nb_6h2(:,5)   = (/  y       ,  -x+y   , 3.0_dp/4 /)
    Nb_6h2(:,6)   = (/ x-y      ,   x     , 3.0_dp/4 /)
    !
    !-----------------------------------
    x => x1
    y => y1
    z => z1
    Se_12i => Se_12i1
    !
    Se_12i(:,1)  = (/  x      ,   y   ,   z  /)
    Se_12i(:,2)  = (/  -y    ,   x-y    ,   z  /)
    Se_12i(:,3)  = (/ -x+y   ,  -x   ,   z  /)
    Se_12i(:,4)  = (/  -x      ,   -y  ,   z+0.5_dp  /)
    Se_12i(:,5)  = (/  y     ,   -x+y  ,   z+0.5_dp  /)
    Se_12i(:,6)  = (/  x-y   ,    x    ,   z+0.5_dp  /)
    Se_12i(:,7)  = (/  -x    ,   -y    ,  -z  /)
    Se_12i(:,8)  = (/  y     ,   -x+y  ,  -z  /)
    Se_12i(:,9)  = (/  x-y   ,    x    ,  -z  /)
    Se_12i(:,10) = (/  x     ,    y    ,  -z+0.5_dp  /)
    Se_12i(:,11) = (/  -y    ,   x-y   ,  -z+0.5_dp  /)
    Se_12i(:,12) = (/  -x+y  ,   -x    ,  -z+0.5_dp  /)
    !
    !-----------------------------------
    x => x2
    y => y2
    z => z2
    Se_12i => Se_12i2
    !
    Se_12i(:,1)  = (/  x      ,   y   ,   z  /)
    Se_12i(:,2)  = (/  -y    ,   x-y    ,   z  /)
    Se_12i(:,3)  = (/ -x+y   ,  -x   ,   z  /)
    Se_12i(:,4)  = (/  -x      ,   -y  ,   z+0.5_dp  /)
    Se_12i(:,5)  = (/  y     ,   -x+y  ,   z+0.5_dp  /)
    Se_12i(:,6)  = (/  x-y   ,    x    ,   z+0.5_dp  /)
    Se_12i(:,7)  = (/  -x    ,   -y    ,  -z  /)
    Se_12i(:,8)  = (/  y     ,   -x+y  ,  -z  /)
    Se_12i(:,9)  = (/  x-y   ,    x    ,  -z  /)
    Se_12i(:,10) = (/  x     ,    y    ,  -z+0.5_dp  /)
    Se_12i(:,11) = (/  -y    ,   x-y   ,  -z+0.5_dp  /)
    Se_12i(:,12) = (/  -x+y  ,   -x    ,  -z+0.5_dp  /)
    !
    !-----------------------------------
    x => x3
    y => y3
    z => z3
    Se_12i => Se_12i3
    !
    Se_12i(:,1)  = (/  x      ,   y   ,   z  /)
    Se_12i(:,2)  = (/  -y    ,   x-y    ,   z  /)
    Se_12i(:,3)  = (/ -x+y   ,  -x   ,   z  /)
    Se_12i(:,4)  = (/  -x      ,   -y  ,   z+0.5_dp  /)
    Se_12i(:,5)  = (/  y     ,   -x+y  ,   z+0.5_dp  /)
    Se_12i(:,6)  = (/  x-y   ,    x    ,   z+0.5_dp  /)
    Se_12i(:,7)  = (/  -x    ,   -y    ,  -z  /)
    Se_12i(:,8)  = (/  y     ,   -x+y  ,  -z  /)
    Se_12i(:,9)  = (/  x-y   ,    x    ,  -z  /)
    Se_12i(:,10) = (/  x     ,    y    ,  -z+0.5_dp  /)
    Se_12i(:,11) = (/  -y    ,   x-y   ,  -z+0.5_dp  /)
    Se_12i(:,12) = (/  -x+y  ,   -x    ,  -z+0.5_dp  /)
    !
