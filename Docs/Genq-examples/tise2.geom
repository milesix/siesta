!---------------------------
!
!   Geometry specification
!   TiSe2, non-distorted structure

    ! Convenient symbolic names for generalized coordinates
    !
    real(p)                   :: xTi, xSe, ySe, zSe
    real(p), pointer          ::  delta1, delta2

    ! Convenient symbolic names for Wyckoff orbits
    !
    real(p), pointer, dimension(:,:)  :: Ti_1a
    real(p), pointer, dimension(:,:)  :: Ti_3e
    real(p), pointer, dimension(:,:)  :: Se_2d    
    real(p), pointer, dimension(:,:)  :: Se_6g

    !  This has to be specified here
    natoms = 12
    if (size(xa,dim=2) /= natoms) call die("natoms, na mismatch")

    if (present(spec_no)) then
       spec_no(1:4)   = 1
       spec_no(5:12)  = 2
    endif

    ! Wyckoff-orbit blocks are convenient, but care should be taken
    ! to double-check the numbering!!

    Ti_1a => xa(:,1:1)
    Ti_3e => xa(:,2:4)
    Se_2d => xa(:,5:6)
    Se_6g => xa(:,7:12)

    !
    !  Now, the coordinates
    !  We can use all of Fortran's operators, so this is a more
    !  readable format. 

    !  For convenience, we could group the blocks by species or
    !  Wyckoff-position symbol, and/or use block indices or pointers
    !
    !  Some atoms might not use any q's
    !

    Ti_1a(:,1)   = (/  0.0_dp , 0.0_dp, 0.0_dp /)  

    !-----------------------------------

    delta1 => q(1)
    delta2 => q(2)

    xTi  = 1.0_dp/2.0_dp
 
    xSe  = 1.0_dp/6.0_dp
    ySe  = 5.0_dp/6.0_dp
    zSe  = 0.030402381_dp

    Ti_3e(:,1)  = (/   xTi + delta1     ,   0.0_dp             ,   0.0_dp  /)
    Ti_3e(:,2)  = (/  0.0_dp            ,    xTi + delta1      ,   0.0_dp  /)
    Ti_3e(:,3)  = (/  -xTi - delta1     ,   -xTi - delta1      ,   0.0_dp  /)
    
    Se_2d(:,1)  = (/  1.0_dp/3.0_dp , 2.0_dp/3.0_dp ,   zSe  /)
    Se_2d(:,2)  = (/  2.0_dp/3.0_dp , 1.0_dp/3.0_dp ,  -zSe  /)
    
    Se_6g(:,1)  = (/  xSe  - delta2    ,   ySe  - delta2  ,   zSe  /)
    Se_6g(:,2)  = (/ -ySe  + delta2    ,   xSe  - ySe     ,   zSe  /)
    Se_6g(:,3)  = (/ -xSe  + ySe       ,  -xSe  + delta2  ,   zSe  /)
    Se_6g(:,4)  = (/  ySe  - delta2    ,   xSe  - delta2  ,  -zSe  /)
    Se_6g(:,5)  = (/  xSe  - ySe       ,  -ySe  + delta2  ,  -zSe  /)
    Se_6g(:,6)  = (/ -xSe  + delta2    ,  -xSe  + ySe     ,  -zSe  /)
