!---------------------------
!
!   Geometry specification
!   Water molecule

    ! Convenient symbolic names for generalized coordinates
    !
    real(p), pointer        :: doh, alpha

    ! Convenient symbolic names for Wyckoff orbits
    !
    real(p), pointer, dimension(:,:)  :: H_at, O_at

    !  This has to be specified here

    natoms = 3
    if (size(xa,dim=2) /= natoms) call die("natoms, na mismatch")

    if (present(spec_no)) then
       spec_no(1)  = 1          ! O
       spec_no(2:3)  = 2        ! H
    endif

    ! Wyckoff-orbit blocks are convenient, but care should be taken
    ! to double-check the numbering!!

    O_at => xa(:,1:1)
    H_at => xa(:,2:3)

    !
    !  Now, the coordinates
    !  We can use all of Fortran's operators, so this is a more
    !  readable format. 

    !  For convenience, we could group the blocks by species or
    !  Wyckoff-position symbol, and/or use block indices or pointers
    !
    !  Some atoms might not use any q's
    !
    doh => q(1)
    alpha => q(2)
    !
    H_at(:,1)   = (/  doh , 0.0_dp, 0.0_dp /)
    H_at(:,2)   = (/  doh*cos(alpha) , doh*sin(alpha), 0.0_dp /)
    O_at(:,1)   = (/  0.0_dp , 0.0_dp, 0.0_dp /)
    !
    !-----------------------------------

