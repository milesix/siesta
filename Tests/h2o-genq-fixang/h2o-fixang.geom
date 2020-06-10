!---------------------------
!
!   Geometry specification
!   Water molecule with a fixed (wrong) angle

    ! Convenient symbolic names for generalized coordinates
    !
    real(p)                 :: alpha = 2.00
    real(p), pointer        :: doh

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
    !  Some atoms might not use any q's
    !
    doh => q(1)
    !
    H_at(:,1)   = (/  doh , 0.0_dp, 0.0_dp /)
    H_at(:,2)   = (/  doh*cos(alpha) , doh*sin(alpha), 0.0_dp /)
    O_at(:,1)   = (/  0.0_dp , 0.0_dp, 0.0_dp /)
    !
    !-----------------------------------

