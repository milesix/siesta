module spher_harm2
use atmfuncs, only : func_t
use spher_harm, only: lofilm, gauleg, rlylm
use precision
Contains
SUBROUTINE YLMYLM( ILM1,func, ILM2, R, YY, DYYDR )
    INTEGER, intent(in)       ::       ILM1, ILM2
    type(func_t), intent(in)  :: func  ! Function label
    REAL(DP), intent(in)      ::       R(3)
    REAL(DP), intent(out)     ::       YY, DYYDR(3)

! Returns the product of two real spherical harmonics (SH),
!  times r**l each, and its gradient.
! *********** INPUT ***************************************************
! INTEGER ILM1, ILM2 : Combined SH indexes. ILM = L*L+L+M+1
! REAL*8  R(3)       : Vector towards the (theta,phi) direction.
! *********** OUTPUT **************************************************
! REAL*8  YY       : Product of the two SH.
! REAL*8  DYYDR(3) : Derivative (gradient) of YY with respect to R
! *********************************************************************
! Written by J.M.Soler. Feb' 96.
! *********************************************************************


    INTEGER MAXLM
    INTEGER           I, L
    
    real(dp), ALLOCATABLE, SAVE ::   Y(:)
    real(dp), ALLOCATABLE, SAVE ::   DYDR(:,:)

    EXTERNAL          MEMORY

    L = MAX( LOFILM(ILM1), LOFILM(ILM2) )
    MAXLM = (L+1)*(L+1)

    allocate(Y(0:MAXLM))
    call memory('A','D',MAXLM+1,'ylmylm')
    allocate(DYDR(3,0:MAXLM))
    call memory('A','D',3*MAXLM+3,'ylmylm')

    CALL RLYLM( L, R, Y, DYDR )

    YY = Y(ILM1-1) * Y(ILM2-1)
    do I = 1,3
       DYYDR(I) = DYDR(I,ILM1-1)*Y(ILM2-1) + Y(ILM1-1)*DYDR(I,ILM2-1)
    enddo

    call memory('D','D',size(DYDR),'ylmylm')
    deallocate(DYDR)
    call memory('D','D',size(Y),'ylmylm')
    deallocate(Y)

  END subroutine ylmylm

  SUBROUTINE YLMEXP( LMAX, RLYLM_F, FUNC, IS, Ftype,IO, IR1, NR, &
       RMAX, NY, ILM, FLM )
    ! Makes a radial times spherical-harmonic expansion of a function.

    integer, intent(in)          :: lmax
    interface
       subroutine rlylm_f(lmax,rvec,y,grady)
         integer, intent(in)   :: lmax
         real(selected_real_kind(14)), intent(in)  :: rvec(3)
         real(selected_real_kind(14)), intent(out) :: y(0:)
         real(selected_real_kind(14)), intent(out) :: grady(1:,0:)
       end subroutine rlylm_f

       subroutine func(i1,ftype,i2,rvec,y,grady)
         use atmfuncs, only : func_t
         integer, intent(in)   :: i1, i2
         type(func_t), intent(in) :: ftype  ! Function label
         real(selected_real_kind(14)), intent(in)  :: rvec(3)
         real(selected_real_kind(14)), intent(out) :: y, grady(3)
       end subroutine func
    end interface

    integer, intent(in)        :: is, io
    type(func_t), intent(in) :: ftype  ! Function label
    integer, intent(in)        :: ir1, nr
    real(dp), intent(in)       :: rmax

    integer, intent(out)       :: ny
    integer, intent(out)       :: ilm(:)
    real(dp), dimension(ir1:nr,*), intent(out)  :: flm

    ! Written by J.M.Soler. September 1995.
    ! ************************* INPUT ***********************************
    ! INTEGER  LMAX                     : Max. ang. momentum quantum number
    ! EXTERNAL RLYLM_F(Lmax,Rvec,Y,GradY) : Returns spherical harmonics,
    !                                     times r**l
    ! EXTERNAL FUNC(IS,IO,Rvec,F,GradF) : Function to be expanded.
    ! INTEGER  IS, IO                   : Indexes to call FUNC.
    ! INTEGER  IR1                      : First radial index.
    ! INTEGER  NR                       : Last radial index.
    ! REAL*8   RMAX                     : Maximum radius: R(IR)=RMAX*IR/NR
    ! ************************* OUTPUT **********************************
    ! INTEGER  NY            : Number of spherical harmonics required.
    ! INTEGER  ILM(NY)       : Spherical harmonics indexes: ILM=L*L+L+M+1.
    ! REAL*8   FLM(IR1:NR,NY): Radial expansion for each spherical harmonic:
    !   F(Rvec) = Sum_iy( FLM(Rmod,iy) * Y(Rdir,ILM(iy)) ),
    !   with   Rmod = RMAX*IR/NR
    ! ************************* UNITS ***********************************
    ! RMAX must be in the same units of the argument Rvec in FUN!.
    ! FLM is in the same units of the argument F in FUNC.
    ! ************************* BEHAVIOUR *******************************
    ! If function FUNC contains angular components with L higher than LMAX,
    !   they will corrupt those computed with lower L.
    ! *******************************************************************


    ! Tolerance for FLM -------------------------------------------------
    REAL(DP) FTOL
    PARAMETER ( FTOL = 1.e-12_dp )
    ! -------------------------------------------------------------------

    ! Dimension parameters for internal variables -----------------------
    INTEGER MAXL, MAXLM, MAXSP
    PARAMETER ( MAXL  = 8 )
    PARAMETER ( MAXLM = (MAXL+1) * (MAXL+1) )
    PARAMETER ( MAXSP = (MAXL+1) * (2*MAXL+1) )
    ! -------------------------------------------------------------------

    ! Declare internal variables ----------------------------------------
    INTEGER IM, IR, ISP, IZ, JLM, JR, JY, NLM, NSP
    REAL(DP) DFDX(3), DDOT, DYDX(3,MAXLM), F(MAXSP), FY, PHI, PI, R, &
         RX(3), THETA, W(MAXL+1), WSP, X(3,MAXSP), Y(MAXLM), &
         YSP(MAXSP,MAXLM), Z(MAXL+1)
    EXTERNAL CHKDIM !, DDOT
    !     EXTERNAL TIMER
    ! -------------------------------------------------------------------

    ! Start time counter ------------------------------------------------
    !    CALL TIMER( 'YLMEXP', 1 )
    ! -------------------------------------------------------------------

    ! Check maximum angular momentum ------------------------------------
    CALL CHKDIM( 'YLMEXP', 'MAXL', MAXL, LMAX, 1 )
    ! -------------------------------------------------------------------

    ! Find special points and weights for gaussian quadrature -----------
    CALL GAULEG( -1._dp, 1._dp, Z, W, LMAX+1 )
    ! -------------------------------------------------------------------

    ! Find weighted spherical harmonics at special points ---------------
    PI = 4._dp * ATAN(1._dp)
    NLM = (LMAX+1)**2
    NSP = 0
    DO  IZ = 1,LMAX+1
       THETA = ACOS( Z(IZ) )
       DO  IM = 0,2*LMAX
          NSP = NSP + 1
          PHI = IM * 2._dp * PI / (2*LMAX+1)
          X(1,NSP) = SIN(THETA) * COS(PHI)
          X(2,NSP) = SIN(THETA) * SIN(PHI)
          X(3,NSP) = COS(THETA)
          CALL RLYLM_F( LMAX, X(1,NSP), Y, DYDX )
          WSP = W(IZ) * (2._dp*PI) / (2*LMAX+1)
          DO JLM = 1,NLM
             YSP(NSP,JLM) = WSP * Y(JLM)
          ENDDO
       ENDDO
    ENDDO
    ! -------------------------------------------------------------------
    ILM = 0

    ! Expand FUNC in spherical harmonics at each radius -----------------
    NY = 0

    !     Loop on radial points
    DO IR = IR1,NR
       R = IR * RMAX / NR

       !       Find function at points on a sphere of radius R
       DO ISP = 1,NSP
          RX(1) = R * X(1,ISP)
          RX(2) = R * X(2,ISP)
          RX(3) = R * X(3,ISP)
          CALL FUNC( IS, FTYPE, IO, RX, F(ISP), DFDX )
       ENDDO

       !       Expand F(R) in spherical harmonics
       DO  JLM = 1,NLM
          FY = DDOT(NSP,F,1,YSP(1,JLM),1)
          IF ( ABS(FY) .GT. FTOL ) THEN
             
             !           Find JY corresponding to JLM
             DO JY = 1,NY
                IF ( ILM(JY) .EQ. JLM ) THEN
                   FLM(IR,JY) = FY
                   GOTO 70
                ENDIF
             ENDDO

             !           New spherical harmonic required.
             NY = NY + 1
             ILM(NY) = JLM
             DO  JR = IR1,NR
                FLM(JR,NY) = 0._dp
             ENDDO
             FLM(IR,NY) = FY
          ENDIF
70     ENDDo

    ENDDO
    ! -------------------------------------------------------------------
       
    ! Special case for zero function ------------------------------------
    IF (NY .EQ. 0) THEN
       NY = 1
       ILM(1) = 1
       DO  IR = IR1,NR
          FLM(IR,1) = 0._dp
       ENDDO
    ENDIF
    ! -------------------------------------------------------------------
    
    ! Stop time counter -------------------------------------------------
    !     CALL TIMER( 'YLMEXP', 2 )
    ! -------------------------------------------------------------------

  end subroutine ylmexp
end module spher_harm2
