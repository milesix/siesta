      SUBROUTINE RADFFT( L, NR, RMAX, F, G )
C *********************************************************************
C Makes a fast Fourier transform of a radial function.
C If function f is of the form
C   f(r_vec) = F(r_mod) * Ylm(theta,phi)
C where Ylm is a spherical harmonic with l = argument L, and
C argument F contains on input the real function F(r_mod), in a uniform
C radial grid:
C   r_mod = ir * RMAX / NR,  ir = 0,1,...,NR,
C and if g is the 3-dimensional Fourier transform of f:
C   g(k_vec) = 1/(2*pi)**(3/2) *
C              Integral( d3_r_vec * exp(-i * k_vec * r_vec) * f(r_vec) )
C then g has the form
C   g(k_vec) = (-i)**L * G(k_mod) * Ylm(theta,phi)
C where argument G contains on output the real function G(k_mod) in
C a uniform radial grid:
C   k_mod = ik * k_max / NR, ik = 0,1,...,NR,  k_max = NR*pi/RMAX
C Ref: J.M.Soler notes of 16/08/95.
C *************** INPUT ***********************************************
C INTEGER L       : Angular momentum quantum number
C INTEGER NR      : Number of radial intervals.
C                   2*NR must be an acceptable number of points for the
C                   FFT routine used.
C REAL*8  RMAX    : Maximum radius
C REAL*8  F(0:NR) : Function to be tranformed, in a radial mesh
C *************** OUTPUT **********************************************
C REAL*8  G(0:NR) : Fourier transform of F
C *************** UNITS ***********************************************
C Units of RMAX and F are arbitrary.
C Units of k_max and G are related with those of RMAX and F in the
C   obvious way (see above).
C *************** BEHAVIOUR *******************************************
C 1) F and G may be the same physical array, i.e. it is allowed:
C      CALL RADFFT( L, NR, RMAX, F, F )
C 2) It also works in the opposite direction, but then the factor
C    multiplying the output is (+i)**L. Thus, the following two calls
C      CALL RADFFT( L, NR, RMAX, F, G )
C      CALL RADFFT( L, NR, NR*PI/RMAX, G, H )
C    make H = F
C 3) If you will divide the output by q**l, truncation errors may be
C    quite large for small k's if L and NR are large. Therefore, these
C    components are calculated by direct integration rather than FFT.
C    Parameter ERRFFT is the typical truncation error in the FFT, and
C    controls which k's are integrated directly. A good value is 1e-8.
C    If you will not divide by k**l, make ERRFFT=1.e-30.
C *********************************************************************
C Written by J.M.Soler. August 1996.
C *********************************************************************

C Next line is non-standard and may be suppressed -------------------
      IMPLICIT NONE
C -------------------------------------------------------------------

C Declare argument types and dimensions -----------------------------
      INTEGER           L, NR
      DOUBLE PRECISION  F(0:NR), G(0:NR), RMAX
C -------------------------------------------------------------------

C ERRFFT is the typical truncation error in the FFT routine ---------
      DOUBLE PRECISION  ERRFFT
      PARAMETER ( ERRFFT = 1.D-8 )
C -------------------------------------------------------------------

C Internal variable types and dimensions ----------------------------

      INTEGER
     .  I, IQ, IR, JR, M, MQ, N, NQ
      DOUBLE PRECISION
     .  BESSPH, C, DQ, DR, FR, PI, R, RN, Q

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE ::
     .  GG

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE ::
     .  FN

      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, SAVE ::
     .  P

      EXTERNAL
*    .  BESSPH, CFT, CHKDIM, TIMER
     .  BESSPH, CHKDIM, TIMER, MEMORY
C -------------------------------------------------------------------

C Start time counter ------------------------------------------------
*     CALL TIMER( 'RADFFT', 1 )
C -------------------------------------------------------------------

C Allocate local memory ---------------------------------------------
      allocate(GG(0:2*NR))
      call memory('A','D',2*NR+1,'radfft')
      allocate(FN(2,0:2*NR))
      call memory('A','D',2*(2*NR+1),'radfft')
      allocate(P(2,0:l,0:l))
      call memory('A','D',2*(l+1)*(l+1),'radfft')

C Find some constants -----------------------------------------------
      PI = 4.D0 * ATAN( 1.D0 )
      NQ = NR
      DR = RMAX / NR
      DQ = PI / RMAX
      C = DR / SQRT( 2.D0*PI )
C -------------------------------------------------------------------


C Set up a complex polynomial such that the spherical Bessel function:
C   j_l(x) = Real( Sum_n( P(n,l) * x**n ) * exp(i*x) ) / x**(l+1)
      P(1,0,0) =  0.D0
      P(2,0,0) = -1.D0
      if (l.gt.0) then
        P(1,0,1) =  0.D0
        P(2,0,1) = -1.D0
        P(1,1,1) = -1.D0
        P(2,1,1) =  0.D0
        if (l.gt.1) then
          DO M = 2,L
          DO N = 0,M
          DO I = 1,2
            P(I,N,M) = 0.D0
            IF (N .LT. M) P(I,N,M) = P(I,N,M) + (2*M-1) * P(I,N,M-1)
            IF (N .GE. 2) P(I,N,M) = P(I,N,M) - P(I,N-2,M-2)
          ENDDO
          ENDDO
          ENDDO
        endif
      endif
C -------------------------------------------------------------------

C Initialize accumulation array -------------------------------------
      DO IQ = 0,NQ
        GG(IQ) = 0.D0
      ENDDO
C -------------------------------------------------------------------

C Iterate on terms of the j_l(q*r) polynomial -----------------------
      DO N = 0,L

C       Set up function to be fast fourier transformed
        FN(1,0) = 0.D0
        FN(2,0) = 0.D0
        DO JR = 1, 2*NR-1

          IF (JR .LE. NR) THEN
            IR = JR
            R = IR * DR
            FR = F(IR)
          ELSE
            IR = 2*NR - JR
            R = - (IR * DR)
            FR = F(IR) * (-1.D0)**L
          ENDIF

C         Find  r**2 * r**n / r**(l+1)
          RN = R**(N-L+1)

          FN(1,JR) = C * FR * RN * P(1,N,L)
          FN(2,JR) = C * FR * RN * P(2,N,L)
        ENDDO

C       Perform one-dimensional complex FFT
        CALL FOUR1( FN, 2*NR, +1 )

C       Accumulate contribution
        DO IQ = 1,NQ
          Q = IQ * DQ
          GG(IQ) = ( GG(IQ) + FN(1,IQ) ) / Q
        ENDDO

      ENDDO
C -------------------------------------------------------------------

C Special case for Q=0 ---------------------------------------------
      GG(0) = 0.D0
      IF ( L .EQ. 0 ) THEN
        DO IR = 1,NR
          R = IR * DR
          GG(0) = GG(0) + R*R * F(IR)
        ENDDO
        GG(0) = GG(0) * 2.D0 * C
      ENDIF
C -------------------------------------------------------------------

C Direct integration for the smallest Q's ---------------------------
      IF (L.EQ.0) THEN
        MQ = 0
      ELSE
        MQ = NQ * ERRFFT**(1.D0/L)
      ENDIF
      DO IQ = 1,MQ
        Q = IQ * DQ
        GG(IQ) = 0.D0
        DO IR = 1,NR
          R = IR * DR
          GG(IQ) = GG(IQ) + R*R * F(IR) * BESSPH(L,Q*R)
        ENDDO
        GG(IQ) = GG(IQ) * 2.D0 * C
      ENDDO
C -------------------------------------------------------------------

C Copy from local to output array -----------------------------------
      DO IQ = 0,NQ
        G(IQ) = GG(IQ)
      ENDDO
C -------------------------------------------------------------------

C Deallocate local memory -------------------------------------------
      call memory('D','D',size(GG),'radfft')
      deallocate(GG)
      call memory('D','D',size(FN),'radfft')
      deallocate(FN)
      call memory('D','D',size(P),'radfft')
      deallocate(P)

C Stop time counter ------------------------------------------------
*     CALL TIMER( 'RADFFT', 2 )
C -------------------------------------------------------------------

      END

