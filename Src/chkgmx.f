      SUBROUTINE CHKGMX (K,BG,MESH,G2MAX)
C
C  Modules
C
      use precision, only: dp
      use sys,       only: die
      use parallel,  only : Node

      integer :: i, j, i1, i2, i3

      real(dp), parameter :: ZERO=0.0_dp,HALF=0.5_dp,
     $                       TINY=1.0e-8_dp,BIG=1.0e20_dp

      real(dp), intent(inout) :: g2max
      real(dp) K(3),BG(3,3),BM(3,3),G(3)
      real(dp) :: r, gmod, gmax, g2msh
      INTEGER MESH(3)

      real(dp) :: ddot
      external :: ddot

      DO I=1,3
        DO J=1,3
          BM(J,I)=BG(J,I)*MESH(I)
        ENDDO
      ENDDO
      CALL MINVEC (BM,BM)
      GMAX=BIG
      DO I3=-1,1
        DO I2=-1,1
          DO I1=-1,1
            IF (I1.EQ.0.AND.I2.EQ.0.AND.I3.EQ.0) GO TO 20
              G(1)=BM(1,1)*I1+BM(1,2)*I2+BM(1,3)*I3
              G(2)=BM(2,1)*I1+BM(2,2)*I2+BM(2,3)*I3
              G(3)=BM(3,1)*I1+BM(3,2)*I2+BM(3,3)*I3
              GMOD=SQRT(ddot(3,G,1,G,1))
              R=HALF*GMOD-ddot(3,K,1,G,1)/GMOD
              GMAX=MIN(GMAX,R)
  20        CONTINUE
          ENDDO
        ENDDO
      ENDDO
      IF (GMAX.LT.ZERO) call die('CHKGMX: K NOT IN FIRST BZ')
      G2MSH=GMAX*GMAX-TINY
      IF (G2MSH.LT.G2MAX) THEN
        if (Node.eq.0) then
*        WRITE(6,*) 'CHKGMX: MESH TOO SPARSE. GMAX HAS BEEN REDUCED'
*        WRITE(6,*) 'CHKGMX: OLD G2MAX =',G2MAX
*        WRITE(6,*) 'CHKGMX: NEW G2MAX =',G2MSH
        endif
        G2MAX=G2MSH
      ENDIF
      END
