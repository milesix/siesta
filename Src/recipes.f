      module m_recipes
!**********************************************************************
!    This file contains routines adapted from 'Numerical Recipes, 
!    The Art of Scientific Computing' by W.H. Press, S.A. Teukolsky, 
!    W.T. Veterling and B.P. Flannery, Cambridge U.P. 1987 and 1992.
!**********************************************************************
! The routines contained in this file are:
!     SUBROUTINE POLINT
!     SUBROUTINE SPLINE
!     SUBROUTINE SPLINT
!     SUBROUTINE SORT
!**********************************************************************

      use sys, only: die

      public :: polint, spline, splint
      public :: sort

      contains

      SUBROUTINE POLINT(XA,YA,N,X,Y,DY) 
!*****************************************************************
! Polinomic interpolation. Modified and adapted to double 
! precision from same routine of Numerical Recipes.
! D. Sanchez-Portal, Oct. 1996
!*****************************************************************
! Input:
!   real*8  XA(N) : x values of the function y(x) to interpolate
!   real*8  YA(N) : y values of the function y(x) to interpolate
!   integer N     : Number of data points
!   real*8  X     : x value at which the interpolation is desired
! Output:
!   real*8  Y     : interpolated value of y(x) at X
!   real*8  DY    : accuracy estimate
!*****************************************************************

      IMPLICIT NONE
      INTEGER          :: N
      DOUBLE PRECISION :: XA(N),YA(N), X, Y, DY

      INTEGER          :: I, M, NS
      DOUBLE PRECISION :: C(N), D(N), DEN, DIF, DIFT, HO, HP, W
      DOUBLE PRECISION, PARAMETER :: ZERO=0.D0

      NS=1
      DIF=ABS(X-XA(1))
      DO I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
      END DO ! I
      Y=YA(NS)
      NS=NS-1
      DO M=1,N-1
        DO I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF (DEN.EQ.ZERO) call die('polint: ERROR. Two XAs are equal')
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
        END DO ! I
        IF (2*NS.LT.N-M) THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
      END DO ! M

      END SUBROUTINE POLINT



      SUBROUTINE SPLINE(DX,Y,N,YP1,YPN,Y2) 
!*********************************************************** 
! Cubic Spline Interpolation. Adapted from Numerical Recipes 
! routine of same name for a uniform grid and double precision
! D. Sanchez-Portal, Oct. 1996.
! Input:
!   real*8  DX   : x interval between data points
!   real*8  Y(N) : value of y(x) at data points
!   integer N    : number of data points
!   real*8  YP1  : value of dy/dx at X1 (first point)
!   real*8  YPN  : value of dy/dx at XN (last point)
! Output:
!   real*8  Y2(N): array to be used by routine SPLINT
! Behavior:
! - If YP1 or YPN are larger than 1E30, the natural spline
!   condition (d2y/dx2=0) at the corresponding edge point.
!************************************************************

      IMPLICIT NONE
      INTEGER          :: N
      DOUBLE PRECISION :: DX, Y(N), YP1, YPN, Y2(N)

      INTEGER          :: I, K
      DOUBLE PRECISION :: QN, P, SIG, U(N), UN
      DOUBLE PRECISION, PARAMETER :: YPMAX=0.99D30, 
     .  HALF=0.5D0, ONE=1.D0, THREE=3.D0, TWO=2.D0, ZERO=0.D0
    
      IF (YP1.GT.YPMAX) THEN
        Y2(1)=ZERO
        U(1)=ZERO
      ELSE
        Y2(1)=-HALF
        U(1)=(THREE/DX)*((Y(2)-Y(1))/DX-YP1)
      ENDIF
      DO I=2,N-1
        SIG=HALF
        P=SIG*Y2(I-1)+TWO
        Y2(I)=(SIG-ONE)/P
        U(I)=(THREE*( Y(I+1)+Y(I-1)-TWO*Y(I) )/(DX*DX)
     .       -SIG*U(I-1))/P
      END DO ! I
      IF (YPN.GT.YPMAX) THEN
        QN=ZERO
        UN=ZERO
      ELSE
        QN=HALF
        UN=(THREE/DX)*(YPN-(Y(N)-Y(N-1))/DX)
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+ONE)
      DO K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
      END DO ! K

      END SUBROUTINE SPLINE



      SUBROUTINE SPLINT(DX,YA,Y2A,N,X,Y,DYDX) 
!***************************************************************
! Cubic Spline Interpolation. Adapted from Numerical Recipes 
! routine of same name for a uniform grid, double precision,
! and to return the function derivative in addition to its value
! D. Sanchez-Portal, Oct. 1996.
! Input:
!   real*8  DX    : x interval between data points
!   real*8  YA(N) : value of y(x) at data points
!   real*8  Y2A(N): array returned by routine SPLINE
!   integer N     : number of data points
!   real*8  X     : point at which interpolation is desired
!   real*8  Y     : interpolated value of y(x) at point X
!   real*8  DYDX  : interpolated value of dy/dx at point X
!***************************************************************

      IMPLICIT NONE
      INTEGER          :: N
      DOUBLE PRECISION :: DX, YA(N), Y2A(N), X, Y, DYDX

      INTEGER          :: NHI, NLO
      DOUBLE PRECISION :: A, B
      DOUBLE PRECISION, PARAMETER ::
     .    ONE=1.D0, THREE=3.D0, SIX=6.D0, ZERO=0.D0

      IF (DX.EQ.ZERO) call die('splint: ERROR: DX=0')
      NLO=INT(X/DX)+1
      NHI=NLO+1
      A=NHI-X/DX-1
      B=ONE-A
      Y=A*YA(NLO)+B*YA(NHI)+
     .  ((A**3-A)*Y2A(NLO)+(B**3-B)*Y2A(NHI))*(DX**2)/SIX
      DYDX=(YA(NHI)-YA(NLO))/DX+
     .     (-((THREE*(A**2)-ONE)*Y2A(NLO))+
     .     (THREE*(B**2)-ONE)*Y2A(NHI))*DX/SIX

      END SUBROUTINE SPLINT
c
      subroutine sort(n,arrin,indx)
c     sorts an array by the heapsort method
c     w. h. preuss et al. numerical recipes  (called indexx there)
C     .. Scalar Arguments ..
      integer n
C     ..
C     .. Array Arguments ..
      double precision arrin(n)
      integer indx(n)
C     ..
C     .. Local Scalars ..
      double precision q
      integer i, indxt, ir, j, l
C     ..

      if (n < 2) call die("Sort called for n<2")
      do 10 j = 1, n
         indx(j) = j
 10       continue
      l = n/2 + 1
      ir = n
 20    continue
      if (l .gt. 1) then
         l = l - 1
         indxt = indx(l)
         q = arrin(indxt)
      else
         indxt = indx(ir)
         q = arrin(indxt)
         indx(ir) = indx(1)
         ir = ir - 1
         if (ir .eq. 1) then
            indx(1) = indxt
c
            return
c
         end if
      end if
      i = l
      j = l + l
 30    continue
      if (j .le. ir) then
         if (j .lt. ir) then
            if (arrin(indx(j)) .lt. arrin(indx(j+1))) j = j + 1
         end if
         if (q .lt. arrin(indx(j))) then
            indx(i) = indx(j)
            i = j
            j = j + j
         else
            j = ir + 1
         end if
c
         go to 30
c
      end if
      indx(i) = indxt
c
      go to 20
c
      end subroutine sort

      end module m_recipes
