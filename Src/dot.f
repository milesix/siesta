      DOUBLE PRECISION FUNCTION DOT(A,B,N)
C
C     RETURNS REAL*8 SCALAR PRODUCT OF TWO REAL*8 VECTORS
C
      integer :: n, i
      DOUBLE PRECISION A(N),B(N),SUM
      SUM=0.D0
      DO I=1,N
         SUM = SUM + A(I) * B(I)
      ENDDO
      DOT=SUM
      RETURN
      END
