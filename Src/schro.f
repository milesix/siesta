      module schro
!
!     Routines related to the integration of the Schrodinger Equation
!     in atoms

      use precision, only : dp
      implicit none
      private
      public schro_eq, energ_deriv, rc_vs_e, rphi_vs_e, vhrtre,
     $       polarization

      CONTAINS

      subroutine schro_eq(Zval,rofi,vps,ve,s,drdi,
     .     nrc,l,a,b,nnodes,nprin,e,g) 
      
      implicit none   

      integer, intent(in)  :: nrc,l,nnodes,nprin
      real(dp), intent(in) :: Zval, rofi(1:nrc),vps(1:nrc),ve(1:nrc),
     $     s(1:nrc),drdi(1:nrc),a,b
      real(dp), intent(out) :: e, g(1:nrc)
      
!
!       Automatic arrays
!
      real(dp) :: h(nrc), y(nrc)

      real(dp) ::  a2b4, r2, vtot, rmax, dr, dnrm, phi, dsq
      integer  ::  ir

      a2b4=a*a*0.25d0

      do ir=2,nrc 
         g(ir)=0.0d0
         r2=(rofi(ir)**2)
         vtot=vps(ir)+ve(ir)+dble(l*(l+1))/r2
         h(ir)=vtot*s(ir)+a2b4
      enddo
      h(1)=h(2)
      g(1)=0.0d0 

      e=-((zval/dble(nprin))**2)
      dr=-1.0d6    
      rmax=rofi(nrc)
      
      call egofv(h,s,nrc,e,g,y,l,zval,a,b,rmax,
     .     nprin,nnodes,dr)

      do ir=2,nrc
         phi=g(ir)
         dsq=sqrt(drdi(ir))
         phi=phi*dsq
         g(ir)=phi
      enddo 
      g(1)=0.0d0 

      dnrm=0.0d0
      do ir=2,nrc
         phi=g(ir)
         dnrm=dnrm+phi*phi*drdi(ir) 
      enddo
      dnrm=sqrt(dnrm)

      do ir=2,nrc
         g(ir)=g(ir)/dnrm
      enddo 

      end subroutine schro_eq

!---------------------------------------------------------------
      subroutine rc_vs_e(a,b,r,vps,ve,nrval,l,el,nnode,rnode)

C     Calculate the position, rnode, of the node number nnode of the 
C     radial wavefunction of the pseudopotential Vps, with angular 
C     momentum  l, and energy el.
C     D. Sanchez-Portal, July 1997.
C     Modify by DSP, July 1999

      implicit none

      integer, intent(in) ::  nrval, nnode, l
      double precision, intent(in) :: a, b, el
      double precision, intent(in) ::  r(:), vps(:), ve(:)
      double precision, intent(out) :: rnode

      real*8 drdi(nrval), g(nrval), h(nrval)  !  Automatic arrays

      integer  nn, ir
      real*8 dexpa, ab, hi, gold, r0, g0, r1, g1

      dexpa=exp(a)
      ab=a*b
      do ir=1,nrval
         drdi(ir)=ab
         ab=dexpa*ab
      enddo            

      
      do ir=2,nrval
         hi=vps(ir)+ve(ir)+dble(l*(l+1))/r(ir)**2-el
         hi=hi*(drdi(ir)**2)
         hi=hi+0.25d0*a**2
         h(ir)=hi
        
      enddo 
      h(1)=h(2)

      
      g(1)=0.0d0
      g(2)=1.0d0
      gold=1.0d0
      rnode=r(nrval)
      nn=0
      do ir=3,nrval

         hi=(10.0d0*h(ir-1)*g(ir-1)+h(ir-2)*g(ir-2))/12.0d0
         hi=hi+2.0d0*g(ir-1)-g(ir-2)
         g(ir)=hi/(1.0d0-h(ir)/12.0d0)
  
         if((g(ir).eq.0.0d0).or.(g(ir)*gold.lt.0.0d0)) then
            r0=r(ir-1)
            g0=gold
            r1=r(ir)
            g1=g(ir)
            rnode=r0-g0*(r1-r0)/(g1-g0)
            nn=nn+1
            if(nn.eq.nnode) goto 50
         endif
         gold=g(ir)

      enddo 
 50   continue 

      end subroutine rc_vs_e
!-----------------------------------------------------------------
!
      subroutine rphi_vs_e(a,b,r,vps,ve,nrval,l,el,rphi,rmax)

C   Calculate the atomic 
C   radial wavefunction of the pseudopotential Vps, with angular
C   momentum  l, and energy el, inside r<Rmax
C   The Schrodinger equation is solved using a simple Numerov 
C   scheme. Rmax should not be taken too big. 
C   D. Sanchez-Portal, July 1999.

      implicit none

      real*8, intent(in) ::  a, b, el, rmax
      integer, intent(in) ::  nrval
      real*8, intent(in) ::   r(:), vps(:), ve(:)
      real*8, intent(out) ::   rphi(:)

      real*8  g(nrval),drdi(nrval),h(nrval)  ! Automatic

      real*8 dexpa, ab, hi, dnrm
      real*8, parameter ::  big=1.0d6
      integer  l, nrc, jr, ir


      dexpa=exp(a)
      ab=a*b
      do ir=1,nrval
         drdi(ir)=ab
         ab=dexpa*ab
      enddo

      do ir=2,nrval
         hi=vps(ir)+ve(ir)+dble(l*(l+1))/r(ir)**2-el
         hi=hi*(drdi(ir)**2)
         hi=hi+0.25d0*a**2
         h(ir)=hi
      enddo
      h(1)=h(2)

      g(1)=0.0d0
      g(2)=1.0d0
      nrc=nint(log(rmax/b+1.0d0)/a)+1
      nrc=min(nrc,nrval)
      nrc=nrc+1-mod(nrc,2)
      do ir=3,nrc

         hi=(10.0d0*h(ir-1)*g(ir-1)+h(ir-2)*g(ir-2))/12.0d0
         hi=hi+2.0d0*g(ir-1)-g(ir-2)
         g(ir)=hi/(1.0d0-h(ir)/12.0d0)
         
         if(abs(g(ir)).gt.big) then 
            dnrm=0.0d0
            do jr=1,ir
               dnrm=dnrm+drdi(jr)*(g(jr)*sqrt(drdi(jr)))**2
            enddo 
            dnrm=sqrt(dnrm)
            do jr=1,ir
               g(jr)=g(jr)/dnrm
            enddo 
         endif 
      enddo

C     Normalize the wavefunction
      dnrm=0.0d0
      do ir=1, nrc
         g(ir)=g(ir)*sqrt(drdi(ir))
         dnrm=dnrm+drdi(ir)*(g(ir)**2)
      enddo
      dnrm=sqrt(dnrm)
      do ir=1, nrc
         rphi(ir)=g(ir)/dnrm
      enddo
      
      end subroutine rphi_vs_e

      subroutine energ_deriv(a,r,psi,vps,
     .     ve,drdi,nrc,l,el,psidev,nrval)
      implicit none
C     This routine calculate the energy derivative of 
C     a given wavefunction.
C     The routine solve and inhomogeneus version of 
C     Schrodinger eqn.  
C     It is not an optimized algorithm!!!!!!!!!!!!!!!!!
C     Written by Daniel Sanchez-Portal, July 1999
C     
      integer, intent(in)     ::   l, nrval
      integer, intent(inout)  ::   nrc   
      real*8, intent(in)      ::   el

      real*8, intent(out)     ::   psidev(:)
      real*8, intent(in)      ::   r(:),psi(:), vps(:), ve(:)
      real*8, intent(in)      ::   drdi(:)

      real*8 g(nrval), h(nrval)
      real*8 hi, dnrm, cons, a, ortog, dnrm2
      
      integer, parameter      :: nrmin=1

      integer ir
      nrc=min(nrc,nrval)
      
C     Solving the inhomogeneus Schrodinger equation
      do ir=2,nrc
         hi=vps(ir)+ve(ir)+l*(l+1)/r(ir)**2-el
         hi=hi*(drdi(ir)**2)
         hi=hi+0.25d0*a**2
         h(ir)=hi
      enddo 
      h(1)=h(2)
      
      cons=psi(nrmin+1)/(vps(nrmin+1)+ve(nrmin+1)-el)
      cons=cons/r(nrmin+1)**(l+1) 
      g(1)=0.0d0
      do ir=1,nrmin+1
         g(ir)=cons*(r(ir)**(l+1))/sqrt(drdi(ir))
      enddo 

      do ir=nrmin+2,nrc
         hi=-((psi(ir)+10.0d0*psi(ir-1)
     .        +psi(ir-2))/12.0d0)

         hi=hi+(10.0d0*h(ir-1)*g(ir-1)+h(ir-2)*g(ir-2))/12.0d0
         
         hi=hi+2.0d0*g(ir-1)-g(ir-2)

         g(ir)=hi/(1.0d0-h(ir)/12.0d0)

      enddo 

C     Orthogonalize the energy derivative to the original wavefunction
C     and normalize
      dnrm2=0.0d0
      ortog=0.0d0
      do ir=1, nrc
         g(ir)=g(ir)*sqrt(drdi(ir))
         dnrm2=dnrm2+drdi(ir)*(psi(ir)**2)
         ortog=ortog+drdi(ir)*g(ir)*psi(ir)
      enddo
      dnrm=0.0d0
      do ir=1, nrc
         g(ir)=g(ir)-ortog*psi(ir)/dnrm2
         dnrm=dnrm+drdi(ir)*(g(ir)**2)
      enddo 
      dnrm=sqrt(dnrm)
      do ir=1,nrc
         psidev(ir)=g(ir)/dnrm
      enddo

      end subroutine energ_deriv

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    This file contains routines written by A.R. Williams
C    around 1985 (possibly based on previous work).
C    They are used for the solutions of the radial Schrodinger's
C    and Poisson's equations in the atomic program.
C    It also contains some routines written by J.M.Soler
C    in collaboration with A.R.Williams and based on 
C    algorithms developed by ARW.
C    J. M. Soler, Nov. 1998.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C * The routines contained in this file are:
C
c     SUBROUTINE EGOFV
c     SUBROUTINE YOFE
c     SUBROUTINE NRMLZG
c     SUBROUTINE BCORGN
c     SUBROUTINE BCRMAX
c     SUBROUTINE NUMIN
c     SUBROUTINE NUMOUT
c     SUBROUTINE VHRTRE
c     SUBROUTINE QVLOFZ  ---> moved to periodic_table.f
c     SUBROUTINE LMXOFZ  ---> moved to periodic_table.f
c     SUBROUTINE CNFIG   ---> moved to periodic_table.f
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




      SUBROUTINE EGOFV(H,S,N,E,G,Y,L,Z,A,B,RMAX,NPRIN,NNODE,DR)
      !implicit none
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  EIOFV DETERMINES THE EIGENENERGY AND WAVEFUNCTION CORRESPONDING
C  TO A PARTICULAR L, PRINCIPAL QUANTUM NUMBER AND BOUNDARY CONDITION.
C
C  TWO FUNDAMENTAL TECHNIQUES ARE USED TO LOCATE THE SOLUTION:
C       1) NODE COUNTING AND BISECTION
C       2) VARIATIONAL ESTIMATE BASED ON A SLOPE DISCONTINUITY IN PSI
C  THE ARGUMENTS ARE DEFINED AS FOLLOWS:
C       H,S: G" = (H-E*S)*G
C       NR: MAXIMUM ALLOWED NUMBER OF RADIAL POINTS
C       E: E(I) IS THE I-TH ENERGY FOUND
C       NE: NUMBER OF ENERGIES FOUND
C       L: THE ANGULAR MOMENTUM
C       NCOR: THE NUMBER OF LOWER-ENERGY STATES
C
C  THE INDIVIDUAL ENERGIES ARE RESOLVED BY PERFORMING A FIXED NUMBER
C  OF BISECTIONS AFTER A GIVEN EIGENVALUE HAS BEEN ISOLATED
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
!      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      Integer, intent(in) :: N
      real(dp), intent(in) :: H(1:N),S(1:N)
      real(dp), intent(out) :: Y(1:N)
      real(dp), intent(out) :: G(1:N)
      real(dp), intent(out) :: E
      integer, intent(in) :: L
      real(dp),intent(in) :: Z
      real(dp),intent(in) :: A,B
      real(dp),intent(in) :: rmax, dr
      integer, intent(in) :: nprin,nnode
!      DIMENSION H(N),S(N),G(N),Y(*)

      !Internal vars
      real(dp) :: de, del, et, e1, e2, t
      integer  :: nt,n1,n2, niter,i,ncor
      real(dp), parameter :: TOL=1.D-5

C Added by JDG as NT was used before being initialised
      NT = 0 
      NCOR=NPRIN-L-1
      N1=NNODE
      N2=NNODE-1
      E1=E
      E2=E
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  THE LABELS 1 AND 2 REFER TO THE BISECTION PROCESS, DEFINING THE
C  RANGE IN WHICH THE DESIRED SOLUTION IS LOCATED.  THE INITIAL
C  SETTINGS OF N1, N2, E1 AND E2 ARE NOT CONSISTENT WITH THE BISECTION
C  ALGORITHM; THEY ARE SET TO CONSISTENT VALUES WHEN THE DESIRED
C  ENERGY INTERVAL HAS BEEN LOCATED.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DEL=5.D-1
      DE=0.D0
      NITER = 0
 1    NITER = NITER + 1
      IF(NITER.GT.40) GO TO 3
      ET=E+DE
C  THE FOLLOWING LINE IS THE FUNDAMENTAL "BISECTION"

      E=0.5*(E1+E2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  THE FOLLOWING CONCATENATION OF LOGICAL ORS ENSURES THAT NODE
C  COUNTING IS USED UNLESS THE PREVIOUS INTEGRATION OF THE RADIAL
C  EQ PRODUCED BOTH THE CORRECT NUMBER OF NODES AND A SENSIBLE
C  PREDICTION FOR THE ENERGY.
C
C     SENSIBLE MEANS THAT ET MUST BE GREATER THAN E1 AND LESS THAN E2
C     CORRECT NUMBER OF NODES MEANS THAT NT = NNODE OR NNODE-1.
C
C     LEAVING E SET TO ITS BISECTION VALUE, AND TRANSFERING TO
C     THE CALL TO YOFE MEANS THAT WE ARE PERFORMING BISECTION,
C     WHEREAS SETTING E TO ET IS USE OF THE VARIATIONAL ESTIMATE.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (ET.LE.E1 .OR. ET.GE.E2 .OR.
     1       NT.LT.NNODE-1 .OR. NT.GT.NNODE) GO TO 2
      E=ET
      IF(DABS(DE).LT.TOL) GO TO 6
 2    CALL YOFE(E,DE,DR,RMAX,H,S,Y,N,L,NCOR,NT,Z,A,B)
C     WRITE(6,101) L,DR,N1,NT,NNODE,N2,E1,E,E2,DE
C101  FORMAT('  L     DR     N1  NT   N  N2       E1           E',
C    1       '          E2          DE'/I3,D10.3,4I4,4F12.5)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  YOFE INTEGRATES THE SCHRO EQ.; NOW THE BISECTION LOGIC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(NT.GE.NNODE) GO TO 5
C  TOO FEW NODES; SET E1 AND N1
      E1=E
      N1=NT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  AT THIS POINT, WE HAVE JUST SET THE BOTTOM OF THE BISECTION RANGE;
C  IF THE TOP IS ALSO SET, WE PROCEDE.  IF THE TOP OF THE RANGE HAS NOT
C  BEEN SET, IT MEANS THAT WE HAVE YET TO FIND AN E GREATER THAN THE
C  DESIRED ENERGY.  THE UPPER END OF THE RANGE IS EXTENDED.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(N2.GE.NNODE) GO TO 1
      DEL=DEL*2.D0
      E2=E1+DEL
      GO TO 1
C  TOO MANY NODES; SET E2 AND N2
 5    E2=E
      N2=NT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  AT THIS POINT, WE HAVE JUST SET THE TOP OF THE BISECTION RANGE;
C  IF THE TOP IS ALSO SET, WE PROCEDE.  IF THE TOP OF THE RANGE HAS
C  NOT BEEN SET, IT MEANS THAT WE HAVE YET TO FIND AN E LESS THAN THE
C  DESIRED ENERGY.  THE LOWER END OF THE RANGE IS EXTENDED.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(N1.LT.NNODE) GO TO 1
      DEL=DEL*2.D0
      E1=E2-DEL
      GO TO 1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  THE NUMEROV METHOD USES A TRANSFORMATION OF THE RADIAL WAVE FCN.
C  THAT WE CALL "Y".  HAVING LOCATED THE EIGENENERGY, WE TRANSFORM
C  Y TO "G", FROM WHICH THE DENSITY IS EASILY CONSTRUCTED.
C  FINALLY, THE CALL TO "NRMLZG" NORMALIZES G TO ONE ELECTRON.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 6    G(1) = 0.D0
      DO 7 I=2,N
           T=H(I)-E*S(I)
           G(I)=Y(I)/(1.D0-T/12.D0)
 7    CONTINUE
      CALL NRMLZG(G,S,N)
      RETURN
 3    WRITE(6,4) Z,L,NNODE,E,DE
 4    FORMAT(' EGOFV: TOO MANY ITERATIONS; EXECUTION STOPPING'/
     1       ' Z=',F3.0,'  L=',I2,'  NNODE=',I2,'  E=',F12.5,
     2       '  DE=',F12.5)
      STOP 8
      END subroutine egofv



      SUBROUTINE YOFE(E,DE,DR,RMAX,H,S,Y,NMAX,L,NCOR,NNODE,Z,A,B)
      !implicit none
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   YOFE INTEGRATES THE RADIAL SCHRODINGER EQN USING THE NUMEROV
C   METHOD.
C
C   THE ARGUMENTS ARE DEFINED AS FOLLOWS:
C       E IS THE OLD ENERGY(OVERWRITTEN) BY THE NEW ENERGY
C       DE IS THE E CHANGE PREDICTED TO ELIM THE KINK IN PSI
C       DR IS THE LOG DERIV (THE BOUNDARY CONDITION)
C       G" = (H-ES)G (ALL DIAGONAL IN I (RADIUS) )
C       Y IS THE NUMEROV INDEPENDENT VARIABLE Y = G - G"/12
C       N IS THE NUMBER OF RADIAL MESH POINTS
C       L IS THE ANGULAR MOMENTUM
C       NCOR IS THE NUMBER OF STATES OF LOWER ENERGY
C       NNODE IS 1 + THE NUMBER OF INTERIOR NODES IN PSI
C       Z IS THE ATOMIC NUMBER
C       A AND B SPECIFY THE RADIAL MESH R(I)=(EXP(A*(I-1))-1)*B
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real(dp), intent(inout) :: E

      real(dp), intent(out) :: de
      real(dp), intent(in) :: rmax,dr
      integer, intent(in) :: nmax, l,ncor
      integer, intent(inout) :: nnode
      real(dp), intent(in) :: z,a,b
      real(dp), intent(in) :: H(1:NMAX),S(1:NMAX)
      real(dp), intent(inout) :: Y(1:NMAX)

      !Internal vars
      real(dp) :: zdr,yn
      integer  :: n, knk,i,nndin

      ZDR = Z*A*B
      N=NMAX

 8    IF( H(N)-E*S(N) .LT. 1.D0 ) GO TO 9

      Y(N)=0.D0
      N=N-1
      GO TO 8
 9    CONTINUE

      
      CALL BCORGN(E,H,S,N,L,ZDR,Y2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  BCORGN COMPUTES Y2, WHICH EMBODIES THE BOUNDARY CONDITION
C  SATISFIED BY THE RADIAL WAVE FUNCTION AT THE ORIGIN
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=
      KNK=N
      CALL NUMOUT(E,H,S,Y,NCOR,KNK,NNODE,Y2,G,GSG,X)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  THE OUTWARD INTEGRATION IS NOW COMPLETE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     WE FIRST DECIDE IF THE KINETIC ENERGY IS SUFFICIENTLY NON
C     NEGATIVE TO PERMIT USE OF THE NUMEROV EQ AT RMAX.  IF
C     IT IS NOT, THEN ZERO-VALUE BOUNDARY CONDITION IS USED
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      YN=0.D0
      IF(N.LT.NMAX .OR. DABS(DR).GT.1.D3) GO TO 7
      CALL BCRMAX(E,DR,RMAX,H,S,N,YN,A,B)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  BCRMAX COMPUTES YN, WHICH EMBODIES THE BOUNDARY CONDITION
C  SATISFIED BY THE RADIAL WAVE FUNCTION AT RMAX
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=
 7    CALL NUMIN(E,H,S,Y,N,NNDIN,YN,GIN,GSGIN,XIN,KNK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  NUMIN PERFORMS THE INWARD INTEGRATION BY THE NUMEROV METHOD
C
C  THE ENERGY INCREMENT IS NOW EVALUATED FROM THE KINK IN PSI
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=
      RATIO = G/GIN
      XIN=XIN*RATIO
      GSG=GSG+GSGIN*RATIO*RATIO
      T=H(KNK)-E*S(KNK)
      DE=G*(X+XIN+T*G)/GSG
      NNODE=NNODE+NNDIN
      IF(DE.LT.0.D0) NNODE=NNODE+1
      DO 6 I=KNK,N
         Y(I) = Y(I)*RATIO
 6    CONTINUE
      RETURN
      END subroutine yofe


      SUBROUTINE NRMLZG(G,S,N)
!***********************************************************************
! Normalizes the radial wavefunction, using Simpson's rule for the norm
! Input:
!   real*8  G(N) : Wavefunction to be normalized. G is related to the
!                  true radial wavefunction psi by
!                     G(j) = (dr/dj)^(3/2) * r(j) * psi(r(j))
!   real*8  S(N) : Metric function defined for a logarithmic mesh
!                    r(j) = B*(exp(A*(j-1)) - 1),   j=1,2,...,N
!                  as  S(j) = (dr/dj)^2 = (A*r(j))^2
!   integer N    : Number of radial points (including r(1)=0)
! Output:
!   real*8  G(N) : Normalized wavefunction
!***********************************************************************

      IMPLICIT NONE
      INTEGER          :: N
      DOUBLE PRECISION :: G(N), S(N)

      INTEGER          :: I
      DOUBLE PRECISION :: NORM, SRNRM

      call integrator(g(1:n)*g(1:n),s,n,norm)

      ! Normalize wavefunction
      SRNRM = SQRT(NORM)
      DO I=1,N
         G(I) = G(I)/SRNRM
      END DO

      END SUBROUTINE NRMLZG

      SUBROUTINE BCORGN(E,H,S,N,L,ZDR,Y2)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   YOFE INTEGRATES THE RADIAL SCHRODINGER EQN USING THE NUMEROV
C   METHOD.
C
C   THE ARGUMENTS ARE DEFINED AS FOLLOWS:
C       E IS THE OLD ENERGY(OVERWRITTEN) BY THE NEW ENERGY
C       DE IS THE E CHANGE PREDICTED TO ELIM THE KINK IN PSI
C       DR IS THE LOG DERIV (THE BOUNDARY CONDITION)
C       G" = (H-ES)G (ALL DIAGONAL IN I (RADIUS) )
C       Y IS THE NUMEROV INDEPENDENT VARIABLE Y = G - G"/12
C       N IS THE NUMBER OF RADIAL MESH POINTS
C       L IS THE ANGULAR MOMENTUM
C       NNODE IS 1 + THE NUMBER OF INTERIOR NODES IN PSI
C       Z IS THE ATOMIC NUMBER
C       A AND B SPECIFY THE RADIAL MESH R(I)=(EXP(A*(I-1))-1)*B
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
      real(dp), intent(in) :: E
      integer, intent(in)  :: N,L
      real(dp), intent(in) :: H(1:N),S(1:N)
      real(dp), intent(in) :: zdr
      real(dp), intent(out) :: Y2

      !Internal vars
      real(dp) :: T2,T3,D2,C0,C1,C2
C
C   THE QUANTITY CALLED D(I) IN THE PROGRAM IS ACTUALLY THE INVERSE
C   OF THE DIAGONAL OF THE TRI-DIAGONAL NUMEROV MATRIX
C
      T2=H(2)-E*S(2)
      D2=-((24.D0+10.D0*T2)/(12.D0-T2))
C
C=================================================================
C  THE FOLLOWING SECTION DEALS WITH THE FACT THAT THE INDEPENDENT
C  VARIABLE "Y" IN THE NUMEROV EQUATION IS NOT ZERO AT THE ORIGIN
C  FOR L LESS THAN 2
C  THE L=0 SOLUTION G VANISHES, BUT THE FIRST AND SECOND
C  DERIVATIVES ARE FINITE, MAKING THE NUMEROV VARIABLE Y FINITE
C  THE L=1 SOLUTION G VANISHES, AND G' ALSO VANISHES, BUT
C  THE SECOND DERIVATIVE G" IS FINITE MAKING Y FINITE.  FOR L > 1,
C  G AND ITS FIRST TWO DERIVATIVES VANISH, MAKING Y ZERO.
C=================================================================
      IF(L.GE.2) GOTO 3
      IF(L.GT.0) GOTO 1
      C0=ZDR/6.D0
      C0=C0/(1.D0-0.75*ZDR)
      GO TO 2
 1    C0=1.D0/12.D0
      C0=(-C0)*8.D0/3.D0
 2    C1=C0*(12.D0+13.D0*T2)/(12.D0-T2)
      T3=H(3)-E*S(3)
      C2=(-5.D-1)*C0*(24.D0-T3)/(12.D0-T3)
      D2=(D2-C1)/(1.D0-C2)
 3    Y2=(-1.D0)/D2
      RETURN
      END SUBROUTINE BCORGN

      SUBROUTINE BCRMAX(E,DR,RMAX,H,S,N,YN,A,B)
      implicit none
C
C 22.7.85
C
      !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer, intent(in) :: N
      real(dp), intent(in) :: H(1:N),S(1:N)
      real(dp),intent(in) :: E,DR,RMAX,A,B
      real(dp), intent(out) :: YN

      !Internal vars
      real(dp) ::TNM1,TN,TNP1,BETA,DG,C1,C2,C3,DN
C
C     WRITE(6,*) 'BCRMAX:',DR
      TNM1=H(N-1)-E*S(N-1)
      TN  =H(N  )-E*S(N  )
      TNP1=H(N+1)-E*S(N+1)
      BETA=1.D0+B/RMAX
      DG=A*BETA*(DR+1.D0-5.D-1/BETA)
C
C
      C2=24.D0*DG/(12.D0-TN)
      DN=-((24.D0+10.D0*TN)/(12.D0-TN))
C
      C1= (1.D0-TNM1/6.D0)/(1.D0-TNM1/12.D0)
      C3=-((1.D0-TNP1/6.D0)/(1.D0-TNP1/12.D0))
      YN=-((1.D0-C1/C3)/(DN-C2/C3))
C
C
      RETURN
      END       SUBROUTINE BCRMAX



      SUBROUTINE NUMIN(E,H,S,Y,N,NNODE,YN,G,GSG,X,KNK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   YOFE INTEGRATES THE RADIAL SCHRODINGER EQN USING THE NUMEROV
C   METHOD.
C
C   THE ARGUMENTS ARE DEFINED AS FOLLOWS:
C       E IS THE OLD ENERGY(OVERWRITTEN) BY THE NEW ENERGY
C       DE IS THE E CHANGE PREDICTED TO ELIM THE KINK IN PSI
C       DR IS THE LOG DERIV (THE BOUNDARY CONDITION)
C       G" = (H-ES)G (ALL DIAGONAL IN I (RADIUS) )
C       Y IS THE NUMEROV INDEPENDENT VARIABLE Y = G - G"/12
C       N IS THE NUMBER OF RADIAL MESH POINTS
C       L IS THE ANGULAR MOMENTUM
C       NNODE IS 1 + THE NUMBER OF INTERIOR NODES IN PSI
C       Z IS THE ATOMIC NUMBER
C       A AND B SPECIFY THE RADIAL MESH R(I)=(EXP(A*(I-1))-1)*B
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
      real(dp), intent(in) :: E
      integer, intent(in)  :: N
      real(dp), intent(in) :: H(1:N),S(1:N)
      real(dp), intent(out) :: Y(1:N)
      integer, intent(out) :: NNODE
      real(dp), intent(in) :: YN
      real(dp), intent(out) :: G
      real(dp), intent(out) :: GSG
      real(dp), intent(out) :: X
      integer, intent(out)   :: KNK
      
      !Internal vars
      real(dp) :: T
      integer  :: I

      Y(N)=YN
      T=H(N)-E*S(N)
      G=Y(N)/(1.D0-T/12.D0)
      GSG=G*S(N)*G
      I=N-1
      Y(I)=1.D0
      T=H(I)-E*S(I)
      G=Y(I)/(1.D0-T/12.D0)
      GSG=GSG+G*S(I)*G
      X=Y(I)-Y(N)
      NNODE=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  BEGIN THE INWARD INTEGRATIONBY THE NUMEROV METHOD
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=
 1    X=X+T*G
      I=I-1
      Y(I)=Y(I+1)+X
      IF( Y(I)*Y(I+1) .LT. 0.D0) NNODE=NNODE+1
      T=H(I)-E*S(I)
      G=Y(I)/(1.D0-T/12.D0)
      GSG=GSG+G*S(I)*G
      IF(I.GT.KNK) GO TO 1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  THE LAST STATEMENT DEFINES THE KINK RADIUS AS THE POINT WHERE
C  PSI FIRST TURNS DOWNWARD.  THIS USUALLY MEANS AT THE OUTERMOST
C  MAXIMUM
C
C  THE INWARD INTEGRATION IS NOW COMPLETE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=
      KNK=I
      RETURN
      END       SUBROUTINE NUMIN


      SUBROUTINE NUMOUT(E,H,S,Y,NCOR,KNK,NNODE,Y2,G,GSG,X)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   YOFE INTEGRATES THE RADIAL SCHRODINGER EQN USING THE NUMEROV
C   METHOD.
C
C   THE ARGUMENTS ARE DEFINED AS FOLLOWS:
C       E IS THE OLD ENERGY(OVERWRITTEN) BY THE NEW ENERGY
C       DE IS THE E CHANGE PREDICTED TO ELIM THE KINK IN PSI
C       DR IS THE LOG DERIV (THE BOUNDARY CONDITION)
C       G" = (H-ES)G (ALL DIAGONAL IN I (RADIUS) )
C       Y IS THE NUMEROV INDEPENDENT VARIABLE Y = G - G"/12
C       N IS THE NUMBER OF RADIAL MESH POINTS
C       L IS THE ANGULAR MOMENTUM
C       NNODE IS 1 + THE NUMBER OF INTERIOR NODES IN PSI
C       Z IS THE ATOMIC NUMBER
C       A AND B SPECIFY THE RADIAL MESH R(I)=(EXP(A*(I-1))-1)*B
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
      real(dp), intent(in) :: E
      real(dp), intent(in) :: H(1:KNK), S(1:KNK)
      real(dp), intent(out) :: Y(1:KNK)
      integer, intent(in) :: ncor
      integer, intent(inout) :: knk
      integer, intent(out) :: nnode
      real(dp), intent(in) :: y2
      real(dp), intent(out) :: g,gsg,x
      
      !Internal vars
      real(dp) :: T, XL
      integer  :: I,NM4
      
      !DOUBLE PRECISION H(KNK),S(KNK),Y(KNK)
      Y(1)=0.D0
      Y(2)=Y2
      T=H(2)-E*S(2)
      G=Y(2)/(1.D0-T/12.D0)
      GSG=G*S(2)*G
      Y(3)=1.D0
      T=H(3)-E*S(3)
      G=Y(3)/(1.D0-T/12.D0)
      GSG=GSG+G*S(3)*G
      X=Y(3)-Y(2)
      I=3
      NNODE=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  BEGIN THE OUTWARD INTEGRATIONBY THE NUMEROV METHOD
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=
      NM4=KNK-4
 1    XL=X
      X=X+T*G
      I=I+1
      Y(I)=Y(I-1)+X
C     WRITE(6,300) I,Y(I),X,T,H(I),S(I)
C300  FORMAT(I5,5D14.5)
      IF( Y(I)*Y(I-1) .LT. 0.D0) NNODE=NNODE+1
      T=H(I)-E*S(I)
      G=Y(I)/(1.D0-T/12.D0)
      GSG=GSG+G*S(I)*G
      IF(I.EQ.NM4) GO TO 2
      IF(NNODE.LT.NCOR) GO TO 1
      IF(XL*X.GT.0.D0) GO TO 1
 2    KNK=I
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  THE OUTWARD INTEGRATION IS NOW COMPLETE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=
      RETURN
      END       SUBROUTINE NUMOUT


      SUBROUTINE VHRTRE(RHO,V,R,DRDI,SRDRDI,NR,A)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   VHRTRE CONSTRUCTS THE ELECTROSTATIC POTENTIAL DUE TO A SUPPLIED
C   ELECTRON DENSITY.  THE NUMEROV METHOD IS USED TO INTEGRATE
C   POISSON'S EQN.
C
C   DESCRIPTION OF ARGUMENTS:
C      RHO....4*PI*R**2 * THE ELECTRON DENSITY FOR WHICH WE CALCULATING
C             THE ELECTROSTATIC POTENTIAL
C      V......THE ELECTROSTATIC POTENTIAL DUE TO THE ELECTRON DENSITY
C             RHO.  THE CONSTANTS OF INTEGRATION ARE FIXED SO THAT THE
C             POTENTIAL TENDS TO A CONSTANT AT THE ORIGIN AND TO
C             2*Q/R AT R=R(NR), WHERE Q IS THE INTEGRATED CHARGE
C             CONTAINED IN RHO(R)
C      R......THE RADIAL MESH R(I) = B*(EXP(A(I-1))-1)
C      NR.....THE NUMBER OF RADIAL MESH POINTS
C      DRDI...DR(I)/DI
C      SRDRDI.SQRT(DR/DI)
C      A......THE PARAMETER APPEARING IN R(I) = B*(EXP(A(I-1))-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer, intent(in)  :: NR
      real(dp), intent(in) :: RHO(1:NR),R(1:NR),DRDI(1:NR),SRDRDI(1:NR)
      real(dp), intent(out) :: V(1:NR)
      real(dp), intent(in)  :: A

      !Internal vars
      integer :: nrm1, nrm2,ir
      real(dp) :: ybyq,qbyy,v0,q,qt,dz,t,beta,x,y,dv,qpartc,a2by4

      NRM1=NR-1
      NRM2=NR-2
      A2BY4=A*A/4.D0
      YBYQ=1.D0-A*A/48.D0
      QBYY=1.D0/YBYQ
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  SIMPSON'S RULE IS USED TO PERFORM TWO INTEGRALS OVER THE ELECTRON
C  DENSITY.  THE TOTAL CHARGE QT IS USED TO FIX THE POTENTIAL AT R=R(NR)
C  AND V0 (THE INTEGRAL OF THE ELECTRON DENSITY DIVIDED BY R) FIXES
C  THE ELECTROSTATIC POTENTIAL AT THE ORIGIN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      V0=0.D0
      QT=0.D0
      DO IR=2,NRM1,2
        DZ=DRDI(IR)*RHO(IR)
        QT=QT+DZ
        V0=V0+DZ/R(IR)
      ENDDO
      V0=V0+V0
      QT=QT+QT
      DO IR=3,NRM2,2
        DZ=DRDI(IR)*RHO(IR)
        QT=QT+DZ
        V0=V0+DZ/R(IR)
      ENDDO
      DZ=DRDI(NR)*RHO(NR)
      QT=(QT+QT+DZ)/3.D0
      V0=(V0+V0+DZ/R(NR))/3.D0
      V(1)=2.D0*V0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  THE ELECTROSTATIC POTENTIAL AT R=0 IS SET EQUAL TO
C                       THE AVERAGE VALUE OF RHO(R)/R
C  BEGIN CONSTRUCTION OF THE POTENTIAL AT FINITE R
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IR=2
      T=SRDRDI(IR)/R(IR)
      BETA=DRDI(IR)*T*RHO(IR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  THE NEXT 4 STATEMENTS INDICATE THAT WE FIRST FIND THE PARTICULAR
C  SOLUTION TO THE INHOMOGENEOUS EQN. FOR WHICH Q(2)=0, WE THEN
C  ADD TO THIS PARTICULAR SOLUTION A SOLUTION OF THE HOMOGENEOUS EQN.
C  (A CONSTANT IN V OR A Q PROPORTIONAL TO R)
C  WHICH WHEN DIVIDED BY R IN GOING FROM Q TO V GIVES
C  THE POTENTIAL THE DESIRED COULOMB TAIL OUTSIDE THE ELECTRON DENSITY.
C  THE SIGNIFICANCE OF THE SOLUTION VANISHING AT THE SECOND RADIAL
C  MESH POINT IS THAT, SINCE ALL REGULAR SOLUTIONS OF THE EQUATION
C  FOR Q=R*V VANISH AT THE ORIGIN, THE KNOWLEDGE OF THE SOLUTION
C  VALUE AT THE SECOND MESH POINT PROVIDES THE TWO SOLUTION VALUES
C  REQUIRED TO START THE NUMEROV PROCEDURE.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      X=0.D0
      Y=0.D0
      Q=(Y-BETA/12.D0)*QBYY
      V(IR)=2.D0*T*Q
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  BEGINNING OF THE NUMEROV ALGORITHM
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 3    X=X+A2BY4*Q-BETA
      Y=Y+X
      IR=IR+1
      T=SRDRDI(IR)/R(IR)
      BETA=T*DRDI(IR)*RHO(IR)
      Q=(Y-BETA/12.D0)*QBYY
      V(IR)=2.D0*T*Q
      IF(IR.LT.NR) GO TO 3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  END OF THE NUMEROV ALGORITHM
C
C  WE HAVE NOW FOUND A PARTICULAR SOLUTION TO THE INHOMOGENEOUS EQN.
C  FOR WHICH Q(R) AT THE SECOND RADIAL MESH POINT EQUALS ZERO.
C  NOTE THAT ALL REGULAR SOLUTIONS TO THE EQUATION FOR Q=R*V
C  VANISH AT THE ORIGIN.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      QPARTC = R(NR)*V(NR)/2.D0
      DZ=QT-QPARTC
      DV=2.D0*DZ/R(NR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  THE LOOP FOLLOWING ADDS THE CONSTANT SOLUTION OF THE HOMOGENEOUS
C  EQN TO THE PARTICULAR SOLUTION OF THE INHOMOGENEOUS EQN.
C  NOTE THAT V(1) IS CONSTRUCTED INDEPENDENTLY
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO IR=2,NR
         V(IR)=V(IR)+DV
      ENDDO
      RETURN
      END       SUBROUTINE VHRTRE


      subroutine polarization(a,r,psi,vps,
     .      ve,drdi,nrc,l,el,psipol,nrval)


        integer, intent(in) :: nrval, nrc, l
        real*8, intent(in)  :: a, el
        real*8, intent(in) :: psi(:) !! **AG E.Anglada:use psi_copy
        real*8, intent(in)  :: r(:),vps(:),ve(:),drdi(:)
        real*8, intent(out) :: psipol(:)
        
C      This routine calculates the polarization (unoccupied) 
C      orbitals with angular momentum l from the atomic 
C      orbitals with l-1 angular momentum, using a perturbative
C      approach.
C      The routine solves an inhomogeneus version of 
C      Schrodinger eqn.  
C      It is not an optimized algorithm!!!!!!!!!!!!!!!!!
C
C
C       cons1 is a big number. If it is choosen too 
C       big, too many iterations will be needed.
C       If cons1 is too small it may happen that 
C       the routine never converges to the solution.
C
C       If Rc is too big the very simple (unoptimized)
C       algorithm used here cannot converge. That is 
C       why Rc's are required to be smaller than Rint,
C       where Rint should not be greater than aprox. 15 Bohr
C       Written by Daniel Sanchez-Portal, July 1997
C 

        integer:: nrmin, niter
        real(dp):: cons1, rint
        parameter(nrmin=1,niter=1000,Cons1=1.0d5,rint=15.0d0)


        real(dp):: g(nrval), h(nrval), psi_copy(nrval)   !  Automatic

        real(dp):: rmax, reduc, dl, hi, rnd1, c1, c2, rnodo, cons, gold
        real(dp):: gmax, r0, g0, r1, g1, grmx, dff1, dff2, savecons,dnrm
        integer :: index, nnodes, iter, nnd, ir

        psi_copy=0.0d0
        psi_copy(1:nrc)=psi(1:nrc)
        rmax=r(nrval)

        if(rmax.gt.rint) then 
          write(6,*) 'POLARIZATION: Rc for the polarization orbitals'
          write(6,*) 'must be smaller than ',rint,' Bohr'
          STOP
        endif

        do ir=nrc+1,nrval 
          psi_copy(ir)=0.0d0 
        enddo 
 
          reduc=-0.5d0
C**** We calculate the polarization function with angular** 
C momentum l=l+dl
          dl=1
C****
          do ir=2,nrval
            hi=vps(ir)+ve(ir)+(l+dl)*(l+dl+1)/r(ir)**2-el
            hi=hi*(drdi(ir)**2)
            hi=hi+0.25d0*a**2
            h(ir)=hi
          enddo 
          h(1)=h(2)

          rnd1=0.0d0
          index=1
          nnodes=1
CInitialized c1 and c2 to arbitrary values.............
          c1=0.0d0
          c2=0.0d0
          do iter=1,niter
           rnodo=0.0d0
           if(index.eq.1) then 
              cons=cons1
              index=2
           else
              cons=c2
           endif
          
          g(1)=0.0d0
          do ir=1,nrmin+1
            g(ir)=cons*(r(ir)**(l+dl+1))/sqrt(drdi(ir))
          enddo 
          gold=g(nrmin+1)

          nnd=0
          gmax=0.0d0
          do ir=nrmin+2,nrval
            hi=-((r(ir)*psi_copy(ir)+10.0d0*r(ir-1)*psi_copy(ir-1)
     .         +r(ir-2)*psi_copy(ir-2))/12.0d0)

            hi=hi+(10.0d0*h(ir-1)*g(ir-1)+h(ir-2)*g(ir-2))/12.0d0
 
            hi=hi+2.0d0*g(ir-1)-g(ir-2)

            g(ir)=hi/(1.0d0-h(ir)/12.0d0)
            gmax=max(gmax,abs(g(ir)))
            if((g(ir).eq.0.0d0).or.(g(ir)*gold.lt.0.0d0)) then
              nnd=nnd+1
              if (nnd.eq.nnodes) then 
                  r0=r(ir-1)
                  g0=gold
                  r1=r(ir)
                  g1=g(ir)
                  rnodo=r0-g0*(r1-r0)/(g1-g0)
              endif 
            endif
           gold=g(ir)
          enddo 

          grmx=g(nrval)/gmax

          if(((abs(rnodo-rmax).lt.1.0d-3).and.
     .      (abs(grmx).lt.1.0d-7) )
     .        .or. 
     .        ((rnodo.eq.0.0d0).and.
     .        (abs(grmx).lt.1.0d-7) ) ) goto 100

*We begin by finding a node!!!!**

          if((rnd1.eq.0.0d0).and.(rnodo.eq.0.0d0)) then  
             c2=(-reduc)*cons
             if(abs(c2).le.1.0d0/abs(cons1)) then 
               index=1
               rnd1=0.0d0
               reduc=(1.0d0+reduc)/2.0d0
             endif  
          elseif((rnd1.eq.0.0d0).and.(rnodo.ne.0.0d0)) then
              rnd1=rnodo
              c1=cons
              c2=2.0d0*cons
          endif  
       
****

      
*Now we lead this node to Rc*
          if((rnd1.ne.0.0d0).and.(rnodo.eq.0.0d0)) then 
              c2=0.50d0*(c1+c2)
          elseif((rnd1.ne.0.0d0).and.(rnodo.ne.0.0d0)) then 
              if(abs(rnd1-rnodo).gt.1.0d-6)then 
                 dff1=abs(rnd1-rmax) 
                 dff2=abs(rnodo-rmax) 
                 if(dff1.gt.dff2) then 
                   savecons=c2
                   c2=(rmax-rnd1)*(c1-c2)/(rnd1-rnodo)+c1
                   c1=savecons
                   rnd1=rnodo
                 else
                   c2=1.10d0*c2
                 endif 
              else

               if(abs(cons).gt.1.0d15) then 
                  nnodes=nnodes+1
                  index=1
                  rnd1=0.0d0
               else
                 c2=1.1d0*c2
               endif  

              endif 
           endif 

          enddo 
            write(6,*)'POLARIZATION: Iteration to find the polarization'
            write(6,*)'orbital has failed !!!!!!!!!'
            write(6,*)'Please try with a Rc no bigger than ',rnd1,
     .        ' Bohr'
            STOP
                          
100       continue
          dnrm=0.0d0
          do ir=1,nrval
             g(ir)=g(ir)*sqrt(drdi(ir))
             dnrm=dnrm+drdi(ir)*(g(ir)**2)
          enddo 
          dnrm=sqrt(dnrm)
          do ir=1,nrval
               psipol(ir) = g(ir)/dnrm
          enddo 

      end subroutine polarization


      SUBROUTINE INTEGRATOR(F,S,NP,VAL)
!***********************************************************************
! Integrates a radial function tabulated on a logarithmic grid,
! using a generalized Simpson's rule valid for both even and odd
! number of points. Note that the "h" is 1 as the reparametrization
! involves a mapping of integers to reals.
!
! Alberto Garcia, Dec. 2006, based on code by Art Williams.
!
! Input:
!   real*8  F(NP) : Function to be integrated.
!   real*8  S(NP) : Metric function defined for a logarithmic mesh
!                    r(j) = B*(exp(A*(j-1)) - 1),   j=1,2,...,NP
!                  as  S(j) = (dr/dj)^2 = (A*r(j))^2
!   integer NP    : Number of radial points (including r(1)=0)
! Output:
!   real*8  VAL   : Value of the integral
!***********************************************************************

      IMPLICIT NONE
      INTEGER          :: NP
      DOUBLE PRECISION :: F(NP), S(NP)

      INTEGER          :: I, N
      DOUBLE PRECISION :: VAL

      IF (MOD(NP,2).EQ.1) THEN        
         N = NP               ! ODD
      ELSE
         IF (NP .EQ. 2) THEN
          ! Special case of trapezoidal rule
            VAL = 0.5D0 * (F(1)*S(1) + F(2)*S(2))
            RETURN
         ENDIF
         N = NP - 3           ! EVEN: TAKE A FINAL FOUR-POINT INTERVAL
      ENDIF
!
!     STANDARD EXTENDED SIMPSON RULE WITH ALTERNATING 4/3 AND 2/3 FACTORS
!     FOR THE SECTION MADE UP OF PAIRS OF INTERVALS (THREE-POINT SEGMENTS)
!
      VAL = 0.D0
      DO I = 2,N-1,2
         VAL=VAL + F(I)*S(I)
      END DO
      VAL = VAL * 2.D0
      DO I = 3,N-2,2
         VAL=VAL + F(I)*S(I)
      END DO
      VAL = VAL * 2.D0
      DO I = 1,N,N-1                ! first and last points at 1/3
         VAL=VAL + F(I)*S(I)
      END DO
      VAL = VAL/3.D0

      IF (MOD(NP,2).EQ.0) THEN           ! EVEN
!        ADD THE CONTRIBUTION OF THE 
!        FINAL FOUR-POINT SEGMENT USING SIMPSON'S 3/8 RULE
!        (SEE NUMERICAL RECIPES (FORTRAN), P. 105)
         I = NP - 3
         VAL = VAL +
     $      (3.D0/8.D0) * ( (F(I)*S(I) + F(I+3)*S(I+3)) +
     $         3.D0 * (F(I+1)*S(I+1) + F(I+2)*S(I+2)) )
      ENDIF

      END SUBROUTINE INTEGRATOR

      end module schro
