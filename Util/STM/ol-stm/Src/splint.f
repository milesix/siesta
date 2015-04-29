! ---
! Copyright (C) 1996-2014	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      SUBROUTINE SPLINT(DELT,YA,Y2A,N,X,Y,DYDX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION YA(N),Y2A(N)
      
      NLO=INT(X/DELT)+1
      NHI=NLO+1

 
      IF (DELT.EQ.0.) PAUSE 'Bad DELT input.'

      A=NHI-X/DELT-1
      B=1.0D0-A
      Y=A*YA(NLO)+B*YA(NHI)+
     *      ((A**3-A)*Y2A(NLO)+(B**3-B)*Y2A(NHI))*(DELT**2)/6.

      DYDX=(YA(NHI)-YA(NLO))/DELT +
     * (-((3*(A**2)-1.0)*Y2A(NLO))+(3*(B**2)-1.0)*Y2A(NHI))*DELT/6.
      RETURN
      END
