! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!

      CHARACTER*(*) FUNCTION PASTE( STR1, STR2 )

C CONCATENATES THE STRINGS STR1 AND STR2 REMOVING BLANKS IN BETWEEN
C Written by J. Soler

      CHARACTER*(*) STR1, STR2
      integer :: l

      DO 10 L = LEN( STR1 ), 1, -1
         IF (STR1(L:L) .NE. ' ') GOTO 20
   10 CONTINUE
   20 PASTE = STR1(1:L)//STR2
      END
c Concatenates the integer nitter and STR2 - Rafi Ullah, July 2014
      CHARACTER*(*) FUNCTION NPASTE( nitter, STR2 )
      CHARACTER*(*)  STR2
      character(len=8) :: fmt, xx11
      integer :: l, nitter
      fmt = '(I4.4)'
      write (xx11,fmt) nitter
      NPASTE=trim(xx11)//STR2
      END


      CHARACTER*(*) FUNCTION PASTEB( STR1, STR2 )

C CONCATENATES THE STRINGS STR1 AND STR2 LEAVING ONLY ONE BLANK IN BETWEEN
C Written by J. Soler

      CHARACTER*(*) STR1, STR2 
      integer :: l
      CHARACTER*1 BLANK
      DATA BLANK /' '/
      DO 10 L = LEN( STR1 ), 1, -1
         IF (STR1(L:L) .NE. ' ') GOTO 20
   10 CONTINUE
   20 PASTEB = STR1(1:L)//BLANK
      PASTEB = PASTEB(1:L+1)//STR2
      END




