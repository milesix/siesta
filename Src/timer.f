      SUBROUTINE TIMER (PROG,IOPT)

C  FINDS AND PRINTS THE CPU TIME SPENT IN DIFFERENT ROUTINES AND/OR
C   PROGRAM SECTIONS. IT MUST BE CALLED WITH IOPT=0 AT THE BEGINNING
C   OF EACH ROUTINE AND WITH IOPT=1 AT THE END OF IT.
C  ARGUMENTS:
C    PROG: INPUT ARBITRARY NAME FOR THE ROUTINE AND/OR PROGRAM SECTION
C    IOPT: INPUT OPTION PARAMETER:
C      IOPT = 0  => SET ALL TIMES TO ZERO (OPTIONAL)
C      IOPT = 1  => BEGIN COUNTING TIME FOR A ROUTINE
C      IOPT = 2  => STOP  COUNTING TIME FOR A ROUTINE
C      IOPT = 3  => PRINT TIME FOR A ROUTINE OR FOR ALL (IF PROG='ALL')
C  ROUTINE TIMES INCLUDE THAT SPENT IN THE ROUTINES THEY CALL
C  WRITTEN BY J.SOLER (JSOLER AT EMDUAM11) DEC/90

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=500,ZERO=0.D0,HUNDRD=100.D0,TIMMIN=1.D-6)
      DIMENSION TIME1(NMAX),TIMET(NMAX)
      INTEGER NCALLS(NMAX)
      CHARACTER*10 PROGS(NMAX),PROG*(*)
      SAVE PROGS,NPROGS,TIME0,TIME1,TIMET,NCALLS
      DATA NPROGS,TIME0 /0,ZERO/

      CALL CPUTIM (TIME)
      IF (IOPT.EQ.0) THEN
         NPROGS=0
         TIME0=TIME
      ELSEIF (IOPT.EQ.1 .OR. IOPT.EQ.2) THEN
         DO 10 IPROG=1,NPROGS
            IF (PROGS(IPROG).EQ.PROG) GO TO 20
   10    CONTINUE
            NPROGS=NPROGS+1
            IF (NPROGS.GT.NMAX) THEN
               WRITE (6,*) 'timer: NMAX IS SATURATED. PROG = ',PROG
               RETURN
            ENDIF
            IPROG=NPROGS
            PROGS(IPROG)=PROG
            NCALLS(IPROG)=0
            TIMET(IPROG)=ZERO
   20    CONTINUE
         IF (IOPT.EQ.1) THEN
            NCALLS(IPROG)=NCALLS(IPROG)+1
            TIME1(IPROG)=TIME
         ELSE
            TIMET(IPROG)=TIMET(IPROG)+TIME-TIME1(IPROG)
         ENDIF
      ELSEIF (IOPT.EQ.3) THEN
         TIMTOT=TIME-TIME0
         IF (TIMTOT.LT.TIMMIN) RETURN
         IF (PROG.EQ.'ALL' .OR. PROG.EQ.'all') THEN
           WRITE (6,'(/,A)') 'timer: CPU execution times:'
           WRITE (6,'(A,2X,A10,A9,2A12,A9)') 'timer:',
     .       'Routine   ', 'Calls', 'Time/call', 'Tot.time', '%'
           DO 40 IPROG=1,NPROGS
             AVGTME=TIMET(IPROG)/NCALLS(IPROG)
             FRACTN=TIMET(IPROG)/TIMTOT*HUNDRD
             WRITE(6,'(A,2X,A10,I9,2F12.3,F9.2)') 'timer:',
     .         PROGS(IPROG),NCALLS(IPROG),AVGTME,TIMET(IPROG),FRACTN
   40      CONTINUE
           WRITE(6,*) ' '
         ELSE
           DO 50 IPROG=1,NPROGS
             IF (PROGS(IPROG).NE.PROG) GOTO 50
               FRACTN=TIMET(IPROG)/TIMTOT*HUNDRD
               WRITE(6,'(A,A10,I6,F12.3,F7.2)')
     .          'timer: Routine,Calls,Time,% = ',
     .           PROGS(IPROG),NCALLS(IPROG),TIMET(IPROG),FRACTN
   50      CONTINUE
         ENDIF
      ELSE
         WRITE(6,*) 'timer: INVALID OPTION IOPT =',IOPT
      ENDIF
      END
