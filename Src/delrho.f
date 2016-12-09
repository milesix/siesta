      SUBROUTINE DELRHO( no, na, nspin, maxnh, nuo, maxo, GAMMA,  
     &                   indxuo,xijo,maxspn,EO,tol, eigtol, nk, kpoint,
     &                   wk, Qo, H, S, Hper, Oper, numh, listh, 
     &                   listhptr, ef, temp, Rhoper, Erhoper, iscf)
C **********************************************************************
C Linres: Subroutine to calculate the change in the elements of the Density Matrix
C Written by J. Junquera and J. M. Alonso Pruneda, Novembre 1999.
C adapted to Siesta 3.x by LR, summer '15
C *********************************************************************

      use precision,  only: dp,sp
      use densematrix,    only: psi
      use alloc

      implicit none
C     EO : Eigenvalues
C     QO : Occpuaptions of Eigenstates
C     NK:  Number of K- points
      integer :: no, na, nspin, maxnh, nuo, maxuo, maxspn, nk, 
     &           ispin, ix, listh(maxnh), numh(*), listhptr(*),
     &           i,j,ind, indxuo(no),iscf, maxo
      real(dp) :: tol,EO(nuo,maxspn,NK), QO(nuo,maxspn,NK),
     &            Hper(maxnh,3,nspin), Oper(maxnh,3),ef, temp, eigtol,
     &            Rhoper(maxnh,maxspn,3), Erhoper(maxnh,maxspn,3),
     &            H(maxnh,nspin), S(maxnh)
      real(dp) :: kpoint(3,nk), xijo(3,maxnh), wk(nk)

!      real(dp), pointer :: Haux(:), Saux(:)

      logical :: gamma, frstme

      data frstme /.TRUE./
C Start time counter ---------------------------------------------------
      print *, 'DEBUG TRACK: In DELRHO'

!      if (nspin.le.2 .and. gamma) then
!        nhs  = nuo * nuo
!      elseif (nspin.le.2 .and. .not.gamma) then
!        nhs  = 2 * maxo * nuo
!      endif
      
!      call re_alloc(Haux,1,nhs, 'Haux', 'delrho')
!      call re_alloc(Saux,1,nhs, 'Saux', 'delrho')
     
!	print*,'Opr en delrho',Oper(:,1)

        print*,'iscf',iscf
!        print*,'Hper en delrho (antes de delrhok-delrhog)'
!        print*,'Hper(:,1,1)',Hper(:,1,1)
!        print*,'Hper(:,2,1)',Hper(:,2,1)
!        print*,'Hper(:,3,1)',Hper(:,3,1)

	if(nspin.le.2 .and. GAMMA) then
        print *, 'GAMMA POINT CALCULATION'
        do ix = 1,3
         do ispin = 1,nspin

          call delrhog(nuo, nspin, no, eo(:,ispin,1),tol,
     &                 eigtol, Qo, Hper(:,ix,ispin),
     &                 Oper(:,ix), maxnh, numh, listh, listhptr,
     &                 ix,ef,temp,Rhoper(:,ispin,ix),
     &                 Erhoper(:,ispin,ix),psi)
         enddo
        enddo  
      endif

      if(.not. GAMMA) then
        print *, 'K-point SAMPLING CALCULATION'
        call delrhok(no, nuo, maxo, maxspn, nspin, eo, eigtol, indxuo,  
     &               nk, kpoint, wk, eo,xijo, Qo, H, S, Hper, Oper,
     &               maxnh, numh, listh, listhptr, ef, temp, Rhoper, 
     &               Erhoper, psi)
        print *, 'DEBUG TRACK: just done with delrhok'

      endif 

!        print*,'iscf',iscf
!        print*,'rhoper despues de delrhok'
!        print*,'rhoper(:,1,1)',rhoper(:,1,1)
!        print*,'rhoper(:,1,2)',rhoper(:,1,2)
!        print*,'rhoper(:,1,3)',rhoper(:,1,3)




!	print*,'DB: stop en delrho'
!	stop

      frstme = .false.
!      call de_alloc( Haux, 'Haux', 'delrho')
!      call de_alloc( Saux, 'Saux', 'delrho')

      return
      END 
