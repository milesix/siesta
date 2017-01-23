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

      logical :: gamma

C Start time counter ---------------------------------------------------
      call timer( 'delrho', 1)

      if(nspin.le.2 .and. GAMMA) then
        do ix = 1,3
         do ispin = 1,nspin
          call delrhog(nuo, nspin, no, eo(:,ispin,1),tol,
     &                 eigtol, Qo, Hper(:,ix,ispin),
     &                 Oper(:,ix), maxnh, numh, listh, listhptr,
     &                 ix,ef,temp,Rhoper(:,ispin,ix),
     &                 Erhoper(:,ispin,ix),psi,iscf)
         enddo
        enddo  
      elseif ((.not. GAMMA).and.(nspin.le.2)) then
        call delrhok(no, nuo, maxo, maxspn, nspin, eo, eigtol, indxuo,  
     &               nk, kpoint, wk, eo,xijo, Qo, H, S, Hper, Oper,
     &               maxnh, numh, listh, listhptr, ef, temp, Rhoper, 
     &               Erhoper, psi)
      endif 
      call timer( 'delrho', 2)

      return
      END 
