      SUBROUTINE DELRHO( no, na, nspin, maxnh, no_l, no_u, Gamma,  
     &                   indxuo,xijo,maxspn,EO,tol, eigtol, nk, kpoint,
     &                   wk, Qo, H, S, Hper, Oper, numh, listh, 
     &                   listhptr, ef, temp, Rhoper, Erhoper, iscf)
C **********************************************************************
C Linres: Subroutine to calculate the change in the elements of the Density Matrix
C Written by J. Junquera and J. M. Alonso Pruneda, Novembre 1999.
C adapted to Siesta 3.x by LR, summer '15
C *********************************************************************

      use precision,  only: dp
      use densematrix, only: Haux, Saux, psi
      use alloc

      implicit none
C     EO : Eigenvalues
C     QO : Occpuaptions of Eigenstates
C     NK:  Number of K- points
      integer :: no, na, nspin, maxnh, no_l, maxuo, maxspn, nk, 
     &           ispin, ix, listh(maxnh), numh(*), listhptr(*),
     &           i,j,ind, indxuo(no),iscf, no_u
      real(dp) :: tol,EO(no_l,maxspn,NK), QO(no_l,maxspn,NK),
     &            Hper(maxnh,3,nspin), Oper(maxnh,3),ef, temp, eigtol,
     &            Rhoper(maxnh,maxspn,3), Erhoper(maxnh,maxspn,3),
     &            H(maxnh,nspin), S(maxnh)
      real(dp) :: kpoint(3,nk), xijo(3,maxnh), wk(nk)

      integer :: nhs, npsi
      logical :: gamma

C Start time counter ---------------------------------------------------
      call timer( 'delrho', 1)

      if ( nspin > 2 ) then
         call die('delrho: Only non-polarized or polarized
     &calculations are implemented.')
      end if
      
      if ( Gamma ) then
         nhs  = 1 ! currently the algorithm is too complex to change
         npsi = no_l * no_u
      else
         nhs  = 2 * no_l * no_u
         npsi = 2 * no_l * no_u
      end if

      call re_alloc(Haux, 1, nhs, 'Haux','densematrix',shrink=.false.)
      call re_alloc(Saux, 1, nhs, 'Saux','densematrix',shrink=.false.)
      call re_alloc(psi, 1, npsi, 'psi','densematrix',shrink=.false.)
      
      if ( Gamma ) then
         do ix = 1,3
            do ispin = 1,nspin
               call delrhog(no_l, nspin, no, eo(:,ispin,1),tol,
     &              eigtol, Qo, Hper(:,ix,ispin),
     &              Oper(:,ix), maxnh, numh, listh, listhptr,
     &              ef,temp,Rhoper(:,ispin,ix),
     &              Erhoper(:,ispin,ix),psi,iscf)
            enddo
         enddo  
      else
         call delrhok(no, no_l, no_u, maxspn, nspin, eo, eigtol, indxuo,  
     &        nk, kpoint, wk, eo,xijo, Qo, H, S, Hper, Oper,
     &        maxnh, numh, listh, listhptr, ef, temp,
     &        Haux, Saux,
     &        Rhoper, Erhoper, psi)
      end if 
      call timer( 'delrho', 2)
      
      return
      END 
