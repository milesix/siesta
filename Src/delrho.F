      SUBROUTINE DELRHO( no, na, nspin, maxnh, nuotot, maxo, maxuo,
     &                   Gamma, indxuo, xijo, maxspn, eo,
     &                   tol, eigtol, nk, kpoint,
     &                   wk, Qo, H, S, Hper, Oper, numh, listh,
     &                   listhptr, ef, temp, Rhoper, Erhoper, 
     &                   iscf, usePsi)
C **********************************************************************
C Linres: Subroutine to calculate the change in the elements of the Density Matrix
C Written by J. Junquera and J. M. Alonso Pruneda, Novembre 1999.
C adapted to Siesta 3.x by LR, summer '15
C *********************************************************************

      use precision,  only: dp
      use parallel,     only : Node, Nodes
      use parallelsubs, only : GlobalToLocalOrb, GetNodeOrbs
      use densematrix, only: Haux, Saux, psi
      use alloc
#ifdef MPI
      use mpi_siesta
      use m_diag_option, only: ParallelOverK
#endif

      implicit none
C     EO : Eigenvalues
C     QO : Occpuaptions of Eigenstates
C     NK:  Number of K- points
      integer :: no, na, nspin, maxnh, nuotot, maxo, maxuo, maxspn, nk, 
     &           iscf, listh(maxnh), numh(*), listhptr(*), indxuo(no)
      real(dp) :: tol, ef, temp, eigtol
      real(dp) :: kpoint(3,nk), xijo(3,maxnh), wk(nk),
     &            EO(maxo,maxspn,NK), QO(maxo,maxspn,NK),
     &            Hper(maxnh,3,nspin), Oper(maxnh,3),
     &            Rhoper(maxnh,maxspn,3), Erhoper(maxnh,maxspn,3),
     &            H(maxnh,nspin), S(maxnh)

      logical :: Gamma, UsePsi

      integer :: ix, ispin, io, iuo, nhs, npsi, nuo
      integer,pointer :: muo(:)

      external delrhog, delrhok
#ifdef MPI
      external delrhokp
#endif

C Start time counter ---------------------------------------------------
      call timer( 'delrho', 1)

C Stop if non-collinear spins
      if ( nspin > 2 ) then
         call die('delrho: Only non-polarized or polarized
     &calculations are implemented.')
      end if

C Get Node number and calculate local orbital range
#ifdef MPI
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)
#else
      nuo = nuotot
#endif

      nullify(muo)
      call re_alloc( muo,  1, nuo,  'muo',  'delrho' )
C Check indxuo (as for diagon subrotines)
      do iuo = 1,nuo
        muo(iuo) = 0
      enddo
      do io = 1,no
        iuo = indxuo(io)
        if (iuo.le.0 .or. iuo.gt.nuotot) then
          if (Node.eq.0)
     $        write(6,*) 'delrho: ERROR: invalid index: io, indxuo =',
     .                     io, indxuo(io)
          call die()
        endif
        call GlobalToLocalOrb(indxuo(io),Node,Nodes,iuo)
        if (iuo.gt.0) then
          muo(iuo) = muo(iuo) + 1
        endif
      enddo
      do iuo = 1,nuo
        if (muo(iuo) .ne. muo(1)) then
          if (Node.eq.0)
     $      write(6,'(/,2a,3i6)') 'delrho: ERROR: inconsistent indxuo',
     .             '. iuo, muo(iuo), muo(1) =', iuo, muo(iuo), muo(1)
          call die()
        endif
      enddo

C Check internal dimensions ..........................................
      if ( Gamma ) then
         nhs  = 1 ! currently the algorithm is too complex to change
         npsi = nuotot * maxuo * nspin
      else
         nhs  = 2 * nuotot * nuo
         npsi = 2 * nuotot * nuo
      end if
#ifdef MPI
      nhs  = 2 * nuotot * nuotot
      npsi = 2 * nuotot * nuotot
#endif
      
      call re_alloc( psi,  1, npsi, 'psi',  'densematrix' )
      call re_alloc( Haux, 1, nhs, 'Haux', 'densematrix' )
      call re_alloc( Saux, 1, nhs, 'Saux', 'densematrix' )

! Call DFPT solvers...

      if ( Gamma ) then
         do ix = 1,3
            do ispin = 1,nspin
               call delrhog(nuotot, nuo, nspin, no,
     &                 eo(:,ispin,1), tol,
     &                 eigtol, Qo, Hper(:,ix,ispin),
     &                 Oper(:,ix), maxnh, numh, listh, listhptr,
     &                 ef,temp,Rhoper(:,ispin,ix),
     &                 Erhoper(:,ispin,ix),psi,iscf)
            enddo
         enddo  

#ifdef MPI
      else if (Nodes .gt. 1) then
        if (.not. ParallelOverK) then
          call die('DELRHOK: for MPI calculations and k-points
     & activate the ParallelOverK flag in the fdf file.')
        endif

        call delrhokp(no, nuo, maxo, maxspn, nspin, eo, eigtol, indxuo,
     &               nk, kpoint, wk, eo,xijo, Qo, H, S, Hper, Oper,
     &               maxnh, numh, listh, listhptr, ef, temp, Rhoper,
     &               Erhoper, Haux, Saux, psi, nuotot, iscf, UsePsi)
#endif
      else
         call delrhok(no, nuo, maxo, maxspn, nspin, eo, eigtol, indxuo, 
     &        nk, kpoint, wk, eo,xijo, Qo, H, S, Hper, Oper,
     &        maxnh, numh, listh, listhptr, ef, temp,
     &        Rhoper, Erhoper, Haux, Saux, psi, iscf, UsePsi)
      end if 
      call timer( 'delrho', 2)
      
      return
      END 
