      subroutine delrhok(no, no_l, no_u, maxspn, nspin, eval, eigtol,
     &                   indxuo, nk, kpoint, wk, eo,  
     &                   xij, Qo, H, S, Hper, Oper, maxnh, numh, 
     &                   listh, listhptr, ef, T,
     &                   Rhoper, Erhoper, 
     &                   Haux, Saux, psi)
C **********************************************************************
C FINDS THE CHANGE IN DENSITY MATRIX ELEMENTS DUE TO DISPLACEMENTS OF
C THE ATOMS
C K - SAMPLING VERSION.
C CODED BY J. JUNQUERA AND J. M. ALONSO PRUNEDA. Dec '98
C FOR SIESTA 3.X L. Riches, SUMMER '15
C CHECKED AND CORRECTED BY S. ILLERA APRIL '16
C **********************INPUT*******************************************
C INTEGER NO		     :Number of orbitals in the supercell
C INTEGER NO_L                :Number of basis orbitals in unit cell
C INTEGER NO_U		     :Number of orbitals
C INTEGER MAXSPN	     :Maximum number of differents spin polarizations
C INTEGER NSPIN              :Spin polarization
C REAL*8 EVAL(NO_U)           :Eigenvalues of non-perturbated Hamiltonian
C REAL*8  EIGTOL             :Tolerance to assume degenerate energy levels
C INTEGER INDXUO(NO)	     :Index of equivalent orbital in unit cell
C INTEGER NK		     :Number of kpoints
C REAL*8 KPOINT(3,NK)        :k point vectors 
C REAL*8 WK(NK)		     :k points weight		
C REAL*8 EO(NO_U,MAXSPN,NK)  :Eigenvalues
C REAL*8 XIJ(3,MAXNH)        :Vectors between orbital centers
C REAL*8 QO(NO_U)             :Occupations of unpertubed eigenstates
C REAL*8 H(MAXNH,NSPIN)	     :Hamiltonian in sparse format
C REAL*8 S(MAXNH)            : Overlap in sparse format
C REAL*8  HPER(MAXNH)        :Matrix elements of the perturbated 
C                             Hamiltonian
C REAL*8  OPER(MAXNH,3)      :Matrix elements of the perturbated 
C                             Overlap
C INTEGER MAXNH              :First dimension of listh
C INTEGER NUMH(NO_L)          :Number of nonzero density matrix elements
C                             for each matrix row
C INTEGER LISTH(MAXNH)       :Nonzero density matrix element column
C                             indexes
C INTEGER LISTHPTR(NO_L)      :Pointer to each row (-1) of the
C                             density matrix
C REAL*8  EF                  :Fermi level
C REAL*8  T                   :Temperature
C REAL*8 Haux(2,NO_U,NO_L)    :Auxiliary space for the Hamiltonian
C REAL*8 Saux(2,NO_U,NO_L)    :Auxiliary space for the overlap matrix
C REAL*8 PSI(2,NO_U,NO_L)     :Auxiliary space for the eigenvectors
C ******************  OUTPUT  ******************************************
C REAL*8  RHOPER(MAXNH,3,MAXSPN)  :Matrix elements of the perturbated DM
C REAL*8  ERHOPER(MAXNH,3,MAXSPN) :Matrix elements of the perturbated 
C                                  Energy Density Matrix
C **********************************************************************


      use precision,      only : dp
      use alloc
#ifdef _OPENMP
      use omp_lib
#endif

      implicit none

      integer ::  no, no_l, iscf, 
     &            ispin, maxnh, nspin, numh(*),
     &            listh(maxnh), listhptr(*), nk,indxuo(no) 
      real(dp) :: eval(no_u,nspin,nk), Qo(no_u,maxspn,nk),
     &            Hper(maxnh,3,nspin), eigtol, S(maxnh),
     &            Oper(maxnh,3), ef, T,Rhoper(maxnh,nspin,3),
     &            Erhoper(maxnh,nspin,3), H(maxnh,nspin),
     &            psi(2,no_u,no_l), kpoint(3,nk), xij(3,maxnh),
     &     eo(no_u,maxspn,nk),wk(nk)
      real(dp), target :: Haux(2,no_u,no_l)
      real(dp) :: Saux(2,no_u,no_l)

C     Internal Variables
      integer :: deg, N, io, numb(no_u), k,j, ik, ix,
     &           jden, kden, jo, mu, nu, ind, indmn, i, indbetas,
     &           ialpha, ibeta, jbeta,indi,indj,ierror,indden,
     &           iuo, juo, ind2, nbands, no_u, maxspn 
      real(dp) :: def(3), A(3), B(3), 
     &            dQo(no_u,3), aux(no), qio, 
     &            psiper(2,no_u,no_u), prod1(2), prod2(2), kXij, 
     &            prod3(2), prod4(2), eio, ejo, dStepF, evper(2,3,no_u),
     &            pipj1, pipj2,dpipj1,dpipj2,qde

      real(dp), pointer :: auxcoef(:,:), rotaux(:,:)
      real(dp), pointer :: Psiden(:,:,:) => null()
      real(dp), pointer :: Hauxden(:,:,:) => null()
      real(dp), pointer :: Sauxden(:,:,:) => null()
      real(dp) :: ckXij, skXij, pert

      call timer('delrhok',1)

#ifdef _OPENMP
!$OMP parallel default(shared)

      ! Check that we can use the appropiate partition in the
      ! coefficient change routine
      ! This does put a constraint on the number of OMP threads
      ! we can use. But for now I think this is good.
      ! Else, we force this to be sequential
!$OMP single
      if ( omp_get_num_threads() > no_l ) then
         write(*,*) 'Reduce OMP_NUM_THREADS to (maximum): ',no_l
         call die('delrhok: too many OMP_NUM_THREADS, compared to no_l')
      end if
!$OMP end single

!$OMP end parallel
#endif
      
C     Initialize variables --- 
      A(1:3) = 0.0_dp
      B(1:3) = 0.0_dp

C     Find eigenvectors (where stored for only one k-point) ---------
      do ik = 1,nk     ! for each k-point 
       do ispin = 1, nspin !for each spin component

!$OMP parallel default(shared), private(iuo,j,ind,jo,juo),
!$OMP&private(kXij,ckXij,skXij)
          
!$OMP do
        do iuo = 1, no_l 
         Haux(:,:,iuo) = 0._dp
         Saux(:,:,iuo) = 0._dp
         do j = 1, numh(iuo)
          ind = listhptr(iuo) + j
          jo = listh(ind)
          juo = indxuo(jo)
C         calculates the phases k*r_ij -------------------------------
          kXij = kpoint(1,ik) * Xij(1,ind) +
     .          kpoint(2,ik) * Xij(2,ind) +
     .          kpoint(3,ik) * Xij(3,ind) 
          ckXij = cos(kXij)  
          skXij = sin(kXij)
C         Calculates the hamiltonian and the overlap in k space ------
C         H(k) = Sum(R) exp(i*k*R) * H(R) ----------------------------
          Saux(1,juo,iuo) = Saux(1,juo,iuo) + S(ind)*ckxij
          Saux(2,juo,iuo) = Saux(2,juo,iuo) - S(ind)*skxij
          Haux(1,juo,iuo) = Haux(1,juo,iuo) + H(ind,ispin)*ckxij
          Haux(2,juo,iuo) = Haux(2,juo,iuo) - H(ind,ispin)*skxij
         enddo
        enddo
!$OMP end do
        
C     Symmetrize H & S ---------------------------------------------
!$OMP do schedule(static,1)
        do iuo = 1, no_u
         do juo = 1, iuo-1
          Saux(1,juo,iuo) = 0.5_dp * ( Saux(1,juo,iuo) +
     .                               Saux(1,iuo,juo) )
          Saux(1,iuo,juo) = Saux(1,juo,iuo)
          Saux(2,juo,iuo) = 0.5_dp * ( Saux(2,juo,iuo) -
     .                               Saux(2,iuo,juo) )
          Saux(2,iuo,juo) = - Saux(2,juo,iuo)
          Haux(1,juo,iuo) = 0.5_dp * ( Haux(1,juo,iuo) +
     .                               Haux(1,iuo,juo) )
          Haux(1,iuo,juo) = Haux(1,juo,iuo)
          Haux(2,juo,iuo) = 0.5_dp * ( Haux(2,juo,iuo) -
     .                               Haux(2,iuo,juo) )
          Haux(2,iuo,juo) = - Haux(2,juo,iuo)
         enddo
         Saux(2,iuo,iuo) = 0.0_dp
         Haux(2,iuo,iuo) = 0.0_dp
        enddo
!$OMP end do nowait

!$OMP end parallel
C       Diagonalize for each k-point --------------------------------
        call cdiag( Haux, Saux, no_u, no_l, no_u, eo(1,ispin,ik),
     &                psi, no_u, 1, ierror )

        if (ierror.ne.0) then
          call die('DELRHOK: Terminating due to failed 
     &                                        k-diagonalisation')
        endif
C       loop on spatial xyz -----------------------------------------
        do ix = 1,3
C        Find out the degenerated subespaces ------------------------
         nbands = no_u
         numb(1:no_u) = 0
         N = 0
         do  io = 1, nbands
          eio = eo(io,ispin,ik)
          if(io .lt. nbands) then
           ejo = eo(io+1,ispin,ik)
          else
           ejo = 1.0e7_dp
          endif 
          ! building the dH_nn' = <psi_in|dH|psi_in'> 
          ! where n refers to the degenerate space.
          if(abs(eio-ejo) .lt. eigtol) then
            N = N + 1
          else 
            numb(io) = N + 1
            N = 0
            if(numb(io).gt.1) then !last degenerated state
              call re_alloc(Psiden, 1, 2, 1, numb(io), 1, numb(io),
     &                       'Psiden', 'delrhok')
              call re_alloc(Hauxden, 1, 2, 1, numb(io), 1, numb(io),
     &                       'Hauxden', 'delrhok')
              call re_alloc(Sauxden, 1, 2, 1, numb(io), 1, numb(io),
     &                       'Sauxden', 'delrhok')
              ! Re-use memory for rotation.
              rotaux => Sauxden(:,:,1)
              
              ! Initialize to 0
              Sauxden = 0._dp

              ! TODO transpose this construction, i.e. 
              ! the access is sub-optimal
              do j = io-numb(io)+1, io
                jden = j - io + numb(io)
                numb(j) = numb(io)
                do k = io-numb(io)+1, io
                  kden = k - io + numb(io)
                  if ( jden == kden ) Sauxden(1,jden,kden) = 1.0_dp
                  do mu=1,no_l
                    do nu = 1,numh(mu)
                      indmn=listhptr(mu) + nu
                      jo = listh(indmn)
                      juo = indxuo(jo)
                      kXij = kpoint(1,ik) * Xij(1,indmn) +
     .                       kpoint(2,ik) * Xij(2,indmn) +
     .                       kpoint(3,ik) * Xij(3,indmn) 
                      ckXij = cos(kXij)  
                      skXij = sin(kXij) 
                      pipj1 = psi(1,mu,j) * psi(1,juo,k) +
     &                        psi(2,mu,j) * psi(2,juo,k)
                      pipj2 = psi(1,mu,j) * psi(2,juo,k) -
     &                        psi(2,mu,j) * psi(1,juo,k)
                      pert = Hper(indmn,ix,ispin) -
     &                       eio * Oper(indmn,ix)
                      Hauxden(1,jden,kden) = Hauxden(1,jden,kden) +
     &                        (pipj1 * ckXij - pipj2 * skXij) * pert
                      Hauxden(2,jden,kden) = Hauxden(2,jden,kden) +
     &                        (pipj1 * skXij + pipj2 * ckXij) * pert   
                    enddo
                  enddo
                enddo
              enddo !j loop

              ! Haux is a dummy argument, we don't need the eigenvalues
              ! and here Haux is always big enough to hold them all!
              call cdiag( Hauxden, Sauxden, numb(io), 
     &             numb(io), numb(io),
     .             Haux, psiden,numb(io),1,ierror )
              if (ierror.ne.0) then
               call die('DELRHOK: Terminating due to failed
     &                                deg- diagonalisation')
              endif

C             Rotate the coefficients inside the subspace ---------- 
              do mu = 1,no_l
               do ialpha = 1, numb(io)
                 ! This is actually Sauxden(:,:,1)
                 rotaux(:,ialpha) = 0._dp
                 do ibeta = 1, numb(io)
                  rotaux(1,ialpha) = rotaux(1,ialpha) + 
     &                  psiden(1,ibeta,ialpha) *
     &                  psi(1,mu,io-numb(io)+ibeta) -
     &                  psiden(2,ibeta,ialpha) *
     &                  psi(2,mu,io-numb(io)+ibeta)
                  rotaux(2,ialpha) = rotaux(2,ialpha) +
     &                 psiden(1,ibeta,ialpha) *
     &                 psi(2,mu,io-numb(io)+ibeta) +
     &                 psiden(2,ibeta,ialpha) *
     &                 psi(1,mu,io-numb(io)+ibeta)
                end do
               end do
               do ialpha = 1, numb(io)
                  psi(1,mu,io-numb(io)+ialpha) = rotaux(1,ialpha)
                  psi(2,mu,io-numb(io)+ialpha) = rotaux(2,ialpha)
               end do
              enddo !mu loop

              call de_alloc( Psiden, 'Psiden', 'delrhok' )
              call de_alloc( Hauxden, 'Hauxden', 'delrhok' )
              call de_alloc( Sauxden, 'Sauxden', 'delrhok' )     
            endif !numb() last degenerated state
          endif !lt tol (degenerated)
         enddo !bands io loop

!$OMP parallel default(shared),
!$OMP&private(io,qio,eio,mu,juo,iuo,nu,ind,jo),
!$OMP&private(kXij,ckXij,skXij,auxcoef),
!$OMP&private(ialpha,prod1,prod2,prod3,prod4),
!$OMP&private(indmn,ibeta,ejo),
!$OMP&reduction(+:A,B)

! OMP remark:
! evper is accessed without possibility of data-race
! psipero is accessed without possibility of data-race
! dQo is accessed without possibility of data-race

#ifdef _OPENMP
         io = omp_get_thread_num( ) + 1
         auxcoef => Haux(:,:,io)
#else
         auxcoef => Haux(:,:,1)
#endif

!$OMP workshare
         evper(1:2,ix,1:no_u) = 0.0_dp
!$OMP end workshare

C        Compute the change in the coefficients ---------------------
!$OMP do schedule(static,1)
         do io = 1, nbands
          qio = Qo(io,ispin,ik)
          if (qio > 1.0e-6_dp) then !is it occupied?

          eio = eo(io,ispin,ik)

          auxcoef(:,:) = 0.0_dp
          
          do ialpha = 1, no_l
           prod1(1:2) = 0.0_dp
           prod2(1:2) = 0.0_dp
            do jbeta = 1, numh(ialpha)
            indmn = listhptr(ialpha) + jbeta
            ibeta = indxuo( listh(indmn) )
            kXij = kpoint(1,ik) * Xij(1,indmn) +
     &             kpoint(2,ik) * Xij(2,indmn) +
     &             kpoint(3,ik) * Xij(3,indmn)
            ckXij = cos(kXij)
            skXij = sin(kXij)
            prod1(1) = prod1(1) + Hper(indmn,ix,ispin) *
     &                (psi(1,ibeta,io) * ckXij - 
     &                 psi(2,ibeta,io) * skXij)
            prod1(2) = prod1(2) + Hper(indmn,ix,ispin) *
     &                (psi(1,ibeta,io) * skXij + 
     &                 psi(2,ibeta,io) * ckXij)
            prod2(1) = prod2(1) + Oper(indmn,ix) *
     &                (psi(1,ibeta,io) * ckXij - 
     &                psi(2,ibeta,io) * skXij)
            prod2(2) = prod2(2) + Oper(indmn,ix) *
     .                (psi(1,ibeta,io) * skXij + 
     &                 psi(2,ibeta,io) * ckXij) 
            enddo !jbeta

           prod3(1) = prod1(1) - eio * prod2(1)
           prod3(2) = prod1(2) - eio * prod2(2)
           prod4(1) = prod2(1) / 2.0_dp
           prod4(2) = prod2(2) / 2.0_dp
            do j = 1,nbands
            ejo = eo(j,ispin,ik)
             if(abs(eio-ejo) .gt. eigtol) then
             auxcoef(1,j) = auxcoef(1,j) + (1.0_dp/(eio-ejo)) *
     .            (psi(1,ialpha,j)*prod3(1) + psi(2,ialpha,j)*prod3(2))
             auxcoef(2,j) = auxcoef(2,j) + (1.0_dp/(eio-ejo)) *
     .            (psi(1,ialpha,j)*prod3(2) - psi(2,ialpha,j)*prod3(1))
             else
             auxcoef(1,j) = auxcoef(1,j) -
     .             (psi(1,ialpha,j)*prod4(1) + psi(2,ialpha,j)*prod4(2))
             auxcoef(2,j) = auxcoef(2,J) -
     .             (psi(1,ialpha,j)*prod4(2) - psi(2,ialpha,j)*prod4(1))
             endif
            enddo ! jbands
           evper(1,ix,io) = evper(1,ix,io) +
     .            psi(1,ialpha,io)*prod3(1) + psi(2,ialpha,io)*prod3(2)
           evper(2,ix,io) = evper(2,ix,io) +
     .            psi(1,ialpha,io)*prod3(2) - psi(2,ialpha,io)*prod3(1)
           enddo !ialpha      
          
           do ialpha = 1,no_u
           psiper(1:2,ialpha,io) = 0.0_dp
            do j = 1, nbands
            psiper(1,ialpha,io) = psiper(1,ialpha,io) +
     .                           auxcoef(1,j)*psi(1,ialpha,j) -
     .                           auxcoef(2,j)*psi(2,ialpha,j)
            psiper(2,ialpha,io) = psiper(2,ialpha,io) +
     .                           auxcoef(1,j)*psi(2,ialpha,j) +
     .                           auxcoef(2,j)*psi(1,ialpha,j)
           enddo
          enddo
         
***DEF***
          dQo(io,ix) = wk(ik)*evper(1,ix,io)*dstepf((eio-ef)/T)/
     .                 (T*nspin+1.0e-12)
          A(ix) = A(ix) + wk(ik)*evper(1,ix,io)*dstepf((eio-ef)/T)
          B(ix) = B(ix) + wk(ik)*dstepf((eio-ef)/T)
*********
         endif ! occupancy
         enddo !io bands loop
!$OMP end do
         
C Compute the change in density matrix elements -------------

!$OMP workshare
         Saux = 0.0_dp
         Haux = 0.0_dp
!$OMP end workshare

         do io = 1, nbands
          qio = qO(io,ispin,ik)
          if ( qio < 1.0e-6_dp ) cycle
          eio = eo(io,ispin,ik) * qio
          qde = qio * evper(1,ix,io) + eio * dQo(io,ix)
!$OMP do collapse(2)
          do mu = 1, no_u
           do juo = 1, no_u
            dpipj1 = (psiper(1,mu,io) * psi(1,juo,io) +
     .              psiper(2,mu,io) * psi(2,juo,io) +
     .              psi(1,mu,io) * psiper(1,juo,io) +
     .              psi(2,mu,io) * psiper(2,juo,io))
            dpipj2 = (psiper(1,mu,io) * psi(2,juo,io) -
     .              psiper(2,mu,io) * psi(1,juo,io) +
     .              psi(1,mu,io) * psiper(2,juo,io) -
     .              psi(2,mu,io) * psiper(1,juo,io))
            pipj1 = (psi(1,mu,io)*psi(1,juo,io) +
     .              psi(2,mu,io)*psi(2,juo,io))
            pipj2 = (psi(1,mu,io)*psi(2,juo,io) -
     .              psi(2,mu,io)*psi(1,juo,io))
            Saux(1,juo,mu) = Saux(1,juo,mu) + dpipj1 * qio
     .                     + pipj1 * dQo(io,ix)
            Saux(2,juo,mu) = Saux(2,juo,mu) + dpipj2 * qio
     .                     + pipj2 * dQo(IO,ix)
            Haux(1,juo,mu) = Haux(1,juo,mu) + dpipj1 * eio
     .                     + qde * pipj1
            Haux(2,juo,mu) = Haux(2,juo,mu) + dpipj2 * eio
     .                     + qde * pipj2
!
!            Saux(1,juo,mu) = Saux(1,juo,mu) +
!     .             (psiper(1,mu,io) * psi(1,juo,io) +
!     .              psiper(2,mu,io) * psi(2,juo,io) +
!     .              psi(1,mu,io) * psiper(1,juo,io) +
!     .              psi(2,mu,io) * psiper(2,juo,io)) * qio
!     .              + (psi(1,mu,io)*psi(1,juo,io) +
!     .              psi(2,mu,io)*psi(2,juo,io))*dQo(io,ix)
!            Saux(2,juo,mu) = Saux(2,juo,mu) +
!     .             (psiper(1,mu,io) * psi(2,juo,io) -
!     .              psiper(2,mu,io) * psi(1,juo,io) +
!     .              psi(1,mu,io) * psiper(2,juo,io) -
!     .              psi(2,mu,io) * psiper(1,juo,io)) * qio
!     .              + (psi(1,mu,io)*psi(2,juo,io) -
!     .              psi(2,mu,io)*psi(1,juo,io))*dQo(IO,ix)
!             Haux(1,juo,mu) = Haux(1,juo,mu) +
!     .             (psiper(1,mu,io) * psi(1,juo,io) +
!     .              psiper(2,mu,io) * psi(2,juo,io) +
!     .              psi(1,mu,io) * psiper(1,juo,io) +
!     .              psi(2,mu,io) * psiper(2,juo,io)) * qio * eio +
!     .              qio * evper(1,ix,io) *
!     .              (psi(1,mu,io) * psi(1,juo,io) +
!     .               psi(2,mu,io) * psi(2,juo,io))
!     .            + (psi(1,mu,io)*psi(1,juo,io) +
!     .              psi(2,mu,io)*psi(2,juo,io))*eio*dQo(io,ix)
!             Haux(2,juo,mu) = Haux(2,juo,mu) +
!     .             (psiper(1,mu,io) * psi(2,juo,io) -
!     .              psiper(2,mu,io) * psi(1,juo,io) +
!     .              psi(1,mu,io) * psiper(2,juo,io) -
!     .              psi(2,mu,io) * psiper(1,juo,io)) * qio * eio +
!     .              qio * evper(1,ix,io) *
!     .              (psi(1,mu,io) * psi(2,juo,io) -
!     .               psi(2,mu,io) * psi(1,juo,io) )
!     .       + (psi(1,mu,io)*psi(2,juo,io) -
!     .          psi(2,mu,io)*psi(1,juo,io))*eio*dQo(io,ix)
           enddo
          enddo
!$OMP end do ! retain barrier
       enddo

!$OMP do
         do iuo = 1, no_l
          do nu = 1, numh(iuo)
           ind = listhptr(iuo) + nu
           jo = listh(ind)
           juo = indxuo(jo)
           kXij = kpoint(1,ik) * Xij(1,ind) +
     &            kpoint(2,ik) * Xij(2,ind) +
     &            kpoint(3,ik) * Xij(3,ind)
           ckXij = cos(kXij)
           skXij = sin(kXij)
           Rhoper(ind,ispin,ix) = Rhoper(ind,ispin,ix) +
     .                            Saux(1,juo,iuo)*ckxij -
     .                            Saux(2,juo,iuo)*skxij
           erhoper(ind,ispin,ix) = erhoper(ind,ispin,ix) +
     .                             Haux(1,juo,iuo)*ckxij -
     .                             Haux(2,juo,iuo)*skxij
          enddo
         enddo  
!$OMP end do nowait
         
!$OMP end parallel
         
        enddo !spatial coordinates
       enddo !nspin loop

      enddo ! nk loop 

      do ix = 1,3
        def(ix) = A(ix) / (B(ix) + 1.0e-12_dp)
      enddo

c    Variation of the energy level ocupation--------------------------
      if ( (def(1).gt.1.0e-6_dp) .or. (def(2).gt.1.0e-6_dp) .or.
     .                             (def(3).gt.1.0e-6_dp) ) then
        print *, 'DEBUG TRACK: delrhok, change in the fermi level!!'
        do 800 ik=1,nk
          do 770 ispin=1,nspin
            Haux=0.0_dp
            Saux=0.0_dp
            do iuo = 1, no_l
              do j = 1, numh(iuo)
                ind = listhptr(iuo) + j
                jo = listh(ind)
                juo = indxuo(jo)
C         calculates the phases k*r_ij -------------------------------
                kXij = kpoint(1,ik) * Xij(1,ind) +
     .          kpoint(2,ik) * Xij(2,ind) +
     .          kpoint(3,ik) * Xij(3,ind)
                ckXij = cos(kXij)
                skXij = sin(kXij)
C         Calculates the hamiltonian and the overlap in k space ------
C         H(k) = Sum(R) exp(i*k*R) * H(R) ----------------------------
                Saux(1,juo,iuo) = Saux(1,juo,iuo) + S(ind)*ckxij
                Saux(2,juo,iuo) = Saux(2,juo,iuo) - S(ind)*skxij!-
                Haux(1,juo,iuo) = Haux(1,juo,iuo) + H(ind,ispin)*ckxij
                Haux(2,juo,iuo) = Haux(2,juo,iuo) - H(ind,ispin)*skxij!-
              enddo
            enddo
C       Symmetrize H & S ---------------------------------------------
            do iuo = 1, no_u
              do juo = 1, iuo-1
                Saux(1,juo,iuo) = 0.5_dp * ( Saux(1,juo,iuo) +
     .                               Saux(1,iuo,juo) )
                Saux(1,iuo,juo) = Saux(1,juo,iuo)
                Saux(2,juo,iuo) = 0.5_dp * ( Saux(2,juo,iuo) -
     .                               Saux(2,iuo,juo) )
                Saux(2,iuo,juo) = - Saux(2,juo,iuo)
                Haux(1,juo,iuo) = 0.5_dp * ( Haux(1,juo,iuo) +
     .                               Haux(1,iuo,juo) )
                Haux(1,iuo,juo) = Haux(1,juo,iuo)
                Haux(2,juo,iuo) = 0.5_dp * ( Haux(2,juo,iuo) -
     .                               Haux(2,iuo,juo) )
                Haux(2,iuo,juo) = - Haux(2,juo,iuo)
              enddo
              Saux(2,iuo,iuo) = 0.0_dp
              Haux(2,iuo,iuo) = 0.0_dp
            enddo
C       Diagonalize for each k-point --------------------------------
            call cdiag( Haux, Saux, no_u, no_l, no_u, eo(:,ispin,ik),
     &                psi, no_u, 1, ierror )
C        Compute the change in density matrix elements -------------
            Saux = 0.0_dp
            Haux = 0.0_dp
            do 550 io = 1, nbands
              qio = qO(io,ispin,ik)
              if(qio .lt. 1.0e-6_dp) goto 550
                eio = eo(io,ispin,ik)
                dQo(io,1) = -wk(ik)*dstepf((eio-ef)/T)/
     .                 (T*nspin+1.0e-12)
              do mu = 1, no_u
                do juo = 1, no_u
                  Saux(1,juo,mu) = Saux(1,juo,mu) 
     .              + (psi(1,mu,io)*psi(1,juo,io) +
     .              psi(2,mu,io)*psi(2,juo,io))*dQo(io,1)
                  Saux(2,juo,mu) = Saux(2,juo,mu) 
     .              + (psi(1,mu,io)*psi(2,juo,io) -
     .              psi(2,mu,io)*psi(1,juo,io))*dQo(io,1)

                  haux(1,juo,mu) = haux(1,juo,mu) +
     .             (psi(1,mu,io)*psi(1,juo,io) +
     .              psi(2,mu,io)*psi(2,juo,io))*eio*dQo(io,1)
                  haux(2,juo,mu) = haux(2,juo,mu) +
     .             (psi(1,mu,io)*psi(2,juo,io) -
     .              psi(2,mu,io)*psi(1,juo,io))*eio*dQo(io,1)
                enddo
              enddo
550         enddo
            do iuo=1,no_l
              do j=1,numh(iuo)
                ind = listhptr(iuo) + j
                jo = listh(ind)
                juo = indxuo(jo)
                kXij = kpoint(1,ik) * Xij(1,ind) +
     .          kpoint(2,ik) * Xij(2,ind) +
     .          kpoint(3,ik) * Xij(3,ind)
                do ix=1,3
                  ckXij = cos(kXij)*def(ix)
                  skXij = sin(kXij)*def(ix)
                  Rhoper(ind,ispin,ix) = Rhoper(ind,ispin,ix) +
     .                            Saux(1,juo,iuo)*ckxij -
     .                            Saux(2,juo,iuo)*skxij
                  erhoper(ind,ispin,ix) = erhoper(ind,ispin,ix) +
     .                             Haux(1,juo,iuo)*ckxij -
     .                             Haux(2,juo,iuo)*skxij
                enddo !ix loop
              enddo
            enddo
770       enddo
800     enddo
      endif !Variation of the energy level ocupation

      call timer('delrhok',2)
      return
      end

