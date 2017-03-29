      subroutine delrhok(no, nuo, maxo, maxspn, nspin, eval, eigtol,
     &                   indxuo, nk, kpoint, wk, eo,  
     &                   xij, Qo, H, S, Hper, Oper, maxnh, numh, 
     &                   listh, listhptr, ef, T, Rhoper, Erhoper, 
     &                   psi)
C **********************************************************************
C FINDS THE CHANGE IN DENSITY MATRIX ELEMENTS DUE TO DISPLACEMENTS OF
C THE ATOMS
C K - SAMPLING VERSION.
C CODED BY J. JUNQUERA AND J. M. ALONSO PRUNEDA. Dec '98
C FOR SIESTA 3.X L. Riches, SUMMER '15
C CHECKED AND CORRECTED BY S. ILLERA APRIL '16
C **********************INPUT*******************************************
C INTEGER NO		     :Number of orbitals in the supercell
C INTEGER NUO                :Number of basis orbitals in unit cell
C INTEGER MAXO		     :Number of orbitals
C INTEGER MAXSPN	     :Maximum number of differents spin polarizations
C INTEGER NSPIN              :Spin polarization
C INTEGER MAXORB             :Number of orbitals
C REAL*8 EVAL(NUO)           :Eigenvalues of non-perturbated Hamiltonian
C REAL*8  EIGTOL             :Tolerance to assume degenerate energy levels
C INTEGER INDXUO(NO)	     :Index of equivalent orbital in unit cell
C INTEGER NK		     :Number of kpoints
C REAL*8 KPOINT(3,NK)        :k point vectors 
C REAL*8 WK(NK)		     :k points weight		
C REAL*8 EO(MAXO,MAXSPN,NK)  :Eigenvalues
C REAL*8 XIJ(3,MAXNH)        :Vectors between orbital centers
C REAL*8 QO(NUO)             :Occupations of unpertubed eigenstates
C REAL*8 H(MAXNH,NSPIN)	     :Hamiltonian in sparse format
C REAL*8 S(MAXNH)            : Overlap in sparse format
C REAL*8  HPER(MAXNH)        :Matrix elements of the perturbated 
C                             Hamiltonian
C REAL*8  OPER(MAXNH,3)      :Matrix elements of the perturbated 
C                             Overlap
C INTEGER MAXNH              :First dimension of listh
C INTEGER NUMH(NUO)          :Number of nonzero density matrix elements
C                             for each matrix row
C INTEGER LISTH(MAXNH)       :Nonzero density matrix element column
C                             indexes
C INTEGER LISTHPTR(NUO)      :Pointer to each row (-1) of the
C                             density matrix
C REAL*8  EF                 :Fermi level
C REAL*8  T                  :Temperature
C REAL*8 PSI(2,NUO,NUO)      :Auxiliary space for the eigenvectors
C ******************  OUTPUT  ******************************************
C REAL*8  RHOPER(MAXNH,3,MAXSPN)  :Matrix elements of the perturbated DM
C REAL*8  ERHOPER(MAXNH,3,MAXSPN) :Matrix elements of the perturbated 
C                                  Energy Density Matrix
C **********************************************************************


      use precision,      only : dp,sp
      use alloc

      implicit none

      integer ::  no, nuo, iscf, 
     &            ispin, maxnh, nspin, numh(*),
     &            listh(maxnh), listhptr(*), nk,indxuo(no) 
      real(dp) :: eval(nuo,nspin,nk), Qo(maxo,maxspn,nk),
     &            Hper(maxnh,3,nspin), eigtol, S(maxnh),
     &            Oper(maxnh,3), ef, T,Rhoper(maxnh,3,nspin),
     &            Erhoper(maxnh,nspin,3), H(maxnh,nspin),
     &            psi(2,nuo,nuo), kpoint(3,nk), xij(3,maxnh),
     &            eo(maxo,maxspn,nk),wk(nk)

C     Internal Variables
      integer :: deg, N, io, maxorb, numb(maxo), k,j, ik, ix,
     &           jden, kden, jo, mu, nu, ind, indmn, i, indbetas,
     &           ialpha, ibeta, jbeta,indi,indj,ierror,indden,
     &           iuo, juo, ind2, nbands, maxo, maxspn 
      real(dp) :: def(3), A(3), B(3), 
     &            dQo(nuo,3), aux(no),auxcoef(2,maxo), qio, 
     &            psiper(2,maxo,maxo), prod1(2), prod2(2), kXij, 
     &            prod3(2), prod4(2), eio, ejo, dStepF, evper(2,3,maxo),
     &            pipj1, pipj2,Haux(2,maxo,nuo),
     &            Saux(2,maxo,nuo) !maxo = nuotot

      real(dp), pointer :: Psiden(:,:,:),
     &                     eden(:),rotaux(:,:), Hauxden(:,:,:), 
     &                     Sauxden(:,:,:)
      real(dp) :: ckXij, skXij, pert

      call timer('delrhok',1)


C     Initialize variables --- 
      def(1:3) = 0.0_dp
      A(1:3) = 0.0_dp
      B(1:3) = 0.0_dp

C     Find eigenvectors (where stored for only one k-point) ---------
      do ik = 1,nk     ! for each k-point 
       do ispin = 1, nspin !for each spin component
        Haux = 0.0_dp
        Saux = 0.0_dp
        do iuo = 1, nuo 
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
        do iuo = 1, nuo
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
        call cdiag( Haux, Saux, nuo, nuo, nuo, eo(:,ispin,ik),
     &                psi, nuo, 1, ierror )

        if (ierror.ne.0) then
          call die('DELRHOK: Terminating due to failed 
     &                                        k-diagonalisation')
        endif
C       loop on spatial xyz -----------------------------------------
        do 215 ix = 1,3
C        Find out the degenerated subespaces ------------------------
         nbands = nuo
         numb(1:maxo) = 0
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
              nullify(Psiden,Eden,rotaux,Hauxden,Sauxden)
              call re_alloc(Psiden, 1, 2, 1, numb(io), 1, numb(io),
     &                       'Psiden', 'delrhok')
              call re_alloc(eden,1, numb(io),'eden', 'delrhok')
              call re_alloc(rotaux,1,2,1, numb(io),'rotaux', 'delrhok')
              call re_alloc(Hauxden,1,2,1, numb(io),1, numb(io),
     &                       'Hauxden', 'delrhok')
              call re_alloc(Sauxden, 1, 2, 1, numb(io), 1, numb(io),
     &                       'Sauxden', 'delrhok')

              do j = io-numb(io)+1, io
                jden = j - io + numb(io)
                numb(j) = numb(io)
                do k = io-numb(io)+1, io
                  kden = k - io + numb(io)
                  if(jden .eq. kden) then
                    Sauxden(1,jden,kden) = 1.0_dp
                    Sauxden(2,jden,kden) = 1.0_dp
                  endif
                  do mu=1,nuo
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
                      pert = Hper(indmn,ix,ispin) - !looks ok 
     &                       eio * Oper(indmn,ix)
                      Hauxden(1,jden,kden) = Hauxden(1,jden,kden) +
     &                        (pipj1 * ckXij - pipj2 * skXij) * pert
                      Hauxden(2,jden,kden) = Hauxden(2,jden,kden) +
     &                        (pipj1 * skXij + pipj2 * ckXij) * pert   
                    enddo
                  enddo
                enddo
              enddo !j loop
              call cdiag( Hauxden, Sauxden, numb(io), 
     &             numb(io), numb(io),
     .             eden, psiden,numb(io),1,ierror )
              if (ierror.ne.0) then
               call die('DELRHOK: Terminating due to failed
     &                                deg- diagonalisation')
              endif

C             Rotate the coefficients inside the subspace ---------- 
              do mu = 1,nuo
               do ialpha = 1, numb(io)
                do ibeta = 1, numb(io)
                  rotaux(1,ialpha) = rotaux(1,ialpha) +
     &                   psiden(1,ibeta,ialpha) *
     &                   psi(1,mu,io-numb(io)+ibeta) -
     &                   psiden(2,ibeta,ialpha) *
     &                   psi(2,mu,io-numb(io)+ibeta)
                  rotaux(2,ialpha) = rotaux(2,ialpha) +
     &                   psiden(1,ibeta,ialpha) *
     &                   psi(2,mu,io-numb(io)+ibeta) +
     &                   psiden(2,ibeta,ialpha) *
     &                   psi(1,mu,io-numb(io)+ibeta)

                enddo
               enddo
               do ialpha = 1, numb(io)
                psi(1,mu,io-numb(io)+ialpha) = rotaux(1,ialpha)
                psi(2,mu,io-numb(io)+ialpha) = rotaux(2,ialpha)
                rotaux(1,ialpha) = 0.0_dp
                rotaux(2,ialpha) = 0.0_dp
               enddo
              enddo !mu loop
              do jo = 1, numb(io)
               do k = 1, numb(io)
                Hauxden(1:2,jo,k) = 0.D0
                Sauxden(1:2,jo,k) = 0.D0
               enddo
              enddo

              call de_alloc( Psiden, 'Psiden', 'delrhok' )
              call de_alloc( rotaux, 'rotaux', 'delrhok' )
              call de_alloc( eden, 'eden', 'delrhok' )
              call de_alloc( Hauxden, 'Hauxden', 'delrhok' )
              call de_alloc( Sauxden, 'Sauxden', 'delrhok' )     
            endif !numb() last degenerated state
          endif !lt tol (degenerated)
         enddo !bands io loop
         evper(1:2,ix,1:maxo) = 0.0_dp

C        Compute the change in the coefficients ---------------------
         do 300 io = 1, nbands
          qio = Qo(io,ispin,ik)
          if(qio .gt. 1.0e-6_dp) then !is it occupied?

          eio = eo(io,ispin,ik)

          auxcoef(1,1:maxo) = 0.0_dp
          auxcoef(2,1:maxo) = 0.0_dp
          
           do ialpha = 1, nuo
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
            !both psis look good
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
          
           do ialpha = 1,nuo
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
          endif !occupancy
300      enddo !io bands loop

C        Compute the change in density matrix elements -------------
         Saux = 0.0_dp
         Haux = 0.0_dp

         do 370 io = 1, nbands
          qio = qO(io,ispin,ik)
          if(qio .lt. 1.0e-6_dp) goto 370
          eio = eo(io,ispin,ik)
          do mu = 1, nuo
           do juo = 1, nuo
            nu = juo
            Saux(1,juo,mu) = Saux(1,juo,mu) +
     .             (psiper(1,mu,io) * psi(1,juo,io) +
     .              psiper(2,mu,io) * psi(2,juo,io) +
     .              psi(1,mu,io) * psiper(1,juo,io) +
     .              psi(2,mu,io) * psiper(2,juo,io)) * qio
     .              + (psi(1,mu,io)*psi(1,juo,io) +
     .              psi(2,mu,io)*psi(2,juo,io))*dQo(io,ix)
            Saux(2,juo,mu) = Saux(2,juo,mu) +
     .             (psiper(1,mu,io) * psi(2,juo,io) -
     .              psiper(2,mu,io) * psi(1,juo,io) +
     .              psi(1,mu,io) * psiper(2,juo,io) -
     .              psi(2,mu,io) * psiper(1,juo,io)) * qio
     .              + (psi(1,mu,io)*psi(2,juo,io) -
     .              psi(2,mu,io)*psi(1,juo,io))*dQo(IO,ix)
             Haux(1,juo,mu) = Haux(1,juo,mu) +
     .             (psiper(1,mu,io) * psi(1,juo,io) +
     .              psiper(2,mu,io) * psi(2,juo,io) +
     .              psi(1,mu,io) * psiper(1,juo,io) +
     .              psi(2,mu,io) * psiper(2,juo,io)) * qio * eio +
     .              qio * evper(1,ix,io) *
     .              (psi(1,mu,io) * psi(1,juo,io) +
     .               psi(2,mu,io) * psi(2,juo,io))
     .            + (psi(1,mu,io)*psi(1,juo,io) +
     .              psi(2,mu,io)*psi(2,juo,io))*eio*dQo(io,ix)
             Haux(2,juo,mu) = Haux(2,juo,mu) +
     .             (psiper(1,mu,io) * psi(2,juo,io) -
     .              psiper(2,mu,io) * psi(1,juo,io) +
     .              psi(1,mu,io) * psiper(2,juo,io) -
     .              psi(2,mu,io) * psiper(1,juo,io)) * qio * eio +
     .              qio * evper(1,ix,io) *
     .              (psi(1,mu,io) * psi(2,juo,io) -
     .               psi(2,mu,io) * psi(1,juo,io) )
     .       + (psi(1,mu,io)*psi(2,juo,io) -
     .          psi(2,mu,io)*psi(1,juo,io))*eio*dQo(io,ix)
           enddo
          enddo
 370     enddo

         do iuo = 1, nuo
          do nu = 1, numh(iuo)
           ind = listhptr(iuo) + nu
           jo = listh(ind)
           juo = indxuo(jo)
           kXij = kpoint(1,ik) * Xij(1,ind) +
     &            kpoint(2,ik) * Xij(2,ind) +
     &            kpoint(3,ik) * Xij(3,ind)
           ckXij = cos(kXij)
           skXij = sin(kXij)
           Rhoper(ind,ix,ispin) = Rhoper(ind,ix,ispin) +
     .                            Saux(1,juo,iuo)*ckxij -
     .                            Saux(2,juo,iuo)*skxij
           erhoper(ind,ispin,ix) = erhoper(ind,ispin,ix) +
     .                             Haux(1,juo,iuo)*ckxij -
     .                             Haux(2,juo,iuo)*skxij
          enddo
         enddo  

215     enddo !spatial coordinates
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
            do iuo = 1, nuo
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
            do iuo = 1, nuo
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
            call cdiag( Haux, Saux, nuo, nuo, nuo, eo(:,ispin,ik),
     &                psi, nuo, 1, ierror )
C        Compute the change in density matrix elements -------------
            Saux = 0.0_dp
            Haux = 0.0_dp
            do 550 io = 1, nbands
              qio = qO(io,ispin,ik)
              if(qio .lt. 1.0e-6_dp) goto 550
                eio = eo(io,ispin,ik)
                dQo(io,1) = -wk(ik)*dstepf((eio-ef)/T)/
     .                 (T*nspin+1.0e-12)
              do mu = 1, nuo
                do juo = 1, nuo
                  nu = juo
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
            do iuo=1,nuo
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
                  Rhoper(ind,ix,ispin) = Rhoper(ind,ix,ispin) +
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

