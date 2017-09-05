      subroutine delrhokp(no, nuo, maxo, maxspn, nspin, eval, eigtol,
     &                   indxuo, nk, kpoint, wk, eo,  
     &                   Xij, Qo, H, S, Hper, Oper, maxnh, numh, 
     &                   listh, listhptr, ef, T, Rhoper, Erhoper, 
     &                   Haux, Saux, psi, nuotot, iscf)
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
      use sys
      use parallel,     only : Node, Nodes
      use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb
      use mpi_siesta
      use alloc,        only : re_alloc, de_alloc

      implicit none

      integer ::  no, nuo, iscf, maxo, nuotot,
     &            ispin, maxnh, nspin, numh(*),
     &            listh(maxnh), listhptr(*), nk,indxuo(no) 
      real(dp) :: eval(maxo,nspin,nk), Qo(maxo,maxspn,nk),
     &            Hper(maxnh,3,nspin), eigtol, S(maxnh),
     &            Oper(maxnh,3), ef, T,Rhoper(maxnh,nspin,3),
     &            Erhoper(maxnh,nspin,3), H(maxnh,nspin),
     &            Haux(2,nuotot,nuotot), Saux(2,nuotot,nuotot),
     &            psi(2,nuotot,nuotot), kpoint(3,nk), Xij(3,maxnh),
     &            eo(maxo,maxspn,nk),wk(nk)
      external     cdiag

C     Internal Variables
      integer           :: BNode, MPIerror
      integer,  pointer :: numhg(:) => null(),
     &                     listhptrg(:) => null(),
     &                     listhg(:) => null()
      real(dp), pointer :: Snew(:) => null(),
     &                     Hnew(:,:) => null(),
     &                     Xijloc(:,:) => null(),
     &                     dSnew(:,:) => null(),
     &                     dHnew(:,:,:) => null(),
     &                     Rhoperloc(:,:,:) => null(),
     &                     Erhoperloc(:,:,:) => null(),
     &                     Dloc(:) => null()

      integer :: deg, N, io, iio, nuog, maxnhg, maxnumh, maxorb, 
     &           numb(maxo), k,j, ik, ix,
     &           jden, kden, jo, mu, nu, ind, indmn, i, indbetas,
     &           ialpha, ibeta, jbeta,indi,indj,ierror,indden,
     &           iuo, juo, ind2, nbands, maxspn 
      real(dp) :: def(3), A(3), B(3), 
     &            dQo(maxo,3), aux(no),auxcoef(2,maxo), qio, qde,
     &            psiper(2,maxo,maxo), prod1(2), prod2(2), kXij, 
     &            prod3(2), prod4(2), eio, ejo, dStepF, 
     &            evper(maxo,3),
     &            pipj1, pipj2, dpipj1, dpipj2
*     &            , Haux(2,maxo,nuo),
*     &            Saux(2,maxo,nuo) !maxo = nuotot


      real(dp), pointer :: Psiden(:,:,:),
     &                     eden(:),rotaux(:,:), Hauxden(:,:,:), 
     &                     Sauxden(:,:,:)
      real(dp) :: ckXij, skXij, pert

      call timer('delrhok',1)

C Globalise list arrays - assumes listh and listd are the same

      nullify(numhg, listhptrg)
C Allocate local memory for global list arrays
      call re_alloc( numhg, 1, nuotot, name='numhg',
     &               routine= 'delrhokp' )
      call re_alloc( listhptrg, 1, nuotot, name='listhptrg',
     &               routine= 'delrhokp' )

C Find maximum value in numh and create local storage
      maxnumh = 0
      do io = 1,nuo
        maxnumh = max(maxnumh,numh(io))
      enddo
      nullify(Dloc)
      call re_alloc( Dloc, 1, maxnumh, name='Dloc',
     &               routine= 'delrhokp' )

C Globalise numh
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          numhg(io) = numh(iio)
        endif
        call MPI_Bcast(numhg(io),1,MPI_integer,BNode,
     &    MPI_Comm_World,MPIerror)
      enddo

C Build global listhptr
      listhptrg(1) = 0
      do io = 2,nuotot
        listhptrg(io) = listhptrg(io-1) + numhg(io-1)
      enddo

C Globalse listh
      maxnhg = listhptrg(nuotot) + numhg(nuotot)
      nullify(listhg)
      call re_alloc( listhg, 1, maxnhg, name='listhg',
     &               routine= 'delrhokp' )
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do jo = 1,numhg(io)
            listhg(listhptrg(io)+1:listhptrg(io)+numhg(io)) =
     &        listh(listhptr(iio)+1:listhptr(iio)+numh(iio))
          enddo
        endif
        call MPI_Bcast(listhg(listhptrg(io)+1),numhg(io),MPI_integer,
     &    BNode,MPI_Comm_World,MPIerror)
      enddo

C     Create new distribution of H and S and the likes
      nuog = nuotot
      nullify(Snew,Hnew,dSnew,dHnew,xijloc)
      call re_alloc( Snew, 1, maxnhg, name='Snew',
     &               routine= 'delrhokp' )
      call re_alloc( Hnew, 1, maxnhg, 1, nspin, name='Hnew',
     &               routine= 'delrhokp' )
      call re_alloc( dSnew, 1, maxnhg, 1, 3, name='dSnew',
     &               routine= 'delrhokp' )
      call re_alloc( dHnew,1,maxnhg,1,3,1,nspin,name='dHnew',
     &               routine= 'delrhokp' )
      call re_alloc( xijloc, 1, 3, 1, maxnhg, name='xijloc',
     &               routine= 'delrhokp' )

      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do ispin = 1,nspin
            do jo = 1,numh(iio)
              Hnew(listhptrg(io)+jo,ispin) = H(listhptr(iio)+jo,ispin)
            enddo
            do ix = 1, 3
              do jo = 1,numh(iio)
                dHnew(listhptrg(io)+jo,ix,ispin) = 
     &                             Hper(listhptr(iio)+jo,ix,ispin)
              enddo
            enddo
          enddo
          do jo = 1,numh(iio)
            Snew(listhptrg(io)+jo) = S(listhptr(iio)+jo)
          enddo
          do ix = 1, 3
            do jo = 1,numh(iio)
              dSnew(listhptrg(io)+jo,ix) = Oper(listhptr(iio)+jo,ix)
            enddo
          enddo
          do jo = 1,numh(iio)
            xijloc(1:3,listhptrg(io)+jo) = xij(1:3,listhptr(iio)+jo)
          enddo
        endif
        do ispin = 1,nspin
          call MPI_Bcast(Hnew(listhptrg(io)+1,ispin),numhg(io),
     &      MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
          do ix = 1, 3
            call MPI_Bcast(dHnew(listhptrg(io)+1,ix,ispin),numhg(io),
     &        MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
          enddo
        enddo
        call MPI_Bcast(Snew(listhptrg(io)+1),numhg(io),
     &    MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
        do ix = 1, 3
          call MPI_Bcast(dSnew(listhptrg(io)+1,ix),numhg(io),
     &      MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
        enddo
        call MPI_Bcast(xijloc(1,listhptrg(io)+1),3*numhg(io),
     &    MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
      enddo


C     Initialize variables --- 

      nullify(Rhoperloc,Erhoperloc)
C Allocate local copies of Dk and Ek for building matrices
      call re_alloc( Rhoperloc, 1, maxnhg, 1,2, 1,3, name='Rhoperloc',
     &               routine= 'delrhokp' )
      call re_alloc( Erhoperloc, 1, maxnhg, 1,2, 1,3, name='ERhoperloc',
     &               routine= 'delrhokp' )


      Rhoperloc(1:maxnhg,1:nspin,1:3) = 0.0_dp
      Erhoperloc(1:maxnhg,1:nspin,1:3) = 0.0_dp

      def(1:3) = 0.0_dp
      A(1:3) = 0.0_dp
      B(1:3) = 0.0_dp

      nbands = nuotot

C     Find eigenvectors (where stored for only one k-point) ---------
      do ik = 1+Node,nk,Nodes     ! for each k-point 
       do ispin = 1, nspin !for each spin component
!$OMP parallel default(shared),
!$OMP&private(iuo,juo,io,j,ind,jo,kxij,ckxij,skxij)
!$OMP do collapse(2)
         do iuo = 1,nuog
          do juo = 1,nuotot
            Saux(1,juo,iuo) = 0.0d0
            Saux(2,juo,iuo) = 0.0d0
            Haux(1,juo,iuo) = 0.0d0
            Haux(2,juo,iuo) = 0.0d0
          enddo
        enddo
!$OMP end do
!$OMP do
        do io = 1, nuog 
         iuo = indxuo(io)
         do j = 1, numhg(io)
          ind = listhptrg(io) + j
          jo = listhg(ind)
          juo = indxuo(jo)
C         calculates the phases k*r_ij -------------------------------
          kXij = kpoint(1,ik) * Xijloc(1,ind) +
     &           kpoint(2,ik) * Xijloc(2,ind) +
     &           kpoint(3,ik) * Xijloc(3,ind) 
          ckXij = cos(kXij)  
          skXij = sin(kXij)
C         Calculates the hamiltonian and the overlap in k space ------
C         H(k) = Sum(R) exp(i*k*R) * H(R) ----------------------------
          Saux(1,juo,iuo) = Saux(1,juo,iuo) + Snew(ind)*ckxij
          Saux(2,juo,iuo) = Saux(2,juo,iuo) - Snew(ind)*skxij!-
          Haux(1,juo,iuo) = Haux(1,juo,iuo) + Hnew(ind,ispin)*ckxij
          Haux(2,juo,iuo) = Haux(2,juo,iuo) - Hnew(ind,ispin)*skxij!-
         enddo
        enddo
!$OMP end do nowait
!$OMP end parallel
C       Symmetrize H & S ---- (MP not completly sure why) ------------
!$OMP do
        do iuo = 1, nuog
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
!$OMP end do

C       Diagonalize for each k-point --------------------------------
        call cdiag( Haux, Saux, nuotot, nuotot, nuog, eo(1,ispin,ik),
     &                psi, nuotot, iscf, ierror, -1 ) !dummy blocksize

        if (ierror.ne.0) then
          call die('DELRHOK: Terminating due to failed 
     &                                        k-diagonalisation')
        endif


C       loop on spatial xyz -----------------------------------------
        do 215 ix = 1,3
C ------ Find out the degenerated subespaces ------------------------
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
                    Sauxden(2,jden,kden) = 0.0_dp
                  endif
                  do mu=1,nuog
                    do nu = 1,numhg(mu)
                      indmn=listhptrg(mu) + nu
                      jo = listhg(indmn)
                      juo = indxuo(jo)
                      kXij = kpoint(1,ik) * Xijloc(1,indmn) +
     .                       kpoint(2,ik) * Xijloc(2,indmn) +
     .                       kpoint(3,ik) * Xijloc(3,indmn) 
                      ckXij = cos(kXij)  
                      skXij = sin(kXij) 
                      pipj1 = psi(1,mu,j) * psi(1,juo,k) +
     &                        psi(2,mu,j) * psi(2,juo,k)
                      pipj2 = psi(1,mu,j) * psi(2,juo,k) -
     &                        psi(2,mu,j) * psi(1,juo,k)
                      pert = dHnew(indmn,ix,ispin) - !looks ok 
     &                       eio * dSnew(indmn,ix)
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
     .             eden, psiden,numb(io),1,ierror,-1 )
              if (ierror.ne.0) then
               call die('DELRHOK: Terminating due to failed
     &                                deg- diagonalisation')
              endif

C             Rotate the coefficients inside the subspace ---------- 
              do mu = 1,nuog
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
*C --- End of work on degenerated subespaces ------------------------
         evper(1:maxo,ix) = 0.0_dp

C ------ Compute the change in the coefficients ---------------------
         do 300 io = 1, nbands
          qio = Qo(io,ispin,ik)
          if(abs(qio) .lt. 1.0e-6_dp) cycle

          eio = eo(io,ispin,ik)

          auxcoef(1,1:maxo) = 0.0_dp
          auxcoef(2,1:maxo) = 0.0_dp
*          auxcoef(1,1:nuotot) = 0.0_dp
*          auxcoef(2,1:nuotot) = 0.0_dp
          
           do ialpha = 1, nuog
           prod1(1:2) = 0.0_dp
           prod2(1:2) = 0.0_dp
            do jbeta = 1, numhg(ialpha)
            indmn = listhptrg(ialpha) + jbeta
            ibeta = indxuo( listhg(indmn) )
            kXij = kpoint(1,ik) * Xijloc(1,indmn) +
     &             kpoint(2,ik) * Xijloc(2,indmn) +
     &             kpoint(3,ik) * Xijloc(3,indmn)
            ckXij = cos(kXij)
            skXij = sin(kXij)
            prod1(1) = prod1(1) + dHnew(indmn,ix,ispin) *
     &                (psi(1,ibeta,io) * ckXij - 
     &                 psi(2,ibeta,io) * skXij)
            !both psis look good
            prod1(2) = prod1(2) + dHnew(indmn,ix,ispin) *
     &                (psi(1,ibeta,io) * skXij + 
     &                 psi(2,ibeta,io) * ckXij)
            prod2(1) = prod2(1) + dSnew(indmn,ix) *
     &                (psi(1,ibeta,io) * ckXij - 
     &                psi(2,ibeta,io) * skXij)
            prod2(2) = prod2(2) + dSnew(indmn,ix) *
     .                (psi(1,ibeta,io) * skXij + 
     &                 psi(2,ibeta,io) * ckXij) 
            enddo !jbeta

           prod3(1) = prod1(1) - eio * prod2(1)
           prod3(2) = prod1(2) - eio * prod2(2)
           prod4(1) = prod2(1) * 0.5_dp
           prod4(2) = prod2(2) * 0.5_dp
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
           evper(io,ix) = evper(io,ix) +
     .            psi(1,ialpha,io)*prod3(1) + psi(2,ialpha,io)*prod3(2)
           enddo !ialpha      
          
           do ialpha = 1,nuog
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
          dQo(io,ix) = wk(ik)*evper(io,ix)*dstepf((eio-ef)/T)/
     .                 (T*nspin+1.0e-15)
          A(ix) = A(ix) + wk(ik)*evper(io,ix)*dstepf((eio-ef)/T)
          B(ix) = B(ix) + wk(ik)*dstepf((eio-ef)/T)
*********
!          endif !occupancy
300      enddo !io bands loop
! ----- End change in coefficients -----------

C        Compute the change in density matrix elements -------------
         Saux = 0.0_dp
         Haux = 0.0_dp

         do 370 io = 1, nbands
          qio = qO(io,ispin,ik)
          if(qio .lt. 1.0e-6_dp) cycle
          eio = eo(io,ispin,ik) * qio
          qde = qio * evper(io,ix) + eio * dQo(io,ix)
          do mu = 1, nuog
           do juo = 1, nuotot
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
            nu = juo
            Saux(1,juo,mu) = Saux(1,juo,mu) + dpipj1 * qio
     .                     + pipj1 * dQo(io,ix)
            Saux(2,juo,mu) = Saux(2,juo,mu) + dpipj2 * qio
     .                     + pipj2 * dQo(IO,ix)
            Haux(1,juo,mu) = Haux(1,juo,mu) + dpipj1 * eio 
     .                     + qde * pipj1
            Haux(2,juo,mu) = Haux(2,juo,mu) + dpipj2 * eio 
     .                     + qde * pipj2 
           enddo
          enddo
 370     enddo

         do io = 1, nuog
          iuo = indxuo(io)
          do nu = 1, numhg(iuo)
           ind = listhptrg(iuo) + nu
           jo = listhg(ind)
           juo = indxuo(jo)
           kXij = kpoint(1,ik) * Xijloc(1,ind) +
     &            kpoint(2,ik) * Xijloc(2,ind) +
     &            kpoint(3,ik) * Xijloc(3,ind)
           ckXij = cos(kXij)
           skXij = sin(kXij)
           Rhoperloc(ind,ispin,ix) = Rhoperloc(ind,ispin,ix) 
     .                             + Saux(1,juo,iuo)*ckxij 
     .                             - Saux(2,juo,iuo)*skxij
           Erhoperloc(ind,ispin,ix) = Erhoperloc(ind,ispin,ix) 
     .                              + Haux(1,juo,iuo)*ckxij 
     .                              - Haux(2,juo,iuo)*skxij
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
        do 800 ik=1+Node,nk,Nodes
          do 770 ispin=1,nspin
            Haux=0.0_dp
            Saux=0.0_dp
            do io = 1, nuog
              iuo = indxuo(io)
              do j = 1, numhg(iuo)
                ind = listhptrg(iuo) + j
                jo = listhg(ind)
                juo = indxuo(jo)
C         calculates the phases k*r_ij -------------------------------
                kXij = kpoint(1,ik) * Xijloc(1,ind) +
     &                 kpoint(2,ik) * Xijloc(2,ind) +
     &                 kpoint(3,ik) * Xijloc(3,ind)
                ckXij = cos(kXij)
                skXij = sin(kXij)
C         Calculates the hamiltonian and the overlap in k space ------
C         H(k) = Sum(R) exp(i*k*R) * H(R) ----------------------------
               Saux(1,juo,iuo) = Saux(1,juo,iuo) + Snew(ind)*ckxij
               Saux(2,juo,iuo) = Saux(2,juo,iuo) - Snew(ind)*skxij!-
               Haux(1,juo,iuo) = Haux(1,juo,iuo) + Hnew(ind,ispin)*ckxij
               Haux(2,juo,iuo) = Haux(2,juo,iuo) - Hnew(ind,ispin)*skxij!-
              enddo
            enddo
C       Symmetrize H & S ---------------------------------------------
            do iuo = 1, nuog
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
            call cdiag( Haux, Saux, nuotot, nuotot, nuog,
     &                  eo(1,ispin,ik), psi, nbands, 1, ierror, -1 )
C        Compute the change in density matrix elements -------------
            Saux = 0.0_dp
            Haux = 0.0_dp
            do 550 io = 1, nbands
              qio = qO(io,ispin,ik)
              if(qio .lt. 1.0e-6_dp) cycle
              dQo(io,1) = -wk(ik)*dstepf((eio-ef)/T)/
     .                    (T*nspin+1.0e-12)
              eio = eo(io,ispin,ik)*dQo(io,1)

              do mu = 1, nuog
                do juo = 1, nuotot
                  nu = juo
                  pipj1 = (psi(1,mu,io)*psi(1,juo,io) +
     .              psi(2,mu,io)*psi(2,juo,io))
                  pipj2 = (psi(1,mu,io)*psi(2,juo,io) -
     .              psi(2,mu,io)*psi(1,juo,io))

                  Saux(1,juo,mu) = Saux(1,juo,mu) + pipj1*dQo(io,1)
                  Saux(2,juo,mu) = Saux(2,juo,mu) + pipj2*dQo(io,1)

                  haux(1,juo,mu) = haux(1,juo,mu) + pipj1 * eio
                  haux(2,juo,mu) = haux(2,juo,mu) + pipj2 * eio
                enddo
              enddo
550         enddo
            do io=1,nuog
              iuo = indxuo(io)
              do j=1,numhg(iuo)
                ind = listhptrg(iuo) + j
                jo = listhg(ind)
                juo = indxuo(jo)
                kXij = kpoint(1,ik) * Xijloc(1,ind) +
     &                 kpoint(2,ik) * Xijloc(2,ind) +
     &                 kpoint(3,ik) * Xijloc(3,ind)
                do ix=1,3
                  ckXij = cos(kXij)*def(ix)
                  skXij = sin(kXij)*def(ix)
                  Rhoperloc(ind,ispin,ix) = Rhoperloc(ind,ispin,ix) 
     &                                    + Saux(1,juo,iuo)*ckxij 
     &                                    - Saux(2,juo,iuo)*skxij
                  Erhoperloc(ind,ispin,ix) = Erhoperloc(ind,ispin,ix) 
     &                                     + Haux(1,juo,iuo)*ckxij 
     &                                     - Haux(2,juo,iuo)*skxij
                enddo !ix loop
              enddo
            enddo
770       enddo
800     enddo
      endif !Variation of the energy level ocupation
C

      do io = 1,nuog
        call WhichNodeOrb(io,Nodes,BNode)
        call GlobalToLocalOrb(io,BNode,Nodes,iio)
        do ix = 1,3
          do ispin = 1,nspin
            call MPI_Reduce(Rhoperloc(listhptrg(io)+1,ispin,ix),
     &        Dloc(1),numhg(io),MPI_double_precision,
     &        MPI_sum,BNode,MPI_Comm_World,MPIerror)
            if (Node.eq.BNode) then
              do j = 1,numh(iio)
                Rhoper(listhptr(iio)+j,ispin,ix) = Dloc(j)
              enddo
            endif
            call MPI_Reduce(Erhoperloc(listhptrg(io)+1,ispin,ix),
     &        Dloc(1),numhg(io),MPI_double_precision,
     &        MPI_sum,BNode,MPI_Comm_World,MPIerror)
            if (Node.eq.BNode) then
              do j = 1,numh(iio)
                Erhoper(listhptr(iio)+j,ispin,ix) = Dloc(j)
              enddo
            endif
          enddo
        enddo
      enddo

C Free local memory
      call de_alloc( xijloc, name='xijloc', routine= 'delrhokp' )
      call de_alloc( Hnew, name='Hnew', routine= 'delrhokp' )
      call de_alloc( Snew, name='Snew', routine= 'delrhokp' )
      call de_alloc( dHnew, name='dHnew', routine= 'delrhokp' )
      call de_alloc( dSnew, name='dSnew', routine= 'delrhokp' )
      call de_alloc( Dloc, name='Dloc', routine= 'delrhokp' )
      call de_alloc( listhg, name='listhg', routine= 'delrhokp' )
      call de_alloc( listhptrg, name='listhptrg', routine= 'delrhokp' )
      call de_alloc( numhg, name='numhg', routine= 'delrhokp' )
      call de_alloc( Rhoperloc, name='Rhoperloc', routine= 'delrhokp' )
      call de_alloc( Erhoperloc, name='Erhoperloc', routine= 'delrhokp')

      call timer('delrhok',2)
      return
      end

