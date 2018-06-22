      subroutine delrhok(no, nuo, maxo, maxspn, nspin, eval, eigtol,
     &                   indxuo, nk, kpoint, wk, eo,  
     &                   xij, Qo, H, S, Hper, Oper, maxnh, numh, 
     &                   listh, listhptr, ef, T,
     &                   Rhoper, Erhoper, Haux, Saux, psi, iscf, UsePsi)
C **********************************************************************
C FINDS THE CHANGE IN DENSITY MATRIX ELEMENTS DUE TO DISPLACEMENTS OF
C THE ATOMS
C K - SAMPLING VERSION (SERIAL AND PARALEL OVER ORBITALS).
C CODED BY J. JUNQUERA AND J. M. ALONSO PRUNEDA. Dec '98
C FOR SIESTA 3.X L. Riches, SUMMER '15
C CHECKED AND CORRECTED BY S. ILLERA APRIL '16
C **********************INPUT*******************************************
C INTEGER NO		     :Number of orbitals in the supercell
C INTEGER NUO                :Number of basis orbitals in unit cell (local)
C INTEGER MAXO		     :Number of orbitals
C INTEGER MAXSPN	     :Maximum number of differents spin polarizations
C INTEGER NSPIN              :Spin polarization
C REAL*8 EVAL(MAXO)          :Eigenvalues of non-perturbated Hamiltonian
C REAL*8  EIGTOL             :Tolerance to assume degenerate energy levels
C INTEGER INDXUO(NO)	     :Index of equivalent orbital in unit cell
C INTEGER NK		     :Number of kpoints
C REAL*8 KPOINT(3,NK)        :k point vectors 
C REAL*8 WK(NK)		     :k points weight		
C REAL*8 EO(MAXO,MAXSPN,NK)  :Eigenvalues
C REAL*8 XIJ(3,MAXNH)        :Vectors between orbital centers
C REAL*8 QO(MAXO)             :Occupations of unpertubed eigenstates
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
C INTEGER LISTHPTR(NIO)      :Pointer to each row (-1) of the
C                             density matrix
C REAL*8  EF                  :Fermi level
C REAL*8  T                   :Temperature
C REAL*8 Haux(2,MAXO,NUO)    :Auxiliary space for the Hamiltonian
C REAL*8 Saux(2,MAXO,NUO)    :Auxiliary space for the overlap matrix
C REAL*8 PSI(2,MAXO,NUO)     :Auxiliary space for the eigenvectors
C ******************  OUTPUT  ******************************************
C REAL*8  RHOPER(MAXNH,MAXSPN,3)  :Matrix elements of the perturbated DM
C REAL*8  ERHOPER(MAXNH,MAXSPN,3) :Matrix elements of the perturbated 
C                                  Energy Density Matrix
C **********************************************************************


      use precision,      only : dp
      use alloc
      use parallel,     only : Node, Nodes, BlockSize, IONode
      use parallelsubs,  only: GlobalToLocalOrb, LocalToGlobalOrb,
     &                         WhichNodeOrb, GetNodeOrbs
      use m_diag_option, only: diag_set_serial, diag_set_parallel

#ifdef MPI
      use m_mpi_utils,  only : broadcast, Globalize_sum
      use mpi_siesta
#endif

#ifdef _OPENMP
      use omp_lib
#endif

      implicit none

      integer ::  no, nuo, maxo, iscf, 
     &            ispin, maxnh, nspin, numh(*),
     &            listh(maxnh), listhptr(*), nk,indxuo(no) 
      real(dp) :: eval(maxo,nspin,nk), Qo(maxo,maxspn,nk),
     &            Hper(maxnh,3,nspin), eigtol, S(maxnh),
     &            Oper(maxnh,3), ef, T,Rhoper(maxnh,nspin,3),
     &            Erhoper(maxnh,nspin,3), H(maxnh,nspin),
     &            psi(2,maxo,nuo), kpoint(3,nk), xij(3,maxnh),
     &            eo(maxo,maxspn,nk),wk(nk)
      real(dp), target :: Haux(2,maxo,nuo)
      real(dp) :: Saux(2,maxo,nuo)
      logical :: UsePsi


C     Internal Variables
      integer :: deg, N, io, numb(maxo), k, j, ik, ix, iio, 
     &           jden, kden, jo, mu, mug, nu, ind, indmn, i, indbetas,
     &           ialpha, ibeta, jbeta,indi,indj,ierror,indden,
     &           iuo, juo, ind2, nbands, maxspn, BNode, maxnhg, MPIerror
 
      real(dp) :: def(3), A(3), B(3), Ag(3), Bg(3),
     &            dQo(maxo,3), aux(no), qio, 
     &            psiper(2,maxo,maxo), prod1(2), prod2(2), kXij, 
     &            prod3(2), prod4(2), eio, ejo, dStepF, evper(2,3,maxo),
     &            pipj1, pipj2,dpipj1,dpipj2,qde

      real(dp) :: ckXij, skXij, pert


      real(dp), pointer :: auxcoef(:,:), rotaux(:,:)
      real(dp), pointer :: Psiden(:,:,:) => null()
      real(dp), pointer :: Hauxden(:,:,:) => null()
      real(dp), pointer :: HauxdenL(:,:,:) => null()
      real(dp), pointer :: Hauxdentmp(:) => null()
      real(dp), pointer :: Sauxden(:,:,:) => null()
      real(dp), pointer :: Hperg(:) => null()
      real(dp), pointer :: Operg(:) => null()
      real(dp), pointer :: xijg(:,:) => null()

      real(dp),pointer :: RhoperL(:)=>null(),
     &                    erhoperL(:)=>null(),
     &                    psig(:,:,:)=>null(),
     &                    psitmp(:,:,:)=>null()

      integer, pointer :: numhg(:)=>null(),
     &                    listhg(:)=>null(),
     &                    listhptrg(:)=>null()

#ifdef MPI
      real(dp), pointer :: psigtmp(:,:)=>null(),
     &                     Rhoperg(:)=>null(),
     &                     erhoperg(:)=>null(),
     &                     xijgtmp(:,:)=>null()
#endif

      call timer('delrhok',1)

#ifdef _OPENMP
!$OMP parallel default(shared)

      ! Check that we can use the appropiate partition in the
      ! coefficient change routine
      ! This does put a constraint on the number of OMP threads
      ! we can use. But for now I think this is good.
      ! Else, we force this to be sequential
!$OMP single
      if ( omp_get_num_threads() > nuo ) then
         write(*,*) 'Reduce OMP_NUM_THREADS to (maximum): ',nuo
         call die('delrhok: too many OMP_NUM_THREADS, compared to nuo')
      end if
!$OMP end single

!$OMP end parallel
#endif


C The MPI implementation plays with some variables in local and other variables
C in global (shared over all nodes).
C Thus, lets start "globalizing" some variables--usual notation ***g

! globalize numh to numhg
      nullify(numhg)
      call re_alloc( numhg, 1, maxo, 'numhg','delrhok')
      do io = 1,maxo
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          numhg(io) = numh(iio)
        endif
#ifdef MPI
        if (Nodes .gt. 1) then
          call MPI_Bcast(numhg(io),1,MPI_integer,BNode,
     &    MPI_Comm_World,MPIerror)
        endif
#endif
      enddo

! globalize listhptr to listhptrg
      nullify(listhptrg)
      call re_alloc( listhptrg, 1, maxo, 'listhptrg', 'delrhok' )
      listhptrg(1) = 0
      do io = 2,maxo
        listhptrg(io) = listhptrg(io-1) + numhg(io-1)
      enddo

! globalize listh to listhg
      maxnhg = listhptrg(maxo) + numhg(maxo)
      nullify(listhg)
      call re_alloc( listhg, 1, maxnhg, 'listhg', 'delrhok' )
      do io = 1,maxo
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do jo = 1,numhg(io)
            listhg(listhptrg(io)+1:listhptrg(io)+numhg(io)) =
     &        listh(listhptr(iio)+1:listhptr(iio)+numh(iio))
          enddo
        endif
#ifdef MPI
        if (Nodes .gt. 1) then
          call MPI_Bcast(listhg(listhptrg(io)+1),numhg(io),
     &  MPI_integer,BNode,MPI_Comm_World,MPIerror)
        endif
#endif
      enddo

! globalize xij(3,maxnh) to xijg(3,maxnhg)
      nullify(xijg)
      call re_alloc(xijg, 1, 3, 1, maxnhg, 'xijg', 'delrhok')
      xijg(:,:) = 0.0_dp
      do io = 1,maxo
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
           call GlobalToLocalOrb(io,Node,Nodes,iio)
          do ix=1,3
            do jo = 1,numh(iio)
              xijg(ix,listhptrg(io)+jo)=xij(ix,listhptr(iio)+jo)
            enddo
          enddo
        endif
      enddo
#ifdef MPI 
        if (Nodes .gt.1 ) then !share xijg over all nodes
          nullify(xijgtmp)
          call re_alloc(xijgtmp,1,3,1,maxnhg,'xijgtmp','delrhok')
          call globalize_sum(xijg(1:3,1:maxnhg),xijgtmp(1:3,1:maxnhg))
          xijg(:,:)=xijgtmp(:,:)
          call de_alloc(xijgtmp)
        endif
#endif
C******************************************************************************

C     Initialize variables --- 
      A(1:3) = 0.0_dp
      B(1:3) = 0.0_dp
      rhoper(:,:,:)=0.0_dp
      Erhoper(:,:,:)=0.0_dp


C     Find eigenvectors (where stored for only one k-point) ---------
      do ik = 1,nk     ! for each k-point 
       do ispin = 1, nspin !for each spin component

C     Global psis
        nullify(psig) ! will contain the global psis
        call re_alloc(psig  ,1,2,1,maxo,1,maxo,'psig','delrhok')
        nullify(psitmp) !buffer
        call re_alloc(psitmp,1,2,1,maxo,1,maxo,'psitmp','delrhok')
C     Initialize psis-auxiliar if MPI
#ifdef MPI
        nullify(psigtmp)
        call re_alloc(psigtmp,1,maxo,1,maxo,'psig','delrhok')
#endif

        if (UsePsi .and. (iscf.ne.1)) go to 999 !skip the psi calculation

        call timer( 'delrhok-psik', 1 )

! The following part looks like diagk: obtain the psis of H0 for 
! a given k-point and spin

!$OMP parallel do default(shared), 
!$OMP&privateprivate(iuo,j,ind,jo,juo,kXij,ckXij,skXij)
          
        do iuo = 1, nuo 
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
!$OMP end parallel do

C       Diagonalize for each k-point --------------------------------
        call cdiag( Haux, Saux, maxo, nuo, maxo, eo(1,ispin,ik),
     &                psi, maxo, iscf, ierror, BlockSize )

        if (ierror.ne.0) then
          call die('DELRHOK: Terminating due to failed 
     &                                        k-diagonalisation')
        endif

        if (UsePsi) then !For iscf==1, save the psi into a file
          call save_psi(psi) ! for this kpoint and spin
        endif
        call timer( 'delrhok-psik', 2 )

  999 continue !skipped the psi calculation

        if (UsePsi .and. (iscf.gt.1)) then !For iscf!=1, read the psi from file 
          call read_psi(ik,ispin,nspin,psi) ! for this kpoint and spin
        endif

C       For a given k-point, globalize the psi (wavefunction)
C       Up to here, the psi are local-->change to global variable (psig)
C       Note that eo is already global
        do io=1,nuo
          call WhichNodeOrb(io,Nodes,BNode)
          call LocalToGlobalOrb(io,Node,Nodes,iio)
          do jo=1,maxo
            do j=1,2
              psitmp(j,jo,iio)=psi(j,jo,io)
            enddo
          enddo
        enddo
        psig(:,:,:)=psitmp(:,:,:) !direct copy in serial

#ifdef MPI
        if (Nodes .gt.1 ) then
          psig(:,:,:)=0._dp
          do j=1,2
            psigtmp(:,:)=0._dp
            do io=1,maxo
              call globalize_sum(psitmp(j,io,1:maxo),
     &                         psigtmp(io,1:maxo))
            enddo
            psig(j,:,:)=psigtmp(:,:)
          enddo
        endif
        call de_alloc(psigtmp)
#endif
        call de_alloc(psitmp)

! Here, the psis are global for a given k-pint and spin
! Now, the "real" linres procedure start

C       loop on spatial xyz -----------------------------------------
        do ix = 1,3

! Lets start doing global the Hpert and Opert for a given ix and spin
         nullify( Hperg, Operg)
         call re_alloc( Operg, 1, maxnhg, 'Operg', 'delrhok' )
         call re_alloc( Hperg, 1, maxnhg, 'Hperg', 'delrhok' )

         do io = 1,maxo
           call WhichNodeOrb(io,Nodes,BNode)
           if (Node.eq.BNode) then
             call GlobalToLocalOrb(io,Node,Nodes,iio)
             do jo = 1,numh(iio)
               Hperg(listhptrg(io)+jo) =
     &                            Hper(listhptr(iio)+jo,ix,ispin)
               Operg(listhptrg(io)+jo) =
     &                            Oper(listhptr(iio)+jo,ix)
             enddo
           endif
#ifdef MPI
           if (Nodes .gt.1 ) then
             call MPI_Bcast(Hperg(listhptrg(io)+1),numhg(io),
     &         MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
             call MPI_Bcast(Operg(listhptrg(io)+1),numhg(io),
     &         MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
           endif
#endif
         enddo

C        Find out the degenerated subespaces ------------------------
         nbands = maxo
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
              nullify(Psiden,Hauxden,HauxdenL,Sauxden)
              call re_alloc(Psiden, 1, 2, 1, numb(io), 1, numb(io),
     &                       'Psiden', 'delrhok')
              call re_alloc(Hauxden, 1, 2, 1, numb(io), 1, numb(io),
     &                       'Hauxden', 'delrhok')
              call re_alloc(HauxdenL, 1, 2, 1, numb(io), 1, numb(io),
     &                       'HauxdenL', 'delrhok')
              call re_alloc(Sauxden, 1, 2, 1, numb(io), 1, numb(io),
     &                       'Sauxden', 'delrhok')

              ! Re-use memory for rotation.
              rotaux => Sauxden(:,:,1)

!$OMP parallel default(shared),
!$OMP&private(j,jden,k,kden,mu,mug,nu,indmn,jo,juo,kXij,ckXij,skXij),
!$OMP&private(pipj1,pipj2,pert)
              
              ! Initialize to 0
              Sauxden = 0._dp
              HauxdenL(:,:,:) = 0._dp
!$OMP do
              do j = 1,numb(io)
                Sauxden(1,j,j)=1.0_dp
              enddo
!$OMP end do

              ! TODO transpose this construction, i.e. 
              ! the access is sub-optimal
              ! Construction of the LOCAL deg Hamiltonian
!$OMP do
              do j = io-numb(io)+1, io
                jden = j - io + numb(io)
                numb(j) = numb(io)
                do k = io-numb(io)+1, io
                  kden = k - io + numb(io)
                  do mu=1,nuo
                    call LocalToGlobalOrb(mu,Node,Nodes,mug)
                    do nu = 1,numh(mu)
                      indmn=listhptr(mu) + nu
                      jo = listh(indmn)
                      juo = indxuo(jo)
                      kXij = kpoint(1,ik) * Xij(1,indmn) +
     .                       kpoint(2,ik) * Xij(2,indmn) +
     .                       kpoint(3,ik) * Xij(3,indmn) 
                      ckXij = cos(kXij)  
                      skXij = sin(kXij) 
                      pipj1 = psig(1,mug,j) * psig(1,juo,k) +
     &                        psig(2,mug,j) * psig(2,juo,k)
                      pipj2 = psig(1,mug,j) * psig(2,juo,k) -
     &                        psig(2,mug,j) * psig(1,juo,k)
                      pert = Hper(indmn,ix,ispin) -
     &                       eio * Oper(indmn,ix)
                      HauxdenL(1,jden,kden) = HauxdenL(1,jden,kden) +
     &                        (pipj1 * ckXij - pipj2 * skXij) * pert
                      HauxdenL(2,jden,kden) = HauxdenL(2,jden,kden) +
     &                        (pipj1 * skXij + pipj2 * ckXij) * pert   
                    enddo
                  enddo
                enddo
              enddo !j loop
!$OMP end do 
!$OMP end parallel

! The construction of the perturbed H (Hauxden) was done locally. It s necessary to
! sum over all nodes
#ifdef MPI    
              if (Nodes .gt. 1) then !From Local to global degenerated Hamiltonian
                Hauxden(:,:,:)=0._dp

                nullify(Hauxdentmp)
                call re_alloc(Hauxdentmp,1,numb(io),
     &                        'Hauxdentmp','delrhok')
                do j=1,2
                  do i = 1, numb(io)
                    Hauxdentmp(:)=0._dp
                    call globalize_sum(HauxdenL(j,i,:) , Hauxdentmp(:))
                    Hauxden(j,i,:)=Hauxdentmp(:)
                  enddo
                enddo
                call de_alloc(Hauxdentmp)

              else
                Hauxden(:,:,:)=0._dp
                Hauxden(:,:,:)=HauxdenL(:,:,:) !just copy MPI-one node
              endif
#else
              Hauxden(:,:,:)=0._dp
              Hauxden(:,:,:)=HauxdenL(:,:,:) !just copy (serial case)
#endif

              ! Haux is a dummy argument, we don't need the eigenvalues
              ! and here Haux is always big enough to hold them all!
              call diag_set_serial()
              call cdiag( Hauxden, Sauxden, numb(io),numb(io), numb(io),
     .             Haux, psiden,numb(io),iscf,ierror,BlockSize )
              call diag_set_parallel()

              if (ierror.ne.0) then
               call die('DELRHOK: Terminating due to failed
     &                                deg- diagonalisation')
              endif

#ifdef MPI
              ! psiden returned by cdiag only in ONE node, share psiden
              if (Nodes.gt.1) call broadcast(psiden)
#endif

C             Rotate the coefficients inside the subspace ---------- 

!$OMP parallel do default(shared), private(mu,ialpha,ibeta)
              do mu = 1,maxo
               do ialpha = 1, numb(io)
                 ! This is actually Sauxden(:,:,1)
                 rotaux(:,ialpha) = 0._dp
                 do ibeta = 1, numb(io)
                  rotaux(1,ialpha) = rotaux(1,ialpha) + 
     &                  psiden(1,ibeta,ialpha) *
     &                  psig(1,mu,io-numb(io)+ibeta) -
     &                  psiden(2,ibeta,ialpha) *
     &                  psig(2,mu,io-numb(io)+ibeta)
                  rotaux(2,ialpha) = rotaux(2,ialpha) +
     &                 psiden(1,ibeta,ialpha) *
     &                 psig(2,mu,io-numb(io)+ibeta) +
     &                 psiden(2,ibeta,ialpha) *
     &                 psig(1,mu,io-numb(io)+ibeta)
                end do
               end do

               do ialpha = 1, numb(io)
                  psig(1,mu,io-numb(io)+ialpha) = rotaux(1,ialpha)
                  psig(2,mu,io-numb(io)+ialpha) = rotaux(2,ialpha)
               end do
              enddo !mu loop
!$OMP end parallel do
              call de_alloc( Psiden, 'Psiden', 'delrhok' )
              call de_alloc( Hauxden, 'Hauxden', 'delrhok' )
              call de_alloc( HauxdenL, 'HauxdenL', 'delrhok' )
              call de_alloc( Sauxden, 'Sauxden', 'delrhok' )     
            endif !numb() last degenerated state
          endif !lt tol (degenerated)
         enddo !bands io loop

C        End degenerate case--------------------------------------------

!$OMP parallel default(shared),
!$OMP&private(i,io,qio,eio,mu,j,juo,iuo,nu,ind,jo),
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
         evper(1:2,ix,1:maxo) = 0.0_dp
         psiper(1:2,:,:) = 0.0_dp
!$OMP end workshare

C        Compute the change in the coefficients ---------------------
!$OMP do schedule(static,1)
         do i = 1, nuo
          call LocalToGlobalOrb(i,Node,Nodes,io)
          qio = Qo(io,ispin,ik)
          if (qio > 1.0e-6_dp) then !is it occupied?

          eio = eo(io,ispin,ik)

          auxcoef(:,:) = 0.0_dp
          
          do ialpha = 1, maxo
           prod1(1:2) = 0.0_dp
           prod2(1:2) = 0.0_dp
            do jbeta = 1, numhg(ialpha)
            indmn = listhptrg(ialpha) + jbeta
            ibeta = indxuo( listhg(indmn) )
            kXij = kpoint(1,ik) * Xijg(1,indmn) +
     &             kpoint(2,ik) * Xijg(2,indmn) +
     &             kpoint(3,ik) * Xijg(3,indmn)
            ckXij = cos(kXij)
            skXij = sin(kXij)
            prod1(1) = prod1(1) + Hperg(indmn) *
     &                (psig(1,ibeta,io) * ckXij - 
     &                 psig(2,ibeta,io) * skXij)
            prod1(2) = prod1(2) + Hperg(indmn) *
     &                (psig(1,ibeta,io) * skXij + 
     &                 psig(2,ibeta,io) * ckXij)
            prod2(1) = prod2(1) + Operg(indmn) *
     &                (psig(1,ibeta,io) * ckXij - 
     &                 psig(2,ibeta,io) * skXij)
            prod2(2) = prod2(2) + Operg(indmn) *
     .                (psig(1,ibeta,io) * skXij + 
     &                 psig(2,ibeta,io) * ckXij) 
            enddo !jbeta

            prod3(1) = prod1(1) - eio * prod2(1)
            prod3(2) = prod1(2) - eio * prod2(2)
            prod4(1) = prod2(1) / 2.0_dp
            prod4(2) = prod2(2) / 2.0_dp
            do j = 1,nbands
             ejo = eo(j,ispin,ik)
             if(abs(eio-ejo) .gt. eigtol) then
               auxcoef(1,j) = auxcoef(1,j) + (1.0_dp/(eio-ejo)) *
     .           (psig(1,ialpha,j)*prod3(1) + psig(2,ialpha,j)*prod3(2))
               auxcoef(2,j) = auxcoef(2,j) + (1.0_dp/(eio-ejo)) *
     .           (psig(1,ialpha,j)*prod3(2) - psig(2,ialpha,j)*prod3(1))
             else
               auxcoef(1,j) = auxcoef(1,j) -
     .           (psig(1,ialpha,j)*prod4(1) + psig(2,ialpha,j)*prod4(2))
               auxcoef(2,j) = auxcoef(2,J) -
     .           (psig(1,ialpha,j)*prod4(2) - psig(2,ialpha,j)*prod4(1))
             endif
            enddo ! jbands
            evper(1,ix,io) = evper(1,ix,io) +
     .           psig(1,ialpha,io)*prod3(1) + psig(2,ialpha,io)*prod3(2)
            evper(2,ix,io) = evper(2,ix,io) +
     .           psig(1,ialpha,io)*prod3(2) - psig(2,ialpha,io)*prod3(1)
           enddo !ialpha      
          
           do ialpha = 1,maxo
            psiper(1:2,ialpha,io) = 0.0_dp
            do j = 1, nbands
              psiper(1,ialpha,io) = psiper(1,ialpha,io) +
     .                           auxcoef(1,j)*psig(1,ialpha,j) -
     .                           auxcoef(2,j)*psig(2,ialpha,j)
              psiper(2,ialpha,io) = psiper(2,ialpha,io) +
     .                           auxcoef(1,j)*psig(2,ialpha,j) +
     .                           auxcoef(2,j)*psig(1,ialpha,j)
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
!$OMP end parallel
         
C dQo, Evper, psiper have dimension of total number of orbitals.
C But they are distributed on the nodes being different to zero 
C only the elements that correspond to its node.

C Compute the change in density matrix elements -------------

         Saux = 0.0_dp
         Haux = 0.0_dp

! The construction of the perturbed density/energy matrix is
! paralelized over bands

         nullify(Sauxden,Hauxden)
         call re_alloc(Sauxden, 1, 2, 1, maxo, 1, maxo,
     &                       'Sauxden', 'delrhok')
         call re_alloc(Hauxden, 1, 2, 1, maxo, 1, maxo,
     &                       'Hauxden', 'delrhok')

         Sauxden = 0._dp
         Hauxden = 0._dp

!$OMP parallel default(shared),
!$OMP&private(i,io,qio,eio,qde,mu,juo),
!$OMP&private(dpipj1,dpipj2,pipj1,pipj2)


         do i = 1, nuo
          call LocalToGlobalOrb(i,Node,Nodes,io)
          qio = qO(io,ispin,ik)
          if ( qio < 1.0e-6_dp ) cycle
          eio = eo(io,ispin,ik) * qio
          qde = qio * evper(1,ix,io) + eio * dQo(io,ix)
!$OMP do collapse(2)
          do mu = 1, maxo
           do juo = 1, maxo
            dpipj1 = (psiper(1,mu,io) * psig(1,juo,io) +
     &                psiper(2,mu,io) * psig(2,juo,io) +
     &                psig(1,mu,io) * psiper(1,juo,io) +
     &                psig(2,mu,io) * psiper(2,juo,io))
            dpipj2 = (psiper(1,mu,io) * psig(2,juo,io) -
     &                psiper(2,mu,io) * psig(1,juo,io) +
     &                psig(1,mu,io) * psiper(2,juo,io) -
     &                psig(2,mu,io) * psiper(1,juo,io))
            pipj1 = (psi(1,mu,i)*psi(1,juo,i) +
     &               psi(2,mu,i)*psi(2,juo,i))
            pipj2 = (psi(1,mu,i)*psi(2,juo,i) -
     &               psi(2,mu,i)*psi(1,juo,i))
            Sauxden(1,juo,mu) = Sauxden(1,juo,mu) + dpipj1 * qio
     &                     + pipj1 * dQo(io,ix)
            Sauxden(2,juo,mu) = Sauxden(2,juo,mu) + dpipj2 * qio
     &                     + pipj2 * dQo(IO,ix)
            Hauxden(1,juo,mu) = Hauxden(1,juo,mu) + dpipj1 * eio
     &                     + qde * pipj1
            Hauxden(2,juo,mu) = Hauxden(2,juo,mu) + dpipj2 * eio
     &                     + qde * pipj2
           enddo
          enddo
!$OMP end do ! retain barrier
       enddo
!$OMP end parallel

       nullify(RhoperL,erhoperL)
       call re_alloc(RhoperL,1,maxnhg,'RhoperL','delrhok') !Local part
       call re_alloc(erhoperL,1,maxnhg,'erhoperL','delrhok') !Local "
       RhoperL(:) = 0._dp
       erhoperL(:) = 0._dp

!$OMP parallel default(shared),
!$OMP&private(iuo,nu,ind,jo,juo,kXij,ckXij,skXij)
!$OMP do
         do iuo = 1, maxo
          do nu = 1, numhg(iuo)
           ind = listhptrg(iuo) + nu
           jo = listhg(ind)
           juo = indxuo(jo)
           kXij = kpoint(1,ik) * Xijg(1,ind) +
     &            kpoint(2,ik) * Xijg(2,ind) +
     &            kpoint(3,ik) * Xijg(3,ind)
           ckXij = cos(kXij)
           skXij = sin(kXij)
           RhoperL(ind) = RhoperL(ind) +
     .                    Sauxden(1,juo,iuo)*ckxij -
     .                    Sauxden(2,juo,iuo)*skxij
           erhoperL(ind) = erhoperL(ind) +
     .                     Hauxden(1,juo,iuo)*ckxij -
     .                     Hauxden(2,juo,iuo)*skxij
          enddo
         enddo 
!$OMP end do nowait
!$OMP end parallel 
 
         call de_alloc(Sauxden)
         call de_alloc(Hauxden)

! Sum all the local contributions to Rhoper for a given ix,spin
#ifdef MPI
        if (Nodes .gt.1 ) then
          nullify(Rhoperg,erhoperg)
          call re_alloc(Rhoperg,1,maxnhg,'Rhoperg','delrhok') !global
          call re_alloc(erhoperg,1,maxnhg,'erhoperg','delrhok') !global

          Rhoperg(:)=0._dp
          erhoperg(:)=0._dp
! Done Locally, sum over all nodes 
          call globalize_sum(RhoperL(1:maxnhg),
     &                       Rhoperg(1:maxnhg))
          call globalize_sum(erhoperL(1:maxnhg),
     &                       erhoperg(1:maxnhg)) !shared in all nodes

! Final step, distribute the density/energy matrix over all nodes
          do io = 1, maxo
            call WhichNodeOrb(io,Nodes,BNode)
            if (Node .eq. BNode) then
              call GlobalToLocalOrb(io,BNode,Nodes,iio)
              indmn=listhptrg(io) !just copy some global parts into local ones
              do j = 1,numh(iio)
                rhoper(listhptr(iio)+j,ispin,ix) =
     &              rhoper(listhptr(iio)+j,ispin,ix) + Rhoperg(indmn+j)
                Erhoper(listhptr(iio)+j,ispin,ix) =
     &              Erhoper(listhptr(iio)+j,ispin,ix)+erhoperg(indmn+j)
              enddo
            endif
          enddo

          call de_alloc(Rhoperg)
          call de_alloc(erhoperg)
        else
! MPI one node, just copy
          rhoper(:,ispin,ix)  = rhoper(:,ispin,ix)  + RhoperL(:)
          Erhoper(:,ispin,ix) = Erhoper(:,ispin,ix) + erhoperL(:)
        endif
#else
! Serial case, just copy
        rhoper(:,ispin,ix)  = rhoper(:,ispin,ix)  + RhoperL(:)
        Erhoper(:,ispin,ix) = Erhoper(:,ispin,ix) + erhoperL(:)
#endif

        call de_alloc(RhoperL)
        call de_alloc(erhoperL)
        call de_alloc(Hperg)
        call de_alloc(Operg)
         
        enddo !spatial coordinates
        call de_alloc(psig)
       enddo !nspin loop
      enddo ! nk loop 

! A and B are local
      Ag(:) = 0.0_dp
      Bg(:) = 0.0_dp

#ifdef MPI
      call globalize_sum(A(:),Ag(:))
      call globalize_sum(A(:),Bg(:))
#else
      Ag(:) = A(:)
      Bg(:) = B(:)
#endif

      do ix = 1,3
        def(ix) = Ag(ix) / (Bg(ix) + 1.0e-12_dp)
      enddo

C **********************************************************************
c    Variation of the energy level ocupation--------------------------
C    if the fermi level changes, repeat the previous process again
C *********************************************************************

c    Variation of the energy level ocupation--------------------------
      if ( (def(1).gt.1.0e-6_dp) .or. (def(2).gt.1.0e-6_dp) .or.
     .                             (def(3).gt.1.0e-6_dp) ) then

        call timer( 'delrhok-psik-fermi', 1 )


        if (IONode) then
        print *, 'DEBUG TRACK: delrhok, change in the fermi level!!'
        endif

        do 800 ik=1,nk
          do 770 ispin=1,nspin

C     Global psis
            nullify(psig)
            call re_alloc(psig  ,1,2,1,maxo,1,maxo,'psig','delrhok')
            nullify(psitmp)
            call re_alloc(psitmp,1,2,1,maxo,1,maxo,'psitmp','delrhok')
C     Initialize psis-auxiliar if MPI
#ifdef MPI
            nullify(psigtmp)
            call re_alloc(psigtmp,1,maxo,1,maxo,'psig','delrhok')
#endif

C Init some variables-------
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

C       Diagonalize for each k-point --------------------------------
            call cdiag( Haux, Saux, maxo, nuo, maxo, eo(1,ispin,ik),
     &                psi, maxo, iscf, ierror, BlockSize )

            if (ierror.ne.0) then
              call die('DELRHOK: Terminating due to failed 
     &              k-diagonalisation in the change of fermi level')
            endif

C       Share the psis---------------------------------------------
            do io=1,nuo
              call WhichNodeOrb(io,Nodes,BNode)
              call LocalToGlobalOrb(io,Node,Nodes,iio)
              do jo=1,maxo
                do j=1,2
                  psitmp(j,jo,iio)=psi(j,jo,io)
                enddo
              enddo
            enddo
            psig(:,:,:)=psitmp(:,:,:) !direct copy in serial

#ifdef MPI
            if (Nodes .gt.1 ) then
              psig(:,:,:)=0._dp
              do j=1,2
                psigtmp(:,:)=0._dp
                do io=1,maxo
                  call globalize_sum(psitmp(j,io,1:maxo),
     &                         psigtmp(io,1:maxo))
                enddo
                psig(j,:,:)=psigtmp(:,:)
              enddo
            endif
            call de_alloc(psigtmp)
#endif
            call de_alloc(psitmp)

C        Compute the change in density matrix elements -------------

            nullify(Sauxden,Hauxden) !re-use names....
            call re_alloc(Sauxden, 1, 2, 1, maxo, 1, maxo,
     &                       'Sauxden', 'delrhok')
            call re_alloc(Hauxden, 1, 2, 1, maxo, 1, maxo,
     &                       'Hauxden', 'delrhok')

            Sauxden = 0.0_dp
            Hauxden = 0.0_dp
            do 550 i = 1, nuo
              call LocalToGlobalOrb(i,Node,Nodes,io)
              qio = qO(io,ispin,ik)
              if(qio .lt. 1.0e-6_dp) cycle
                eio = eo(io,ispin,ik)
                dQo(io,1) = -wk(ik)*dstepf((eio-ef)/T)/
     .                 (T*nspin+1.0e-12)
              do mu = 1, maxo
                do juo = 1, maxo
                  dpipj1 = (psig(1,mu,io)*psig(1,juo,io) +
     .                      psig(2,mu,io)*psig(2,juo,io))
                  dpipj2 = (psig(1,mu,io)*psig(2,juo,io) -
     .                      psig(2,mu,io)*psig(1,juo,io))

                  Sauxden(1,juo,mu) = Sauxden(1,juo,mu)
     .                           + dpipj1 * dQo(io,1)
                  Sauxden(2,juo,mu) = Sauxden(2,juo,mu)
     .                           + dpipj2 * dQo(io,1)

                  hauxden(1,juo,mu) = hauxden(1,juo,mu) +
     .                             dpipj1 * eio * dQo(io,1)
                  hauxden(2,juo,mu) = hauxden(2,juo,mu) +
     .                             dpipj2 * eio*dQo(io,1)
                enddo
              enddo
550         enddo

            call de_alloc(psig)

! hauxden and Sauxden are local
! Lets add all the parts
! sum_Nodes (hauxden) --> Psiden
! sum_Nodes (Sauxden) --> HauxdenL (to avoid define new variables)
            nullify(Psiden,HauxdenL)
            call re_alloc(Psiden,1,2,1,maxo,1,maxo,
     .                                  'Psiden','delrhok')
            call re_alloc(HauxdenL,1,2,1,maxo,1,maxo,
     .                                  'HauxdenL','delrhok')
#ifdef MPI
            if (Nodes.gt.1) then
              do j=1,2
                do io=1,maxo
              call globalize_sum(hauxden(j,io,1:maxo),
     &                         psiden(j,io,1:maxo))
              call globalize_sum(Sauxden(j,io,1:maxo),
     &                         HauxdenL(j,io,1:maxo))
                enddo
              enddo
            else
              Psiden(:,:,:) = hauxden(:,:,:) !just copy
              HauxdenL(:,:,:) = Sauxden(:,:,:)
            endif
#else
              Psiden(:,:,:) = hauxden(:,:,:) !just copy
              HauxdenL(:,:,:) = Sauxden(:,:,:)
#endif
            call de_alloc(Hauxden)
            call de_alloc(Sauxden)

!Psiden and HauxdenL are already global

            do i=1,nuo
              call LocalToGlobalOrb(i,Node,Nodes,iuo)
              do j=1,numh(i)
                ind = listhptr(i) + j
                jo = listh(ind)
                juo = indxuo(jo)
                kXij = kpoint(1,ik) * Xij(1,ind) +
     .          kpoint(2,ik) * Xij(2,ind) +
     .          kpoint(3,ik) * Xij(3,ind)
                do ix=1,3
                  ckXij = cos(kXij)*def(ix)
                  skXij = sin(kXij)*def(ix)
                  Rhoper(ind,ispin,ix) = Rhoper(ind,ispin,ix) +
     .                            HauxdenL(1,juo,iuo)*ckxij -
     .                            HauxdenL(2,juo,iuo)*skxij
                  erhoper(ind,ispin,ix) = erhoper(ind,ispin,ix) +
     .                             psiden(1,juo,iuo)*ckxij -
     .                             psiden(2,juo,iuo)*skxij
                enddo !ix loop
              enddo
            enddo

            call de_alloc(HauxdenL)
            call de_alloc(psiden)

770       enddo !nspin
800     enddo! kpoint

      call timer( 'delrhok-psik-fermi', 2 )

      endif !Variation of the energy level ocupation


      call de_alloc(numhg)
      call de_alloc(listhg)
      call de_alloc(listhptrg)

      call timer('delrhok',2)
      return

      contains !---------------------------------------------------------------

      subroutine save_psi(psi)
! This subrotuine saves for the first time (iscf=1) the psi (solved the degeneracy problem)
! for each spin for the distributed k-points in each node
!               Node 1                                           Node 2            ...
!       kpoin_1              Kpoint_2     ...   || kpoin_1  Kpoint_2   Kpoint_3 || ....   
!psi(:,1) psi(:,spin2)  psi(:,1) psi(:,spin2) ... 
! in one file per node called #.LRPSI+NODE

! *****************************INPUTS***************************************
! REAL PSI (2*nuotot*nuotot) : wavefunction for a k-point and spin
! ****************************OUTPUTS***************************************
! NONE

        use files,        only: slabel, label_length

        real(dp),intent (in) :: psi(2*maxo*nuo)

! Internal variables
        character(len=label_length+11) :: fname
        character(len=5) :: currNode
        character(len=*), parameter :: intfmt = '(I11)'
        integer :: unit1, m, init


        write(currNode,'(i5)') Node
        fname = trim(slabel)//'.LRPSI'//adjustl(currNode)

        call io_assign(unit1)

        open( unit1, file=fname, form='unformatted', action='write',
     .        position='append', status='unknown' )

        do m=1, 2*maxo*nuo
          write(unit1) psi(m)
        enddo

        call io_close(unit1)

      end subroutine save_psi

      subroutine read_psi(ik,ispin,nspin,psi)
! This subroutine reads from the #.LRPSI+NODE file the wavefunctions for a specific
! k-point and spin
! *****************************INPUTS***************************************
! REAL PSI (2*nuotot*nuotot) : wavefunction for a k-point and spin (set to 0)
! INTEGER IK                 : counter of k points in the node
! INTEGER ISPIN              : current spin of the wavefunction
! INTEGER NSPIN              : number of spin components
! ****************************OUTPUTS***************************************
! REAL PSI (2*nuotot*nuotot) : read  wavefunction for a k-point and spin

        use files,        only: slabel, label_length

        real(dp),intent (inout) :: psi(2*maxo*nuo)
        integer, intent (in)    :: ik
        integer, intent (in)    :: ispin
        integer, intent (in)    :: nspin

! Internal variables
        character(len=label_length+11) :: fname
        character(len=5) :: currNode
        character(len=*), parameter :: intfmt = '(I11)'
        integer :: unit1, m, j, init

        psi(:)=0.0_dp

        init=(ik-1)*nspin*2*maxo*nuo +
     .               (ispin-1)*2*maxo*nuo

        write(currNode,'(i5)') Node
        fname = trim(slabel)//'.LRPSI'//adjustl(currNode)
 
        call io_assign(unit1)
        open( unit1, file=fname, form='unformatted', position='rewind',
     .        action='read', status='old' )

        do m=1,init
          read(unit1)  !skip psis that not correspond to this spin or kpoint
        enddo

        do j=1,2*maxo*nuo
          read(unit1) psi(j)
        enddo

        call io_close(unit1)

      end subroutine read_psi

      end !delrhok end

