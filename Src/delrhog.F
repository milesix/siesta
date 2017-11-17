      SUBROUTINE DELRHOG(nuotot, nuo, nspin, maxorb, eval, tol, eigtol,  
     &                   occ, Hper, Oper, maxnh, numh, listh,listhptr,
     &                   ef,T,Rhoper,Erhoper,psi,iscf)


C **********************************************************************
C FINDS THE CHANGE IN DENSITY MATRIX ELEMENTS DUE TO DISPLACEMENTS OF
C THE ATOMS
C ONLY GAMMA POINT PER PHONON AND ELECTRONS.
C CODED BY J. JUNQUERA AND J. M. ALONSO PRUNEDA. Dec '98
c For Siesta 3.x LR, summer '15
C CHECKED AND CORRECTED BY S. ILLERA APRIL '16
C **********************INPUT*******************************************
C INTEGER NUO                :Number of basis orbitals in unit cell
C INTEGER ISPIN              :Index of spin (not used)
C INTEGER NSPIN              :Spin polarization
C INTEGER MAXORB             :Number of orbitals
C REAL*8 EVAL(NUO)           :Eigenvalues of non-perturbated Hamiltonian
C REAL*8  EIGTOL             :Tolerance to assume degenerate energy levels
C REAL*8 OCC(NUO)            :Occupations of unpertubed eigenstates.
C REAL*8  HPER(MAXNH)        :Matrix elements of the perturbated 
C                             Hamiltonian
C REAL*8  OPER(MAXNH)        :Matrix elements of the perturbated 
C                             Overlap
C INTEGER MAXNH              :First dimension of listh
C INTEGER NUMH(NUO)          :Number of nonzero density matrix elements
C                             for each matrix row
C INTEGER LISTH(MAXNH)       :Nonzero density matrix element column
C                             indexes
C INTEGER LISTHPTR(NUO)      :Pointer to each row (-1) of the
C                             density matrix
C INTEGER IX                 :Spatial coordinate
C REAL*8  EF                 :Fermi level
C REAL*8  T                  :Temperature
C REAL*8 PSI(NUO,NUO)        :Coeficients of the wavefunctions
C INTEGER ISCF		     :Counter of main Linres SCF loop
C ******************  OUTPUT  ******************************************
C REAL*8  RHOPER(MAXNH)  :Matrix elements of the perturbated 
C REAL*8  ERHOPER(MAXNH) :Matrix elements of the perturbated 
C                             Energy Density Matrix
C **********************************************************************

      use precision,      only : dp
      use alloc
      use parallel,      only: Node, Nodes, BlockSize
      use parallelsubs,  only: GlobalToLocalOrb, LocalToGlobalOrb,
     &                         WhichNodeOrb, GetNodeOrbs
      use m_mpi_utils,  only : broadcast, Globalize_sum

#ifdef MPI
      use mpi_siesta
#endif


      implicit none

      integer :: nuotot, nuo, nspin, maxorb, maxnh, numh(*),
     &            listh(maxnh), listhptr(*), ix, iscf
      real(dp) :: eval(maxorb), tol, occ(maxorb), 
     &            Hper(maxnh), eigtol, 
     &            Oper(maxnh), ef, T,Rhoper(maxnh),Erhoper(maxnh)
      real(dp) :: psi(nuotot,nuo)


C     Internal Variables
      integer :: deg, N, io, numb(maxorb), k,j, maxden,
     &           jden, kden, jo, mu, nu, ind, indmn, i, indbetas,
     &           ialpha, ibeta, jbeta,indi,indj,indden, ierror,
     &           indk, nbands, maxnhg 

      integer :: BNode, MPIerror,iio,jjo,noccloc

      real(dp) :: def, A, B, 
     &            dQo(maxorb), Qo,aux(maxorb), 
     &            psiper(maxorb,maxorb), prod1, prod2, 
     &            prod3, prod4, ei0, ej0, dStepF, evper(maxorb)

      real(dp), pointer :: Haux(:,:), Saux(:,:), Psiden(:,:),
     &                     rotaux(:)

      real(dp), pointer :: Hperg(:)=>null(),
     &                     Operg(:)=>null(),
     &                     psig(:,:)=>null(),
     &                     psitmp(:,:)=>null(),
     &                     Daux(:)=>null(),
     &                     DEaux(:)=>null()


      integer, pointer :: numhg(:)=>null(),
     &                    listhg(:)=>null(),
     &                    listhptrg(:)=>null()


#ifdef MPI
      real(dp),pointer :: rhoperg(:)=>null(),
     &                    Erhoperg(:)=>null()
#endif

      save :: maxden

      call timer('delrhog',1)

! For MPI process, all teh Nodes do the same job

! First, globalize everything that are distributed over nodes

!  eval is already in all the nodes

!  globalize psi to psig
      call re_alloc(psig,1,nuotot,1,nuotot,'psig','delrhog')
      call re_alloc( psitmp, 1, nuotot, 1,nuotot, 'psitmp','delrhog')

      do io=1,nuo
        call WhichNodeOrb(io,Nodes,BNode)
          call LocalToGlobalOrb(io,Node,Nodes,iio)
          do jo=1,nuotot
            psitmp(jo,iio)=psi(jo,io)
          enddo
      enddo
      psig(:,:)=psitmp(:,:) !direct copy in serial
#ifdef MPI
        if (Nodes .gt.1 ) then
           call globalize_sum(psitmp,psig)
        endif
#endif
      call de_alloc(psitmp)

! globalize numh to numhg
      nullify(numhg)
      call re_alloc( numhg, 1, nuotot, 'numhg','delrhog')
      do io = 1,nuotot
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
      call re_alloc( listhptrg, 1, nuotot, 'listhptrg', 'delrhog' )
      listhptrg(1) = 0
      do io = 2,nuotot
        listhptrg(io) = listhptrg(io-1) + numhg(io-1)
      enddo

! globalize listh to listhg
      maxnhg = listhptrg(nuotot) + numhg(nuotot)
      nullify(listhg)
      call re_alloc( listhg, 1, maxnhg, 'listhg', 'delrhog' )
      do io = 1,nuotot
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

! globalize Hper/Oper to Hperg/Operg
      nullify( Hperg, Operg)
      call re_alloc( Operg, 1, maxnhg, 'Operg', 'delrhog' )
      call re_alloc( Hperg, 1, maxnhg, 'Hperg', 'delrhog' )

      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do jo = 1,numh(iio)
            Hperg(listhptrg(io)+jo) =
     &                            Hper(listhptr(iio)+jo)
            Operg(listhptrg(io)+jo) =
     &                            Oper(listhptr(iio)+jo)
          enddo
        endif
#ifdef MPI
        if (Nodes .gt.1 ) then
          call MPI_Bcast(Hperg(listhptrg(io)+1),numhg(io),
     &        MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
          call MPI_Bcast(Operg(listhptrg(io)+1),numhg(io),
     &        MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
        endif
#endif
      enddo

      nbands = nuotot
      if(iscf.eq.1) then
       deg = 1
       N = 1
       maxden = 1
        do io = 1,nbands-1
         ei0 = eval(io)
         ej0 = eval(io+1)
         if(abs(ei0-ej0) .lt. eigtol) then
           N = N + 1 !Degenerate Eig
           if(N.gt.MAXDEN) MAXDEN = N
         else
           DEG = max(DEG,N)
           N = 1 !Non deg Eig
         endif
        enddo
      endif      
 
      def = 0.0_dp
      A = 0.0_dp
      B = 0.0_dp
      dQo(1:nbands) = 0.0_dp

      N = 0
      do io = 1, nbands
        qo = occ (io)
        if(qo .ge. 1.0e-5_dp) then
          ei0 = eval(io) 
          if(io .lt. nbands) then
            ej0 = eval(io+1)
          else
            ej0 = 1.0e7_dp
          endif
          ! building the dH_nn' = <psi_in|dH|psi_in'> 
          ! where n refers to the degenerate space.
          if(abs(ei0-ej0) .lt. eigtol) then
            N = N + 1
          else
            numb(io) = N + 1
            N = 0
            if((numb(io).gt.1) .and. (iscf.gt.1)) then
               ! initialize with size of deg subspace
               nullify(psiden,Haux,Saux)
               call re_alloc(Haux, 1,numb(io), 1,numb(io),
     &                       'Haux','delrhog')
               call re_alloc(Saux, 1,numb(io), 1,numb(io),
     &                       'Saux','delrhog')
               call re_alloc(Psiden, 1,numb(io), 1,numb(io),
     &                       'Psiden','delrhog')
               rotaux => Saux(:,1)

               ! Initialize to 0
               Saux = 0._dp
               
               do j = io-numb(io)+1, io
                jden = j - io + numb(io)
                numb(j) = numb(io)
                do k = io-numb(io)+1, io
                  kden = k - io + numb(io)
                  if ( jden == kden ) Saux(jden,kden) = 1.0_dp
                  do mu=1,nuotot 
                    do nu = 1,numhg(mu)
                      indmn=listhptrg(mu) + nu
!                      jo = ucorb(listh(indmn),nuo)
                      jo = listhg(indmn)    ! LET ME TRY THIS ALTERNATIVE
                       Haux(jden,kden) = Haux(jden,kden)  
     .                     + psig(mu,j) * psig(jo,k) * 
     .                     (Hperg(indmn)-eval(io)*Operg(indmn))
                    enddo
                  enddo
                enddo
              enddo
             
              ! Compute eigenvalues of dH_nn' to give dE_in
              call rdiag( Haux, Saux, numb(io), numb(io),numb(io),
     .             evper(io-numb(io)+1), psiden,numb(io),iscf,ierror,
     .             BlockSize)
              if (ierror.ne.0) then
                 call die('DELRHOG: Terminating due to failed
     &                     deg- diagonalisation')
              end if

#ifdef MPI
! psiden returned by rdiag only in ONE node
              if (Nodes.gt.1) call broadcast(psiden)
#endif
              ! Rotate eigenvectors in degenerate subspace
              do mu = 1,nuotot
               do ialpha = 1, numb(io)
                rotaux(ialpha) = 0._dp
                do ibeta = 1, numb(io)
                  rotaux(ialpha) = rotaux(ialpha) +
     &                   psiden(ibeta,ialpha) * 
     &                   psig(mu,io-numb(io)+ibeta)
                enddo
               enddo
               do ialpha = 1, numb(io)
                  psig(mu,io-numb(io)+ialpha)=rotaux(ialpha)
               enddo
              enddo
                    
             call de_alloc(psiden,'psiden','delrhog')
             call de_alloc(Haux,'Haux','delrhog')
             call de_alloc(Saux,'Saux','delrhog')
             
            endif !numb(io)>1 
         endif !abs
        endif !qo.ge.
      enddo

C Initialize the perturbed coefficients
      do j = 1,maxorb
       evper(j) = 0.0_dp 
       do i = 1,maxorb
        psiper(i,j) = 0.0_dp
       enddo
      enddo

C Calculate perturbed coefficients
      do 110 i=1,nuo
        call LocalToGlobalOrb(i,Node,Nodes,io)
        QO = OCC(io)
        if(QO .lt. 1.0e-5_dp) goto 110
        ei0 = eval(io)
        do j=1,nbands
         aux(j) = 0.0_dp
        enddo
        do ialpha = 1, nuotot
          prod1 = 0.0_dp
          prod2 = 0.0_dp
          do jbeta = 1, numhg(ialpha)
            indbetas = listhptrg(ialpha) + jbeta
            ibeta = ucorb(listhg(indbetas),nuotot)
            prod1 = prod1 + psig(ibeta,io) * 
     &                      Hperg(indbetas)
            prod2 = prod2 + psig(ibeta,io) *
     &                       Operg(indbetas)
          enddo
          prod3 = prod1 - ei0 * prod2
          prod4 = prod2 / 2.0_dp
          do j = 1, nbands
            ej0 = eval(j)
            if(io.ne.j) then 
             aux(j) = aux(j) + prod3 * psig(ialpha,j) 
     &                        / (ei0-ej0) 
            else
             aux(j) = aux(j) - prod4 * psig(ialpha,j)
            endif
          enddo
         evper(io) = evper(io) + prod3 * psig(ialpha,io)
        enddo

        do ialpha = 1,nuotot
          psiper(io,ialpha) = 0.0_dp
          do j = 1,nbands
           psiper(io,ialpha) = psiper(io,ialpha) +  
     .                        aux(j) * psig(ialpha,j) 
          enddo
        enddo
        dQo(io) = evper(io) * dStepF((ei0-ef)/T)/(T*nspin+1.0e-12_dp)
        A = A + evper(io) * dStepF((ei0-ef)/T)
        B = B + dStepF((ei0-ef)/T)
 110  enddo
      def = A / (B + 1.0e-12_dp)
      do i = 1, nuo
        call LocalToGlobalOrb(i,Node,Nodes,io)
        Qo = occ(io)
        if(Qo .gt. 1.0e-5_dp) then
          ei0 = eval(io)
          dQo(io) = dQo(io) - def * dStepF((ei0-ef)/T) / 
     &                      (T*nspin + 1.0e-12_dp)
        endif
      enddo
C dQo, Evper, psiper are distributed on the node in the form (...0... nuo)

C Up to here we have distributed over nodes the nbands, this means that:
C To compute the perturbed density matrix, we have to perform orbital
C multiplication in the SAME band and at the end, sum over all the bands
C Lets paralelize over bands: each node will contain some part of the
C final sum

C Calculate matrix elements of the perturbed Density Matrix
      call re_alloc(Daux,1,maxnhg,'Daux','delrhog')
      call re_alloc(DEaux,1,maxnhg,'DEaux','delrhog')

      Daux(:)=0.0_dp
      DEaux(:)=0.0_dp

C Calculate matrix elements of the perturbed Density Matrix
      do 1200 i = 1, nuo
        call LocalToGlobalOrb(i,Node,Nodes,io) !since occ is global
        Qo = occ(io)

        if(Qo .lt. 1.0e-5_dp) goto 1200
        ei0 = eval(io)
        do 1400 mu = 1, nuotot
          do 1300 nu = 1, numhg(mu)
            indmn = listhptrg(mu) + nu
            j = ucorb(listhg(indmn),nuotot)
            Daux(indmn) = Daux(indmn) + Qo *  
     &                  (psiper(io,mu) * psig(j,io) +
     &                  psig(mu,io) * psiper(io,j)) +
     &                   dQo(io) *(psig(mu,io) * psig(j,io))
            DEaux(indmn) = DEaux(indmn) +
     .                           qo* ei0*
     .                           ( psiper(io,mu) * psig(j,io) +
     .                             psig(mu,io) * psiper(io,j))+
     .               evper(io)*qo * ( psig(mu,io)*psig(j,io)) 
     .             + dQo(i) * ei0 * ( psig(mu,io) * psig(j,io) )
 1300     enddo
 1400   enddo
 1200 enddo 
      
! Sum over all nodes (the bands)
! The final sum is already distributed, so, is GLOBAL.

! Do it in serial (just copy)
        rhoper(:) = Daux(:)
        erhoper(:) = DEaux(:)

#ifdef MPI
        if (Nodes .gt.1 ) then
          call re_alloc(rhoperg,1,maxnhg,'rhoperg','delrhog')
          call re_alloc(Erhoperg,1,maxnhg,'Erhoperg','delrhog')

          call Globalize_sum(Daux,rhoperg) !shared in all nodes
          call Globalize_sum(DEaux,Erhoperg) !shared in all nodes

! Final step, distribute the density/energy matrix over all nodes
          do io = 1, nuotot
            call WhichNodeOrb(io,Nodes,BNode)
            if (Node .eq. BNode) then
              call GlobalToLocalOrb(io,BNode,Nodes,iio)
              indmn=listhptrg(io) !just copy some global parts into local ones
              do j = 1,numh(iio)
                rhoper(listhptr(iio)+j) = rhoperg(indmn+j)
                Erhoper(listhptr(iio)+j) = Erhoperg(indmn+j)
              enddo
            endif
          enddo
          call de_alloc(rhoperg)
          call de_alloc(Erhoperg)

        else ! MPI with 1 node, just copy
          rhoper(:) = Daux(:)
          erhoper(:) = DEaux(:)
        endif
#endif

      call de_alloc(Daux)
      call de_alloc(DEaux)
      call de_alloc(numhg)
      call de_alloc(listhptrg)
      call de_alloc(listhg)
      call de_alloc(Hperg)
      call de_alloc(Operg)

      call timer('delrhog',2)
      return
      contains

      elemental function ucorb(a,p)
      integer, intent(in) :: a,p
      integer :: ucorb
      if ( a > p ) then
         ucorb = MOD(a,p)
         if ( ucorb == 0 ) ucorb = p
      else
         ucorb = a
      end if
      end function
      END

c ===================================================================
      DOUBLE PRECISION FUNCTION DSTEPF(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     Complementary error function. Ref: Fu & Ho, PRB 28, 5480 (1983)
*     STEPF=DERFC(X)

C     Improved step function. Ref: Methfessel & Paxton PRB40 (15/Aug/89)
*     PARAMETER (C=0.5641895835D0)
*     STEPF=DERFC(X)-C*X*DEXP(-X*X)

C     Fermi-Dirac distribution
      IF (X.GT.100.0D0) THEN
        DSTEPF = 0.0D0
      ELSEIF (X.LT.-100.0D0) THEN
        DSTEPF = 0.0D0
      ELSE
        DSTEPF = - 2.0D0*EXP(X) / ( (1.0D0 + EXP(X))*(1.0D0 +EXP(X)))
      ENDIF

      END

