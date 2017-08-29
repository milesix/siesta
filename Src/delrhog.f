      SUBROUTINE DELRHOG(nuo, nspin, maxorb, eval, tol, eigtol,  
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
      implicit none

      integer :: nuo, nspin, maxorb, maxnh, numh(*),
     &            listh(maxnh), listhptr(*), ix, iscf
      real(dp) :: eval(maxorb), tol, occ(maxorb), 
     &            Hper(maxnh), eigtol, 
     &            Oper(maxnh), ef, T,Rhoper(maxnh),Erhoper(maxnh)
      real(dp) :: psi(nuo,nuo)
C     Internal Variables
      integer :: deg, N, io, numb(maxorb), k,j, maxden,
     &           jden, kden, jo, mu, nu, ind, indmn, i, indbetas,
     &           ialpha, ibeta, jbeta,indi,indj,indden, ierror,
     &           indk, nbands 
      real(dp) :: def, A, B, 
     &            dQo(maxorb), Qo,aux(maxorb), 
     &            psiper(maxorb,maxorb), prod1, prod2, 
     &            prod3, prod4, ei0, ej0, dStepF, evper(maxorb)
      real(dp), pointer :: Haux(:,:), Saux(:,:), Psiden(:,:),
     &                     rotaux(:), eden(:)
      save :: maxden

      call timer('delrhog',1)

      nbands = nuo
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
      ! initialize with size of deg subspace
!      call re_alloc(Haux, 1,maxden, 1,maxden,'Haux', 'delrhog')
!      call re_alloc(Saux, 1,maxden, 1,maxden,'Saux', 'delrhog')
!      call re_alloc(Psiden, 1,maxden, 1,maxden,'Psiden', 'delrhog')
!      call re_alloc(rotaux, 1,maxden,'rotaux', 'delrhog')
!      call re_alloc(eden, 1,maxden,'eden', 'delrhog')
 
      def = 0.0_dp
      A = 0.0_dp
      B = 0.0_dp
!      rotaux(1:maxden) = 0.0_dp
!      Haux(1:maxden,1:maxden) = 0.0_dp
!      Saux(1:maxden,1:maxden) = 0.0_dp
!      numb(1:nbands) = 0
      dQo(1:nbands) = 0.0_dp

      N = 0
      do 520 io = 1, nbands
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
      nullify(psiden,eden,rotaux,Haux,Saux)
      call re_alloc(Haux, 1,numb(io), 1,numb(io),'Haux', 'delrhog')
      call re_alloc(Saux, 1,numb(io), 1,numb(io),'Saux', 'delrhog')
      call re_alloc(Psiden, 1,numb(io), 1,numb(io),'Psiden', 'delrhog')
      call re_alloc(rotaux, 1,numb(io),'rotaux', 'delrhog')
      call re_alloc(eden, 1,numb(io),'eden', 'delrhog')
              do j = io-numb(io)+1, io
                jden = j - io + numb(io)
                numb(j) = numb(io)
                do k = io-numb(io)+1, io
                  kden = k - io + numb(io)
                  if(jden .eq. kden) Saux(jden,kden) = 1.0_dp
                  do mu=1,nuo 
                    do nu = 1,numh(mu)
                      indmn=listhptr(mu) + nu
!                      jo = ucorb(listh(indmn),nuo)
                      jo = listh(indmn)    ! LET ME TRY THIS ALTERNATIVE
                       Haux(jden,kden) = Haux(jden,kden)  
     .                     + psi(mu,j) * psi(jo,k) * 
     .                     (Hper(indmn)-eval(io)*Oper(indmn))
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
        ! Compute eigenvalues of dH_nn' to give dE_in
!              call rdiag( Haux, Saux, maxden, numb(io),maxden,
!     .             eden, psiden,numb(io),iscf,ierror )
              call rdiag( Haux, Saux, numb(io), numb(io),numb(io),
     .             eden, psiden,numb(io),iscf,ierror )

! Rotate eigenvectors in degenerate subspace
              do mu = 1,nuo
               do ialpha = 1, numb(io)
                do ibeta = 1, numb(io)
                  rotaux(ialpha) = rotaux(ialpha) +
     &                   psiden(ibeta,ialpha) * 
     &                   psi(mu,io-numb(io)+ibeta)
                enddo
               enddo
               do ibeta = 1, numb(io)
                  psi(mu,io-numb(io)+ibeta)=rotaux(ibeta)
                  rotaux(ibeta) = 0.D0
               enddo
              enddo
              do jo = 1, numb(io)
                evper(io-numb(io)+jo) = eden(jo)
              enddo
                    
!              do jo = 1, numb(io)
!               do k = 1, numb(io)
!                Haux(k,jo) = 0.D0
!                Saux(k,jo) = 0.D0
!               enddo
!              enddo
             call de_alloc(psiden,'psiden','delrhog')
             call de_alloc(rotaux,'rotaux','delrhog')
             call de_alloc(eden,'eden','delrhog')
             call de_alloc(Haux,'Haux','delrhog')
             call de_alloc(Saux,'Saux','delrhog')
            endif !numb(io)>1 
         endif !abs
        endif !qo.ge.
 520  enddo

C Initialize the perturbed coefficients
      do j = 1,maxorb
       evper(j) = 0.0_dp 
       do i = 1,maxorb
        psiper(i,j) = 0.0_dp
       enddo
      enddo 
      
C Calculate perturbed coefficients
      do 110 i=1,nbands
        QO = OCC(i)
        if(QO .lt. 1.0e-5_dp) goto 110
        ei0 = eval(i)
        do j=1,nbands
         aux(j) = 0.0_dp
        enddo
        do ialpha = 1, nuo
          prod1 = 0.0_dp
          prod2 = 0.0_dp
          do jbeta = 1, numh(ialpha)
            indbetas = listhptr(ialpha) + jbeta
            ibeta = ucorb(listh(indbetas),nuo)
            prod1 = prod1 + psi(ibeta,i) * 
     &                      Hper(indbetas)
            prod2 = prod2 + psi(ibeta,i) *
     &                       Oper(indbetas)
          enddo
          prod3 = prod1 - ei0 * prod2
          prod4 = prod2 / 2.0_dp
          do j = 1, nbands
            ej0 = eval(j)
            if(i.ne.j) then 
             aux(j) = aux(j) + prod3 * psi(ialpha,j) 
     &                        / (ei0-ej0) 
            else
             aux(j) = aux(j) - prod4 * psi(ialpha,j)
            endif
          enddo
         evper(i) = evper(i) + prod3 * psi(ialpha,i)
        enddo
        do ialpha = 1,nuo
          psiper(i,ialpha) = 0.0_dp
          do j = 1,nbands
           psiper(i,ialpha) = psiper(i,ialpha) +  
     .                        aux(j) * psi(ialpha,j) 
          enddo
        enddo
        dQo(i) = evper(i) * dStepF((ei0-ef)/T)/(T*nspin+1.0e-12_dp)
        A = A + evper(i) * dStepF((ei0-ef)/T)
        B = B + dStepF((ei0-ef)/T)
 110  enddo
      def = A / (B + 1.0e-12_dp)
      do i = 1, nbands
        Qo = occ(i)
        if(Qo .gt. 1.0e-5_dp) then
          ei0 = eval(i)
          dQo(i) = dQo(i) - def * dStepF((ei0-ef)/T) / 
     &                      (T*nspin + 1.0e-12_dp)
        endif
      enddo

C Calculate matrix elements of the perturbed Density Matrix
      do 1200 i = 1, nbands
        Qo = occ(i)
        if(Qo .lt. 1.0e-5_dp) goto 1200
        ei0 = eval(i)
        do 1400 mu = 1, nuo
          do 1300 nu = 1, numh(mu)
            indmn = listhptr(mu) + nu
            j = ucorb(listh(indmn),nuo)
            rhoper(indmn) = rhoper(indmn) + Qo *  
     &                  (psiper(i,mu) * psi(j,i) +
     &                  psi(mu,i) * psiper(i,j)) +
     &                   dQo(i) *(psi(mu,i) * psi(j,i))
            erhoper(indmn) = erhoper(indmn) +
     .                           qo* ei0*
     .                           ( psiper(i,mu) * psi(j,i) +
     .                             psi(mu,i) * psiper(i,j))+
     .               evper(i)*qo * ( psi(mu,i)*psi(j,i)) 
     .             + dQo(i) * ei0 * ( psi(mu,i) * psi(j,i) )
 1300     enddo
 1400   enddo
 1200 enddo 
      
      call de_alloc( Haux, 'Haux', 'delrhog' )
      call de_alloc( Saux, 'Saux', 'delrhog' )
      call de_alloc( Psiden, 'Psiden', 'delrhog' )
      call de_alloc( rotaux, 'rotaux', 'delrhog' )
      call de_alloc( eden, 'eden', 'delrhog' )

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

