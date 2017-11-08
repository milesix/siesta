! 
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!      


      module m_dvnloc

      implicit none

      public :: dvnloc

      private

      CONTAINS

      subroutine dvnloc( scell, nua, na, isa, xa, indxua,Dscf,
     .                   maxnh, maxnd, lasto, lastkb, iphorb, 
     .                   iphKB, numd, listdptr, listd, numh, 
     .                   listhptr, listh, nspin, IALR, norb, iaorb, 
     .                   dDscf, dHmat, dynmat, FIRST) 
C *********************************************************************
C Energies in Ry. Lengths in Bohr.
c Module created for LinRes. Adapted version from M. Pruneda. 
C Based on NLEFSM (Writen by J.Soler and P.Ordejon, June 1997) 
C with first and second derivatives contributions of non-local potential
C contribution needed for the calculation of the perturbed hamiltonian 
C as well as elements of the dynamical matrix. 
C Adapted for Linres by L.Riches and S.Illera (2016).
C **************************** INPUT **********************************
C real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
C integer nua              : Number of atoms in unit cell
C integer na               : Number of atoms in supercell
C integer isa(na)          : Species index of each atom
C real*8  xa(3,na)         : Atomic positions in cartesian coordinates
C integer indxua(na)       : Index of equivalent atom in unit cell
C integer maxnh            : First dimension of H and listh
C integer maxnd            : Maximum number of elements of the
C                            density matrix
C integer lasto(0:na)      : Position of last orbital of each atom
C integer lastkb(0:na)     : Position of last KB projector of each atom
C integer iphorb(no)       : Orbital index of each orbital in its atom,
C                            where no=lasto(na)
C integer iphKB(nokb)      : Index of each KB projector in its atom,
C                            where nokb=lastkb(na)
C integer numd(nuo)        : Number of nonzero elements of each row of the
C                            density matrix
C integer listdptr(nuo)    : Pointer to the start of each row (-1) of the
C                            density matrix
C integer listd(maxnd)     : Nonzero hamiltonian-matrix element column 
C                            indexes for each matrix row
C integer numh(nuo)        : Number of nonzero elements of each row of the
C                            hamiltonian matrix
C integer listhptr(nuo)    : Pointer to the start of each row (-1) of the
C                            hamiltonian matrix
C integer listh(maxnh)     : Nonzero hamiltonian-matrix element column 
C                            indexes for each matrix row
C integer nspin            : Number of spin components of Dscf and H
C                            If computing only matrix elements, it
C                            can be set to 1
C real*8 Dscf(maxnd,nspin) : Ground state Density matrix. Only used to 
C			     compute dynamical matrix terms
C real*8 dDscf(maxnd,nspin): Perturbed denisty matrix. Only used to compute 
C 			     dynamical matrix terms
C integer IALR		   : atom which is perturbed
C logical FIRST		   : flag
C			TRUE= compute elements of the perturbed Hamiltonian
C			FALSE= compute dynamical matrix terms (non-SCF and SCF terms)
C ******************* INPUT and OUTPUT *********************************
c dynmat(3,nua,3,nua)	: Dynamical Matrix
C dHMAT(maxnh,3, nspin) : Perturbed Hamiltonian Matrix Elements

C  Modules
      use precision,     only : dp
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : GetNodeOrbs, LocalToGlobalOrb
      use parallelsubs,  only : GlobalToLocalOrb
      use atmfuncs,      only : rcut, epskb, orb_gindex,kbproj_gindex,
     &                          nofis, nkbfis
      use atm_types,     only : nspecies 
      use neighbour,     only : iana=>jan, r2ki=>r2ij, xki=>xij
      use neighbour,     only : mneighb, reset_neighbour_arrays
      use alloc,         only : re_alloc, de_alloc
      use m_new_matel,   only : new_matel
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif


      integer, intent(in)     ::
     .   maxnh, na, maxnd, nspin, nua, IALR

      integer, intent(in)     ::
     .  indxua(na), iphKB(*), iphorb(*), isa(na),  
     .  lasto(0:na), lastkb(0:na), listd(maxnd), listh(maxnh),
     .  numd(*), numh(*), listdptr(*), listhptr(*)

      real(dp), intent(in)   :: scell(3,3), Dscf(maxnd,nspin),
     .                        xa(3,na), dDscf(maxnh,nspin,3)
      logical, Intent(in)    :: first
      real(dp), intent(inout) ::dynmat(3,nua,3,nua),
     &                         dHMAT(maxnh,3, nspin)

      real(dp) ::   volcel
      external ::   timer, volcel

C Internal variables ................................................

      integer, save ::  maxno = 2000
  
      integer
     .  ia, ikb, ina, ind, ino, iua, ja,
     .  io, iio, ioa, is, ispin, ix, jua, 
     .  j, jno, jo, jx, ka, ko, koa, ks, kua,
     .  nkb, nna, nno, no, nuo, nuotot, maxkba,
     .  norb, iaorb(norb), ig, jg

      integer, dimension(:), pointer :: iano, iono

      real(dp)
     .  Cijk, epsk, fik, rki, rmax, rmaxkb, rmaxo,prot, 
     .  Sik, Sjk, volume, dCijk(3), Htest(maxnh,nspin)

      real(dp), dimension(:), pointer :: Di, Vi
      real(dp), dimension(:,:), pointer :: Ski, xno
      real(dp), dimension(:,:,:), pointer :: grSki
      real(dp), dimension(:,:), pointer :: dDi, dVi
      real(dp), dimension(:,:,:,:), pointer :: gr2Ski, d2Vi
      real(dp), dimension(:,:,:,:), pointer ::  dynmatl
      real(dp), dimension(:,:,:,:), pointer :: t_DYL

      logical ::   within
      logical, dimension(:), pointer ::  listed, listedall

      real(dp), allocatable :: rorbmax(:), rkbmax(:)

C ------------------------------------------------------------------

C Linres ------------------------------------------------------------
      INTERFACE
        SUBROUTINE MATEL(A, ioa, joa, rij,
     &                   Sij, grSij, gr2Sij)
          use precision,     only : dp
             integer,   intent(inout) ::  ioa, joa
             real(dp),  intent(in) :: rij, Sij
             real(dp),  intent(out) :: grSij(3)
             real(dp),  intent(out), optional ::  gr2Sij(3,3)
             character, intent(in) :: A
        END SUBROUTINE MATEL
      END INTERFACE
C -------------------------------------------------------------------

      call timer('dvnloc',1)
      
C Find unit cell volume
      volume = volcel( scell ) * nua / na

C Find maximum range
      maxkba = 0

      allocate(rorbmax(nspecies),rkbmax(nspecies))
      do is = 1, nspecies

         ! Species orbital range
         rorbmax(is) = 0.0_dp
         do io = 1, nofis(is)
            rorbmax(is) = max(rorbmax(is), rcut(is,io))
         enddo

         ! Species KB range
         io = nkbfis(is)
         rkbmax(is) = 0.0_dp
         do ikb = 1, io
            rkbmax(is) = max(rkbmax(is), rcut(is,-ikb))
         enddo
         maxkba = max(maxkba,io)

      enddo
      rmaxo = maxval(rorbmax(1:nspecies))
      rmaxkb = maxval(rkbmax(1:nspecies))
      ! Calculate max extend
      rmax = rmaxo + rmaxkb

C Initialize arrays Di and Vi only once
      no = lasto(na)
      nuotot = lasto(nua)
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)

C Allocate local memory

      nullify( Vi )
      call re_alloc( Vi, 1, no, 'Vi', 'dvnloc' )
      Vi(1:no) = 0.0_dp
      nullify( listed )
      call re_alloc( listed, 1, no, 'listed', 'dvnloc' )
      listed(1:no) = .false.
      nullify( listedall )
      call re_alloc( listedall, 1, no, 'listedall', 'dvnloc' )
      listedall(1:no) = .false.
      nullify( Di )
      call re_alloc( Di, 1, no, 'Di', 'dvnloc' )
      Di(1:no) = 0.0_dp

      nullify(dVi)
      call re_alloc( dVi, 1, no, 1, 3, 'dVi', 'dvnloc' )
      dVi(1:no,1:3) = 0.0_dp

      nullify(dDi)
      call re_alloc( dDi, 1, no, 1, 3, 'dDi', 'dvnloc' )
      dDi(1:no,1:3) = 0.0_dp

      nullify(d2Vi)
      call re_alloc( d2Vi, 1, no, 1, 3, 1, 3, 1, nua, 
     &              'd2Vi', 'dvnloc' )
      d2Vi(1:no,1:3,1:3,1:nua) = 0.0_dp

      if (.not. first) then
        nullify(dynmatl) !stores local dvnloc dynamat element
        call re_alloc( dynmatl, 1, 3, 1, nua, 1, 3, 1, nua,
     &                                  'dynmatl','dvnloc')
        dynmatl(:,:,:,:)=0.0_dp
#ifdef MPI 
! For parallel purpose: will store the global sum
          nullify(t_DYL)
          if (Nodes.gt.1) then
            call re_alloc( t_DYL, 1, 3, 1, nua, 1, 3, 1, nua,
     &                                  't_DYL','dvnloc')
            t_DYL(:,:,:,:)=0.0_dp
          endif
#endif
      endif

C Make list of all orbitals needed for this node
      do io = 1,nuo
        call LocalToGlobalOrb(io,Node,Nodes,iio)
        listedall(iio) = .true.
        do j = 1,numh(io)
          jo = listh(listhptr(io)+j)
          listedall(jo) = .true.
        enddo
      enddo

C Allocate local arrays that depend on saved parameters
      nullify( iano )
      call re_alloc( iano, 1, maxno, 'iano', 'dvnloc' )
      nullify( iono )
      call re_alloc( iono, 1, maxno, 'iono', 'dvnloc' )
      nullify( xno )
      call re_alloc( xno, 1, 3, 1, maxno, 'xno',  'dvnloc' )
      nullify( Ski )
      call re_alloc( Ski, 1, maxkba, 1, maxno, 'Ski', 'dvnloc' )
      nullify( grSki )
      call re_alloc( grSki, 1, 3, 1, maxkba, 1, maxno, 'grSki',
     &               'dvnloc' )
      nullify( gr2Ski )
      call re_alloc( gr2Ski, 1, 3, 1, 3, 1, maxkba, 1, maxno, 
     &               'gr2Ski', 'dvnloc' )

C     Initialize neighb subroutine
      call mneighb( scell, rmax, na, xa, 0, 0, nna )

C     Loop on atoms with KB projectors      
      do ka = 1,na !all the atoms
        kua = indxua(ka) !index atom in the u. cell
        ks = isa(ka)  !atom specie
        nkb = lastkb(ka) - lastkb(ka-1) !indices kb projectors of atom

C       Find neighbour atoms
        call mneighb( scell, rmax, na, xa, ka, 0, nna )
C       Find neighbour orbitals
        nno = 0
        do ina = 1,nna
          ia = iana(ina)
          is = isa(ia)
          rki = sqrt(r2ki(ina))
          if (rki - rkbmax(ks) - rorbmax(is) > 0.d0) CYCLE

          ! Loop over orbitals close enough to overlap
          do io = lasto(ia-1)+1,lasto(ia)

C           Only calculate if needed locally
            if (.not. listedall(io)) CYCLE

            ioa = iphorb(io)

C           Find if orbital is within range
            within = (rki-rkbmax(ks)) < rcut(is,ioa)
            if (.not. within) CYCLE

C          Find overlap between neighbour orbitals and KB projectors
C          Check maxno - if too small then increase array sizes
            if (nno.eq.maxno) call increase_maxno()

            nno = nno + 1 !index orbital
            iono(nno) = io
            iano(nno) = ia
            do ix = 1,3
              xno(ix,nno) = xki(ix,ina)
            enddo
                
            ikb = 0
            ig=orb_gindex(is,ioa)
            do ko = lastkb(ka-1)+1,lastkb(ka)
              ikb = ikb + 1 !index orbital
              koa = iphKB(ko)
              jg=kbproj_gindex (ks,koa)
              call new_MATEL( 'S',jg,ig, xki(1:3,ina),
     &                  Ski(ikb,nno), grSki(1:3,ikb,nno),
     &                  gr2Ski(1:3,1:3,ikb,nno) )
            enddo
          enddo
        enddo !ina do (all neighbour atoms)

C       Loop on neighbour orbitals
        do ino = 1,nno
          io = iono(ino)
          ia = iano(ino)
          iua = indxua(ia)
          if (ia>nua) CYCLE

          call GlobalToLocalOrb(io,Node,Nodes,iio)
          if (iio==0) CYCLE
C         Valid orbital
C         Scatter density matrix row of orbital io
          do j = 1,numd(iio)
            ind = listdptr(iio)+j
            jo = listd(ind)
            do ispin = 1,nspin
              do ix =1,3
                dDi(jo,ix) = dDi(jo,ix) + dDscf(ind,ispin,ix)
              enddo
            enddo
          enddo

          do j = 1,numh(iio)
            jo = listh(listhptr(iio)+j)
            listed(jo) = .true.
          enddo
        
C     Find matrix elements with other neighbour orbitals
          do jno = 1,nno
            jo = iono(jno)
            ja = iano(jno)
            jua = indxua(ja)
            if ( .not. listed(jo)) CYCLE

C           Loop on KB projectors
            ikb = 0
            do ko = lastkb(ka-1)+1,lastkb(ka)
              ikb = ikb + 1
              koa = iphKB(ko)
              epsk = epskb(ks,koa)
              Sik = Ski(ikb,ino)
              Sjk = Ski(ikb,jno)
C Three cases:
C Case A ...        
              if( kua .eq. ialr ) then
                do ix = 1,3
                  dVi(jo,ix) = dVi(jo,ix) 
     &                       - epsk*Sik*grSki(ix,ikb,jno) 
     &                       - epsk*Sjk*grSki(ix,ikb,ino)
                  do jx = 1,3
                    if(ia .ne. ka) then
                      d2Vi(jo,ix,jx,ia) = d2Vi(jo,ix,jx,ia) -
     &                    epsk*gr2Ski(jx,ix,ikb,ino)*Sjk
                    endif

                    if((ia .ne. ka) .and. (ja .ne. ka)) then
                      d2Vi(jo,ix,jx,ia) = d2Vi(jo,ix,jx,ia) -
     &                   epsk*grSki(jx,ikb,ino)*grSki(ix,ikb,jno)
                    endif
                           
                    if(ja .ne. ka) then
                      d2Vi(jo,ix,jx,jua) = d2Vi(jo,ix,jx,jua)
     &                    - epsk*Sik*gr2Ski(ix,jx,ikb,jno)
                    endif
                          
                    if((ja .ne. ka) .and. (ia .ne. ka)) then
                      d2Vi(jo,ix,jx,jua) = d2Vi(jo,ix,jx,jua)
     &                   - epsk*grSki(ix,ikb,ino)*grSki(jx,ikb,jno)
                    endif

                    if(ia .ne. ka) then
                      d2Vi(jo,ix,jx,kua) = d2Vi(jo,ix,jx,kua)
     &                   + epsk*gr2Ski(ix,jx,ikb,ino)*Sjk
                    endif

                    if(ja.ne.ka) then
                      d2Vi(jo,ix,jx,kua) = d2Vi(jo,ix,jx,kua)
     &                   + epsk*gr2Ski(ix,jx,ikb,jno)*Sik
                    endif
 
                    if((ia.ne.ka) .and. (ja .ne. ka)) then
                      d2Vi(jo,ix,jx,kua) = d2Vi(jo,ix,jx,kua)
     &                   + epsk*grSki(ix,ikb,ino)*grSki(jx,ikb,jno)
     &                   + epsk*grSki(jx,ikb,ino)*grSki(ix,ikb,jno)
                    endif
                  enddo !jx loop
                enddo ! ix loop
              endif !kua==ialr
C Case B ...
              if(iua .eq. ialr) then
                do ix = 1,3
                  dVi(jo,ix) = dVi(jo,ix) +
     &                 epsk*Sjk*grSki(ix,ikb,ino)
                       
                  do jx = 1,3
                    if(ka.ne.ia) then
                      d2Vi(jo,ix,jx,ia) = d2Vi(jo,ix,jx,ia) +
     &                      epsk*gr2Ski(ix,jx,ikb,ino)*Sjk
                    endif

                    if((ka .ne. ia) .and. (ka .ne. ja)) then
                      d2Vi(jo,ix,jx,jua) = d2Vi(jo,ix,jx,jua) 
     &                      + epsk*grSki(ix,ikb,ino)*
     &                      grSki(jx,ikb,jno) 
                    endif
                           
                    if(ka .ne. ia) then
                      d2Vi(jo,ix,jx,kua) = d2Vi(jo,ix,jx,kua)
     &                    - epsk*gr2Ski(ix,jx,ikb,ino)*Sjk
                    endif
                           
                    if((ia .ne. ka) .and. (ja .ne. ka)) then
                      d2Vi(jo,ix,jx,kua) = d2Vi(jo,ix,jx,kua)
     &                     - epsk*grSki(ix,ikb,ino)*
     &                       grSki(jx,ikb,jno)
                    endif
                  enddo !jx loop
                enddo !ix loop
              endif !iua==ialr
C Case C ...
              if( jua .eq. ialr) then
                do ix = 1, 3
                  dVi(jo,ix) = dVi(jo,ix) +
     &                         epsk*Sik*grSki(ix,ikb,jno)

                  do jx = 1, 3
                    if((ia.ne.ka).and.(ja.ne.ka)) then
                      d2Vi(jo,ix,jx,ia) = d2Vi(jo,ix,jx,ia)
     &                  + epsk*grSki(jx,ikb,ino)*grSki(ix,ikb,jno)
                    endif

                    if((ja.ne.ka)) then
                      d2Vi(jo,ix,jx,jua) = d2Vi(jo,ix,jx,jua)
     &                   + epsk*gr2Ski(ix,jx,ikb,jno)*Sik
                    endif

                    if( (ja.ne.ka) ) then
                      d2Vi(jo,ix,jx,kua) = d2Vi(jo,ix,jx,kua)
     &                   - epsk*gr2Ski(jx,ix,ikb,jno)*Sik
                    endif

                    if((ka.ne.ia).and.(ka.ne.ja))then
                      d2Vi(jo,ix,jx,kua) = d2Vi(jo,ix,jx,kua)
     &                  - epsk*grSki(jx,ikb,ino)*grSki(ix,ikb,jno)
                    endif
                  enddo !jx loop
                enddo !ix loop
              endif !jua ==ialr
C ... end cases
C Computing dynmat terms ------------------------------------------------
              if (.not. first) then 
                do ix = 1,3
                  dCijk(ix) = dDi(jo,ix) * epsk
                  do jx = 1,3
                    dynmatl(jx,ia,ix,ialr) = dynmatl(jx,ia,ix,ialr)
     &                - 2.0_dp * dCijk(ix) * Sjk * grSki(jx,ikb,ino)
                    dynmatl(jx,kua,ix,ialr) = dynmatl(jx,kua,ix,ialr)
     &                + 2.0_dp * dCijk(ix) * Sjk * grSki(jx,ikb,ino)
                       ! this is ok 
                  enddo
                enddo
              endif 
C------------------------------------------------------------------------
            enddo !ko loop
          enddo !jno loop

C         Pick up contributions to H and restore Di and Vi
          do j = 1,numh(iio)
            ind = listhptr(iio)+j
            jo = listh(ind)
            ja = iaorb(jo)
            jua = indxua(ja)
            do ispin = 1,nspin
              do ix = 1,3
                dHMAT(ind,ix,ispin) = dHMAT(ind,ix,ispin)
     &                                    + dVi(jo,ix)
C Computing dynmat terms ------------------------------------------------
                if(.not.first) then
                  do jx = 1, 3
                    dynmatl(jx,ia,ix,ialr) = dynmatl(jx,ia,ix,ialr)
     &                            - d2Vi(jo,ix,jx,ia)*dscf(ind,ispin)
                    if((ja.ne.ia).and.(ja.ne.ka)) then
                       dynmatl(jx,jua,ix,ialr) = dynmatl(jx,jua,ix,ialr)
     &                            - d2Vi(jo,ix,jx,jua)*dscf(ind,ispin)
                    endif
                    if( (ka.ne.ia) ) then
                       dynmatl(jx,kua,ix,ialr) = dynmatl(jx,kua,ix,ialr)
     .                             - d2Vi(jo,ix,jx,kua)*dscf(ind,ispin)
                        !this is ok
                    endif
                  enddo
                endif
C------------------------------------------------------------------------
              enddo !ix loop
            enddo !ispin loop
            listed(jo) = .false.
            dVi(jo,:) = 0.d0
            dDi(jo,:) = 0.d0
            d2Vi(jo,:,:,ia) = 0.d0
            d2Vi(jo,:,:,jua) = 0.d0
            d2Vi(jo,:,:,kua) = 0.d0
          enddo !j loop
          do j = 1,numd(iio)
            jo = listd(listdptr(iio)+j)
            Di(jo) = 0.0d0
          enddo
        enddo !ino loop
      enddo !ka loop
    
C Add the calculated dynmatl to the input dynmat
      if (.not. first) then
        if (Nodes .gt. 1) then
#ifdef MPI
! dynmatl is calculated locally in each node, globalize and add to the
! dynmat input
            do ix=1,3
              do jx=1,3
                do j=1,nua
                  call globalize_sum(dynmatl(ix,1:nua,jx,j),
     &                                   t_DYL(ix,1:nua,jx,j))
                enddo
              enddo
            enddo
#endif
            dynmat(:,:,:,:)=dynmat(:,:,:,:)+t_DYL(:,:,:,:)
          else
            dynmat(:,:,:,:)=dynmat(:,:,:,:)+dynmatl(:,:,:,:) !add dvnloc
!dynmat to input 
          endif
      endif

C     Deallocate local memory
      call reset_neighbour_arrays( )
      call de_alloc( grSki, 'grSki', 'dvnloc' )
      call de_alloc( Ski, 'Ski', 'dvnloc' )
      call de_alloc( xno, 'xno', 'dvnloc' )
      call de_alloc( iono, 'iono', 'dvnloc' )
      call de_alloc( iano, 'iano', 'dvnloc' )
      call de_alloc( listedall, 'listedall', 'dvnloc' )
      call de_alloc( listed, 'listed', 'dvnloc' )
      call de_alloc( Vi, 'Vi', 'dvnloc' )
      call de_alloc( Di, 'Di', 'dvnloc' )
      call de_alloc( dVi, 'dVi', 'dvnloc' )
      call de_alloc( dDi, 'dDi', 'dvnloc' )
      call de_alloc( d2Vi, 'd2Vi', 'dvnloc' )
      call de_alloc( gr2Ski, 'gr2Ski', 'dvnloc' )

      if (.not.first) then
        call de_alloc( dynmatl, 'dynmatl', 'dvnloc' )
#ifdef MPI
        if (Nodes.gt.1) then
          call de_alloc( t_DYL, 't_DYL', 'dvnloc' )
        endif
#endif
      endif

      call timer('dvnloc',2)

      CONTAINS

         subroutine increase_maxno()

! if too small then increase array sizes
            maxno = maxno + 100
            call re_alloc( iano, 1, maxno, 'iano', 'dvnloc',
     &                           .true. )
            call re_alloc( iono, 1, maxno, 'iono', 'dvnloc',
     &                           .true. )
            call re_alloc( xno, 1, 3, 1, maxno, 'xno',  'dvnloc',
     &                           .true. )
            call re_alloc( Ski, 1, maxkba, 1, maxno, 'Ski',
     &                          'dvnloc', .true. )
            call re_alloc( grSki, 1, 3, 1, maxkba, 1, maxno,
     &                          'grSki', 'dvnloc', .true. )
            call re_alloc( gr2Ski, 1, 3, 1, 3, 1, maxkba,
     &                    1, maxno,'gr2Ski', 'dvnloc', .true. )
         end subroutine increase_maxno

      end subroutine dvnloc

      end module m_dvnloc
