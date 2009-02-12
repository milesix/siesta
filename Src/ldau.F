! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module m_ldau


      use m_fdf_global
      use precision, only: dp
#ifdef MPI
            use m_mpi_utils, only: globalize_sum
#endif

      implicit none
      
      real(dp), save :: dDmax_current
      logical, save  :: ldau
      logical, save  :: conv_ldaupop
      data ldau /.false./
      data conv_ldaupop /.false./
         
      public :: ldau
      public :: conv_ldaupop
      public :: dDmax_current
      public :: hubbard_term
    
      private  
      real(dp), save  :: dtol_ldaupop
      integer, save   :: iter_ldaupop
      logical, save   :: ldau_init
      real(dp), save  :: dDmax_threshold
      real(dp), parameter  :: dtol_ldaupop_default =1.0e-3_dp
      real(dp), parameter  ::  dDmax_threshold_default= 1.0e-2_dp
      data iter_ldaupop / 0 /
       
      CONTAINS

      subroutine hubbard_term( scell, nua, na, isa, xa, indxua,
     .                   maxnh, maxnd, lasto, iphorb, 
     .                   numd, listdptr, listd, numh, 
     .                   listhptr, listh, nspin, Dscf, Eldau, DEldau,
     .                   fa, stress, H, first, last )
C *********************************************************************
C Calculates Hubbard-like U term contribution to total 
C energy, atomic forces, stress and hamiltonian matrix elements.
C Energies in Ry. Lengths in Bohr.
C Written by D. Sanchez-Portal, October 2008, 
C after subroutine nlefsm by J.Soler and P.Ordejon (June 1997).
C Based on a previous version by S. Riikonen and D. Sanchez-Portal (2005)
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
C integer iphorb(no)       : Orbital index of each orbital in its atom,
C                            where no=lasto(na)
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
C integer nspin            : Number of spin components
C real*8  Dscf(maxnd,nspin): Density matrix
C ******************* INPUT and OUTPUT *********************************
C real*8 fa(3,na)          : NL forces (added to input fa)
C real*8 stress(3,3)       : NL stress (added to input stress)
C real*8 H(maxnh,nspin)    : NL Hamiltonian (added to input H)
C **************************** OUTPUT *********************************
C real*8 Eldau               : U-Hubbard energy 1 
C real*8 DEldau              : U-hubbard energy 2 
C *********************************************************************
C
C  Modules
C
      use precision,     only : dp
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : GetNodeOrbs, LocalToGlobalOrb
      use parallelsubs,  only : GlobalToLocalOrb
      use atmfuncs,      only : rcut, epskb,  orb_f, kbpj_f, ldaupj_f,
     &              nldaupfis, switch_ldau, lofio, Ufio, Jfio, 
     &              cnfig_ldaupjfio
      use neighbour    , only : iana=>jan, r2ki=>r2ij, xki=>xij
      use neighbour    , only : mneighb
      use alloc,         only : re_alloc, de_alloc
      use siesta_options, only : dDtol
      integer, intent(in) ::
     .   maxnh, na, maxnd, nspin, nua

      integer, intent(in)  ::
     .  indxua(na), iphorb(*), isa(na),  
     .  lasto(0:na), listd(maxnd), listh(maxnh),
     .  numd(*), numh(*), listdptr(*), listhptr(*)

      real(dp), intent(in) :: scell(3,3), Dscf(maxnd,nspin),
     .                        xa(3,na)
      logical,  intent(in) :: first, last
      real(dp), intent(inout) :: fa(3,nua), stress(3,3)
      real(dp), intent(inout) :: H(maxnh,nspin)
      real(dp), intent(out)   :: Eldau, DEldau

      real(dp) ::   volcel
      external ::   timer, volcel

C Internal variables ................................................
C maxno  = maximum number of basis orbitals overlapping a KB projector

      integer, save ::  maxno = 2000
      integer, save ::  maxldau

      integer
     .  ia, ikb, ina, ind, ino,
     .  io, iio, ioa, is, ispin, ix, 
     .  j, jno, jo, jx, ka, ko, koa, ks, kua,
     .  nkb, nna, nno, no, nuo, nuotot, ja,
     .  lko, lkb, nprin_ko, nprin_kb

      integer, dimension(:), pointer :: iano, iono

      real(dp)
     .  Cijk, epsk, fik, rki, rmax, rmaxldau, rmaxo, 
     .  Sik, Sjk, volume,  oc(2), Ueff, dn, dnmax,
     .  Dij, stress_ldau(3,3)

      real(dp), dimension(:,:), pointer :: Vi, Di
      real(dp), dimension(:,:), pointer :: Ski, xno
      real(dp), dimension(:,:), pointer, save :: Hldau(:,:)
      real(dp), dimension(:,:,:), pointer :: grSki
      real(dp), dimension(:,:,:,:), pointer, save :: occu, occu_old
#ifdef MPI
      real(dp), dimension(:,:,:,:), pointer, save :: buffer1
#endif


      logical ::   within
      logical, dimension(:), pointer ::  listed, listedall
      logical, save :: frstme
      data frstme /.true./
C ......................

C Start time counter
      call timer( 'hubbard', 1 )

C Check if there is a U term
      if(frstme) then 
         do ia=1,na
             is=isa(ia) 
             if(switch_ldau(is)) ldau=.true.
             if(ldau) exit
         enddo
       endif

       Eldau=0._dp
       if(.not.ldau) return

       if(frstme) then 
C  Find maximum number of LdaU projectors of one atom = maxldau
        maxldau = 0
        do ka = 1,na
          is=isa(ka)
          maxldau = max(maxldau,nldaupfis(is))
        enddo
C Allocate local array to store the occupations of the 
C LDAU projectors
        allocate(occu(maxldau,maxldau,nua,nspin))
        call memory('A','D',size(occu),'hubbard_term')
        allocate(occu_old(maxldau,maxldau,nua,nspin))
        call memory('A','D',size(occu_old),'hubbard_term')
        occu_old=0.0_dp
#ifdef  MPI
	allocate(buffer1(maxldau,maxldau,nua,nspin))
	call memory('A','D',size(buffer1),'hubbard_term')
#endif
        allocate(Hldau(maxnh,nspin))
        call memory('A','D',size(Hldau),'hubbard_term')
        Hldau=0.0_dp 
       endif
       
       call fdf_global_get(dtol_ldaupop,'LDAU.PopTol',
     &           dtol_ldaupop_default)
       call fdf_global_get(ldau_init, 'LDAU.FirstIteration',
     &           .false.)
       call fdf_global_get(dDmax_threshold,'LDAU.ThresholdTol',
     &           dDmax_threshold_default)
         
       frstme=.false.
C End initialization


       if(first.and.(.not.ldau_init) )return
      
       occupations: if((first.and.ldau_init).or.
     &            (dDmax_current.lt.dDmax_threshold)) then 

       

C Find unit cell volume
      volume = volcel( scell ) * nua / na

C Find maximum range
      rmaxo = 0.0_dp
      rmaxldau = 0.0_dp
      do ia = 1,nua
        is = isa(ia)
        do ikb=1,nldaupfis(is)
          rmaxldau = max( rmaxldau, rcut(is,ldaupj_f,ikb) )
c         write(6,*) is,ikb, rcut(is,ldaupj_f,ikb)
        enddo
        do io = lasto(ia-1)+1,lasto(ia)
          ioa = iphorb(io)
          rmaxo = max( rmaxo, rcut(is,orb_f,ioa) )
        enddo
      enddo
      rmax = rmaxo + rmaxldau

C Initialize arrays Di and Vi only once
      no = lasto(na)
      nuotot = lasto(nua)
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)

C Allocate local memory
      nullify( Di )
      call re_alloc( Di, 1, no, 1, nspin, 
     &                  name='Di', routine='hubbard_term' )
      nullify( Vi )
      call re_alloc( Vi, 1, no, 1, nspin,
     & name='Vi', routine='hubbard_term' )
      nullify( listed )
      call re_alloc( listed, 1, no, name='listed', 
     &  routine='hubbard_term')
      nullify( listedall )
      call re_alloc( listedall, 1, no, name='listedall',
     &               routine='hubbard_term' )
  
      do jo = 1,no
        do ispin=1,nspin
           Di(jo,ispin) = 0.0d0   !AG: superfluous after initial re_alloc
           Vi(jo,ispin) = 0.0d0
        enddo 
        listed(jo) = .false.
        listedall(jo) = .false.
      enddo

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
      call re_alloc( iano, 1, maxno, name='iano', 
     & routine='hubbard_term' )
      nullify( iono )
      call re_alloc( iono, 1, maxno, name='iono', 
     & routine='hubbard_term' )
      nullify( xno )
      call re_alloc( xno, 1, 3, 1, maxno, name='xno', 
     &               routine='hubbard_term' )
      nullify( Ski )
      call re_alloc( Ski, 1, maxldau, 1, maxno, name='Ski', 
     &               routine='hubbard_term' )
      nullify( grSki )
      call re_alloc( grSki, 1, 3, 1, maxldau, 1, maxno, name='grSki',
     &               routine='hubbard_term' )
      
       iter_ldaupop=iter_ldaupop+1
       if (Node.eq.0) write(6,'(a,i4)') 
     &   'LDAU: recalculating local occupations ', iter_ldaupop

C Initialize occupations
      occu=0._dp
#ifdef  MPI
      buffer1=0._dp
#endif

C Initialize neighb subroutine
      call mneighb( scell, rmax, na, xa, 0, 0, nna )

C Loop on atoms with LdaU projectors      
      do ka = 1,na
        kua = indxua(ka)
        ks = isa(ka)
        if(.not. switch_ldau(ks)) cycle
        nkb = nldaupfis(ks)

C Find neighbour atoms
        call mneighb( scell, rmax, na, xa, ka, 0, nna )

C Find neighbour orbitals
          nno = 0
          do ina = 1,nna
            ia = iana(ina)
            is = isa(ia)
            rki = sqrt(r2ki(ina))
            do io = lasto(ia-1)+1,lasto(ia)

C Only calculate if needed locally
              if (listedall(io)) then

                ioa = iphorb(io)
          
C Find if orbital is within range
                within = .false.
                do ko = 1, nkb
                 if (rki .lt. rcut(is,orb_f,ioa)+rcut(ks,ldaupj_f,ko)) 
     .              within = .true.
                enddo

C Find overlap between neighbour orbitals and LDAU projectors
                if (within) then
C Check maxno - if too small then increase array sizes
                  if (nno.eq.maxno) then
                    maxno = maxno + 100
                    call re_alloc( iano, 1, maxno, name='iano',
     $                   copy=.true., routine='hubbard_term' )
                    call re_alloc( iono, 1, maxno, name='iono',
     $                   copy=.true., routine='hubbard_term' )
                    call re_alloc( xno, 1, 3, 1, maxno, name='xno', 
     &                   copy=.true., routine='hubbard_term' )
                    call re_alloc( Ski,1, maxldau, 1, maxno, name='Ski',
     &                   copy=.true., routine='hubbard_term' )
                    call re_alloc( grSki, 1, 3, 1, maxldau, 1, maxno,
     $               name='grSki', routine='hubbard_term', copy=.true. )
                  endif
                  nno = nno + 1
                  iono(nno) = io
                  iano(nno) = ia
                  do ix = 1,3
                    xno(ix,nno) = xki(ix,ina)
                  enddo
                  do ko = 1, nkb
                    ioa = iphorb(io)
                    call matel( 'S', ks, is, ldaupj_f, orb_f, ko, ioa,
     .                     xki(1,ina),Ski(ko,nno), grSki(1,ko,nno) )
                  enddo

                endif

              endif

            enddo

          enddo

C Loop on neighbour orbitals
          do ino = 1,nno
            io = iono(ino)
            ia = iano(ino)

            call GlobalToLocalOrb(io,Node,Nodes,iio)
            if (iio.gt.0) then
C Valid orbital

              if (ia .le. nua) then

C Scatter density matrix row of orbital io
                do j = 1,numd(iio)
                  ind = listdptr(iio)+j
                  jo = listd(ind)
                  do ispin = 1,nspin
                    Di(jo,ispin) = Di(jo,ispin) + 
     &               Dscf(ind,ispin)
                  enddo
                enddo
          
C Scatter filter of desired matrix elements
                do j = 1,numh(iio)
                  jo = listh(listhptr(iio)+j)
                  listed(jo) = .true.
                enddo

C Find matrix elements with other neighbour orbitals
                do jno = 1,nno
                  jo = iono(jno)
                  ja = iano(jno)
                  if (listed(jo)) then
C Loop on LdaU projectors
                    do ko = 1, nkb
                      do ikb= 1, nkb
                        Sik = Ski(ko,ino)
                        Sjk = Ski(ikb,jno)
                        do ispin=1,nspin
                          Dij = Di(jo,ispin) 
                          occu(ko,ikb,kua,ispin) = 
     &                    occu(ko,ikb,kua,ispin) + 
     &                    Dij * Sik * Sjk/(3.0_dp-dble(nspin))
                        enddo 
                      enddo
                    enddo
                  endif
                enddo
                
C Restore  Di and listed
                do j = 1,numh(iio)
                  ind = listhptr(iio)+j
                  jo = listh(ind)
                  listed(jo) = .false.
                enddo
                do j = 1,numd(iio)
                  jo = listd(listdptr(iio)+j)
                  do ispin=1,nspin
                     Di(jo,ispin) = 0.0d0
                  enddo
                enddo

              endif

            endif

          enddo

      enddo
     
      dnmax=0.0_dp 
      oc=0.0_dp 
#ifdef MPI
! Global reduction of occupation
      do ka=1,nua
          is=isa(ka)
          if(switch_ldau(is)) then
           nkb=nldaupfis(is)
           do ispin=1,nspin             
             call globalize_sum(occu(1:nkb,1:nkb,ka,ispin),
     &           buffer1(1:nkb,1:nkb,ka,ispin))
           enddo 
         endif
      enddo 
      occu=buffer1
#endif
      do ka=1,nua
          is=isa(ka)
	  if(switch_ldau(is)) then
	  nkb=nldaupfis(is)

            do ko=1,nkb
              do ikb=1,nkb
c                 write(6,'(2i4,2f12.5)') ko, ikb,
c    &          (occu_old(ko,ikb,ka,ispin),ispin=1,nspin)
c               write(6,'(3i4,2f12.5)')  Node, ko, ikb,
c    &          (occu(ko,ikb,ka,ispin),ispin=1,nspin)

                do ispin=1,nspin
                  dn=occu(ko,ikb,ka,ispin)-occu_old(ko,ikb,ka,ispin)
                  dnmax=max(dnmax,dabs(dn))
                enddo 
              enddo 
c             oc(1)=oc(1)+occu(ko,ko,ka,1)
c             oc(2)=oc(2)+occu(ko,ko,ka,2)
            enddo 
          endif
      enddo 
c     write(6,*) ' occup ', oc(1), oc(2)
c     write(6,'(a,f12.6)')'LDAU: maximum change in local occup.', dnmax

      conv_ldaupop=.false.
      if(dnmax.lt.dtol_ldaupop) conv_ldaupop=.true. 
      if(.not.conv_ldaupop) occu_old=occu
       

      recalc_hamilt: if(.not.conv_ldaupop.or. last) then
      if(Node.eq.0) then 
       if(.not.last) then
        write(6,'(a)') 'LDAU: recalculating Hamiltonian'
       else
        write(6,'(a)') 'LDAU: recalculating Hamiltonian and forces'
       endif
      endif
      Hldau=0.0_dp
C Initialize neighb subroutine
      call mneighb( scell, rmax, na, xa, 0, 0, nna )
      
C Loop on atoms with LdaU projectors      
      do ka = 1,na
        kua = indxua(ka)
        ks = isa(ka)
        if(.not. switch_ldau(ks)) cycle
        nkb = nldaupfis(ks)

C Find neighbour atoms
        call mneighb( scell, rmax, na, xa, ka, 0, nna )

C Find neighbour orbitals
          nno = 0
          do ina = 1,nna
            ia = iana(ina)
            is = isa(ia)
            rki = sqrt(r2ki(ina))
            do io = lasto(ia-1)+1,lasto(ia)

C Only calculate if needed locally
              if (listedall(io)) then

                ioa = iphorb(io)
          
C Find if orbital is within range
                within = .false.
                do ko = 1, nkb
                 if (rki .lt. rcut(is,orb_f,ioa)+rcut(ks,ldaupj_f,ko)) 
     .              within = .true.
                enddo

C Find overlap between neighbour orbitals and LDAU projectors
                if (within) then
C Check maxno - if too small then increase array sizes
                  if (nno.eq.maxno) then
                    maxno = maxno + 100
                    call re_alloc( iano, 1, maxno, name='iano',
     $                   copy=.true., routine='hubbard_term' )
                    call re_alloc( iono, 1, maxno, name='iono',
     $                   copy=.true., routine='hubbard_term' )
                    call re_alloc( xno, 1, 3, 1, maxno, name='xno', 
     &                   copy=.true., routine='hubbard_term' )
                    call re_alloc( Ski,1, maxldau, 1, maxno, name='Ski',
     &                   copy=.true., routine='hubbard_term' )
                    call re_alloc( grSki, 1, 3, 1, maxldau, 1, maxno,
     $               name='grSki', routine='hubbard_term', copy=.true. )
                  endif
                  nno = nno + 1
                  iono(nno) = io
                  iano(nno) = ia
                  do ix = 1,3
                    xno(ix,nno) = xki(ix,ina)
                  enddo
                  do ko = 1, nkb
                    ioa = iphorb(io)
                    call matel( 'S', ks, is, ldaupj_f, orb_f, ko, ioa,
     .                     xki(1,ina),Ski(ko,nno), grSki(1,ko,nno) )
                  enddo

                endif

              endif

            enddo

          enddo

C Loop on neighbour orbitals
          do ino = 1,nno
            io = iono(ino)
            ia = iano(ino)

            call GlobalToLocalOrb(io,Node,Nodes,iio)
            if (iio.gt.0) then
C Valid orbital

              if (ia .le. nua) then

C Scatter density matrix row of orbital io
                do j = 1,numd(iio)
                  ind = listdptr(iio)+j
                  jo = listd(ind)
                  do ispin = 1,nspin
                    Di(jo,ispin) = Di(jo,ispin) + 
     &               Dscf(ind,ispin)
                  enddo
                enddo
          
C Scatter filter of desired matrix elements
                do j = 1,numh(iio)
                  jo = listh(listhptr(iio)+j)
                  listed(jo) = .true.
                enddo

C Find matrix elements with other neighbour orbitals
                do jno = 1,nno
                  jo = iono(jno)
                  ja = iano(jno)

                  if (listed(jo)) then
C Loop on LdaU projectors
                    do ko = 1, nkb
                       lko=lofio(ks,ldaupj_f,ko)
                       nprin_ko=cnfig_ldaupjfio(ks,ko)
                       Ueff=Ufio(ks,ko)-Jfio(ks,ko)
                       Sik = Ski(ko,ino)
                       Sjk = Ski(ko,jno)
                       do ispin=1, nspin
                          Vi(jo,ispin) = Vi(jo,ispin) +
     &                    0.5_dp*Sik*Sjk*Ueff
                          Dij=Di(jo,ispin)
                          Cijk=Ueff*Dij*Sjk 
                          if(last) then 
                           do ix=1,3
                             fik=Cijk*grSki(ix,ko,ino)
                             fa(ix,ia)=fa(ix,ia) -fik
                             fa(ix,kua)=fa(ix,kua) +fik
                             do jx=1,3
                                stress(jx,ix)=stress(jx,ix) +
     &                            xno(jx,ino)*fik/volume
                             enddo 
                           enddo 
                          endif
                          do ikb = 1, nkb
                            lkb=lofio(ks,ldaupj_f,ikb)
                            nprin_kb=cnfig_ldaupjfio(ks,ikb)
                            if(lko.eq.lkb.and.
     &                         nprin_ko.eq.nprin_kb) then 
C For the time being we will use the formulation 
C of Dudarev and collaborators
                              Cijk=occu(ko,ikb,kua,ispin)*
     &                            Ueff*Ski(ikb,jno)
                              Vi(jo,ispin) = Vi(jo,ispin) -
     &                  Sik*Cijk
                              if(last) then 
                               do ix=1,3
                                 fik=-2.0_dp*Dij*Cijk*grSki(ix,ko,ino) 
                                 fa(ix,ia)=fa(ix,ia)-fik
                                 fa(ix,kua)=fa(ix,kua)+fik
                                 do jx=1,3
                                    stress(jx,ix)=stress(jx,ix) +
     &                                  xno(jx,ino)*fik/volume
                                 enddo
                               enddo 
                              endif
                            endif 
                          enddo 
                       enddo 
                    enddo
                  endif
                enddo
                
C Restore  Di and listed
                do j = 1,numh(iio)
                  ind = listhptr(iio)+j
                  jo = listh(ind)
                  do ispin = 1,nspin
                     Hldau(ind,ispin) = Hldau(ind,ispin)+
     &                  Vi(jo,ispin)
                     Vi(jo,ispin) = 0._dp
                  enddo 
                  listed(jo) = .false.
                enddo

                do j = 1,numd(iio)
                  jo = listd(listdptr(iio)+j)
                  do ispin=1,nspin
                     Di(jo,ispin) = 0._dp
                  enddo
                enddo

              endif
            endif
          enddo
      enddo

      endif recalc_hamilt

      endif occupations
       
      Eldau=0.0_dp
      DEldau=0.0_dp
      add_hamilt: if(iter_ldaupop.gt.0) then 
      if(Node.eq.0) 
     & write(6,'(a)') 'LDAU: Adding Hamiltonian contribution'

      do ia = 1,nua
        is = isa(ia)
        do io = lasto(ia-1)+1,lasto(ia)
            call GlobalToLocalOrb(io,Node,Nodes,iio)
            if (iio.gt.0) then
               do j = 1,numh(iio)
                  ind = listhptr(iio)+j
                  do ispin = 1,nspin
                     Eldau=Eldau+
     &                0.5_dp*Hldau(ind,ispin)*Dscf(ind,ispin)
                      H(ind,ispin)=H(ind,ispin)+Hldau(ind,ispin)
                  enddo
                enddo
            endif
         enddo
         if(switch_ldau(is)) then
           nkb=nldaupfis(is)
           do ispin=1,nspin
              do ko=1,nkb
                  Ueff=Ufio(is,ko)-Jfio(is,ko)
                  DEldau=DEldau+
     &   0.25_dp*(3.0_dp-dble(nspin))*Ueff*occu(ko,ko,ia,ispin)
              enddo
           enddo   
         endif
       enddo

c      write(6,'(a,f12.6)') 'Eldau(eV)=',Eldau*13.6057_dp
c      write(6,'(a,f12.6)') 'DEldau(eV)=',DEldau*13.6057_dp
       endif add_hamilt

C Deallocate local memory

      call de_alloc( grSki, name='grSki' )
      call de_alloc( Ski, name='Ski' )
      call de_alloc( xno, name='xno' )
      call de_alloc( iono, name='iono' )
      call de_alloc( iano, name='iano' )

      call de_alloc( listedall, name='listedall' )
      call de_alloc( listed, name='listed' )
      call de_alloc( Vi, name='Vi' )
      call de_alloc( Di, name='Di' )

      call timer( 'hubbard', 2 )

      end subroutine hubbard_term

      end module m_ldau
