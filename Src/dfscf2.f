! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine dfscf2( no, np, dvol, nspin, maxnd, 
     &                 numd, listdptr, listd,  
     &                 nuo, nuotot, iaorb, iphorb, isa,
     &                 iter, ialr,datm,Dscf,
     &                 Drhoatm,Drhoscf0 )

C ********************************************************************
C Based on vmat.f. Computes contribution to density: 
c dRhoscf = 2 Dscf dPhi Phi
c dRhoatm = 2 Datm dPhi Phi
C LR, Linres, summer 2015
C on first call this computes equation 2.47 in MP thesis and on 
C second call the non scf contribution to equation 2.25  
C *********************** INPUT **************************************
C integer no              : Number of basis orbitals
C integer np              : Number of columns in C (local)
C real*8  dvol            : Volume per mesh point
C integer nspin           : Number of spin components
C integer iaorb(*)        : Pointer to atom to which orbital belongs
C integer iphorb(*)       : Orbital index within each atom
C integer isa(*)          : Species index of all atoms

C  Modules
      use precision,     only: dp, grid_p,sp
      use atmfuncs,      only: rcut, all_phi
      use atm_types,     only: nsmax=>nspecies
      use atomlist,      only: indxuo,indxua
      use listsc_module, only: LISTSC
      use mesh,          only: dxa, nsp, xdop, xdsp, meshLim
      use meshdscf,      only: matrixMtoO
      use meshdscf,      only: needdscfl, listdl, numdl, nrowsdscfl,
     &                         listdlptr
      use meshphi,       only: directphi, endpht, lstpht, listp2, phi
      use parallel,      only: Nodes, node
      use alloc,         only: re_alloc, de_alloc
      use parallelsubs,  only: GlobalToLocalOrb

      implicit none

C Argument types and dimensions
      integer                  :: no, np, maxnd, nuo, nuotot, iaorb(*),
     &                            nspin, iphorb(*), isa(*), numd(nuo),
     &                            listdptr(nuo), listd(maxnd), jx,
     &                            ialr, iter
      real(dp)                 :: dvol

      real(dp), intent(in)  :: datm(no), Dscf(maxnd,nspin)
      real(grid_p), intent(inout) :: Drhoatm(3,nsp,np),
     &                           Drhoscf0(3,nsp,np,nspin) 
C Internal variables and arrays
      integer,       parameter :: minloc = 1000,  ! Min buffer size
     &                            maxoa  = 100   ! Max # of orb/atom
      integer                  :: i, ia, ic, ii, ijl, il, imp, ind, iop,
     &                            ip, iphi, io, is, isp, ispin, iu, iul,
     &                            ix, j, jc, jl, last, lasta, lastop,jb,
     &                            maxloc, maxloc2, nc, nlocal, nphiloc,
     &                            maxndl, triang, lenx, leny, lenz,
     &                            lenxy,last2,index(no),in,ind2
      logical                  :: ParallelLocal
      real(dp)                 :: r2sp, dxsp(3), Dji, Dij

      integer,         pointer :: ilc(:), ilocal(:), iorb(:)   
      real(dp),        pointer :: DscfL(:,:), 
     &                            Clocal(:,:),dClocal(:,:,:),  
     &                            phia(:,:), 
     &                            r2cut(:), grphia(:,:,:) 
      integer                  :: NTH, TID, ili, ja, ib, ibuff(no),
     &                            maxb 
      integer, pointer, save :: iob(:), ibc(:)
      real(dp), pointer, save :: D(:,:,:) 
C     Start time counter
      call timer('dfscf2',1)
C     Find atomic cutoff radii
      nullify(r2cut)
      call re_alloc( r2cut, 1, nsmax, 'r2cut', 'dfscf2' )
      r2cut = 0.0_dp
      do i = 1,nuotot
        ia = iaorb(i)
        is = isa(ia)
        io = iphorb(i)
        r2cut(is) = max( r2cut(is), rcut(is,io)**2 )
      enddo
C     Set algorithm logical
      lenx  = meshLim(2,1) - meshLim(1,1) + 1
      leny  = meshLim(2,2) - meshLim(1,2) + 1
      lenz  = meshLim(2,3) - meshLim(1,3) + 1
      lenxy = lenx*leny

C     Find value of maxloc
      maxloc2 = maxval(endpht(1:np)-endpht(0:np-1))
      maxloc  = maxloc2 + minloc
      maxloc  = min( maxloc, no )
      triang  = (maxloc+1)*(maxloc+2)/2

C     Allocate local memory
      NTH = 1
      TID = 1

      nullify ( ilocal, ilc, iorb, 
     &          dClocal, phia, grphia,Clocal,D,iob,ibc )
      call re_alloc( ilocal, 1, no, 'ilocal', 'dfscf2' )
      call re_alloc( ilc, 1, maxloc2, 'ilc', 'dfscf2' )
      call re_alloc( iorb, 1, maxloc, 'iorb', 'dfscf2' )

      call re_alloc( dClocal, 1,3, 1, nsp, 1, maxloc2, 
     &               'dClocal', 'dfscf2' )
      call re_alloc( Clocal, 1, nsp, 1, maxloc2,
     &               'Clocal', 'vmat' )

      call re_alloc( phia, 1, maxoa, 1, nsp, 
     &               'phia', 'dfscf2' )
      call re_alloc( grphia, 1, 3, 1, maxoa, 1, nsp, 
     &               'phia', 'dfscf2' )

      call re_alloc( D, 0, maxloc, 0, maxloc, 1, nspin, 'D', 'dfscf2' )
      call re_alloc( iob, 0, maxloc, 'iob', 'dfscf2' )      
      call re_alloc( ibc, 1, maxloc2, 'ibc', 'dfscf2' ) 


C     Full initializations done only once
      ilocal(1:no)             = 0
      iorb(1:maxloc)           = 0
      last                     = 0
      last2                    = 0
C NEW ---
      D(:,:,:) = 0.0_dp
      ibuff(:) = 0
      iob(:) = 0
C NEW END ---
C     Loop over grid points
      do ip = 1,np
C       Find number of nonzero orbitals at this point
        nc = endpht(ip) - endpht(ip-1)
C  iob(ib)>0 means that row ib of D must not be overwritten
C  iob(ib)=0 means that row ib of D is empty
C  iob(ib)<0 means that row ib of D contains a valid row of 
C             Dscf, but which is not required at this point

        do imp = 1+endpht(ip-1), endpht(ip )!BOTH
          i = lstpht(imp)
          ib = ibuff(i)
          if (ib.gt.0) iob(ib) = i
        enddo

        !Look for required rows of Dscf not yet stored in D (dfscf)
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          if(ibuff(i) .eq.0) then

C  Look for an available row in D
            do ib = 1,maxloc
C  last runs circularly over rows of D
              last2 = last2 + 1
              if (last2 .gt. maxloc) last2 = 1
              if (iob(last2) .le. 0) goto 10
            enddo
            call die('rhoofd: no slot available in D')
   10       continue

C  Copy row i of Dscf into row last of D

            j = abs(iob(last2))
            if (j.ne.0) ibuff(j) = 0
            ibuff(i) = last2
            iob(last2) = i
            ib = last2
            iu = indxuo(i)
            call GlobalToLocalOrb( iu, Node, Nodes, iul )
            do ii = 1, numd(iul)
              ind = listdptr(iul)+ii
              j = listd(ind)
              if (i.ne.iu) j = listsc( i, iu, j )
              jb = ibuff(j)
              D(ib,jb,1:nspin) = Dscf(ind,1:nspin)
              D(jb,ib,1:nspin) = Dscf(ind,1:nspin)
            enddo
          endif
          ibc(ic) = ibuff(i)
        enddo ! orbs ic

C  Restore iob for next point (dfscf)
        do imp = 1+endpht(ip-1), endpht(ip)
          i = lstpht(imp)
          ib = ibuff(i)
          iob(ib) = -i
        enddo

C  Look for required orbitals not yet in Vlocal
        if (nlocal .gt. last) then
          do ic = 1,nc
            imp = endpht(ip-1) + ic
            i = lstpht(imp)
            if (ilocal(i) .eq. 0) then
              last = last + 1
              ilocal(i) = last
              iorb(last) = i
            endif
          enddo
        endif

C       Loop on first orbital of mesh point (VMAT)
        lasta = 0
        lastop = 0
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp) !list of nonzero orbs at this point
          il = ilocal(i)
          iu = indxuo(i) !index eqv orbital in unit cell
          ia = iaorb(i) !Atom to which orbitals belong
          is = isa(ia) !atom specie
          iop = listp2(imp) !Maps orbital-mesh point to iop
          ilc(ic) = il
          ib = ibc(ic)
C         Generate or retrieve phi values !respecto al ultimo orbital
          if (ia.ne.lasta .or. iop.ne.lastop) then
            lasta = ia
            lastop = iop
            do isp = 1,nsp
              do ix = 1,3
                dxsp(ix) = xdsp(ix,isp) + xdop(ix,iop) - dxa(ix,ia)
              enddo
              r2sp = sum(dxsp**2)
              if (r2sp.lt.r2cut(is)) then
                call all_phi( is, +1, dxsp, nphiloc, phia(:,isp),
     &                      grphia(:,:,isp))
              else
                phia(:,isp) = 0.0_dp
                grphia(1:3,:,isp) = 0.0_dp
              endif
            enddo
          endif
          iphi = iphorb(i) !Orbital index of each orbital in its atom
          Clocal(1:nsp,ic) = phia(iphi,1:nsp)
          if (iter .eq. 1) then
            dClocal(1:3,1:nsp,ic) = grphia(1:3,iphi,1:nsp)
          elseif (iter .ne. 1) then
            dClocal(1:3,1:nsp,ic) = (-1.0_dp)*grphia(1:3,iphi,1:nsp)
          endif
       enddo !ic loop 
         
       do ic=1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          iu = indxuo(i) !eq orbital unit cell
          ia = iaorb(i) !atom to wich orbital belongs
          ib = ibc(ic)
          do isp = 1,nsp
           do jx = 1, 3
            if ( iter.eq.1 ) then
              !2*Datm_ii*<Ci|dCi> (gradiente) 
              Drhoatm(jx,isp,ip) = Drhoatm(jx,isp,ip) 
     &      + 2.0_dp * Datm(iu) * Clocal(isp,ic) * dClocal(jx,isp,ic)
            elseif(indxua(ia).eq.ialr) then
              !2*Datm_ii*<Ci|dCi>  (derivada orbital)
              Drhoatm(jx,isp,ip) = Drhoatm(jx,isp,ip)
     &         + 2.0_dp * Datm(iu)*Clocal(isp,ic) * dClocal(jx,isp,ic)
            endif
           enddo
          enddo
          do ispin=1,nspin
          do jc = 1,nc
            jb = ibc(jc)
            jl = ilc(jc)
            Dji = D(jb,ib,ispin)
            do isp=1,nsp
             do jx = 1, 3 
              if ( iter.eq.1 ) then
               !Eq 2.47 in MPs thesis (gradiente)
                Drhoscf0(jx,isp,ip,ispin) = Drhoscf0(jx,isp,ip,ispin)
     &          +2.0_dp*Dji*Clocal(isp,jc) * dClocal(jx,isp,ic)
              elseif( indxua(ia) .eq. ialr) then
                !contribution to eq 2.25 in MPs (derivada orbital)
                Drhoscf0(jx,isp,ip,ispin) = Drhoscf0(jx,isp,ip,ispin)
     &          +2.0_dp*Dji*Clocal(isp,jc) * dClocal(jx,isp,ic)
              endif
             enddo
            enddo
          enddo
          enddo !ispin
         enddo !ic loop
          
      enddo !mesh

C     Free memory
      call de_alloc( phia, 'phia', 'dfscf2' )
      call de_alloc( Clocal, 'Clocal', 'dfscf2' )
      call de_alloc( iorb, 'iorb', 'dfscf2' )
      call de_alloc( ilc, 'ilc', 'dfscf2' )
      call de_alloc( ilocal, 'ilocal', 'dfscf2' )
C LINRES
      call de_alloc( grphia, 'phia', 'dfscf2' )
      call de_alloc( dClocal, 'Clocal', 'dfscf2' )
      call de_alloc( D, 'D', 'dfscf2')


      call de_alloc( r2cut, 'r2cut', 'dfscf2' )
      call timer('dfscf2',2)
      return
      end 
