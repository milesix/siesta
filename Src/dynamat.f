      subroutine dynamat(no, nuo, na, nua, nuotot, nspin, ispin, jx, 
     &                  np, ntpl, iphorb, iaorb, ialr, numd, listd, 
     &                  listdptr, indxua, maxnd, isa, iter, Dscf,dDscf, 
     &                  dvol, dRho, dRhoscf, Datm, Vscf, VLR,
     &                  dvnoscf, dynmat) 
C    ----------------------------------------------------------------
C    This module computes terms 1.4, 3, 4, 5, 6 of the Dynamical 
C    Matrix as in MP thesis. It is called twice to avoid storing
C    in memory some matrices on the mesh: on the first LinRes call
C    to DHSCF (with IDYN=0) and on the final call after SCF (IDYN=1)
C    LinRes Riches, Junquera, Ordejon and Pruneda 2015 
C    ----------------------------------------------------------------


C Modules 
      use precision,     only: dp, grid_p
      use meshphi,       only: endpht, lstpht, listp2
      use meshphi,       only: directphi, phi, gradphi
      use atomlist,      only: indxuo
      use listsc_module, only: listsc
      use alloc
      use atmfuncs,      only: rcut, phiatm, all_phi
      use mesh,          only: nsp, dxa, xdop, xdsp, indexp, 
     &                         idop, xdop,cmesh, nmeshg, nsm,ipa
      use atm_types,     only: nsmax=>nspecies
      use parallelsubs,  only: GlobalToLocalOrb
      use parallel,      only: Nodes, node


      implicit none

C Argument types and dimensions
      integer :: no, nuo, na, nua, nspin, ispin, jx, indxua(na), 
     &           np, iphorb(*), iaorb(*), maxnd, listdptr(*), iter,
     &           listd(maxnd), numd(nuo), isa(*), ialr, nuotot, ntpl

      real(dp), intent(in) :: Dscf(maxnd), dDscf(maxnd), dVol,
     &                        Datm(nuotot)
      real(dp)  :: dynmat(nua,3)
     
      real(grid_p), intent(in)  :: dRho(nsp,np), 
     &              dRhoscf(nsp,np), dvnoscf(nsp,np), 
     &              Vscf(ntpl), VLR(ntpl,nspin)

C     Internal variables
      integer, parameter :: minloc  = 1000,  ! Min buffer size 
     &                      maxoa = 100   ! Max # of orb/atom
      integer  :: imp, ip, iu, iphi, ia, ind,  
     &            ii, j, last, maxloc, maxloc2, triang,
     &            isp, ijl, jil, il, io, iul,ix, iop, is, iua,
     &            nc, i, ic, nphiloc, nind, jc, 
     &            lastop, lasta, ja, jua ,iii(3)
      real(dp) :: r2sp, dxsp(3), va, grva(3,nsp),
     &            gr2va(3,3), gr2vna(3,nsp),   
     &            prod1, prod2, prod3, prod4(2),prod5 


      integer,       pointer :: ilc(:), ilocal(:), iorb(:)

      integer,       pointer :: LISTED(:,:)

      real(dp),      pointer :: r2cut(:)
      real(dp),      pointer :: phia(:,:), grphi(:,:,:)
      real(dp),      pointer :: C(:,:), gC(:,:,:)
      real(dp),      pointer :: Dlocal(:), dDlocal(:)

      logical                :: VnaListed

!  Start time counter
      call timer('dynamat', 1)

!  Find atomic cutoff radii
      nullify(r2cut)
      call re_alloc( r2cut, 1, nsmax, 'r2cut', 'dynamat' )
      r2cut(:) = 0.0_dp
      do i = 1,nuotot
        ia = iaorb(i)
        is = isa(ia)
        io = iphorb(i)
        r2cut(is) = max( r2cut(is), rcut(is,io)**2 )
      enddo

!  Find value of maxloc
      maxloc2 = maxval(endpht(1:np)-endpht(0:np-1))
      maxloc = maxloc2 + minloc
      maxloc = min( maxloc, no)
      triang  = (maxloc+1)*(maxloc+2)/2

!  Nullify pointers
      nullify(ilocal, ilc, iorb)
      nullify(C, gC, LISTED, phia, grphi) 

      call re_alloc( ilocal, 1, no, 'ilocal', 'dynamat' )
      call re_alloc( ilc, 1, maxloc2, 'ilc', 'dynamat' )
      call re_alloc( iorb, 1, maxloc, 'iorb', 'dynamat' )


      call re_alloc( Dlocal, 1, triang,'Dlocal','dynamat')
      call re_alloc( dDlocal, 1, triang,'dDlocal','dynamat')
      call re_alloc( C, 1, nsp, 1, maxloc2, 'C', 'dynamat' )
      call re_alloc( gC, 1, 3, 1, nsp, 1, maxloc2, 'gC', 'dynamat' )
      call re_alloc( LISTED, 1,na,1,maxloc2, 'LISTED', 'dynamat' )

      if (DirectPhi) allocate(phia(maxoa,nsp),grphi(3,maxoa,nsp))

C     Full initializations done only once
      ilocal(1:no)  = 0
      iorb(1:maxloc)= 0
      last          = 0

      Dlocal(:) = 0.0_dp
      dDlocal(:) = 0.0_dp
      

      LISTED(:,:)=0
      VnaListed=.false.

C     loop over grid points -----------------------------------------
      do ip = 1,np

C  Find number of nonzero orbitals at this point
        nc = endpht(ip) - endpht(ip-1)
C  iob(ib)>0 means that row ib of D must not be overwritten
C  iob(ib)=0 means that row ib of D is empty
C  iob(ib)<0 means that row ib of D contains a valid row of 
C             Dscf, but which is not required at this point

        do imp = 1+endpht(ip-1), endpht(ip)
          i = lstpht(imp)
          il = ilocal(i)
          if (il.gt.0) iorb(il) = i
        enddo

c      Look for required rows of Dscf not yet stored in D-----------
        do ic = 1, nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp) !list of non zero orbs at mesh p
          if(ilocal(i) .eq. 0 ) then

C  Look for an available row in D
            do il = 1,maxloc
C  last runs circularly over rows of D
              last = last + 1
              if (last .gt. maxloc) last = 1
              if (iorb(last) .le. 0) goto 10
            enddo
            call die('dynamat: no slot available in D or dD')
   10       continue

C  Copy row i of Dscf and dDscf into row last of D and dD

            j = abs(iorb(last))
            if (j.ne.0) ilocal(j) = 0
            ilocal(i) = last
            iorb(last) = i
            il = last
            iu = indxuo(i)
!!!!!!!!!!!!!!!!!
            ! aqui va el if del parallelLocal de rhoofd
!!!!!!!!!!!!!!!!!
            call GlobalToLocalOrb( iu, Node, Nodes, iul )
            if (i .eq. iu) then
              do ii = 1, numd(iul)
                ind = listdptr(iul)+ii
                j = listd(ind)
                ijl=idx_ijl(il,ilocal(j))
                Dlocal(ijl) = Dscf(ind)
                dDlocal(ijl) = dDscf(ind)
              enddo
            else
              do ii = 1, numd(iul)
                ind = listdptr(iul)+ii
                j   = LISTSC( i, iu, listd(ind) )
                ijl=idx_ijl(il,ilocal(j))
                Dlocal(ijl) = Dscf(ind)
                dDlocal(ijl) = dDscf(ind)
              enddo
            endif
          endif
        enddo !orbs ic

C       Loop on first orbital of mesh point
        lasta = 0
        lastop = 0
        do ic = 1, nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          il = ilocal(i)
          iu = indxuo(i)
          ia = iaorb(i)
          is = isa(ia)
          iua = indxua(ia)
          iop = listp2(imp)
          ilc(ic) = il

C  Check condition if we have to compute NA for this atom
C  to avoid problems in gamma points where superecell=u.cell
C  and the neighbor atoms have the same na index
          if (ic.eq.1) then
            VnaListed =.false.
          else
            if (ia.ne.iaorb(lstpht(imp-1))) then
              VnaListed=.false.
            else
              do ix=1,ic-1
                if (LISTED(ia,ix).eq.idop(iop)) VnaListed=.true.
              enddo
            endif
          endif

          if (directphi) then
            if (ia.ne.lasta .or. iop.ne.lastop) then 
              lasta = ia
              lastop = iop
              do isp = 1, nsp
                dxsp(1:3) = xdop(1:3,iop) + xdsp(1:3,isp) - dxa(1:3,ia)     
                r2sp = dxsp(1)**2 + dxsp(2)**2 + dxsp(3)**2
                if (r2sp.lt.r2cut(is)) then
                  call all_phi(is,+1, dxsp, nphiloc, 
     &          phia(:,isp), grphi(:,:,isp))
                else
                  grphi(:,:,isp) = 0.0_dp
                  phia(:,isp) = 0.0_dp
                endif
              enddo
            endif 
            iphi = iphorb(i)
            C(:,ic) = phia(iphi,:)
            gC(1:3,:,ic) = grphi(1:3,iphi,:)
            if (.not.(VnaListed) ) then
              do isp = 1, nsp
                dxsp(1:3) = xdop(1:3,iop) + xdsp(1:3,isp) - dxa(1:3,ia)     
                call phiatm( is, 0, dxsp, va, grva(:,isp), gr2va )
                if(iua .eq. ialr) then
                  gr2vna(1:3,isp) = gr2va(1:3,jx)
                else
                  gr2vna(1:3,isp) = 0.0_dp
                endif
              enddo 
            endif

          else !directphi false

            C(1:nsp,ic) = phi(1:nsp,imp) 
            gC(1:3,1:nsp,ic) = gradphi(1:3,1:nsp,imp)
            if (.not.(VnaListed) ) then
              do isp = 1, nsp
                dxsp(1:3) = xdop(1:3,iop) + xdsp(1:3,isp) - dxa(1:3,ia)     
                call phiatm( is, 0, dxsp, va, grva(:,isp), gr2va )
                if(iua .eq. ialr) then
                  gr2vna(1:3,isp) = gr2va(1:3,jx)
                else
                  gr2vna(1:3,isp) = 0.0_dp
                endif
              enddo
            endif

          endif !directphi

          if( ispin .eq. 1) then
            if(.not.(VnaListed) ) then
              if(iter.eq.1) then
                do isp = 1, nsp  
!            ! d2Vna*Rho (dynmat term 6)
                  do ix = 1,3
                    dynmat(iua,ix) = dynmat(iua,ix) - gr2vna(ix,isp)
     &             *  dRho(isp,ip) * dVol 
                  enddo
                enddo
              elseif(iter.ne.1) then
                do isp = 1, nsp
             ! dVna*dRho (dynmat term 5)
                  do ix = 1,3
                    dynmat(iua,ix) = dynmat(iua,ix) + grva(ix,isp) 
     &                  *  dRhoscf(isp,ip) * dVol 
                  enddo
                enddo
              endif
            endif
            LISTED(ia,ic)=idop(iop)
            VnaListed=.false.

            prod1 = 2.0_dp * Datm(iu) * dvol  
            do isp = 1, nsp
              nind = (ip-1) * nsp + isp
              prod2 = prod1 * C(isp,ic) 
              do ix = 1,3
                prod3 = prod2 * gC(ix,isp,ic)  
             ! 2 Rho*Phi*dPhi*dVna (part of 4 in dynmat)
                if(iter.eq.1 .and. iua .eq. ialr) then
                  dynmat(iua,ix) = dynmat(iua,ix) -
     &                         prod3 * dvnoscf(isp,ip)
                elseif(iter.ne.1) then
                  dynmat(iua,ix) = dynmat(iua,ix) -
     &                         dvnoscf(isp,ip) * prod3
                endif
              enddo
            enddo
          endif !spin = 1
  
          if (iter==1) then
            prod4(:)=0.0_dp
            do jc=1, ic-1
              jua=indxua(iaorb(lstpht((endpht(ip-1) + jc))))
              ja = ilc(jc)
              ijl = idx_ijl(il,ilc(jc))
              jil = idx_ijl(ilc(jc),il)
              prod4(1)=2.0_dp * dVol * Dlocal(ijl) !Dlocal is symmetric
              prod4(2)=2.0_dp * dVol * Dlocal(jil)

              if((jua.eq.ialr).and.(iua.eq.ialr)) then
                do isp=1,nsp
                  nind = (ip-1) * nsp + isp
                  do ix = 1,3
                    dynmat(iua,ix) = dynmat(iua,ix) + 
     &            VLR(nind,ispin)*prod4(1)*gC(ix,isp,ic)*C(isp,jc) +
     &            VLR(nind,ispin)*prod4(2)*C(isp,ic)*gC(ix,isp,jc)

                    dynmat(iua,ix) = dynmat(iua,ix) +
     &            Vscf(nind)*prod4(1)*gC(jx,isp,ic)*gC(ix,isp,jc) +
     &            Vscf(nind)*prod4(2)*gC(ix,isp,ic)*gC(jx,isp,jc)

                     dynmat(jua,ix) = dynmat(jua,ix) -
     &            Vscf(nind)*prod4(1)*gC(jx,isp,ic)*gC(ix,isp,jc) -
     &            Vscf(nind)*prod4(2)*gC(ix,isp,ic)*gC(jx,isp,jc)


                  enddo
                enddo
              elseif(jua.eq.ialr) then
                do isp=1,nsp
                  nind = (ip-1) * nsp + isp
                  do ix = 1,3
                    dynmat(jua,ix) = dynmat(jua,ix) +
     &            VLR(nind,ispin)*prod4(2)*C(isp,ic)*gC(ix,isp,jc)

                    dynmat(jua,ix) = dynmat(jua,ix) +
     &            Vscf(nind)*prod4(2)*gC(jx,isp,jc)*gC(ix,isp,ic)

                    dynmat(iua,ix) = dynmat(iua,ix) -
     &            Vscf(nind)*prod4(2)*gC(jx,isp,jc)*gC(ix,isp,ic)

                  enddo
                enddo
              elseif(iua.eq.ialr) then
                do isp=1,nsp
                  nind = (ip-1) * nsp + isp
                  do ix = 1,3
                    dynmat(iua,ix) = dynmat(iua,ix) + 
     &            VLR(nind,ispin)*prod4(1)*gC(ix,isp,ic)*C(isp,jc)

                    dynmat(iua,ix) = dynmat(iua,ix) +
     &            Vscf(nind)* prod4(1)*gC(jx,isp,ic)*gC(ix,isp,jc)

                    dynmat(jua,ix) = dynmat(jua,ix) -
     &            Vscf(nind)* prod4(1)*gC(jx,isp,ic)*gC(ix,isp,jc)

                  enddo
                enddo
              endif ! ialr condition
            enddo !jc loop
!          mu-mu cases diagonal elements
            if (iua.eq.ialr) then
              ijl = idx_ijl(il,il)
              prod4(1)=2.0_dp * dVol * Dlocal(ijl)
              do isp=1,nsp
                nind = (ip-1) * nsp + isp
                do ix = 1,3
                  dynmat(iua,ix) = dynmat(iua,ix) + 
     &            VLR(nind,ispin)*prod4(1)*gC(ix,isp,ic)*C(isp,ic)

                  dynmat(iua,ix) = dynmat(iua,ix) +
     &            Vscf(nind)*prod4(1)*gC(jx,isp,ic)*gC(ix,isp,ic)

                  dynmat(iua,ix) = dynmat(iua,ix) -
     &            Vscf(nind)*prod4(1)*gC(jx,isp,ic)*gC(ix,isp,ic)

                enddo
              enddo
            endif !diagonal elements
          elseif (iter.ne.1) then
            prod4(:)=0.0_dp
            prod5=0.0_dp
            do jc=1,ic-1
              jua=indxua(iaorb(lstpht((endpht(ip-1) + jc))))
              ja = ilc(jc)
              ijl = idx_ijl(il,ilc(jc))
              jil = idx_ijl(ilc(jc),il)
              prod4(1)=2.0_dp * dVol * Dlocal(ijl)
              prod4(2)=2.0_dp * dVol * Dlocal(jil)
              prod5=2.0_dp * dVol * dDlocal(ijl)
              do isp=1,nsp
                nind = (ip-1) * nsp + isp
                do ix = 1,3
                   dynmat(iua,ix) = dynmat(iua,ix) +
     &                Vscf(nind)*prod5*gC(ix,isp,ic)*C(isp,jc) 
                   dynmat(jua,ix) = dynmat(jua,ix) +
     &                Vscf(nind)*prod5*C(isp,ic)*gC(ix,isp,jc)

                   dynmat(iua,ix) = dynmat(iua,ix) +
     &                 VLR(nind,ispin)*prod4(1)*gC(ix,isp,ic)*C(isp,jc)
                   dynmat(jua,ix) = dynmat(jua,ix) +
     &                VLR(nind,ispin)*prod4(2)*C(isp,ic)*gC(ix,isp,jc)
                enddo
              enddo
            enddo !jc loop
! mu-mu cases diagonal elements
            ijl = idx_ijl(il,il)
            prod4(1)=2.0_dp * dVol * Dlocal(ijl)
            prod5=2.0_dp * dVol * dDlocal(ijl)
            do isp=1,nsp
              nind = (ip-1) * nsp + isp
              do ix = 1,3
                   dynmat(iua,ix) = dynmat(iua,ix) +
     &                Vscf(nind)*prod5*gC(ix,isp,ic)*C(isp,ic) 

                   dynmat(iua,ix) = dynmat(iua,ix) +
     &                VLR(nind,ispin)*prod4(1)*gC(ix,isp,ic)*C(isp,ic)
              enddo
            enddo
          endif !iter if
        enddo !nc/imp loop 


! restore listed for the nexxt point
        LISTED(:,:)=0
        VnaListed=.false.

c      Restore iob for next point
        do imp = 1+endpht(ip-1), endpht(ip)
          i = lstpht(imp)
          il = ilocal(i)
          iorb(il) = -i
        enddo


!	print*,'ip, dynmat',ip, dynmat
!        call timer('dynamat',3)
!	if (ip ==1) stop
      enddo ! grid points loops

!	if (iter .ne.1) then
!      call timer('dynamat',3)
!       print*,'ip,dynmat',ip,dynmat
!        stop
!	endif

      call de_alloc( Dlocal, 'Dlocal','dynamat')
      call de_alloc( dDlocal, 'dDlocal','dynamat')
      call de_alloc( gC, 'gC', 'dynamat' )
      call de_alloc( C, 'C', 'dynamat' )
      call de_alloc( LISTED,'LISTED', 'dynamat' )
      call de_alloc( ilocal, 'ilocal', 'dynamat' )
      call de_alloc( ilc,'ilc', 'dynamat' )
      call de_alloc( iorb, 'iorb', 'dynamat' )
      call de_alloc( r2cut, 'r2cut', 'dynamat' )

      call timer('dynamat',2)


      return
      
      contains

! In any case will the compiler most likely inline this
! small routine. So it should not pose any problem.
      pure function idx_ijl(i,j) result(ij)
        integer, intent(in) :: i,j
        integer :: ij
        if ( i > j ) then
          ij = i * (i + 1)/2 + j + 1
        else
          ij = j * (j + 1)/2 + i + 1
        end if
      end function idx_ijl










      end subroutine dynamat
