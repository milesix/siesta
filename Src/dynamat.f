      subroutine dynamat(no, nuo, na, nua, nuotot, nspin, ispin, jx, 
     &                  np, ntpl, iphorb, iaorb, ialr, numd, listd, 
     &                  listdptr, indxua, maxnd, isa, iter, Dscf,dDscf, 
     &                  dvol, dRho, dRhoscf, Datm, Vscf, VLR,
     &                  dvnoscf, dynmat,idyn) !añadido idyn for check
C    ----------------------------------------------------------------
C    This module computes terms 1.4, 3, 4, 5, 6 of the Dynamical 
C    Matrix as in MP thesis. It is called twice to avoid storing
C    in memory some matrices on the mesh: on the first LinRes call
C    to DHSCF (with IDYN=0) and on the final call after SCF (IDYN=1)
C    LinRes Riches, Junquera, Ordejon and Pruneda 2015 
C    ----------------------------------------------------------------

 
      use precision,     only: dp, grid_p
      use meshphi,       only: endpht, lstpht, listp2
      use atomlist,      only: indxuo
      use listsc_module, only: listsc
      use alloc
      use atmfuncs,      only: rcut, phiatm, all_phi
      use mesh,          only: nsp, dxa, xdop, xdsp, indexp, 
     &                         idop, xdop,cmesh, nmeshg, nsm,ipa
      use atm_types,     only: nsmax=>nspecies
      implicit none

      integer :: no, nuo, na, nua, nspin, ispin, jx, indxua(na), 
     &           np, iphorb(*), iaorb(*), maxnd, listdptr(*), iter,
     &           listd(maxnd), numd(nuo), isa(*), ialr, nuotot, ntpl,
     &           idyn  !añadido para chequear 

      real(dp), intent(in) :: Dscf(maxnd), dDscf(maxnd), dVol,
     &                        Datm(nuotot)
      real(dp)  :: dynmat(nua,3)
     
      real(grid_p), intent(in)  :: dRho(nsp,np), 
     &              dRhoscf(nsp,np), dvnoscf(nsp,np), 
     &              Vscf(ntpl), VLR(ntpl,nspin)
                    ! VLR is the dVscf but with density from dfscf2
                    ! not dependent on scfDens for first run

C integer no              : Number of basis orbitals
C integer np              : Number of columns in C (local)
C integer iphorb(*)       : Orbital index within each atom
C integer iaorb(*)        : Pointer to atom to which orbital belongs


C     Internal variables
      integer, parameter ::
     .   minb  = 100,  ! Min buffer size for local copy of Dscf
     .   maxoa = 100,   ! Max # of orbitals per atom
     .   maxloc = 300,
     .   maxn = 3

      integer :: imp, ip, io, iu, iphi, ia, ind, ib,jb, ibuff(no),
     &           ii, j, last, maxb, maxc, nscmax(3), iscell(3), isp,
     &           ix, iop, is, iua, nc, i, ic, nphiloc, nind, jc, 
     &           imag, iacell(3), iep, lastop, lasta, ja, jo,
     &           jua ,iii(3)
      real(dp) :: r2o, r2sp, dxsp(3), va, grva(3,nsp),
     &            gr2va(3,3), gr2vna(3,nsp), r, ra, phi(maxoa,nsp),
     &            grphi(3,maxoa,nsp),gr2phi(3,3,nsp),
     &            prod1, prod2, prod3, prod4, prod5, prod6, prod7,xr(3), 
     &            prod8, Dij, dDij, r2cut(nsmax), total, xa(3) 

      integer, pointer, save :: iob(:),ibc(:),ilc(:)

      logical      :: flag
      integer, pointer :: LISTED(:,:)

      real(dp), pointer, save  :: Dlocal(:,:), dDlocal(:,:), C(:,:),
     &                             gC(:,:,:)


      call timer('dynamat', 1)
      print *, 'DEBUG TRACK: in DYNAMAT'
      maxc = maxval(endpht(1:np)-endpht(0:np-1))
      maxb = maxc + minb
      maxb = min( maxb, no )
      nullify(ilc,ibc,iob)
      call re_alloc( iob, 0, maxb, 'iob', 'dynamat' )
      call re_alloc( ibc, 1, maxc, 'ibc', 'dynamat' )
      call re_alloc( Dlocal, 0,maxb,0,maxb,'Dlocal','dynamat')
      call re_alloc( dDlocal, 0,maxb,0,maxb,'dDlocal','dynamat')
      call re_alloc( C, 1, nsp, 1, maxc, 'C', 'dynamat' )
      call re_alloc( gC, 1, 3, 1, nsp, 1, maxc, 'gC', 'dynamat' )
      call re_alloc( ilc, 1, maxc, 'ilc', 'dynamat' )
      call re_alloc( LISTED, 1,na,1,maxc, 'LISTED', 'dynamat' )


      Dlocal(:,:) = 0.0_dp
      dDlocal(:,:) = 0.0_dp
      ibuff(:) = 0
      iob(:) = 0
      last = 0
      LISTED(:,:)=0
      flag=.false.

C     Find atomic cutoff radii
      r2cut(:) = 0.0_dp
      do i = 1,nuotot
        ia = iaorb(i)
        is = isa(ia)
        io = iphorb(i)
        r2cut(is) = max( r2cut(is), rcut(is,io)**2 )
      enddo

      nscmax(1:3) = maxn

C     loop over grid points -----------------------------------------
      do ip = 1,np

C  Find number of nonzero orbitals at this point
       nc = endpht(ip) - endpht(ip-1)
       do imp = 1+endpht(ip-1), endpht(ip)
         io = lstpht(imp)
         ib = ibuff(io)
         if (ib.gt.0) iob(ib) = io
       enddo
c      Look for required rows of Dscf not yet stored in D-----------
       do ic = 1, nc
        imp = endpht(ip-1) + ic
        i = lstpht(imp) !list of non zero orbs at mesh p
        if(ibuff(i) .eq. 0 ) then
 
         do ib = 1,maxb
C         last runs circularly over rows of D
          last = last + 1
          if (last .gt. maxb) last = 1
          if (iob(last) .le. 0) goto 10
         enddo
   10    continue
C        Copy row i of Dscf into row last of D
         j = abs(iob(last))
         if (j.ne.0) ibuff(j) = 0
         ibuff(i) = last
         iob(last) = i
         ib = last
         iu = indxuo(i)
         do ii = 1, numd(iu)
          ind = listdptr(iu)+ii
          j = listd(ind)
          if(i.ne.iu) j = listsc( i, iu, j)
          jb = ibuff(j)
          Dlocal(ib,jb) = Dscf(ind)
          Dlocal(jb,ib) = Dscf(ind)
          dDlocal(ib,jb) = dDscf(ind)
          dDlocal(jb,ib) = dDscf(ind)
         enddo  
        endif  
        ibc(ic) = ibuff(i)
       enddo !ic loop

c      Restore iob for next point
       do imp = 1+endpht(ip-1), endpht(ip)
         i = lstpht(imp)
         ib = ibuff(i)
         iob(ib) = -i
       enddo
c      loop over all non-zero orbitals at mesh point ----------------
       lasta = 0
       lastop = 0


        call ipack(-1,3,nmeshg/nsm,iii,ip) 
c	point coordinates (ip)  respect origin in unit cell 	
        do ix = 1,3
          xr(ix) = iii(1) * cmesh(ix,1) + iii(2)*cmesh(ix,2) +
     .             iii(3) * cmesh(ix,3)
        enddo

       do ic = 1, nc
        imp = endpht(ip-1) + ic
        i = lstpht(imp)
        iu = indxuo(i)
        iphi = iphorb(i)! Orbital index of each  orbital in its atom       
        ia = iaorb(i)
        is = isa(ia)
        iua = indxua(ia)
        iop = listp2(imp)
        ilc(ic) = i

        do isp = 1, nsp

         dxsp(1:3) = xdop(1:3,iop) + xdsp(1:3,isp) - dxa(1:3,ia)     
         r2sp = dxsp(1)**2 + dxsp(2)**2 + dxsp(3)**2
          
         if (r2sp.lt.r2cut(is)) then
           call all_phi(is,+1,dxsp,nphiloc, 
     &       phi(:,isp),grphi(:,:,isp))
         else
           grphi(:,:,isp) = 0.0_dp
           phi(:,isp) = 0.0_dp
         endif

         if(ia.ne.lasta .or. iop .ne. lastop ) then 
          call phiatm( is, 0, dxsp, va, grva(:,isp), gr2va )
          if(iua .eq. ialr) gr2vna(1:3,isp) = gr2va(1:3,jx)
         else
           grva(:,isp) = 0.0_dp
           gr2va = 0.0_dp
         endif
         
          C(isp,ic) = phi(iphi,isp)
          gC(1:3,isp,ic) = grphi(1:3,iphi,isp)

        enddo !isp loop
C       compute contribution of neutral atom potential to dynamical
C       matrix ------------------------------------------------------
	if (ic.eq.1) then
	  flag=.false.
	else
          if (ia.ne.iaorb(lstpht(imp-1))) then
             flag=.false.
          else
            do ix=1,ic-1
              if (LISTED(ia,ix).eq.idop(iop)) flag=.true.
            enddo
          endif
        endif

!	if (ip.eq.365) then
!	print*,'ip,ic',ip,ic
!	print*,'ip,ia,flag',ip,ia,idop(iop),flag
!	endif

        if( ispin .eq. 1) then
         if(.not.(flag) ) then
          lasta = ia
          lastop = iop
          if((iter.eq.1) .and. (iua.eq.ialr)) then
           do isp = 1, nsp  
!            ! d2Vna*Rho (dynmat term 6)
            do ix = 1,3
             dynmat(iua,ix) = dynmat(iua,ix) - gr2vna(ix,isp)
     &             *  dRho(isp,ip) * dVol !Ok-ref1
            enddo
           enddo
          elseif(iter.ne.1) then
            do isp = 1, nsp
             ! dVna*dRho (dynmat term 5)
              do ix = 1,3
            dynmat(iua,ix) = dynmat(iua,ix) + grva(ix,isp) 
     &                  *  dRhoscf(isp,ip) * dVol !OK-ref2
              enddo
            enddo
          endif
         endif
        endif !ispin
	LISTED(ia,ic)=idop(iop)
	flag=.false.
C       Calulate the terms with the diagonal density matrix that 
C       contribute to the dynamical matrix ---------------------------
C       Only run this loop for the first component of the spin to 
C       avoid the double counting ------------------------------------
        if(ispin .eq.1) then
         prod1 = 2.0_dp * Datm(iu) * dvol !ok-ref3 
         do isp = 1, nsp
           nind = (ip-1) * nsp + isp
           prod2 = prod1 * C(isp,ic) !ok-ref4
           do ix = 1,3
             prod3 = prod2 * gC(ix,isp,ic) !ni-ref5 
             ! 2 Rho*Phi*dPhi*dVna (part of 4 in dynmat)
             if(iter.eq.1 .and. iua .eq. ialr) then
             !decent
              dynmat(iua,ix) = dynmat(iua,ix) - prod3 * dvnoscf(isp,ip)
             ! ok ref dvno
             elseif(iter.ne.1) then
              !dvnoscf now only contains dVh + dVna, part of dynmat term 4 
             dynmat(iua,ix) = dynmat(iua,ix) - dvnoscf(isp,ip) * prod3
              ! ok ref dvno
             endif
           enddo
         enddo
        endif !spin = 1
       enddo !nc/imp loop 


!	if (ip.eq.365) then
!		print*,'nc',nc
!	endif

	LISTED(:,:)=0
	flag=.false.

!	print*,'ip,dynmat',ip,dynmat 
!	if (ip.eq. 365) stop
 
C      Calculate the integrals for the dynamical matrix -------------
       do ic = 1, nc 
        imp  = endpht(ip-1) + ic
        io   = lstpht(imp) 
        iu   = indxuo(io)
        iphi = iphorb(io)

        ia   = iaorb(io)
        iua  = indxua(ia)
        is   = isa(ia)
        ib = ibc(ic)
        do jc = 1,nc
          jb = ibc(jc)
          Dij = Dlocal(jb,ib) !or jb,ib
          dDij = dDlocal(jb,ib)
          prod4 = 2.0_dp * dVol * Dij !ok into 5
          prod6 = 2.0_dp * dVol !OK
          ja = iaorb(ilc(jc))
          jua = indxua(ja)
          do isp = 1, nsp
            nind = (ip-1) * nsp + isp
            prod5 = prod4 * C(isp,jc) !Ok ref7
            prod7 = prod6 * C(isp,jc) * Vscf(nind)!OK ref6
            prod3 = prod4 * gC(jx,isp,ic) !Ok ref7
            do ix = 1,3
             prod2 = prod5 * gC(ix,isp,ic) !ok
             prod8 = prod7 * gC(ix,isp,ic) !ok 
             if(iter .eq. 1 .and. iua .eq. ialr) then
              dynmat(iua,ix) = dynmat(iua,ix) + prod2
     &                         * VLR(nind,ispin) + prod3 
     &                         * Vscf(nind) * gC(ix,isp,jc)
              dynmat(jua,ix) = dynmat(jua,ix) - prod3 * gC(ix,isp,jc) 
     &                         * Vscf(nind) ! for dynmat term 3
             elseif( iter.ne. 1) then
             dynmat(iua,ix) = dynmat(iua,ix) + !dynmat term 1.3
     &                        prod8  * dDij 
           dynmat(iua,ix) = dynmat(iua,ix) +
     &                        prod2 * VLR(nind,ispin)  
                              !dynmat term 3
             endif
            enddo
          enddo ! isp
        enddo !jc loop
       enddo !ic loop
      enddo ! grid points loops


      call de_alloc( iob, 'iob', 'dynamat' )
      call de_alloc( Dlocal, 'Dlocal','dynamat')
      call de_alloc( dDlocal, 'dDlocal','dynamat')
      call de_alloc( ibc, 'ibc', 'dynamat' )
      call de_alloc( gC, 'gC', 'dynamat' )
      call de_alloc( C, 'C', 'dynamat' )
      call de_alloc( LISTED,'LISTED', 'dynamat' )
      call timer('dynamat',2)

      end subroutine dynamat
