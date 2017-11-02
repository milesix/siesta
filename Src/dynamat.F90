subroutine dynamat(no, nuo, na, nua, nuotot, nspin, ispin, jx, &
     np, ntpl, iphorb, iaorb, ialr, numd, listd, &
     listdptr, indxua, maxnd, isa, iter, Dscf,dDscf, & 
     dvol, dRho, dRhoscf, Datm, Vscf, VLR, &
     dvnoscf, dynmat) 
!    ----------------------------------------------------------------
!    This module computes terms 1.4, 3, 4, 5, 6 of the Dynamical 
!    Matrix as in MP thesis. It is called twice to avoid storing
!    in memory some matrices on the mesh: on the first LinRes call
!    to DHSCF (with IDYN=0) and on the final call after SCF (IDYN=1)
!    LinRes Riches, Junquera, Ordejon and Pruneda 2015 
!    ----------------------------------------------------------------

! Modules 
  use precision,     only: dp, grid_p
  use meshphi,       only: endpht, lstpht, listp2
  use meshphi,       only: DirectPhi, phi, gradphi
  use atomlist,      only: indxuo
  use listsc_module, only: listsc
  use alloc,         only: re_alloc, de_alloc, allocDefaults, alloc_default
  use atmfuncs,      only: rcut, phiatm, all_phi
  use mesh,          only: nsp, dxa, xdop, xdsp, indexp, meshLim 
  use mesh,          only: idop, xdop,cmesh, nmeshg, nsm,ipa
  use meshdscf,      only: nrowsDscfL, listdl, listdlptr, NeedDscfL
  use meshdscf,      only: numdl, matrixOtoM, matrixMtoO
  use atm_types,     only: nsmax=>nspecies
  use parallelsubs,  only: GlobalToLocalOrb
  use parallel,      only: Nodes, node,IOnode
#ifdef MPI
  use mpi_siesta
  use m_mpi_utils, only: globalize_sum
#endif
#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

!   Argument types and dimensions
  integer, intent(in) :: no, nuo, na, nua, nspin, ispin, jx, indxua(na) 
  integer, intent(in) :: np, iphorb(*), iaorb(*), maxnd, listdptr(*), iter
  integer, intent(in) :: listd(maxnd), numd(nuo), isa(*), ialr, nuotot, ntpl
  
  real(dp), intent(in) :: Dscf(maxnd), dDscf(maxnd), dVol
  real(dp), intent(in) :: Datm(nuotot)
  real(dp), target, intent(inout)  :: dynmat(3,nua)
  
  real(grid_p), intent(in) :: dRho(nsp,np) 
  real(grid_p), intent(in) :: dRhoscf(nsp,np), dvnoscf(nsp,np) 
  real(grid_p), intent(in) :: Vscf(ntpl), VLR(ntpl,nspin)
  
  !     Internal variables and arrays
  integer, parameter :: minloc  = 1000  ! Min buffer size 
  integer, parameter :: maxoa = 100   ! Max # of orb/atom
  integer  :: imp, ip, iu, iphi, ia, ind  
  integer  :: ii, j, last, maxloc, maxloc2, triang, maxndl
  integer  :: isp, ijl, jil, il, io, iul, ix, iop, is, iua
  integer  :: nc, i, ic, nphiloc, nind, jc 
  integer  :: lastop, lasta, ja, jua ,iii(3)
  integer  :: lenx, leny, lenz, lenxy
  integer  :: NTH, TID      

  real(dp) :: r2sp, dxsp(3), va, grva(3,nsp)
  real(dp) :: gr2va(3,3), gr2vna(3,nsp)   
  real(dp) :: prod1, prod2, prod3, prod4(2),prod5 

  integer,       pointer :: ilc(:), ilocal(:), iorb(:)
  integer,       pointer :: LISTED(:,:)

  real(dp),      pointer :: r2cut(:)
  real(dp),      pointer :: phia(:,:), grphi(:,:,:)
  real(dp),      pointer :: C(:,:), gC(:,:,:)
  real(dp),      pointer :: Dlocal(:), dDlocal(:)
  real(dp),      pointer :: DscfL(:), dDscfL(:)
  ! parallel:  Will store the sumation of t_dynmat from each node
  real(dp),      pointer :: t_DYL(:,:,:) 
  ! calculated contributions in each node
  real(dp),      pointer :: DY(:,:) !will store common (s/p) dynmat
  real(dp),      pointer :: dynmat_g(:,:)

  logical                :: VnaListed, ParallelLocal

  type(allocDefaults) :: oldDefaults
  
#ifdef _TRACE_
  integer :: MPIerror
#endif
#ifdef _TRACE_
  call MPI_Barrier( MPI_Comm_World, MPIerror )
  call MPItrace_event( 1000, 4 )
#endif
  
  !  Start time counter
  call timer('dynamat', 1)

  call alloc_default( old=oldDefaults, &
       copy=.false., shrink=.false., &
       imin=1, routine='dynamat' )

  ! Set algorithm logical
  ParallelLocal = (Nodes > 1)
  lenx  = meshLim(2,1) - meshLim(1,1) + 1
  leny  = meshLim(2,2) - meshLim(1,2) + 1
  lenz  = meshLim(2,3) - meshLim(1,3) + 1
  lenxy = lenx*leny
  
  ! Find value of maxloc
  maxloc2 = maxval(endpht(1:np)-endpht(0:np-1))
  maxloc = maxloc2 + minloc
  maxloc = min( maxloc, no)
  triang  = (maxloc+1)*(maxloc+2)/2
  
  if (ParallelLocal) then
     if (nrowsDscfL > 0) then
        maxndl = listdlptr(nrowsDscfL) + numdl(nrowsDscfL)
     else
        maxndl = 1
     end if
     nullify(DscfL)
     nullify(dDscfL)
     call re_alloc( DscfL, 1, maxndl, 'DscfL')
     call re_alloc( dDscfL, 1, maxndl, 'dDscfL')
     ! Redistribute Dscf/dDscf to DscfL/dDscfL form
     call matrixOtoM( maxnd, numd, listdptr, maxndl, nuo, &
          1, Dscf, DscfL )
     call matrixOtoM( maxnd, numd, listdptr, maxndl, nuo, &
          1, dDscf, dDscfL )
  end if

  ! Find atomic cutoff radii
  nullify(r2cut)
  call re_alloc( r2cut, 1, nsmax, 'r2cut')
  r2cut(:) = 0.0_dp
  do i = 1,nuotot
     ia = iaorb(i)
     is = isa(ia)
     io = iphorb(i)
     r2cut(is) = max( r2cut(is), rcut(is,io)**2 )
  enddo
  
!$OMP parallel default(shared), &
!$OMP&private(TID), &
!$OMP&private(ip,nc,imp,i,il,ilocal,iorb,ic,last,j,iu,iul,ii,ind,ijl,Dlocal,dDlocal),&
!$OMP&private(lasta,lastop,ia,is,iua,iop,ilc,Vnalisted,LISTED,dxsp,r2sp,grphi,phia),&
!$OMP&private(iphi,C,gC,va,grva,gr2va,gr2vna,isp,ix,DY,prod1,nind,prod2,prod3,prod4),&
!$OMP&private(jc,jua,ja,jil,prod5)
  
!$OMP single
#ifdef _OPENMP
  NTH = omp_get_num_threads( )
#else
  NTH = 1
#endif
!$OMP end single ! implicit barrier, IMPORTANT

#ifdef _OPENMP
  TID = omp_get_thread_num( ) + 1
#else
  TID = 1
#endif

  ! Nullify pointers
  nullify(ilocal, ilc, iorb)
  nullify(Dlocal, dDlocal)
  nullify(C, gC, LISTED, phia, grphi) 

!$OMP critical
  call re_alloc( C, 1, nsp, 1, maxloc2, 'C')
  call re_alloc( gC, 1, 3, 1, nsp, 1, maxloc2, 'gC')

  call re_alloc( ilocal, 1, no, 'ilocal')
  call re_alloc( ilc, 1, maxloc2, 'ilc')
  call re_alloc( iorb, 1, maxloc, 'iorb')

  call re_alloc( Dlocal, 1, triang,'Dlocal')
  call re_alloc( dDlocal, 1, triang,'dDlocal')
    
  call re_alloc( LISTED, 1,na,1,maxloc2, 'LISTED')

  if ( DirectPhi ) then
     call re_alloc( phia, 1, maxoa, 1, nsp, 'phia')
     call re_alloc( grphi, 1, 3, 1, maxoa, 1, nsp, 'grphi')
  end if
!$OMP end critical

!$OMP single
  if ( ParallelLocal ) then ! Define parallel buffer
     nullify( t_DYL )
     call re_alloc( t_DYL, 1, 3, 1, nua, 1, NTH, 't_DYL')
     nullify(dynmat_g)
     call re_alloc( dynmat_g, 1, 3, 1, nua,'dynmat_g')
     dynmat_g(:,:)=0.0_dp
  else if ( NTH > 1 ) then
     nullify( t_DYL )
     call re_alloc( t_DYL, 1, 3, 1, nua, 2, NTH, 't_DYL')
  end if
!$OMP end single ! implicit barrier

  if ( ParallelLocal ) then
     DY => t_DYL(:,:,TID)
     DY(:,:) = 0._dp
  else
     if ( NTH > 1 ) then
        if ( TID == 1 ) then
           DY => dynmat
        else
           DY => t_DYL(:,:,TID)
           DY(:,:) = 0._dp
        end if
     else
        DY => dynmat
     end if
  end if

  ! Full initializations done only once
  ilocal(1:no) = 0
  iorb(1:maxloc) = 0
  last = 0
  Dlocal(:) = 0.0_dp
  dDlocal(:) = 0.0_dp
  LISTED(:,:) = 0
  VnaListed = .false.

!     loop over grid points -----------------------------------------
!$OMP do schedule(static,1)
  do ip = 1,np

     !  Find number of nonzero orbitals at this point
     nc = endpht(ip) - endpht(ip-1)
     do imp = 1+endpht(ip-1), endpht(ip)
        i = lstpht(imp)
        il = ilocal(i)
        if (il.gt.0) iorb(il) = i
     enddo

!   Look for required rows of DscfL not yet stored in Dlocal-----------
!   Look for required rows of dDscfL not yet stored in dDlocal-----------
        do ic = 1, nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp) 
          if(ilocal(i) == 0 ) then
!           Look for an available row in D
            do il = 1,maxloc
!             last runs circularly over rows of D
              last = last + 1
              if (last > maxloc) last = 1
              if (iorb(last) <= 0) goto 10
            enddo
            call die('dynamat: no slot available in D or dD')
   10       continue
!  Copy row i of Dscf and dDscf into row last of D and dD
            j = abs(iorb(last))
            if (j /= 0) ilocal(j) = 0
            ilocal(i) = last
            iorb(last) = i
            il = last
            iu = indxuo(i)
            if ( ParallelLocal ) then
               iul = NeedDscfL(iu)
               if ( i == iu) then
                  do ii = 1, numdl(iul)
                     ind = listdlptr(iul)+ii
                     j = listdl(ind)
                     ijl=idx_ijl(il,ilocal(j))
                     Dlocal(ijl) = DscfL(ind)
                     dDlocal(ijl) = dDscfL(ind)
                  enddo
               else
                  do ii = 1, numdl(iul)
                     ind = listdlptr(iul)+ii
                     j   = LISTSC( i, iu, listdl(ind) )
                     ijl=idx_ijl(il,ilocal(j))
                     Dlocal(ijl) = DscfL(ind)
                     dDlocal(ijl) = dDscfL(ind)
                  enddo
               endif
            else
               call GlobalToLocalOrb( iu, Node, Nodes, iul )
               if (i == iu) then
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
          endif
        enddo !orbs ic

!       Loop on first orbital of mesh point
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

!  Check condition if we have to compute NA for this atom
!  to avoid problems in gamma points where superecell=u.cell
!  and the neighbor atoms have the same na index
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
!$OMP critical
                  call all_phi(is,+1, dxsp, nphiloc, & 
                       phia(:,isp), grphi(:,:,isp))
!$OMP end critical
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
!$OMP critical
                call phiatm( is, 0, dxsp, va, grva(:,isp), gr2va )
!$OMP end critical
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
!$OMP critical
                call phiatm( is, 0, dxsp, va, grva(:,isp), gr2va )
!$OMP end critical
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
                    DY(ix,iua) = DY(ix,iua) - gr2vna(ix,isp) &
                           * dRho(isp,ip) * dVol 
                  enddo
                enddo
              elseif(iter.ne.1) then
                do isp = 1, nsp
             ! dVna*dRho (dynmat term 5)
                  do ix = 1,3
                    DY(ix,iua) = DY(ix,iua) + grva(ix,isp) &
                         * dRhoscf(isp,ip) * dVol 
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
                ! TODO check these equations??? They are the same
                if(iter.eq.1 .and. iua .eq. ialr) then
                  DY(ix,iua) = DY(ix,iua) - prod3 * dvnoscf(isp,ip)
                elseif(iter.ne.1) then
                  DY(ix,iua) = DY(ix,iua) - dvnoscf(isp,ip) * prod3
                endif
              enddo
            enddo
          endif !spin = 1
  
          if (iter==1) then
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
                  prod1 = VLR(nind,ispin) * prod4(1)
                  prod2 = VLR(nind,ispin) * prod4(2)
                  prod3 = Vscf(nind) * prod4(1) * gC(jx,isp,ic)
                  prod5 = Vscf(nind) * prod4(2) * gC(jx,isp,jc)
                  do ix = 1,3
                     DY(ix,iua) = DY(ix,iua) +   &
                          prod1*gC(ix,isp,ic)*C(isp,jc) + &
                          prod2*C(isp,ic)*gC(ix,isp,jc)
                     
                     DY(ix,iua) = DY(ix,iua) +  &
                          prod3*gC(ix,isp,jc) + &
                          prod5*gC(ix,isp,ic)
                     
                     DY(ix,jua) = DY(ix,jua) - &
                          prod3*gC(ix,isp,jc) - &
                          prod5*gC(ix,isp,ic)
                  enddo
                enddo
              elseif(jua.eq.ialr) then
                do isp=1,nsp
                  nind = (ip-1) * nsp + isp
                  prod1 = VLR(nind,ispin) * prod4(2) * C(isp,ic)
                  prod2 = Vscf(nind) * prod4(2) * gC(jx,isp,jc)
                  do ix = 1,3
                    DY(ix,jua) = DY(ix,jua) + prod1 * gC(ix,isp,jc)
                    DY(ix,jua) = DY(ix,jua) + prod2 * gC(ix,isp,ic)
                    DY(ix,iua) = DY(ix,iua) - prod2 * gC(ix,isp,ic)
                  enddo
                enddo
              elseif(iua.eq.ialr) then
                do isp=1,nsp
                  nind = (ip-1) * nsp + isp
                  prod1 = VLR(nind,ispin) * prod4(1) * C(isp,jc)
                  prod2 = Vscf(nind) * prod4(1) * gC(jx,isp,ic)
                  do ix = 1,3
                    DY(ix,iua) = DY(ix,iua) + prod1 * gC(ix,isp,ic)

                    DY(ix,iua) = DY(ix,iua) + prod2 * gC(ix,isp,jc)

                    DY(ix,jua) = DY(ix,jua) - prod2 * gC(ix,isp,jc)

                  enddo
                enddo
              endif ! ialr condition
            enddo !jc loop
!          mu-mu cases diagonal elements
            if (iua.eq.ialr) then
              ijl = idx_ijl(il,il)
              prod3=2.0_dp * dVol * Dlocal(ijl)
              do isp=1,nsp
                nind = (ip-1) * nsp + isp
                prod1 = VLR(nind,ispin) * prod3 * C(isp,ic)
                prod2 = Vscf(nind) * prod3 * gC(jx,isp,ic)
                do ix = 1,3
                 DY(ix,iua) = DY(ix,iua) + prod1 * gC(ix,isp,ic)

                 DY(ix,iua) = DY(ix,iua) + prod2 * gC(ix,isp,ic)

                 DY(ix,iua) = DY(ix,iua) - prod2 * gC(ix,isp,ic)

                enddo
              enddo
            endif !diagonal elements
          elseif (iter.ne.1) then
            do jc=1,ic-1
              jua=indxua(iaorb(lstpht((endpht(ip-1) + jc))))
              ja = ilc(jc)
              ijl = idx_ijl(il,ilc(jc))
              jil = idx_ijl(ilc(jc),il)
              prod4(1) = 2.0_dp * dVol * Dlocal(ijl)
              prod4(2) = 2.0_dp * dVol * Dlocal(jil)
              prod5 = 2.0_dp * dVol * dDlocal(ijl)
              do isp=1,nsp
                nind = (ip-1) * nsp + isp
                prod1 = Vscf(nind)*prod5
                prod2 = VLR(nind,ispin)*prod4(1)
                prod3 = VLR(nind,ispin)*prod4(2)
                do ix = 1,3
                   DY(ix,iua) = DY(ix,iua) + &
                     prod1*gC(ix,isp,ic)*C(isp,jc) 
                   DY(ix,jua) = DY(ix,jua) + &
                     prod1*C(isp,ic)*gC(ix,isp,jc)

                   DY(ix,iua) = DY(ix,iua) + &
                        prod2*gC(ix,isp,ic)*C(isp,jc)
                   DY(ix,jua) = DY(ix,jua) + &
                        prod3*C(isp,ic)*gC(ix,isp,jc)
                enddo
              enddo
            enddo !jc loop
! mu-mu cases diagonal elements
            ijl = idx_ijl(il,il)
            prod3 = 2.0_dp * dVol * Dlocal(ijl)
            prod5 = 2.0_dp * dVol * dDlocal(ijl)
            do isp=1,nsp
               nind = (ip-1) * nsp + isp
               prod1 = Vscf(nind)*prod5
               prod2 = VLR(nind,ispin)*prod3
               do ix = 1,3
                  DY(ix,iua) = DY(ix,iua) + &
                       prod1*gC(ix,isp,ic)*C(isp,ic)
                  
                  DY(ix,iua) = DY(ix,iua) + &
                       prod2*gC(ix,isp,ic)*C(isp,ic)
              enddo
            enddo
          endif !iter if
        enddo !nc/imp loop 


! restore listed for the nexxt point
        LISTED(:,:)=0
        VnaListed=.false.

!      Restore iob for next point
        do imp = 1+endpht(ip-1), endpht(ip)
          i = lstpht(imp)
          il = ilocal(i)
          iorb(il) = -i
        enddo

      enddo ! grid points loops
!$OMP end do nowait

!$OMP barrier

      if ( ParallelLocal) then
        if (NTH>1) then
!$OMP do collapse(2)
          do ind = 1, nua
            do ix = 1, 3
              do ii = 2, NTH
                t_DYL(ix,ind,1) = t_DYL(ix,ind,1) + &
                     t_DYL(ix,ind,ii)
              end do
            end do
          end do
!$OMP end do
        endif
#ifdef MPI
        call globalize_sum( t_DYL(1:3,1:nua,1),   &
             dynmat_g(1:3,1:nua) )
#endif
        dynmat=dynmat+dynmat_g !add buffer dynmat to the input one
      else
         if (NTH>1) then
            if (TID .ne. 1) then
!$OMP do collapse(2)
              do ind = 1, nua
                do ix = 1, 3
                  do ii = 2, NTH
                    t_DYL(ix,ind,1) = t_DYL(ix,ind,1) + &
                     t_DYL(ix,ind,ii)
                  end do
                end do
              end do
!$OMP end do
            dynmat(:,:) = dynmat(:,:) + t_DYL(:,:,1)
            endif
         endif
       endif

!     Global reduction of dynamical matrix
!$OMP single
    if ( ParallelLocal ) then
       call de_alloc(t_DYL, 't_DYL')
    else if ( NTH > 1 ) then
       call de_alloc(t_DYL, 't_DYL')
    endif
!$OMP end single nowait

!$OMP critical
    call de_alloc( C, 'C')
    call de_alloc( gC, 'gC')
    call de_alloc( ilocal, 'ilocal')
    call de_alloc( ilc,'ilc')
    call de_alloc( iorb, 'iorb')
    call de_alloc( Dlocal, 'Dlocal')
    call de_alloc( dDlocal, 'dDlocal')
    call de_alloc( LISTED, 'LISTED')
    if ( DirectPhi ) then
       call de_alloc( phia, 'phia')
       call de_alloc( grphi, 'grphi')
    end if
!$OMP end critical

!$OMP end parallel

    call de_alloc( r2cut, 'r2cut')

    if (ParallelLocal) then
        call de_alloc( DscfL, 'DscfL')
        call de_alloc( dDscfL, 'dDscfL')
        call de_alloc( dynmat_g,'dynmat_g')
    end if

#ifdef _TRACE_
    call MPI_Barrier( MPI_Comm_World, MPIerror )
    call MPItrace_event( 1000, 0 )
#endif

    call alloc_default( restore=oldDefaults )

    call timer('dynamat',2)

  CONTAINS

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
  
