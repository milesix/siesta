! 
! This file is part of the SIESTA package.
! 
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
  subroutine dfscf2( no, np, dvol, nspin, maxnd,numd, listdptr, listd, &  
     nuo, nuotot, iaorb, iphorb, isa,iter, ialr,datm, &
     Dscf,Drhoatm,Drhoscf0 )

! ********************************************************************
! Based on rhoofd- Computes contribution to perturbed density: 
! dRhoscf = 2 Dscf dPhi Phi (non-scf part of p.density)
! dRhoatm = 2 Datm dPhi Phi (perturbed atomic density)
! LR, Linres, summer 2015
! on first call this computes equation 2.47 in MP thesis
! (based on gradient of perturbed density)  and on 
! second call the non scf contribution to equation 2.25 
! Optimized by S. Illera as rhooda+rhoofd 
! *********************** INPUT **************************************
! integer no              : Number of basis orbitals
! integer np              : Number of columns in C (local)
! real*8  dvol            : Volume per mesh point
! integer nspin           : Number of spin components
! integer iaorb(*)        : Pointer to atom to which orbital belongs
! integer iphorb(*)       : Orbital index within each atom
! integer isa(*)          : Species index of all atoms
! integer iter		  : current linres iteration (1=gradient, 2=non-scf)
! integer ialr		  : index perturbed atom
! real    Dscf		  : ground state density matrix
!*******************INPUT/OUTPUT*************************************
! real  Drhoatm		  : gradient/perturbed atomic density
! real  Droscf0           : gradient/perturbed density
!*********************************************************************

!  Modules
  use precision,     only: dp, grid_p,sp
  use atmfuncs,      only: rcut, all_phi
  use atm_types,     only: nsmax=>nspecies
  use atomlist,      only: indxuo,indxua
  use listsc_module, only: LISTSC
  use mesh,          only: dxa, nsp, xdop, xdsp, meshLim
  use meshdscf,      only: matrixMtoO, matrixOtoM
  use meshdscf,      only: needdscfl, listdl, numdl, nrowsdscfl, &
                           listdlptr
  use meshphi,       only: directphi, endpht, lstpht, listp2, phi
  use meshphi,       only: gradphi
  use parallel,      only: Nodes, node
  use alloc,         only: re_alloc, de_alloc, allocDefaults, alloc_default
  use parallelsubs,  only: GlobalToLocalOrb
#ifdef MPI
  use mpi_siesta
#endif
#ifdef _OPENMP
    use omp_lib
#endif


  implicit none

! Argument types and dimensions
  integer                     :: no, np, maxnd, nuo, nuotot, iaorb(*),&
                                 nspin, iphorb(*), isa(*), numd(nuo),&
                                 listdptr(nuo), listd(maxnd), jx,&
                                 ialr, iter
  real(dp)                    :: dvol

  real(dp),        intent(in) :: datm(no), Dscf(maxnd,nspin)
  real(grid_p), intent(out) :: Drhoatm(3,nsp,np),&
                                 Drhoscf0(3,nsp,np,nspin) 
! Internal variables and arrays
  integer,          parameter :: minloc = 1000,&  ! Min buffer size
                                 maxoa  = 100   ! Max # of orb/atom
  integer                     :: i, ia, ic, ii, ijl, il, imp, ind, iop,&
                                 ip, iphi, io, is, isp, ispin, iu, iul,&
                                 ix, j, jc, last, lasta, lastop,&
                                 maxloc, maxloc2, nc, nphiloc,&
                                 maxndl, triang, lenx, leny, lenz,&
                                 lenxy,last2,jil
  integer                     :: h_spin_dim
  logical                     :: ParallelLocal
  real(dp)                    :: r2sp, dxsp(3)

  integer,           pointer  :: ilc(:), ilocal(:), iorb(:)   
  real(dp),          pointer  :: DscfL(:,:), Clocal(:,:),dClocal(:,:,:),&  
                                 phia(:,:),r2cut(:), grphia(:,:,:)  
  integer                     :: ib, ibuff(no) 
  integer,           pointer  :: iob(:), ibc(:)
  real(dp),          pointer  :: D(:,:) 

  type(allocDefaults) :: oldDefaults
#ifdef _TRACE_
  integer :: MPIerror
#endif

#ifdef DEBUG
  call write_debug( '    PRE dfscf2' )
#endif
#ifdef _TRACE_
  call MPI_Barrier( MPI_Comm_World, MPIerror )
  call MPItrace_event( 1000, 1 )
#endif
!     Start time counter
  call timer('dfscf2',1)

  call alloc_default( old=oldDefaults, &
       copy=.false., shrink=.false., &
       imin=1, routine='dfscf2' )

!     Get spin-size
  h_spin_dim = size(Dscf, 2)

!     Set algorithm logical
  ParallelLocal = (Nodes > 1)
  if (ParallelLocal) then
     if (nrowsDscfL > 0) then
        maxndl = listdlptr(nrowsDscfL) + numdl(nrowsDscfL)
     else
        maxndl = 1
     end if
     nullify(DscfL)
     call re_alloc( DscfL, 1, maxndl, 1, h_spin_dim, 'DscfL')
!     Redistribute Dscf to DscfL form
     call matrixOtoM( maxnd, numd, listdptr, maxndl, nuo, &
          h_spin_dim, Dscf, DscfL )
  end if

!     Find atomic cutoff radii
  nullify(r2cut)
  call re_alloc( r2cut, 1, nsmax, 'r2cut')
  r2cut = 0.0_dp
  do i = 1,nuotot
    ia = iaorb(i)
    is = isa(ia)
    io = iphorb(i)
    r2cut(is) = max( r2cut(is), rcut(is,io)**2 )
  enddo

! Find size of buffers to store partial copies of Dscf and C
  maxloc2 = maxval(endpht(1:np)-endpht(0:np-1))
  maxloc  = maxloc2 + minloc
  maxloc  = min( maxloc, no )
  triang  = (maxloc+1)*(maxloc+2)/2

  lenx  = meshLim(2,1) - meshLim(1,1) + 1
  leny  = meshLim(2,2) - meshLim(1,2) + 1
  lenz  = meshLim(2,3) - meshLim(1,3) + 1
  lenxy = lenx*leny

!$OMP parallel default(shared), &
!$OMP&private(ip,nc,ic,imp,i,il,last2,j,iu,iul,ii,ind,ijl,jil), &
!$OMP&private(lasta,lastop,ia,is,iop,ib,isp,ix,dxsp,r2sp,iphi,ispin), &
!$OMP&private(ilocal,ilc,iorb,dClocal,Clocal,D,iob,ibc,phia,grphia)

! Allocate local memory
  nullify ( ilocal, ilc, iorb, dClocal, Clocal, D, iob, ibc )
  nullify ( phia, grphia )

!$OMP critical
  call re_alloc(Clocal, 1, nsp, 1, maxloc2, 'Clocal')
  call re_alloc(dClocal, 1, 3, 1, nsp, 1, maxloc2, 'dClocal')
  call re_alloc(ilocal, 1, no, 'ilocal')
  call re_alloc(ilc, 1, maxloc2, 'ilc')
  call re_alloc(iorb, 1, maxloc, 'iorb')
  call re_alloc(iob, 0, maxloc, 'iob')
  call re_alloc(ibc, 1, maxloc2, 'ibc')
  call re_alloc(D, 0, triang, 1, nspin, 'D')
  if ( DirectPhi ) then
     call re_alloc(phia, 1, maxoa, 1, nsp, 'phia')
     call re_alloc(grphia, 1, 3, 1, maxoa, 1, nsp, 'grphia')
  end if
!$OMP end critical

! Full initializations done only once
  ilocal(1:no) = 0
  iorb(1:maxloc) = 0
  last = 0
  last2 = 0
  D(:,:) = 0.0_dp
  ibuff(:) = 0
  iob(:) = 0

!$OMP do schedule(static,1)
  do ip = 1,np

! Initializations
   Drhoatm(:,:,ip) = 0.0_grid_p
   Drhoscf0(:,:,ip,:) = 0.0_grid_p

! Find number of nonzero orbitals at this point
    nc = endpht(ip) - endpht(ip-1)
!   iob(ib)>0 means that row ib of D must not be overwritten
!   iob(ib)=0 means that row ib of D is empty
!   iob(ib)<0 means that row ib of D contains a valid row of 
!         Dscf, but which is not required at this point

    do imp = 1+endpht(ip-1), endpht(ip )!BOTH
      i = lstpht(imp)
      il = ilocal(i)
      if (il.gt.0) iorb(il) = i
    enddo

 !Look for required rows of Dscf not yet stored in Dlocal
    do ic = 1,nc
      imp = endpht(ip-1) + ic
      i = lstpht(imp)
      if(ilocal(i) == 0) then
!       Look for an available row in D
        do il = 1,maxloc
!C  last runs circularly over rows of D
          last2 = last2 + 1
          if (last2 .gt. maxloc) last2 = 1
          if (iorb(last2) .le. 0) goto 10
        enddo
        call die('rhoofd: no slot available in D')
10      continue

!  Copy row i of Dscf into row last of D
        j = abs(iorb(last2))
        if (j /= 0) ilocal(j) = 0
        ilocal(i) = last2
        iorb(last2) = i
        il = last2
        iu = indxuo(i)
        if ( ParallelLocal ) then
           iul= NeedDscfL(iu)
          if (i == iu) then
            do ii=1,numdl(iul)
              ind = listdlptr(iul) + ii
              j   = listdl(ind)
              ijl = idx_ijl(il,ilocal(j))
              D(ijl,:) = DscfL(ind,:)
            enddo
          else
            do ii = 1, numdl(iul)
              ind = listdlptr(iul)+ii
              j   = LISTSC( i, iu, listdl(ind) )
              ijl = idx_ijl(il,ilocal(j))
              D(ijl,:) = DscfL(ind,:)
            enddo
          endif
        else !ParallelLocal
          call GlobalToLocalOrb( iu, Node, Nodes, iul )
          if (i == iu) then
            do ii = 1, numd(iul)
              ind = listdptr(iul)+ii
              j = listd(ind)
              ijl=idx_ijl(il,ilocal(j))
              D(ijl,:) = Dscf(ind,:)
            enddo
          else
            do ii = 1, numd(iul)
              ind = listdptr(iul)+ii
              j   = LISTSC( i, iu, listd(ind) )
              ijl=idx_ijl(il,ilocal(j))
              D(ijl,:) = Dscf(ind,:)
            enddo
          endif
        endif !ParallelLocal
      endif !ilocal      
    enddo ! orbs ic

!    Check algorithm
    if ( DirectPhi ) then
      lasta = 0
      lastop = 0
      do ic = 1,nc
        imp = endpht(ip-1) + ic
        i = lstpht(imp) 
        il = ilocal(i)
        iu = indxuo(i) 
        ia = iaorb(i) 
        is = isa(ia)
        iop = listp2(imp) 
        ilc(ic) = il
        ib = ibc(ic)
!       Generate or retrieve phi values 
        if (ia /= lasta .or. iop /= lastop) then
          lasta = ia
          lastop = iop
          do isp = 1,nsp
            do ix = 1,3
              dxsp(ix) = xdsp(ix,isp) + xdop(ix,iop) - dxa(ix,ia)
            enddo
            r2sp = sum(dxsp**2)
            if (r2sp.lt.r2cut(is)) then
!$OMP critical
               call all_phi( is, +1, dxsp, nphiloc, phia(:,isp),&
                           grphia(:,:,isp))
!$OMP end critical
            else
              phia(:,isp) = 0.0_dp
              grphia(1:3,:,isp) = 0.0_dp
            endif
          enddo
        endif
        iphi = iphorb(i) !Orbital index of each orbital in its atom
        Clocal(1:nsp,ic) = phia(iphi,1:nsp)
        if (iter == 1) then
          dClocal(1:3,1:nsp,ic) = grphia(1:3,iphi,1:nsp)
        else
          dClocal(1:3,1:nsp,ic) = - grphia(1:3,iphi,1:nsp)
        endif

! Calculate atomic density =  2* sum_mu Datm*phi_mu * grad phi_mu
        do isp = 1,nsp
          do ix = 1,3
            Drhoatm(ix,isp,ip) = Drhoatm(ix,isp,ip) &
           + 2.0_dp * Datm(iu) * Clocal(isp,ic) * dClocal(ix,isp,ic)
          enddo
        enddo

! Calculate density= 2*sum_mu_nu * phi_mu *  grad phi_nu
        do jc = 1, ic-1 ! Loop on second orbital of mesh point
          ijl = idx_ijl(il,ilc(jc))
          jil= idx_ijl(ilc(jc),il)
          do ispin = 1,nspin
            do isp = 1,nsp
              do ix=1,3
                  Drhoscf0(ix,isp,ip,ispin) = Drhoscf0(ix,isp,ip,ispin) &
               +2.0_dp*D(ijl,ispin)*Clocal(isp,ic) * dClocal(ix,isp,jc) &
               +2.0_dp*D(jil,ispin)*Clocal(isp,jc) * dClocal(ix,isp,ic)
              enddo !ix
            enddo !nsp
          enddo !ispin
        enddo

        ijl = idx_ijl(il,ilc(ic)) ! add the mu-mu case
        do ispin = 1,nspin
          do isp = 1,nsp
            do ix=1,3
                Drhoscf0(ix,isp,ip,ispin) = Drhoscf0(ix,isp,ip,ispin) &
              + 2.0_dp*D(ijl,ispin)*Clocal(isp,ic) * dClocal(ix,isp,ic)
            enddo !ix
          enddo !nsp
        enddo !ispin
      enddo !ic loop 

    else !  DirectPhi==False

      do ic = 1,nc
        imp = endpht(ip-1) + ic
        i = lstpht(imp) !list of nonzero orbs at this point
        il = ilocal(i)
        iu = indxuo(i) !index eqv orbital in unit cell
        ia = iaorb(i) !Atom to which orbitals belong
        ilc(ic) = il
        Clocal(1:nsp,ic) = phi(1:nsp,imp)
        if (iter == 1) then
          dClocal(1:3,1:nsp,ic) = gradphi(1:3,1:nsp,imp) ! for gradients
        else ! for orbital derivatives
          if (indxua(ia) == ialr) then !orb. belongs to perturbed atom
            dClocal(1:3,1:nsp,ic) = - gradphi(1:3,1:nsp,imp)
          else
            dClocal(1:3,1:nsp,ic) = 0.0_dp
          endif
        endif

! Calculate atomic density =  2* sum_mu Datm*phi_mu * grad phi_mu
        do isp=1,nsp
          do ix=1,3
            Drhoatm(ix,isp,ip) = Drhoatm(ix,isp,ip) &
           + 2.0_dp * Datm(iu) * Clocal(isp,ic) * dClocal(ix,isp,ic)
          enddo
        enddo

! Calculate density= 2*sum_mu_nu * phi_mu *  grad phi_nu
        do jc=1, ic-1 ! Loop on second orbital of mesh point
          ijl = idx_ijl(il,ilc(jc))
          jil= idx_ijl(ilc(jc),il)
          do ispin = 1,nspin
            do isp = 1,nsp
              do jx=1,3
                  Drhoscf0(jx,isp,ip,ispin) = Drhoscf0(jx,isp,ip,ispin) &
               +2.0_dp*D(ijl,ispin)*Clocal(isp,ic) * dClocal(jx,isp,jc) &
               +2.0_dp*D(jil,ispin)*Clocal(isp,jc) * dClocal(jx,isp,ic) 
              enddo !ix
            enddo !nsp
          enddo !ispin
        enddo

        ijl = idx_ijl(il,ilc(ic)) ! add the mu-mu case
        do ispin = 1,nspin
          do isp = 1,nsp
            do ix=1,3
                Drhoscf0(ix,isp,ip,ispin) = Drhoscf0(ix,isp,ip,ispin) &
              + 2.0_dp*D(ijl,ispin)*Clocal(isp,ic) * dClocal(ix,isp,ic)
            enddo !ix
          enddo !nsp
        enddo !ispin
      enddo !ic

    endif !DirectPhi end

!    Restore iorb for next point
    do imp = 1+endpht(ip-1), endpht(ip)
      i  = lstpht(imp)
      il = ilocal(i)
      iorb(il) = -i
    end do

  enddo !mesh
!$OMP end do

!$OMP critical
  call de_alloc(Clocal, 'Clocal')
  call de_alloc(dClocal, 'dClocal')
  call de_alloc(ilocal, 'ilocal')
  call de_alloc(ilc, 'ilc')
  call de_alloc(iorb, 'iorb')
  call de_alloc(iob, 'iob')
  call de_alloc(ibc, 'ibc')
  call de_alloc(D, 'D')
  if ( DirectPhi ) then
     call de_alloc(phia, 'phia')
     call de_alloc(grphia, 'grphia')
  end if
!$OMP end critical

!$OMP end parallel
  call de_alloc( r2cut, 'r2cut')

  if (ParallelLocal) then
     call de_alloc( DscfL, 'DscfL')
  end if

#ifdef _TRACE_
  call MPI_Barrier( MPI_Comm_World, MPIerror )
  call MPItrace_event( 1000, 0 )
#endif

  ! Restore old allocation defaults
  call alloc_default( restore=oldDefaults )
        
  call timer('dfscf2',2)

#ifdef DEBUG
  call write_debug( '    POS dfscf2' )
#endif

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

end subroutine  dfscf2
