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
      subroutine dvmat( no, np, dvol, nspin, Vscf, nvmax, 
     &                 numVs, listVsptr, listVs, Vs, 
     &                 nuo, nuotot, iaorb, iphorb, isa,
     &                 ialr )

C ********************************************************************
C LR 2015: Mixture of vmat.f and dfscf.f (Ordejon, Soler and Gale), it 
C computes the matrix elements <dPhi|V|Phi> + <Phi|V|dPhi> also 
C referred to as the Pulay terms that arise due to the localized
C orbitals, necessary for the perturbed hamiltonian matrix elements
C to compute the change in the dRho for DFPT in LinRes. It is non
C SCF so it is called only once on the first scf iteration.
C *********************** INPUT **************************************
C integer no              : Number of basis orbitals
C integer nuo             : Number of orbitals in unit cell (local)
C integer nuotot          : Number of orbitals in unit cell (global)
C integer np              : Number of mesh points (total is nsp*np)
C integer nspin           : Number of spin components
C integer isa(na)         : Species index of each atom
C integer iaorb(no)       : Atom to which orbitals belong
C integer iphorb(no)      : Index of orbital within its atom
C integer nvmax           : First dimension of Vs
C integer numVs(nuo)       : Number of nonzero elemts in each row of Dscf
C integer listVsptr(nuo)   : Pointer to start of row in listd
C integer listVs(dvmax)    : List of nonzero elements of Dscf
C real*4  dVs(nsp,np,nspin): Value of SCF potential at the mesh points
C *********************** INPUT and OUTPUT ****************************
C real*8  Vs    : Perturbed Hamitonian matrix
C *********************************************************************
C    6  10        20        30        40        50        60        7072

C  Modules
      use precision,     only: dp, grid_p
      use atmfuncs,      only: rcut, all_phi
      use atm_types,     only: nsmax=>nspecies
      use atomlist,      only: indxuo, indxua
      use listsc_module, only: LISTSC
      use mesh,          only: dxa, nsp, xdop, xdsp, meshLim
      use meshdscf,      only: matrixMtoO
      use meshdscf,      only: needdscfl, listdl, numdl, nrowsdscfl,
     &                         listdlptr
      use meshphi,       only: directphi, endpht, lstpht, listp2, phi
      use parallel,      only: Nodes, node
      use alloc,         only: re_alloc, de_alloc, alloc_default,
     &                         allocDefaults
      use parallelsubs,  only: GlobalToLocalOrb

      implicit none

C  Argument types and dimensions
      integer                  :: no, nuo, nuotot, np, nspin,  
     &                            isa(*), iaorb(*), iphorb(*), nvmax, 
     &                            numVs(nuo), listVsptr(nuo),
     &                            listVs(nvmax)
      integer,      intent(in) :: ialr
      real(grid_p), intent(in) :: Vscf(nsp,np,nspin) 
      real(dp)                 :: dvol, VolCel
      real(dp),         target :: Vs(nvmax,3,nspin)

C Internal variables
      integer,       parameter :: minloc  = 1000,  ! Min buffer size 
     &                            maxoa = 100   ! Max # of orbitals per atom
      integer                  :: i, ia, ib, ibuff(no), ic, ii, imp,
     &                            ind, iop, ip, iphi, io, is, isp, 
     &                            ispin, iu, iua, iul, ix, ix1, ix2,
     &                            iy, nlocal, j, jb, jc, last, lasta,
     &                            lastop, maxloc, maxloc2, dvmaxl, NTH,
     &                            nc, nphiloc, index(no), jn, k, j0,
     &                            in, il, ij,  jl, triang, ijl
      logical                  :: ParallelLocal

      real(dp)                 :: dxsp(3), grada(3,maxoa,nsp) , Vij(3),
     &                            phia(maxoa,nsp), rvol, r2sp, 
     &                            V(nsp,nspin), Vpart(3,2,nsp)

      integer,         pointer :: ilc(:), ilocal(:), iorb(:)   

      real(dp),        pointer :: C(:,:), gC(:,:,:), Vi(:,:,:)

      real(dp),        pointer :: r2cut(:)
      logical                  :: Parallel_Run
      type(allocDefaults) oldDefaults

C  Start time counter
      call timer('dvmat',1)

C No paralelized subroutine (to merge with vmat in a future...)
      NTH=1     
      ParallelLocal = (Nodes.gt.1) !should be false      

C  Nullify pointers
      nullify( C, gC, Vi, ilc, ilocal, r2cut )

C  Find atomic cutoff radii
      call re_alloc( r2cut, 1, nsmax, 'r2cut', 'dvmat' )
      do i = 1,nuotot
        ia = iaorb(i)
        is = isa(ia)
        io = iphorb(i)
        r2cut(is) = max( r2cut(is), rcut(is,io)**2 )
      enddo

C  Get old allocation defaults and set new ones
      call alloc_default( old=oldDefaults,
     .                    copy=.false., shrink=.false.,
     .                    imin=1, routine='dvmat' )
      
C  Find value of maxloc 
      maxloc2 = maxval(endpht(1:np)-endpht(0:np-1))
      maxloc = maxloc2 + minloc
      maxloc = min( maxloc, no )
      triang  = (maxloc+1)*(maxloc+2)/2

      call re_alloc( C, 1, nsp, 1, maxloc2, 'C', 'dvmat' )
      call re_alloc( gC, 1, 3, 1, nsp, 1, maxloc2, 'gC', 'dvmat' )
      call re_alloc( Vi, 1, triang,1, 3, 1, nspin, 'Vi', 'dvmat' )
      call re_alloc( ilc, 1, maxloc2,'ilc', 'dvmat' )  
      call re_alloc( ilocal, 1, no, 'ilocal', 'dvmat' )
      call re_alloc( iorb, 1, no, 'ilocal', 'dvmat' )

C  Initialise variables
      last                     = 0
      Vi(1:triang,1:3,1:nspin) = 0.0_dp
      ilocal(1:no)             = 0
      iorb(1:no)               = 0
      ilc(:)                   = 0

C  Loop over grid points
      do ip = 1,np
      
C       Find number of nonzero orbitals at this point
        nc = endpht(ip) - endpht(ip-1)
C       new required size of Vlocal 
        nlocal = last
        do ic = 1,nc
         imp = endpht(ip-1) + ic
         i = lstpht(imp)
         if (ilocal(i) .eq. 0) nlocal = nlocal + 1
        enddo

C       If overflooded, add Vlocal to Vs and reinitialize it
        if (nlocal .gt. maxloc) then
          do il = 1,last
            i = iorb(il)
            iu = indxuo(i)
            if (ParallelLocal) then
              !!!TO INCLUDE PARALELIZATION
              print*,'No parallel version,stoppig'
              stop
            else
             call GlobalToLocalOrb( iu, Node, Nodes, iul )
              if (i .eq. iu) then
                do ii = 1, numVs(iul)
                  ind = listVsptr(iul)+ii
                  j = listVs(ind)
                  jl = ilocal(j)
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  do ispin = 1,nspin
                    Vs(ind,1:3,ispin) = Vs(ind,1:3,ispin) + dVol *
     &                Vi(ij,1:3,ispin)
                  enddo
                enddo
              else !i .eq. iu
                do ii = 1, numVs(iul)
                  ind = listVsptr(iul)+ii
                  j = LISTSC( i, iu, listVs(ind) )
                  jl = ilocal(j)
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  do ispin = 1,nspin
                    Vs(ind,1:3,ispin) = Vs(ind,1:3,ispin) + dVol *
     &                Vi(ij,1:3,ispin)
                  enddo
                enddo
              endif !i .eq. iu
            endif !ParallelLocal
          enddo !il
!         Reset local arrays
          do ii= 1, last
            ilocal(iorb(ii)) = 0
          enddo
          iorb(1:last) = 0
          ijl = (last+1)*(last+2)/2
          Vi(1:ijl,1:3,1:nspin) = 0.0_dp
          last = 0
        endif !nlocal .gt. maxloc

C  Look for required orbitals not yet in Vlocal
        if (nlocal .gt. last) then 
          do ic = 1, nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
            if(ilocal(i) .eq. 0) then
              last = last + 1
              ilocal(i) = last
              iorb(last) = i
            endif
          enddo
        endif
 
C  Copy potential to a double precision array
        V(1:nsp,1:nspin) = Vscf(1:nsp,ip,1:nspin)

C  Calculate all phi values and derivatives at all subpoints
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
          if (ia.ne.lasta .or. iop.ne.lastop) then
            lasta = ia
            lastop = iop
            do isp = 1,nsp
              do ix =1,3
              dxsp(ix) = xdsp(ix,isp) + xdop(ix,iop) - dxa(ix,ia)
              enddo
              r2sp = sum(dxsp**2)
              if (r2sp.lt.r2cut(is)) then
C Calculate orbitals value at this ip point
C Calculate gradients only for perturbed atom (unit cell atom)
                grada(1:3,:,isp) = 0.0_dp
                if (indxua(ia).eq.ialr) then
                  call all_phi( is,+1, dxsp, nphiloc,
     .                        phia(:,isp), grada(:,:,isp))
                endif
              else
                phia(:,isp) = 0.0_dp
                grada(1:3,:,isp) = 0.0_dp
              endif !r2sp if
            enddo !isp loop
          endif !lasta lastop if
          iphi = iphorb(i)
          C(1:nsp,ic) = phi(1:nsp,imp)! the value of orbitals is stored in memory
          gC(1:3,1:nsp,ic) = (-1.0_dp)*grada(1:3,iphi,1:nsp)

          do ispin= 1,nspin  !pre-multiplication gradphi(ic)*V and V*phi(ic)
            do isp=1,nsp            
              Vpart(:,1,isp)=gC(:,isp,ic)*V(isp,ispin)
              Vpart(:,2,isp)=V(isp,ispin)*C(isp,ic)
            enddo

C    Loop on the second orbital of mesh point (only for jc.le.ic)
            do jc= 1,ic !tringular matrix form
              jl=ilc(jc)

              Vij(:)=0.0_dp
              do isp= 1,nsp !sum over subpoints
                Vij(:)=Vij(:)+Vpart(:,1,isp)*C(isp,jc) +
     &                    gC(:,isp,jc)*Vpart(:,2,isp)
              enddo
              if (il.gt.jl) then !from V(i,j) format to a single column
                ij = il*(il+1)/2 + jl + 1
              else
                ij = jl*(jl+1)/2 + il + 1
              endif
              Vi(ij,:,ispin)=Vij+Vi(ij,:,ispin)
            enddo !jc loop
          enddo !ispin loop
        enddo !ic loop
      enddo !ip loop
   
      do il= 1,last
        i=iorb(il)      
        iu=indxuo(i)
        if (ParallelLocal) then
          !!!TO INCLUDE PARALELIZATION
          print*,'No parallel version,stoppig'
          stop
        else
          if (i.eq.iu) then
            do ii = 1, numVs(iu)
              ind = listVsptr(iu)+ii
              j = listVs(ind)
              jl=ilocal(j)
              if (il.gt.jl) then
                ij = il*(il+1)/2 + jl + 1
              else
                ij = jl*(jl+1)/2 + il + 1
              endif
              do ispin = 1,nspin
                Vs(ind,1:3,ispin) = Vs(ind,1:3,ispin) + dVol *
     &            Vi(ij,1:3,ispin)
              enddo
            enddo !ii loop
          else
            do ii = 1, numVs(iu)
              ind = listVsptr(iu)+ii
              j = LISTSC( i, iu, listVs(ind) )
              jl = ilocal(j)
              if (il.gt.jl) then
                ij = il*(il+1)/2 + jl + 1
              else
                ij = jl*(jl+1)/2 + il + 1
              endif
              do ispin = 1,nspin
                Vs(ind,1:3,ispin) = Vs(ind,1:3,ispin) + dVol *
     &            Vi(ij,1:3,ispin)
              enddo
            enddo !ii loop
          endif !i .eq. iu
        endif !parallel
      enddo

C  Deallocate local memory
      call de_alloc( gC, 'gC', 'dvmat' )
      call de_alloc( C, 'C', 'dvmat' )
      call de_alloc( ilocal, 'ilocal', 'dvmat' )
      call de_alloc( Vi, 'Vi','dvmat')
      call de_alloc( ilc, 'ilc','dvmat')
      call de_alloc( iorb, 'iorb','dvmat')
      call de_alloc( r2cut, 'r2cut','dvmat')

C  Restore old allocation defaults
      call alloc_default( restore=oldDefaults )
      call timer('dvmat',2)

      end
