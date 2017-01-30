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
      module m_dvmat

      implicit none

      private
      public :: dvmat
      CONTAINS


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
      use meshphi,       only: DirectPhi, endpht, lstpht, listp2, phi
      use meshphi,       only: gradphi
      use parallel,      only: Nodes, node
      use alloc,         only: re_alloc, de_alloc, alloc_default,
     &                         allocDefaults
      use parallelsubs,  only: GlobalToLocalOrb
#ifdef MPI
      use mpi_siesta
#endif
#ifdef _OPENMP
      use omp_lib
#endif

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
      integer                  :: i, ia, ic, ii, imp, 
     &                            ind, iop, ip, iphi, io, is, isp, 
     &                            ispin, iu, iua, iul, ix, 
     &                            nlocal, j, jc, last, lasta,
     &                            lastop, maxloc, maxloc2, NTH,TID,
     &                            nc, nphiloc, index(no), 
     &                            in, il, ij,  jl, triang, ijl
      integer                  :: lenx, leny, lenz, lenxy, nvmaxl

      logical                  :: ParallelLocal

      real(dp)                 :: dxsp(3),  Vij(3),
     &                            r2sp, 
     &                            V(nsp,nspin), Vpart(3,2,nsp)

      integer,         pointer :: ilc(:), ilocal(:), iorb(:)   

      real(dp),        pointer :: C(:,:), gC(:,:,:), phia(:,:) ,
     &                            grada(:,:,:) 
      real(dp),        pointer :: Vi(:,:,:),  r2cut(:)
      real(dp),        pointer :: Vss(:,:,:), t_Vss(:,:,:,:)  !serial potential
      real(dp),        pointer :: Vsp(:,:,:), t_Vsp(:,:,:,:)  !parallel potential

      logical                  :: Parallel_Run
      type(allocDefaults) oldDefaults
#ifdef _TRACE_
      integer :: MPIerror
#endif

#ifdef DEBUG
      call write_debug( '    PRE dvmat' )
#endif

#ifdef _TRACE_
      call MPI_Barrier( MPI_Comm_World, MPIerror )
      call MPItrace_event( 1000, 4 )
#endif

C  Start time counter
      call timer('dvmat',1)

C  Find atomic cutoff radii
      nullify(r2cut)
      call re_alloc( r2cut, 1, nsmax, 'r2cut', 'vmat' )
      r2cut = 0.0_dp
      do i = 1,nuotot
         ia = iaorb(i)
         is = isa(ia)
         io = iphorb(i)
         r2cut(is) = max( r2cut(is), rcut(is,io)**2 )
      end do

!   Set algorithm logical
      ParallelLocal = (Nodes > 1)
      lenx  = meshLim(2,1) - meshLim(1,1) + 1
      leny  = meshLim(2,2) - meshLim(1,2) + 1
      lenz  = meshLim(2,3) - meshLim(1,3) + 1
      lenxy = lenx*leny

C  Find value of maxloc 
      maxloc2 = maxval(endpht(1:np)-endpht(0:np-1))
      maxloc = maxloc2 + minloc
      maxloc = min( maxloc, no )
      triang  = (maxloc+1)*(maxloc+2)/2
      if ( ParallelLocal ) then
         if ( nrowsDscfL > 0 ) then
            nvmaxl = listdlptr(nrowsDscfL) + numdl(nrowsDscfL)
         else
            nvmaxl = 1
         end if
      end if

!   Allocate local memory
!$OMP parallel default(shared), &
!$OMP&shared(NTH,t_Vsp,t_Vss,spin), &
!$OMP&private(TID,last), &
!$OMP&private(ip,nc,nlocal,ic,imp,i,il,iu,iul,ii,ind,j,ijl,ispin), &
!$OMP&private(lasta,lastop,ia,is,iop,isp,dxsp,r2sp,nphiloc,iphi,jc,jl), &
!$OMP&private(Vij,Vpart,Vsp,Vss,ilocal,ilc,iorb,Vi,C,gC,phia,grada)


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

C  Nullify pointers
      nullify( C, gC, Vi, ilc, ilocal, phia, grada)

      call alloc_default( old=oldDefaults,
     .                    copy=.false., shrink=.false.,
     .                    imin=1, routine='dvmat' )

!$OMP critical
      call re_alloc( C, 1, nsp, 1, maxloc2, 'C', 'dvmat' )
      call re_alloc( gC, 1, 3, 1, nsp, 1, maxloc2, 'gC', 'dvmat' )
      call re_alloc( Vi, 1, triang,1, 3, 1, nspin, 'Vi', 'dvmat' )
      call re_alloc( ilc, 1, maxloc2,'ilc', 'dvmat' )
      call re_alloc( ilocal, 1, no, 'ilocal', 'dvmat' )
      call re_alloc( iorb, 1, no, 'ilocal', 'dvmat' )
      if (DirectPhi) allocate(phia(maxoa,nsp),grada(3,maxoa,nsp))
!$OMP end critical

!$OMP single
      if ( ParallelLocal ) then
         nullify( t_Vsp )
         call re_alloc( t_Vsp, 1, nvmaxl, 1, 3, 1 ,nspin, 1, NTH, 
     &         'Vsp',  'dvmat' )
      else
        if ( NTH > 1 ) then
          nullify( t_Vss )
          call re_alloc( t_Vss, 1, nvmax, 1, 3, 1, nspin, 2, NTH, 
     &         'Vss',  'dvmat' )
       end if
      end if
!$OMP end single ! implicit barrier

      if ( ParallelLocal ) then
         Vsp => t_Vsp(1:nvmaxl,:,:,TID)
         Vsp(1:nvmaxl,:,:) = 0._dp
      else
         if ( NTH > 1 ) then
           if ( TID == 1 ) then
             Vss => Vs
           else
             Vss => t_Vss(1:nvmax,:,:,TID)
             Vss(1:nvmax,:,:) = 0._dp
           end if
         else
          Vss => Vs
         end if
      end if

C  Initialise variables
      last                     = 0
      Vi(1:triang,1:3,1:nspin) = 0.0_dp
      ilocal(1:no)             = 0
      iorb(1:no)               = 0
      ilc(:)                   = 0

C  Loop over grid points
!$OMP do
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
        if (nlocal .gt. maxloc .and. last >0) then
           if ( ParallelLocal ) then
             do il = 1,last
                i   = iorb(il)
                iu  = indxuo(i)
                iul = NeedDscfL(iu)
                if ( i == iu ) then
                   do ii = 1, numdl(iul)
                      ind = listdlptr(iul) + ii
                      j   = listdl(ind)
                      ijl = idx_ijl(il,ilocal(j))
                      do ispin = 1, nspin
                         Vsp(ind,:,ispin) = Vsp(ind,:,ispin) + 
     &                                     Vi(ijl,:,ispin) * dVol
                      end do
                   end do
                else
                   do ii = 1, numdl(iul)
                      ind = listdlptr(iul) + ii
                      j   = LISTSC( i, iu, listdl(ind) )
                      ijl = idx_ijl(il,ilocal(j))
                      do ispin = 1, nspin
                         Vsp(ind,:,ispin) = Vsp(ind,:,ispin) + 
     &                                   Vi(ijl,:,ispin) * dVol
                      end do
                   end do
                end if
             end do
           else !Parallel local
             do il = 1,last
                i = iorb(il)
                iu = indxuo(i)
                call GlobalToLocalOrb( iu, Node, Nodes, iul )
                if (i .eq. iu) then
                  do ii = 1, numVs(iul)
                    ind = listVsptr(iul)+ii
                    j = listVs(ind)
                    ijl = idx_ijl(il,ilocal(j))                  
                    do ispin = 1,nspin
                      Vss(ind,1:3,ispin) = Vss(ind,1:3,ispin) + dVol *
     &                Vi(ijl,1:3,ispin)
                    enddo
                  enddo
                else !i .eq. iu
                  do ii = 1, numVs(iul)
                    ind = listVsptr(iul)+ii
                    j = LISTSC( i, iu, listVs(ind) )
                    ijl = idx_ijl(il,ilocal(j))                  
                    do ispin = 1,nspin
                      Vss(ind,1:3,ispin) = Vss(ind,1:3,ispin) + dVol *
     &                Vi(ijl,1:3,ispin)
                    enddo
                  enddo
                endif !i .eq. iu
             enddo !il
          endif ! Parallel local

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
            if (ia.ne.lasta .or. iop.ne.lastop) then
              lasta = ia
              lastop = iop
              do isp = 1,nsp
                do ix =1,3
                  dxsp(ix) = xdsp(ix,isp) + xdop(ix,iop) - dxa(ix,ia)
                enddo
                r2sp = sum(dxsp**2)
                if (r2sp.lt.r2cut(is)) then
                  grada(1:3,:,isp) = 0.0_dp
                  if (indxua(ia).eq.ialr) then
!$OMP critical
                     call all_phi( is,+1, dxsp, nphiloc,
     .                        phia(:,isp), grada(:,:,isp))
!$OMP end critical
                  endif
                else
                  phia(:,isp) = 0.0_dp
                  grada(1:3,:,isp) = 0.0_dp
                endif !r2sp if
              enddo !isp loop
            endif !lasta lastop if
            iphi = iphorb(i)
            C(1:nsp,ic) = phia(iphi,1:nsp)
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
                ijl = idx_ijl(il,ilc(ic))
                Vi(ijl,:,ispin)=Vij+Vi(ijl,:,ispin)
              enddo !jc loop
            enddo !ispin loop
          enddo !ic loop

        else !directphi=false (default)

          do ic = 1,nc
            imp = endpht(ip-1) + ic
            i = lstpht(imp)
            ia = iaorb(i)
            il = ilocal(i)
            ilc(ic) = il
            C(1:nsp,ic) = phi(1:nsp,imp) ! the value of orbitals is stored in memory
            if (indxua(ia).eq.ialr) then
              gC(1:3,1:nsp,ic) = (-1.0_dp)*gradphi(1:3,1:nsp,imp)
            else
              gC(1:3,1:nsp,ic) = 0.0_dp
            endif
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
                ijl = idx_ijl(il,jl)
                Vi(ijl,:,ispin)=Vij+Vi(ijl,:,ispin)
              enddo !jc loop
            enddo !ispin loop
          enddo !ic loop
        endif !directphi if
      enddo !ip loop
!$OMP end do nowait   

! Note that this is already performed in parallel!
!   Add final Vi to Vs
      if ( ParallelLocal .and. last > 0 ) then
        do il= 1,last
          i=iorb(il)      
          iu=indxuo(i)
          iul = NeedDscfL(iu)
          if ( i == iu ) then
             do ii = 1, numdl(iul)
                ind = listdlptr(iul) + ii
                j   = listdl(ind)
                ijl = idx_ijl(il,ilocal(j))
                do ispin = 1, nspin
                   Vsp(ind,1:3,ispin) = Vsp(ind,1:3,ispin) + 
     &                           Vi(ijl,1:3,ispin) * dVol
                enddo
             enddo
          else
             do ii = 1, numdl(iul)
                ind = listdlptr(iul) + ii
                j   = LISTSC( i, iu, listdl(ind) )
                ijl = idx_ijl(il,ilocal(j))
                do ispin = 1, nspin
                   Vsp(ind,1:3,ispin) = Vsp(ind,1:3,ispin) + 
     &                            Vi(ijl,1:3,ispin) * dVol
                end do
             end do
          endif
        enddo
      else if ( last > 0 ) then
        do il = 1 , last
           i  = iorb(il)
           iu = indxuo(i)
           if (i.eq.iu) then
             do ii = 1, numVs(iu)
               ind = listVsptr(iu)+ii
               j = listVs(ind)
               jl=ilocal(j)
               ijl = idx_ijl(il,jl)
               do ispin = 1,nspin
                 Vss(ind,1:3,ispin) = Vss(ind,1:3,ispin) + dVol *
     &            Vi(ijl,1:3,ispin)
               enddo
             enddo !ii loop
           else
             do ii = 1, numVs(iu)
               ind = listVsptr(iu)+ii
               j = LISTSC( i, iu, listVs(ind) )
               jl = ilocal(j)
               ijl = idx_ijl(il,jl)
               do ispin = 1,nspin
                 Vss(ind,1:3,ispin) = Vss(ind,1:3,ispin) + dVol *
     &            Vi(ijl,1:3,ispin)
               enddo
             enddo 
           endif !i .eq. iu
        enddo !il
      endif 

!$OMP barrier

      if ( ParallelLocal .and. NTH > 1 ) then
!$OMP do collapse(2)
        do ispin = 1 , nspin
          do ind = 1, nvmaxl
             do ii = 2, NTH
                t_Vsp(ind,:,ispin,1) = t_Vsp(ind,:,ispin,1) +
     &                 t_Vsp(ind,:,ispin,ii)
             end do
          end do
        end do
!$OMP end do
      else if ( NTH > 1 ) then
!$OMP do collapse(2)
        do ispin = 1 , nspin
          do ind = 1, nvmax
             do ii = 2, NTH
                Vs(ind,:,ispin) = Vs(ind,:,ispin) + 
     &                             t_Vss(ind,:,ispin,ii)
             end do
          end do
        end do
!$OMP end do
      end if

C  Deallocate local memory
      call de_alloc( gC, 'gC', 'dvmat' )
      call de_alloc( C, 'C', 'dvmat' )
      call de_alloc( ilocal, 'ilocal', 'dvmat' )
      call de_alloc( Vi, 'Vi','dvmat')
      call de_alloc( ilc, 'ilc','dvmat')
      call de_alloc( iorb, 'iorb','dvmat')
      if ( DirectPhi ) deallocate( phia, grada )

!$OMP master
      if ( ParallelLocal ) then
!      Redistribute Hamiltonian from mesh to orbital based distribution
         Vsp => t_Vsp(1:nvmaxl,1:3,1:nspin,1)
         do ix=1,3
           call matrixMtoO( nvmaxl, nvmax, numVs, listVsptr, nuo,
     &       nspin, Vsp(:,ix,:), Vs(:,ix,:) )
         enddo
         call de_alloc( t_Vsp, 'Vsp', 'dvmat' )
      else if ( NTH > 1 ) then
         call de_alloc( t_Vss, 'Vss', 'dvmat' )
      end if
!$OMP end master

!$OMP end parallel

      call de_alloc( r2cut, 'r2cut','dvmat')

#ifdef _TRACE_
      call MPI_Barrier( MPI_Comm_World, MPIerror )
      call MPItrace_event( 1000, 0 )
#endif

C  Restore old allocation defaults
      call alloc_default( restore=oldDefaults )
      call timer('dvmat',2)

#ifdef DEBUG
      call write_debug( '    POS dvmat' )
#endif

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

      end subroutine dvmat
      end module m_dvmat
 
