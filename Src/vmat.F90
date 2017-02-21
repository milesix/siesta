! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module m_vmat

  implicit none

  private
  public :: vmat, dvmat

contains

  subroutine vmat( no, np, dvol, spin, V, nvmax, &
       numVs, listVsptr, listVs, Vs, &
       nuo, nuotot, iaorb, iphorb, isa )

! ********************************************************************
! Finds the matrix elements of the potential.
! First version written by P.Ordejon.
! Name and interface modified by J.M.Soler. May'95.
! Re-ordered so that mesh is the outer loop and the orbitals are
! handled as lower-half triangular. J.D.Gale and J.M.Soler, Feb'99
! Version of vmat that use a direct algorithm to save memory.
! Modified by J.D.Gale, November'99
! *********************** INPUT **************************************
! integer no              : Number of basis orbitals
! integer np              : Number of columns in C (local)
! real*8  dvol            : Volume per mesh point
! type(tSpin) spin        : Spin configuration
! real*4  V(np,spin%Grid) : Value of the potential at the mesh points
! integer nvmax           : First dimension of listV and Vs, and maxim.
!                           number of nonzero elements in any row of Vs
! integer numVs(nuo)      : Number of non-zero elements in a row of Vs
! integer listVsptr(nuo)  : Pointer to the start of rows in listVs
! integer listVs(nvmax)   : List of non-zero elements of Vs
! integer iaorb(*)        : Pointer to atom to which orbital belongs
! integer iphorb(*)       : Orbital index within each atom
! integer isa(*)          : Species index of all atoms
! ******************** INPUT AND OUTPUT *******************************
! real*8  Vs(nvmax,spin%H): Value of nonzero elements in each row 
!                           of Vs to which the potential matrix 
!                           elements are summed up
! *********************************************************************

!  Modules
    use precision,     only: dp, grid_p
    use atmfuncs,      only: rcut, all_phi
    use atm_types,     only: nsmax=>nspecies
    use atomlist,      only: indxuo
    use t_spin,        only: tSpin
    use listsc_module, only: LISTSC
    use mesh,          only: dxa, nsp, xdop, xdsp, meshLim
    use meshdscf,      only: matrixMtoO
    use meshdscf,      only: needdscfl, listdl, numdl, nrowsdscfl, listdlptr
    use meshphi,       only: directphi, endpht, lstpht, listp2, phi
    use parallel,      only: Nodes, Node
    use alloc,         only: re_alloc, de_alloc
    use parallelsubs,  only: GlobalToLocalOrb
#ifdef MPI
    use mpi_siesta
#endif
#ifdef _OPENMP
    use omp_lib
#endif

! Argument types and dimensions
    integer                  :: no, np, nvmax, nuo, nuotot, iaorb(*)
    type(tSpin), intent(in)  :: spin
    integer                  :: iphorb(*), isa(*), numVs(nuo)
    integer                  :: listVsptr(nuo), listVs(nvmax)
    real(grid_p), intent(in) :: V(nsp,np,spin%Grid)
    real(dp)                 :: dvol
    real(dp),         target :: Vs(nvmax,spin%H)
! Internal variables and arrays
    integer, parameter :: minloc = 1000 ! Min buffer size
    integer, parameter :: maxoa  = 100  ! Max # of orb/atom
    integer :: i, ia, ic, ii, ijl, il, imp, ind, iop
    integer :: ip, iphi, io, is, isp, ispin, iu, iul
    integer :: j, jc, jl, last, lasta, lastop
    integer :: maxloc, maxloc2, nc, nlocal, nphiloc
    integer :: nvmaxl, triang, lenx, leny, lenz,lenxy

    ! Size of Hamiltonian
    logical :: ParallelLocal
    real(dp) :: Vij, r2sp, dxsp(3), VClocal(nsp)

    integer, pointer :: ilc(:), ilocal(:), iorb(:)
    real(dp), pointer :: DscfL(:,:), t_DscfL(:,:,:)
    real(dp), pointer :: Vss(:,:), t_Vss(:,:,:), Clocal(:,:)
    real(dp), pointer :: Vlocal(:,:), phia(:,:), r2cut(:)
    integer :: NTH, TID
#ifdef _TRACE_
    integer :: MPIerror
#endif

#ifdef DEBUG
    call write_debug( '    PRE vmat' )
#endif
    
#ifdef _TRACE_
    call MPI_Barrier( MPI_Comm_World, MPIerror )
    call MPItrace_event( 1000, 4 )
#endif

!   Start time counter
    call timer('vmat',1)

!   Find atomic cutoff radii
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

!   Find value of maxloc
    maxloc2 = maxval(endpht(1:np)-endpht(0:np-1))
    maxloc  = maxloc2 + minloc
    maxloc  = min( maxloc, no )
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
!$OMP&shared(NTH,t_DscfL,t_Vss,spin), &
!$OMP&private(TID,last), &
!$OMP&private(ip,nc,nlocal,ic,imp,i,il,iu,iul,ii,ind,j,ijl,ispin), &
!$OMP&private(lasta,lastop,ia,is,iop,isp,dxsp,r2sp,nphiloc,iphi,jc,jl), &
!$OMP&private(Vij,VClocal,DscfL,Vss,ilocal,ilc,iorb,Vlocal,Clocal,phia)

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

    nullify(Clocal,phia,ilocal,ilc,iorb,Vlocal)
!$OMP critical
    ! Perhaps the critical section is not needed,
    ! however it "tells" the OS to allocate per
    ! thread, possibly waiting for each thread to
    ! place the memory in the best position.
    allocate( Clocal(nsp,maxloc2) )
    allocate( ilocal(no) , ilc(maxloc2) , iorb(maxloc) )
    allocate( Vlocal(triang,spin%Grid) )
    if ( DirectPhi ) allocate( phia(maxoa,nsp) )
!$OMP end critical

!$OMP single
    if ( ParallelLocal ) then
       nullify( t_DscfL )
       call re_alloc( t_DscfL, 1, nvmaxl, 1, spin%H, 1, NTH, &
            'DscfL',  'vmat' )
    else
       if ( NTH > 1 ) then
          nullify( t_Vss )
          call re_alloc( t_Vss, 1, nvmax, 1, spin%H, 2, NTH, &
               'Vss',  'vmat' )
       end if
    end if
!$OMP end single ! implicit barrier
    
    if ( ParallelLocal ) then
       DscfL => t_DscfL(1:nvmaxl,:,TID)
       DscfL(1:nvmaxl,:) = 0._dp
    else
       if ( NTH > 1 ) then
          if ( TID == 1 ) then
             Vss => Vs
          else
             Vss => t_Vss(1:nvmax,:,TID)
             Vss(1:nvmax,:) = 0._dp
          end if
       else
          Vss => Vs
       end if
    end if

!   Full initializations done only once
    ilocal(1:no)             = 0
    iorb(1:maxloc)           = 0
    Vlocal(1:triang,:) = 0._dp
    last = 0

!   Loop over grid points
!$OMP do
    do ip = 1,np
!      Find number of nonzero orbitals at this point
       nc = endpht(ip) - endpht(ip-1)
!      Find new required size of Vlocal
       nlocal = last
       do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          if (ilocal(i) == 0) nlocal = nlocal + 1
       end do
       
!      If overflooded, add Vlocal to Vs and reinitialize it
       if (nlocal > maxloc .and. last > 0 ) then
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
                      do ispin = 1, spin%Grid
                         DscfL(ind,ispin) = DscfL(ind,ispin) + Vlocal(ijl,ispin) * dVol
                      end do
                      if ( spin%SO ) then
                         DscfL(ind,7:8) = DscfL(ind,7:8) + &
                              Vlocal(ijl,3:4) * dVol
                      end if
                   end do
                else
                   do ii = 1, numdl(iul)
                      ind = listdlptr(iul) + ii
                      j   = LISTSC( i, iu, listdl(ind) )
                      ijl = idx_ijl(il,ilocal(j))
                      do ispin = 1, spin%Grid
                         DscfL(ind,ispin) = DscfL(ind,ispin) + Vlocal(ijl,ispin) * dVol
                      end do
                      if ( spin%SO ) then
                         DscfL(ind,7:8) = DscfL(ind,7:8) + &
                              Vlocal(ijl,3:4) * dVol
                      end if
                   end do
                end if
             end do
          else
             do il = 1,last
                i  = iorb(il)
                iu = indxuo(i)
                call GlobalToLocalOrb( iu, Node, Nodes, iul )
                if ( i == iu ) then
                   do ii = 1, numVs(iul)
                      ind = listVsptr(iul) + ii
                      j   = listVs(ind)
                      ijl = idx_ijl(il,ilocal(j))
                      do ispin = 1, spin%Grid
                         Vss(ind,ispin) = Vss(ind,ispin) + Vlocal(ijl,ispin) * dVol
                      end do
                      if ( spin%SO ) then
                         Vss(ind,7:8) = Vss(ind,7:8) + Vlocal(ijl,3:4) * dVol
                      end if
                   end do
                else
                   do ii = 1, numVs(iul)
                      ind = listVsptr(iul) + ii
                      j   = LISTSC( i, iu, listVs(ind) )
                      ijl = idx_ijl(il,ilocal(j))
                      do ispin = 1, spin%Grid
                         Vss(ind,ispin) = Vss(ind,ispin) + Vlocal(ijl,ispin) * dVol
                      end do
                      if ( spin%SO ) then
                         Vss(ind,7:8) = Vss(ind,7:8) + Vlocal(ijl,3:4) * dVol
                      end if
                   end do
                end if
             end do
          end if
!         Reset local arrays
          do i = 1, last
             ilocal(iorb(i)) = 0
          end do
          iorb(1:last) = 0
          ijl = (last+1)*(last+2)/2
          do ispin = 1 , spin%Grid
             do i = 1 , ijl
                Vlocal(i,ispin) = 0._dp
             end do
          end do
          last = 0
       end if

!      Look for required orbitals not yet in Vlocal
       if ( nlocal > last ) then
          do ic = 1,nc
             imp = endpht(ip-1) + ic
             i   = lstpht(imp)
             if ( ilocal(i) == 0 ) then
                last = last + 1
                ilocal(i) = last
                iorb(last) = i
             end if
          end do
       end if

!      Check algorithm
       if ( DirectPhi ) then
          lasta = 0
          lastop = 0
          do ic = 1 , nc
             imp = endpht(ip-1) + ic
             i   = lstpht(imp)
             il  = ilocal(i)
             ia  = iaorb(i)
             iop = listp2(imp)
             ilc(ic) = il
             
!            Generate or retrieve phi values
             if ( ia /= lasta .or. iop /= lastop ) then
                lasta  = ia
                lastop = iop
                is = isa(ia)
                do isp = 1,nsp
                   dxsp(:) = xdsp(:,isp) + xdop(:,iop) - dxa(:,ia)
                   r2sp = sum(dxsp**2)
                   if ( r2sp < r2cut(is) ) then
!$OMP critical
                      call all_phi( is, +1, dxsp, nphiloc, phia(:,isp) )
!$OMP end critical
                   else
                      phia(:,isp) = 0.0_dp
                   end if
                end do
             end if
             iphi = iphorb(i)
             
             Clocal(:,ic) = phia(iphi,:)
          
             do ispin = 1 , spin%Grid

!               Create VClocal for the first orbital of mesh point
                Vij = 0._dp
                do isp = 1,nsp
                   VClocal(isp) = V(isp,ip,ispin) * Clocal(isp,ic)
!                  This is the jc == ic value
                   Vij = Vij + VClocal(isp) * Clocal(isp,ic)
                end do
                
!               ic == jc, hence we simply do
                ijl = idx_ijl(il,ilc(ic))
                Vlocal(ijl,ispin) = Vlocal(ijl,ispin) + Vij

!               Loop on second orbital of mesh point
                do jc = 1 , ic - 1
                   jl = ilc(jc)

!                  Calculate ic * jc
                   Vij = 0._dp
                   do isp = 1,nsp
                      Vij = Vij + VClocal(isp) * Clocal(isp,jc)
                   end do

                   ijl = idx_ijl(il,jl)
                   if ( il == jl ) then
                      Vlocal(ijl,ispin) = Vlocal(ijl,ispin) + Vij * 2._dp
                   else
                      Vlocal(ijl,ispin) = Vlocal(ijl,ispin) + Vij
                   end if
                   
                end do
             end do

          end do
          
       else

          do ic = 1 , nc
             imp = endpht(ip-1) + ic
             il  = ilocal(lstpht(imp))
             ilc(ic) = il

             Clocal(:,ic) = phi(:,imp)

             do ispin = 1 , spin%Grid

!               Create VClocal for the first orbital of mesh point
                Vij = 0._dp
                do isp = 1 , nsp
                   VClocal(isp) = V(isp,ip,ispin) * Clocal(isp,ic)
!                  This is the jc == ic value
                   Vij = Vij + VClocal(isp) * Clocal(isp,ic)
                end do
                
!               ic == jc, hence we simply do
                ijl = idx_ijl(il,ilc(ic))
                Vlocal(ijl,ispin) = Vlocal(ijl,ispin) + Vij

!               Loop on second orbital of mesh point
                do jc = 1 , ic - 1
                   jl = ilc(jc)

!                  Calculate ic * jc
                   Vij = 0._dp
                   do isp = 1 , nsp
                      Vij = Vij + VClocal(isp) * Clocal(isp,jc)
                   end do

                   ijl = idx_ijl(il,jl)
                   if ( il == jl ) then
                      Vlocal(ijl,ispin) = Vlocal(ijl,ispin) + Vij * 2._dp
                   else
                      Vlocal(ijl,ispin) = Vlocal(ijl,ispin) + Vij
                   end if
                   
                end do
             end do

          end do

       end if
       
    end do
!$OMP end do nowait

! Note that this is already performed in parallel!
!   Add final Vlocal to Vs
    if ( ParallelLocal .and. last > 0 ) then
       do il = 1 , last
          i   = iorb(il)
          iu  = indxuo(i)
          iul = NeedDscfL(iu)
          if ( i == iu ) then
             do ii = 1, numdl(iul)
                ind = listdlptr(iul) + ii
                j   = listdl(ind)
                ijl = idx_ijl(il,ilocal(j))
                do ispin = 1, spin%Grid
                   DscfL(ind,ispin) = DscfL(ind,ispin) + Vlocal(ijl,ispin) * dVol
                end do
                if ( spin%SO ) then
                   DscfL(ind,7:8) = DscfL(ind,7:8) + Vlocal(ijl,3:4) * dVol
                end if
             end do
          else
             do ii = 1, numdl(iul)
                ind = listdlptr(iul) + ii
                j   = LISTSC( i, iu, listdl(ind) )
                ijl = idx_ijl(il,ilocal(j))
                do ispin = 1, spin%Grid
                   DscfL(ind,ispin) = DscfL(ind,ispin) + Vlocal(ijl,ispin) * dVol
                end do
                if ( spin%SO ) then
                   DscfL(ind,7:8) = DscfL(ind,7:8) + Vlocal(ijl,3:4) * dVol
                end if
             end do
          end if
       end do
    else if ( last > 0 ) then
       do il = 1 , last
          i  = iorb(il)
          iu = indxuo(i)
          if ( i == iu ) then
             do ii = 1, numVs(iu)
                ind = listVsptr(iu)+ii
                j   = listVs(ind)
                ijl = idx_ijl(il,ilocal(j))
                do ispin = 1, spin%Grid
                   Vss(ind,ispin) = Vss(ind,ispin) + Vlocal(ijl,ispin) * dVol
                end do
                if ( spin%SO ) then
                   Vss(ind,7:8) = Vss(ind,7:8) + Vlocal(ijl,3:4) * dVol
                end if
             end do
          else
             do ii = 1, numVs(iu)
                ind = listVsptr(iu)+ii
                j   = LISTSC( i, iu, listVs(ind) )
                ijl = idx_ijl(il,ilocal(j))
                do ispin = 1, spin%Grid
                   Vss(ind,ispin) = Vss(ind,ispin) + Vlocal(ijl,ispin) * dVol
                end do
                if ( spin%SO ) then
                   Vss(ind,7:8) = Vss(ind,7:8) + Vlocal(ijl,3:4) * dVol
                end if
             end do
          end if
       end do
    end if

!$OMP barrier

    if ( ParallelLocal .and. NTH > 1 ) then
!$OMP do collapse(2)
       do ispin = 1 , spin%H
          do ind = 1, nvmaxl
             do ii = 2, NTH
                t_DscfL(ind,ispin,1) = t_DscfL(ind,ispin,1) + &
                     t_DscfL(ind,ispin,ii)
             end do
          end do
       end do
!$OMP end do
    else if ( NTH > 1 ) then
!$OMP do collapse(2)
       do ispin = 1 , spin%H
          do ind = 1, nvmax
             do ii = 2, NTH
                Vs(ind,ispin) = Vs(ind,ispin) + t_Vss(ind,ispin,ii)
             end do
          end do
       end do
!$OMP end do
    end if

!   Free memory
    deallocate(Clocal,ilocal,ilc,iorb,Vlocal)
    if ( DirectPhi ) deallocate( phia )

!$OMP master
    if ( ParallelLocal ) then
!      Redistribute Hamiltonian from mesh to orbital based distribution
       DscfL => t_DscfL(1:nvmaxl,1:spin%H,1)
       call matrixMtoO( nvmaxl, nvmax, numVs, listVsptr, nuo, &
            spin%H, DscfL, Vs )
       call de_alloc( t_DscfL, 'DscfL', 'vmat' )
    else if ( NTH > 1 ) then
       call de_alloc( t_Vss, 'Vss', 'vmat' )
    end if
!$OMP end master

!$OMP end parallel

    call de_alloc( r2cut, 'r2cut', 'vmat' )

#ifdef _TRACE_
    call MPI_Barrier( MPI_Comm_World, MPIerror )
    call MPItrace_event( 1000, 0 )
#endif

    call timer('vmat',2)

#ifdef DEBUG
    call write_debug( '    POS vmat' )
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

  end subroutine vmat

  subroutine dvmat( no, np, dvol, nspin, Vscf, nvmax, &
                      numVs, listVsptr, listVs, Vs, &
                      nuo, nuotot, iaorb, iphorb, isa, &
                      ialr )

! ********************************************************************
! LR 2015: Mixture of vmat.f and dfscf.f (Ordejon, Soler and Gale), it 
! computes the matrix elements <dPhi|V|Phi> + <Phi|V|dPhi> also 
! referred to as the Pulay terms that arise due to the localized
! orbitals, necessary for the perturbed hamiltonian matrix elements
! to compute the change in the dRho for DFPT in LinRes. It is non
! SCF so it is called only once on the first scf iteration.
! *********************** INPUT **************************************
! integer no              : Number of basis orbitals
! integer nuo             : Number of orbitals in unit cell (local)
! integer nuotot          : Number of orbitals in unit cell (global)
! integer np              : Number of mesh points (total is nsp*np)
! integer nspin           : Number of spin components
! integer isa(na)         : Species index of each atom
! integer iaorb(no)       : Atom to which orbitals belong
! integer iphorb(no)      : Index of orbital within its atom
! integer nvmax           : First dimension of Vs
! integer numVs(nuo)       : Number of nonzero elemts in each row of Dscf
! integer listVsptr(nuo)   : Pointer to start of row in listd
! integer listVs(dvmax)    : List of nonzero elements of Dscf
! real*4  dVs(nsp,np,nspin): Value of SCF potential at the mesh points
! *********************** INPUT and OUTPUT ****************************
! real*8  Vs    : Perturbed Hamitonian matrix
! *********************************************************************
!    6  10        20        30        40        50        60        7072

!  Modules
      use precision,     only: dp, grid_p
      use atmfuncs,      only: rcut, all_phi
      use atm_types,     only: nsmax=>nspecies
      use atomlist,      only: indxuo, indxua
      use listsc_module, only: LISTSC
      use mesh,          only: dxa, nsp, xdop, xdsp, meshLim
      use meshdscf,      only: matrixMtoO
      use meshdscf,      only: needdscfl, listdl, numdl, nrowsdscfl
      use meshdscf,      only: listdlptr
      use meshphi,       only: DirectPhi, endpht, lstpht, listp2, phi
      use meshphi,       only: gradphi
      use parallel,      only: Nodes, node
      use alloc,         only: re_alloc, de_alloc, alloc_default
      use alloc,         only: allocDefaults
      use parallelsubs,  only: GlobalToLocalOrb
#ifdef MPI
      use mpi_siesta
#endif
#ifdef _OPENMP
      use omp_lib
#endif
!  Argument types and dimensions
      integer                  :: no, nuo, nuotot, np, nspin
      integer                  :: isa(*), iaorb(*), iphorb(*), nvmax
      integer                  :: numVs(nuo), listVsptr(nuo)
      integer                  :: listVs(nvmax)
      integer,      intent(in) :: ialr
      real(grid_p), intent(in) :: Vscf(nsp,np,nspin)
      real(dp)                 :: dvol, VolCel
      real(dp),         target :: Vs(nvmax,3,nspin)

! Internal variables
      integer,       parameter :: minloc  = 1000  ! Min buffer size 
      integer,       parameter :: maxoa = 100   ! Max # of orbitals per atom
      integer                  :: i, ia, ic, ii, imp
      integer                  :: ind, iop, ip, iphi, io, is, isp
      integer                  :: ispin, iu, iua, iul, ix
      integer                  :: nlocal, j, jc, last, lasta
      integer                  :: lastop, maxloc, maxloc2, NTH,TID
      integer                  :: nc, nphiloc, index(no)
      integer                  :: in, il, ij,  jl, triang, ijl
      integer                  :: lenx, leny, lenz, lenxy, nvmaxl

      logical                  :: ParallelLocal

      real(dp)                 :: dxsp(3),  Vij(3)
      real(dp)                 :: r2sp
      real(dp)                 :: V(nsp,nspin), Vpart(3,2,nsp)

      integer,         pointer :: ilc(:), ilocal(:), iorb(:)

      real(dp),        pointer :: C(:,:), gC(:,:,:), phia(:,:) 
      real(dp),        pointer :: grada(:,:,:)
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

!  Start time counter
      call timer('dvmat',1)

!  Find atomic cutoff radii
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

!  Find value of maxloc 
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

!  Nullify pointers
      nullify( C, gC, Vi, ilc, ilocal, phia, grada)

      call alloc_default( old=oldDefaults, &
                         copy=.false., shrink=.false., &
                         imin=1, routine='dvmat' )

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
         call re_alloc( t_Vsp, 1, nvmaxl, 1, 3, 1 ,nspin, 1, NTH, &
              'Vsp',  'dvmat' )
      else
        if ( NTH > 1 ) then
          nullify( t_Vss )
          call re_alloc( t_Vss, 1, nvmax, 1, 3, 1, nspin, 2, NTH, &
              'Vss',  'dvmat' )
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

!  Initialise variables
      last                     = 0
      Vi(1:triang,1:3,1:nspin) = 0.0_dp
      ilocal(1:no)             = 0
      iorb(1:no)               = 0
      ilc(:)                   = 0

!  Loop over grid points
!$OMP do
      do ip = 1,np
!       Find number of nonzero orbitals at this point
        nc = endpht(ip) - endpht(ip-1)
!       new required size of Vlocal 
        nlocal = last
        do ic = 1,nc
         imp = endpht(ip-1) + ic
         i = lstpht(imp)
         if (ilocal(i) .eq. 0) nlocal = nlocal + 1
        enddo

!       If overflooded, add Vlocal to Vs and reinitialize it
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
                         Vsp(ind,:,ispin) = Vsp(ind,:,ispin) + &
                                          Vi(ijl,:,ispin) * dVol
                      end do
                   end do
                else
                   do ii = 1, numdl(iul)
                      ind = listdlptr(iul) + ii
                      j   = LISTSC( i, iu, listdl(ind) )
                      ijl = idx_ijl(il,ilocal(j))
                      do ispin = 1, nspin
                         Vsp(ind,:,ispin) = Vsp(ind,:,ispin) + &
                                        Vi(ijl,:,ispin) * dVol
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
                      Vss(ind,1:3,ispin) = Vss(ind,1:3,ispin) + dVol * &
                     Vi(ijl,1:3,ispin)
                    enddo
                  enddo
                else !i .eq. iu
                  do ii = 1, numVs(iul)
                    ind = listVsptr(iul)+ii
                    j = LISTSC( i, iu, listVs(ind) )
                    ijl = idx_ijl(il,ilocal(j))
                    do ispin = 1,nspin
                      Vss(ind,1:3,ispin) = Vss(ind,1:3,ispin) + dVol * &
                     Vi(ijl,1:3,ispin)
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

!  Look for required orbitals not yet in Vlocal
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

!  Copy potential to a double precision array
        V(1:nsp,1:nspin) = Vscf(1:nsp,ip,1:nspin)

!  Calculate all phi values and derivatives at all subpoints
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
                     call all_phi( is,+1, dxsp, nphiloc, &
                             phia(:,isp), grada(:,:,isp))
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
!    Loop on the second orbital of mesh point (only for jc.<=.ic)
              do jc= 1,ic !tringular matrix form
                jl=ilc(jc)
                Vij(:)=0.0_dp
                do isp= 1,nsp !sum over subpoints
                  Vij(:)=Vij(:)+Vpart(:,1,isp)*C(isp,jc) + &
                         gC(:,isp,jc)*Vpart(:,2,isp)
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
!    Loop on the second orbital of mesh point (only for jc.le.ic)
              do jc= 1,ic !tringular matrix form
                jl=ilc(jc)
                Vij(:)=0.0_dp
                do isp= 1,nsp !sum over subpoints
                Vij(:)=Vij(:)+Vpart(:,1,isp)*C(isp,jc) + &
                         gC(:,isp,jc)*Vpart(:,2,isp)
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
                   Vsp(ind,1:3,ispin) = Vsp(ind,1:3,ispin) + &
                                Vi(ijl,1:3,ispin) * dVol
                enddo
             enddo
          else
             do ii = 1, numdl(iul)
                ind = listdlptr(iul) + ii
                j   = LISTSC( i, iu, listdl(ind) )
                ijl = idx_ijl(il,ilocal(j))
                do ispin = 1, nspin
                   Vsp(ind,1:3,ispin) = Vsp(ind,1:3,ispin) + &
                                 Vi(ijl,1:3,ispin) * dVol
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
                 Vss(ind,1:3,ispin) = Vss(ind,1:3,ispin) + dVol * &
                 Vi(ijl,1:3,ispin)
               enddo
             enddo !ii loop
           else
             do ii = 1, numVs(iu)
               ind = listVsptr(iu)+ii
               j = LISTSC( i, iu, listVs(ind) )
               jl = ilocal(j)
               ijl = idx_ijl(il,jl)
               do ispin = 1,nspin
                 Vss(ind,1:3,ispin) = Vss(ind,1:3,ispin) + dVol * &
                 Vi(ijl,1:3,ispin)
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
                t_Vsp(ind,:,ispin,1) = t_Vsp(ind,:,ispin,1) + &
                      t_Vsp(ind,:,ispin,ii)
             end do
          end do
        end do
!$OMP end do
      else if ( NTH > 1 ) then
!$OMP do collapse(2)
        do ispin = 1 , nspin
          do ind = 1, nvmax
             do ii = 2, NTH
                Vs(ind,:,ispin) = Vs(ind,:,ispin) + &
                                  t_Vss(ind,:,ispin,ii)
             end do
          end do
        end do
!$OMP end do
      end if

!  Deallocate local memory
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
           call matrixMtoO( nvmaxl, nvmax, numVs, listVsptr, nuo, &
            nspin, Vsp(:,ix,:), Vs(:,ix,:) )
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

!  Restore old allocation defaults
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


end module m_vmat
