! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!
!  The old dhscf has been split in two parts: an initialization routine
!  (dhscf_init), which is called after every geometry change but before
!  the main scf loop, and a dhscf proper, which is called at every step
!  of the scf loop
!
!  NOTE that the mesh initialization part is now done *unconditionally*
!  in dhscf_init, i.e, after *every* geometry change, even if the
!  change does not involve a cell change. The reason is to avoid 
!  complexity, since now the mesh parallel distributions will depend on 
!  the detailed atomic positions even if the cell does not change.
!
!  Besides, the relative cost of a "mesh only" initialization is negligible.
!  The only real observable effect would be a printout of "initmesh" data
!  at every geometry iteration.
!
      module m_dhscf

!
!     To facilitate the communication among dhscf_init and dhscf,
!     some arrays that hold data which do not change during the SCF loop
!     have been made into module variables

!     Some others are scratch, such as nmpl, ntpl, etc

      use precision,      only : dp, grid_p
      use m_dfscf,        only : dfscf
      implicit none

      real(grid_p),   pointer :: rhopcc(:), rhoatm(:), Vna(:)
      real(dp)                :: Uharrs     ! Harris energy
      logical                 :: IsDiag, spiral
      character(len=10)       :: shape
      integer                 :: nml(3), ntml(3), npcc,
     &                           nmpl, ntpl
      real(dp)                :: bcell(3,3), cell(3,3),
     &                           dvol, field(3), rmax, scell(3,3)
      real(dp)                :: G2mesh =  0.0_dp
      logical :: debug_dhscf = .false.
      character(len=*), parameter :: debug_fmt =
     &     '('' -- dhscf : Node '',i3,tr1,a,4(tr1,e20.7))'

      public :: dhscf_init, dhscf

      CONTAINS

      subroutine dhscf_init(nspin, norb, iaorb, iphorb,
     &                      nuo, nuotot, nua, na,
     &                      isa, xa, indxua, ucell,
     &                      mscell, G2max, ntm,
     &                      maxnd, numd, listdptr, listd, datm,
     &                      Fal, stressl)

      use precision,      only : dp, grid_p
      use parallel,       only : Node, Nodes
      use atmfuncs,       only : rcut, rcore
      use fdf
      use sys,            only : die
      use mesh,           only : xdsp, nsm, nsp, meshLim
      use parsing
#ifndef BSC_CELLXC
      use siestaXC, only : getXC        ! Returns the XC functional used
#else /* BSC_CELLXC */
      use bsc_xcmod,          only : nXCfunc, XCauth
#endif /* BSC_CELLXC */
      use alloc,          only : re_alloc, de_alloc
      use siesta_options, only : harrisfun
      use meshsubs,       only : PartialCoreOnMesh
      use meshsubs,       only : NeutralAtomOnMesh
      use meshsubs,       only : PhiOnMesh
      use meshsubs,       only : InitMesh
      use meshsubs,       only : InitAtomMesh
      use meshsubs,       only : setupExtMesh
      use meshsubs,       only : distriPhiOnMesh

      use moreMeshSubs,   only : setMeshDistr, distMeshData
      use moreMeshSubs,   only : UNIFORM, QUADRATIC, LINEAR
      use moreMeshSubs,   only : TO_SEQUENTIAL, TO_CLUSTER, KEEP

      use meshdscf,       only : createLocalDscfPointers

      use iogrid_netcdf, only: set_box_limits
#ifdef NCDF_4
      use m_ncdf_io, only : cdf_init_mesh
#endif
#ifdef BSC_CELLXC
      use cellxc_mod,     only : setGGA
#endif /* BSC_CELLXC */
      use m_efield,       only : initialize_efield, acting_efield
      use m_efield,       only : get_field_from_dipole
      use m_efield,       only : dipole_correction
      use m_efield,       only : user_specified_field

      use m_doping_uniform,       only: initialize_doping_uniform
      use m_doping_uniform,       only: compute_doping_structs_uniform,
     $                                  doping_active
      use m_rhog,                 only: rhog, rhog_in
      use m_rhog,                 only: order_rhog
      use siesta_options,         only: mix_charge
#ifdef MPI
      use mpi_siesta
#endif

      use m_mesh_node,   only: init_mesh_node
      use m_charge_add,  only: init_charge_add
      use m_hartree_add, only: init_hartree_add

      use m_ts_global_vars,only: TSmode
      use m_ts_options,   only : IsVolt, N_Elec, Elecs
      use m_ts_voltage,   only : ts_init_voltage
      use m_ts_hartree,   only : ts_init_hartree_fix

      implicit none
      integer, intent(in)     :: nspin, norb, iaorb(norb), iphorb(norb),
     &                           nuo, nuotot, nua, na, isa(na),
     &                           indxua(na), mscell(3,3), maxnd,
     &                           numd(nuo), listdptr(nuo), listd(maxnd)
      real(dp), intent(in)    :: xa(3,na), ucell(3,3), datm(norb)
      real(dp), intent(inout) :: g2max
      integer, intent(inout)  :: ntm(3)
      real(dp), intent(inout) :: Fal(3,nua), stressl(3,3)

      real(dp), parameter     :: tiny  = 1.e-12_dp 
      integer                 :: io, ia, iphi, is, n, i, j
      integer                 :: nsc(3), nbcell, nsd
      real(dp)                :: DStres(3,3), volume
      real(dp), external      :: volcel, ddot
      real(grid_p)            :: dummy_Drho(1,1), dummy_Vaux(1),
     &                           dummy_Vscf(1)
      logical, save           :: frstme = .true.   ! Keeps state
      real(grid_p),   pointer :: Vscf(:,:), rhoatm_par(:)
      integer,        pointer :: numphi(:), numphi_par(:)

      integer                 :: nm(3)   ! For call to initMesh
#ifndef BSC_CELLXC
      integer                 :: nXCfunc
      character(len=20)       :: XCauth(10), XCfunc(10)
#endif /* ! BSC_CELLXC */
      ! Transport direction (unit-cell aligned)
      integer                 :: iE
      real(dp)                :: ortho, field(3), field2(3)
!--------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '    PRE dhscf_init' )
#endif
C ----------------------------------------------------------------------
C     General initialisation
C ----------------------------------------------------------------------
C     Start time counter
      call timer( 'DHSCF_Init', 1 )

      nsd = min(nspin,2)
      nullify(Vscf,rhoatm_par)

      if (frstme) then
        debug_dhscf = fdf_get('Debug.DHSCF', .false.)
         
        nullify( xdsp, rhopcc, Vna, rhoatm )
C       nsm lives in module m_dhscf now    !! AG**
        nsm = fdf_integer( 'MeshSubDivisions', 2 )
        nsm = max( nsm, 1 )

C       Set mesh sub-division variables & perform one off allocation
        nsp = nsm*nsm*nsm
        call re_alloc( xdsp, 1, 3, 1, nsp, 'xdsp', 'dhscf_init' )

C       Check spin-spiral wavevector (if defined)
        if (spiral .and. nspin.lt.4)
     &    call die('dhscf: ERROR: spiral defined but nspin < 4')
      endif   ! First time

#ifndef BSC_CELLXC
C Get functional(s) being used
      call getXC( nXCfunc, XCfunc, XCauth )
#endif /* ! BSC_CELLXC */

      if (harrisfun) then
        do n = 1,nXCfunc
          if (.not.(leqi(XCauth(n),'PZ').or.leqi(XCauth(n),'CA'))) then
            call die("** Harris forces not implemented for non-LDA XC")
          endif
        enddo
      endif

C ----------------------------------------------------------------------
C     Orbital initialisation : part 1
C ----------------------------------------------------------------------

C     Find the maximum orbital radius
      rmax = 0.0_dp
      do io = 1, norb
        ia   = iaorb(io)    ! Atomic index of each orbital
        iphi = iphorb(io)   ! Orbital index of each  orbital in its atom
        is   = isa(ia)      ! Species index of each atom
        rmax = max( rmax, rcut(is,iphi) )
      enddo

C     Start time counter for mesh initialization
      call timer( 'DHSCF1', 1 )

C ----------------------------------------------------------------------
C     Unit cell handling
C ----------------------------------------------------------------------
C     Find diagonal unit cell and supercell
      call digcel( ucell, mscell, cell, scell, nsc, IsDiag )
      if (.not.IsDiag) then
        if (Node.eq.0) then
          write(6,'(/,a,3(/,a,3f12.6,a,i6))')
     &      'DHSCF: WARNING: New shape of unit cell and supercell:',
     &      ('DHSCF:',(cell(i,j),i=1,3),'   x',nsc(j),j=1,3)
        endif
      endif

C     Find the system shape
      call shaper( cell, nua, isa, xa, shape, nbcell, bcell )

C     Find system volume
      volume = volcel( cell )

C ----------------------------------------------------------------------
C     Mesh initialization 
C ----------------------------------------------------------------------
      call InitMesh( na, cell, norb, iaorb, iphorb, isa, rmax, 
     &               G2max, G2mesh, nsc, nmpl, nm,
     &               nml, ntm, ntml, ntpl, dvol )

!     Setup box descriptors for each processor,
!     held in module iogrid_netcdf
      call set_box_limits( ntm, nsm )

      ! Initialize information on local mesh for each node
      call init_mesh_node( cell, ntm , meshLim , nsm )

      ! Setup charge additions in the mesh
      call init_charge_add(cell, ntm)

      ! Setup Hartree additions in the mesh
      call init_hartree_add(cell, ntm)

#ifdef NCDF_4
      ! Initialize the box for each node...
      call cdf_init_mesh( ntm, nsm )
#endif

C     Stop time counter for mesh initialization
      call timer( 'DHSCF1', 2 )
C ----------------------------------------------------------------------
C     End of mesh initialization 
C ----------------------------------------------------------------------

C ----------------------------------------------------------------------
C     Initialize atomic orbitals, density and potential
C ----------------------------------------------------------------------
C     Start time counter for atomic initializations
      call timer( 'DHSCF2', 1 )

      if (nodes.gt.1) then
        call setMeshDistr( UNIFORM, nsm, nsp,
     &                     nml, nmpl, ntml, ntpl )
      endif

C     Initialise quantities relating to the atom-mesh positioning
      call InitAtomMesh( UNIFORM, na, xa )

#ifdef BSC_CELLXC
C     Check if we need extencils in cellxc
      call setGGA( )

#endif /* BSC_CELLXC */
C     Compute the number of orbitals on the mesh and recompute the
C     partions for every processor in order to have a similar load
C     in each of them.
      nullify( numphi )
      call re_alloc( numphi, 1, nmpl, 'numphi', 'dhscf_init' )
!$OMP parallel do default(shared), private(i)
      do i= 1, nmpl
        numphi(i) = 0
      enddo
!$OMP end parallel do

      call distriPhiOnMesh( nm, nmpl, norb, iaorb, iphorb,
     &                      isa, numphi )

C     Find if there are partial-core-corrections for any atom
      npcc = 0
      do ia = 1,na
        if (rcore(isa(ia)) .gt. tiny) npcc = 1
      enddo

C     Find partial-core-correction energy density
C     Vscf and Vaux are not used here
      call re_alloc( rhopcc, 1, ntpl*npcc+1, 'rhopcc', 'dhscf_init' )
      if (npcc .eq. 1) then
        call PartialCoreOnMesh( na, isa, ntpl, rhopcc, indxua,
     &    nsd, dvol, volume, dummy_Vscf, dummy_Vaux, Fal, stressl,
     &    .false., .false. )
        call reord( rhopcc, rhopcc, nml, nsm, TO_SEQUENTIAL )
        if ( debug_dhscf ) then
           write(*,debug_fmt) Node,'rhopcc',sqrt(sum(rhopcc**2))
        end if
      endif

C     Find neutral-atom potential
C     Drho is not used here
      call re_alloc( Vna, 1, ntpl, 'Vna', 'dhscf_init' )
      call NeutralAtomOnMesh( na, isa, ntpl, Vna, indxua, dvol, 
     &                        volume, dummy_DRho, Fal, stressl, 
     &                        .false., .false. )
      call reord( Vna, Vna, nml, nsm, TO_SEQUENTIAL )
      if ( debug_dhscf ) then
         write(*,debug_fmt) Node,'Vna',sqrt(sum(Vna**2))
      end if

      if (nodes.gt.1) then
        if (node .eq. 0) then
           write(6,"(a)") "Setting up quadratic distribution..."
        endif
        call setMeshDistr( QUADRATIC, nsm, nsp,
     &                     nml, nmpl, ntml, ntpl )

C       Create extended mesh arrays for the second data distribution
        call setupExtMesh( QUADRATIC, rmax )

C       Compute atom positions for the second data distribution 
        call InitAtomMesh( QUADRATIC, na, xa )
      endif

C     Calculate orbital values on mesh
!     numphi has already been computed in distriPhiOnMesh
!     in the UNIFORM distribution 
      if (nodes.eq.1) then
        numphi_par => numphi
      else
        nullify(numphi_par)
        call re_alloc( numphi_par, 1, nmpl, 'numphi_par',
     &                 'dhscf_init' )
        call distMeshData( UNIFORM, numphi, QUADRATIC,
     &                     numphi_par, KEEP )
      endif

      call PhiOnMesh( nmpl, norb, iaorb, iphorb, isa, numphi_par )

      if (nodes.gt.1) then
        call de_alloc( numphi_par, 'numphi_par', 'dhscf_init' )
      endif
      call de_alloc( numphi, 'numphi', 'dhscf_init' )

C ----------------------------------------------------------------------
C       Create sparse indexing for Dscf as needed for local mesh
C       Note that this is done in the QUADRATIC distribution
C       since 'endpht' (computed finally in PhiOnMesh and stored in
C       meshphi module) is in that distribution.
C ----------------------------------------------------------------------
      if (Nodes.gt.1) then
        call CreateLocalDscfPointers( nmpl, nuotot, numd, listdptr, 
     &                                listd )
      endif

C ----------------------------------------------------------------------
C     Calculate terms relating to the neutral atoms on the mesh
C ----------------------------------------------------------------------
C     Find Harris (sum of atomic) electron density
      call re_alloc( rhoatm_par, 1, ntpl, 'rhoatm_par', 'dhscf_init' )
      call rhooda( norb, nmpl, datm, rhoatm_par, iaorb, iphorb, isa )
!     rhoatm_par comes out of here in clustered form in QUADRATIC dist

C     Routine Poison should use the uniform data distribution
      if (nodes.gt.1) then
         call setMeshDistr( UNIFORM, nsm, nsp,
     &                      nml, nmpl, ntml, ntpl )
      endif

C     Create Rhoatm using UNIFORM distr, in sequential form
      call re_alloc( rhoatm, 1, ntpl, 'rhoatm', 'dhscf_init' )
      call distMeshData( QUADRATIC, rhoatm_par,
     &     UNIFORM, rhoatm, TO_SEQUENTIAL )

      if ( debug_dhscf ) then
         write(*,debug_fmt) Node,'rhoatm', sqrt(sum(rhoatm**2))
      end if

!
!  AG: The initialization of doping structs could be done here now,
!      in the uniform distribution, and with a simple loop over
!      rhoatm.

        if (frstme) call initialize_doping_uniform()
        if (doping_active) then  
           call compute_doping_structs_uniform(ntpl,rhoatm,nsd)
           ! Will get the global number of hit points
           ! Then, the doping density to be added can be simply computed
        endif

C     Allocate Temporal array
      call re_alloc( Vscf, 1, ntpl, 1, nspin, 'Vscf', 'dhscf_init' )

C     Vscf is filled here but not used later
C     Uharrs is computed (and saved)
C     DStres is computed but not used later
      call poison( cell, ntml(1), ntml(2), ntml(3), ntm, rhoatm,
     &             Uharrs, Vscf, DStres, nsm )
      call de_alloc( Vscf, 'Vscf', 'dhscf_init' )

!     Always deallocate rhoatm_par, as it was used even if nodes=1
      call de_alloc( rhoatm_par, 'rhoatm_par', 'dhscf_init' )

      if (mix_charge) then
         call re_alloc( rhog, 1, 2, 1, ntpl, 1, nspin,
     $                 'Rhog', 'dhscf_init' )
         call re_alloc( rhog_in, 1, 2, 1, ntpl, 1, nspin,
     $                 'Rhog_in', 'dhscf_init' )
         call order_rhog(cell, ntml(1), ntml(2), ntml(3), ntm, nsm)
      endif

C     Stop time counter for atomic initializations
      call timer( 'DHSCF2', 2 )

C ----------------------------------------------------------------------
C     At the end of initializations:
!     We are in the UNIFORM distribution
!     Rhoatm, Rhopcc and Vna are in UNIFORM dist, sequential form
!     The index array endpht is in the QUADRATIC distribution
C ----------------------------------------------------------------------
      if (frstme) then
         call initialize_efield()
      end if

      ! Check if we need to add the potential 
      ! corresponding to the voltage-drop.
      if ( TSmode ) then
        ! These routines are important if there are cell-changes
        call ts_init_hartree_fix( cell, nua, xa, ntm, ntml)
        if ( IsVolt ) then
           call ts_init_voltage( cell, nua, xa, ntm)
        end if

        if ( acting_efield ) then
         ! We do not allow the electric field for 
         ! transiesta runs with V = 0, either.
         ! It does not make sense, only for fields perpendicular
         ! to the applied bias.
           
         ! We need to check that the e-field is perpendicular
         ! to the transport direction, and that the system is
         ! either a chain, or a slab.
         ! However, due to the allowance of a dipole correction
         ! along the transport direction for buffer calculations
         ! we have to allow all shapes. (atom is not transiesta
         ! compatible anyway)
            
         ! check that we do not violate the periodicity
         if ( Node .eq. 0 ) then
            write(*,'(/,2(2a,/))') 'ts-WARNING: ',
     &           'E-field/dipole-correction! ',
     &           'ts-WARNING: ',
     &           'I hope you know what you are doing!'
         end if

         ! This is either dipole or user, or both
         field(:) = user_specified_field(:)
         do iE = 1 , N_Elec
            field2 = Elecs(iE)%cell(:,Elecs(iE)%t_dir)
            ortho = ddot(3,field2,1,field,1)
            if ( abs(ortho) > 1.e-9_dp ) then
               call die('User defined E-field must be
     &perpendicular to semi-infinite directions')
            end if
         end do
        end if ! acting_efield
         
        ! We know that we currently allow people to do more than
        ! they probably should be allowed. However, there are many
        ! corner cases that may require dipole corrections, or 
        ! electric fields to "correct" an intrinsic dipole.
        ! For instance, what should we do with a dipole in a transiesta
        ! calculation?
        ! Should we apply a field to counter act it in a device
        ! calculation?
         
      end if

      frstme = .false.

      call timer( 'DHSCF_Init', 2 )
#ifdef DEBUG
      call write_debug( '    POS dhscf_init' )
#endif
!------------------------------------------------------------------------- END
      end subroutine dhscf_init

      subroutine dhscf( nspin, norb, iaorb, iphorb, nuo, 
     &                  nuotot, nua, na, isa, xa, indxua,
     &                  ntm, ifa, istr, iHmat,
     &                  filesOut, maxnd, numd,
     &                  listdptr, listd, Dscf, datm, maxnh, Hmat,
     &                  Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, 
     &                  Exc, Dxc, dipol, stress, Fal, stressl,
     &                  use_rhog_in, charge_density_only, iai, iaf,
     &                  ialr, lasto, dynmat, dDscf, dHmat,
     &                  LRfirst, dHmat0)

C
C Calculates the self-consistent field contributions to Hamiltonian
C matrix elements, total energy and atomic forces.
C Coded by J.M. Soler, August 1996. July 1997.
C Modified by J.D. Gale, February 2000.
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer nspin         : Number of different spin polarisations
C                         nspin=1 => Unpolarized, nspin=2 => polarized
C                         nspin=4 => Noncollinear spin or spin-orbit
C integer norb          : Total number of basis orbitals in supercell
C integer iaorb(norb)   : Atom to which each orbital belongs
C integer iphorb(norb)  : Orbital index (within atom) of each orbital
C integer nuo           : Number of orbitals in a unit cell in this node
C integer nuotot        : Number of orbitals in a unit cell
C integer nua           : Number of atoms in unit cell
C integer na            : Number of atoms in supercell
C integer isa(na)       : Species index of all atoms in supercell
C real*8  xa(3,na)      : Atomic positions of all atoms in supercell
C integer indxua        : Index of equivalent atom in unit cell
C integer ifa           : Switch which fixes whether the SCF contrib.
C                         to atomic forces is calculated and added to fa
C integer istr          : Switch which fixes whether the SCF contrib.
C                         to stress is calculated and added to stress
C integer iHmat         : Switch which fixes whether the Hmat matrix
C                         elements are calculated or not.
C type(filesOut_t) filesOut : output file names (If blank => not saved)
C integer maxnd             : First dimension of listd and Dscf
C integer numd(nuo)         : Number of nonzero density-matrix
C                             elements for each matrix row
C integer listdptr(nuo)     : Pointer to start of rows of density-matrix
C integer listd(maxnd)      : Nonzero-density-matrix-element column
C                             indexes for each matrix row
C real*8  Dscf(maxnd,h_spin_dim): SCF density-matrix elements
C real*8  datm(norb)        : Harris density-matrix diagonal elements
C                             (atomic occupation charges of orbitals)
C integer maxnh             : First dimension of listh and Hmat
C real*8  Hmat(maxnh,h_spin_dim) : Hamiltonian matrix in sparse form,
C                             to which are added the matrix elements
C                                 <ORB_I | DeltaV | ORB_J>, where
C                             DeltaV = Vna + Vxc(SCF) + 
C                                      Vhartree(RhoSCF-RhoHarris)
C ----------------------------------------------------------------------
C Input/output :
C ----------------------------------------------------------------------
C integer ntm(3) : Number of mesh divisions of each cell
C                  vector, including subgrid.
C ----------------------------------------------------------------------
C Output :
C ----------------------------------------------------------------------
C real*8  Enaatm : Integral of Vna * rhoatm
C real*8  Enascf : Integral of Vna * rhoscf
C real*8  Uatm   : Harris hartree electron-interaction energy
C real*8  Uscf   : SCF hartree electron-interaction energy
C real*8  DUscf  : Electrostatic (Hartree) energy of
C                    (rhoscf - rhoatm) density
C real*8  DUext  : Interaction energy with external electric field
C real*8  Exc    : SCF exchange-correlation energy
C real*8  Dxc    : SCF double-counting correction to Exc
C                    Dxc = integral of ( (epsxc - Vxc) * Rho )
C                    All energies in Rydbergs
C real*8  dipol(3): Electric dipole (in a.u.)
C                   only when the system is a molecule
C ----------------------------------------------------------------------
C Input/output :
C ----------------------------------------------------------------------
C real*8  Fal(3,nua) : Atomic forces, to which the SCF contribution
C                        is added by this routine when ifa=1.
C                        the SCF contribution is minus the derivative
C                        of ( Enascf - Enaatm + DUscf + Exc ) with
C                        respect to atomic positions, in  Ry/Bohr.
C                        contributions local to this node
C real*8 stressl(3,3): Stress tensor, to which the SCF contribution
C                      is added by this routine when ifa=1.
C                      the SCF contribution is minus the derivative of
C                         ( Enascf - Enaatm + DUscf + Exc ) / volume
C                      with respect to the strain tensor, in Ry.
C                        contributions local to this node
C ----------------------------------------------------------------------
C Units :
C ----------------------------------------------------------------------
C Energies in Rydbergs
C Distances in Bohr
C ----------------------------------------------------------------------
C     Modules
      use precision,     only  : dp, grid_p
#ifndef BSC_CELLXC
      use parallel, only : ProcessorY
#endif /* ! BSC_CELLXC */

!     Number of Mesh divisions of each cell vector (global)
!     The status of this variable is confusing
      use parallel,      only  : Node, Nodes
      use atmfuncs,      only  : rcut, rcore
      use units,         only  : Debye, eV, Ang
      use fdf
      use sys,           only  : die, bye
      use mesh,          only  : nsm, nsp
      use parsing
      use m_iorho,       only  : write_rho
      use m_forhar,      only  : forhar
      use alloc,         only  : re_alloc, de_alloc
      use files,         only  : slabel
      use files,         only  : filesOut_t ! derived type for output file names
      use siesta_options, only : harrisfun, save_initial_charge_density
      use siesta_options, only : analyze_charge_density_only
      use meshsubs,       only : LocalChargeOnMesh
      use meshsubs,       only : PartialCoreOnMesh
      use meshsubs,       only : NeutralAtomOnMesh

      use moreMeshSubs,   only : setMeshDistr, distMeshData
      use moreMeshSubs,   only : UNIFORM, QUADRATIC, LINEAR
      use moreMeshSubs,   only : TO_SEQUENTIAL, TO_CLUSTER, KEEP
      use m_partial_charges, only: compute_partial_charges
      use m_partial_charges, only: want_partial_charges
#ifndef BSC_CELLXC

      use siestaXC, only : cellXC       ! Finds xc energy and potential
      use siestaXC, only : myMeshBox    ! Returns my processor mesh box
      use siestaXC, only : jms_setMeshDistr => setMeshDistr
                                        ! Sets a distribution of mesh
                                        ! points over parallel processors
#endif /* BSC_CELLXC */
      use m_vmat,  only  : vmat, dvmat
      use m_rhoofd, only: rhoofd
#ifdef MPI
      use mpi_siesta
#endif
      use iogrid_netcdf, only: write_grid_netcdf
      use iogrid_netcdf, only: read_grid_netcdf
      use siesta_options, only: read_charge_cdf
      use siesta_options, only: savebader
      use siesta_options, only: read_deformation_charge_cdf
      use siesta_options, only: mix_charge

      use m_efield,       only: get_field_from_dipole, dipole_correction
      use m_efield,       only: add_potential_from_field
      use m_efield,       only: user_specified_field, acting_efield
      use m_doping_uniform,       only: doping_active, doping_uniform
      
      use m_charge_add,   only: charge_add
      use m_hartree_add,  only: hartree_add

#ifdef NCDF_4
      use siesta_options, only: write_cdf
      use m_ncdf_siesta, only: cdf_save_grid
#endif
      use m_rhofft,       only: rhofft, FORWARD, BACKWARD
      use m_rhog,         only: rhog_in, rhog
      use m_spin,         only: spin
      use m_spin,         only: Spiral, qSpiral
      use m_iotddft,      only: write_tdrho
      use m_ts_global_vars,only: TSmode, TSrun
      use m_ts_options, only: IsVolt, Elecs, N_elec
      use m_ts_voltage, only: ts_voltage
      use m_ts_hartree, only: ts_hartree_fix


      use linres_matrices   !Linres line

      implicit none

      integer
     &  maxnd, maxnh, nua, na, norb, nspin, nuo, nuotot, 
     &  iaorb(norb), ifa, iHmat,
     &  indxua(na), iphorb(norb), isa(na), 
     &  istr, listd(*), listdptr(nuo), ntm(3), numd(nuo)

      real(dp)
     &  datm(norb), dipol(3), Dscf(:,:),
     &  DUext, DUscf, Dxc, Enaatm, Enascf, Exc,
     &  Hmat(:,:), Uatm, Uscf, xa(3,na),
     &  stress(3,3), Fal(3,nua), stressl(3,3)

      type(filesOut_t) filesOut

      logical, intent(in), optional :: use_rhog_in
      logical, intent(in), optional :: charge_density_only

C Linres optional variables-------------------------------------------------
      integer, intent(in), optional    :: iai,iaf,ialr,lasto(0:na)
      logical, intent(in), optional    :: LRfirst

      real(dp), intent(in), optional :: dDscf(maxnh,nspin,3)

      real(dp), intent(inout),optional :: dynmat(3,nua,3,nua),
     &                                    dHmat0(maxnh,3,nspin),
     &                                    dHmat(maxnh,nspin,3)

C---------------------------------------------------------------------------

C ----------------------------------------------------------------------
C Routines called internally:
C ----------------------------------------------------------------------
C        cellxc(...)    : Finds total exch-corr energy and potential
C        cross(a,b,c)   : Finds the cross product of two vectors
C        dfscf(...)     : Finds SCF contribution to atomic forces
C        dipole(...)    : Finds electric dipole moment
C        doping(...)    : Adds a background charge for doped systems
C        write_rho(...)     : Saves electron density on a file
C        poison(...)    : Solves Poisson equation
C        reord(...)     : Reorders electron density and potential arrays
C        rhooda(...)    : Finds Harris electron density in the mesh
C        rhoofd(...)    : Finds SCF electron density in the mesh
C        rhoofdsp(...)  : Finds SCF electron density in the mesh for
C                         spiral arrangement of spins
C        timer(...)     : Finds CPU times
C        vmat(...)      : Finds matrix elements of SCF potential
C        vmatsp(...)    : Finds matrix elements of SCF potential for
C                         spiral arrangement of spins
C        delk(...)      : Finds matrix elements of exp(i \vec{k} \cdot \vec{r})
C real*8 volcel( cell ) : Returns volume of unit cell
C ----------------------------------------------------------------------
C Internal variables and arrays:
C ----------------------------------------------------------------------
C real*8  bcell(3,3)    : Bulk lattice vectors
C real*8  cell(3,3)     : Auxiliary lattice vectors (same as ucell)
C real*8  const         : Auxiliary variable (constant within a loop)
C real*8  DEc           : Auxiliary variable to call cellxc
C real*8  DEx           : Auxiliary variable to call cellxc
C real*8  dvol          : Mesh-cell volume
C real*8  Ec            : Correlation energy
C real*8  Ex            : Exchange energy
C real*8  field(3)      : External electric field
C integer i             : General-purpose index
C integer ia            : Atom index
C integer io            : Orbital index
C integer ip            : Point index
C integer is            : Species index
C logical IsDiag        : Is supercell diagonal?
C integer ispin         : Spin index
C integer j             : General-purpose index
#ifndef BSC_CELLXC
C integer JDGdistr      : J.D.Gale's parallel distribution of mesh points
C integer myBox(2,3)    : My processor's mesh box
#endif /* ! BSC_CELLXC */
C integer nbcell        : Number of independent bulk lattice vectors
C integer npcc          : Partial core corrections? (0=no, 1=yes)
C integer nsd           : Number of diagonal spin values (1 or 2)
C integer ntpl          : Number of mesh Total Points in unit cell
C                         (including subpoints) locally
C real*4  rhoatm(ntpl)  : Harris electron density
C real*4  rhopcc(ntpl)  : Partial-core-correction density for xc
C real*4  DRho(ntpl)    : Selfconsistent electron density difference
C real*8  rhotot        : Total density at one point
C real*8  rmax          : Maximum orbital radius
C real*8  scell(3,3)    : Supercell vectors
C character shape*10    : Name of system shape
C real*4  Vaux(ntpl)    : Auxiliary potential array
C real*4  Vna(ntpl)     : Sum of neutral-atom potentials
C real*8  volume        : Unit cell volume
C real*4  Vscf(ntpl)    : Hartree potential of selfconsistent density
C real*8  x0(3)         : Center of molecule
C logical harrisfun     : Harris functional or Kohn-Sham?
C ----------------------------------------------------------------------

C     Local variables
      integer  :: i, ia, ip, ispin, nsd, np_vac, jx
#ifndef BSC_CELLXC
!     Interface to JMS's SiestaXC
      integer       :: myBox(2,3)
      integer, save :: JDGdistr=-1
      real(dp) :: stressXC(3,3)
#endif /* ! BSC_CELLXC */
      real(dp) :: b1Xb2(3), const, DEc, DEx, DStres(3,3), 
     &            Ec, Ex, rhotot, x0(3), volume, Vmax_vac, Vmean_vac
#ifdef BSC_CELLXC
!     Dummy arrays for cellxc call
      real(grid_p) :: aux3(3,1)
      real(grid_p) :: dummy_DVxcdn(1,1,1)
#endif /* BSC_CELLXC */

      logical :: use_rhog

      real(dp), external :: volcel, ddot

      external
     &  cross,
     &  dipole, 
     &  poison, 
     &  reord, rhooda, rhoofdsp, 
     &  timer, vmatsp, 
     &  readsp
#ifdef BSC_CELLXC
      external bsc_cellxc
#endif /* BSC_CELLXC */

C     Work arrays
      real(grid_p), pointer :: Vscf(:,:), Vscf_par(:,:),
     &                         DRho(:,:), DRho_par(:,:),
     &                         Vaux(:), Vaux_par(:), Chlocal(:),
     &                         Totchar(:), fsrc(:), fdst(:),
     &                         rhoatm_quad(:) => null(),
     &                         DRho_quad(:,:) => null()
      ! Temporary reciprocal spin quantity
      real(grid_p) :: rnsd
#ifdef BSC_CELLXC
      real(grid_p), pointer :: Vscf_gga(:,:), DRho_gga(:,:)
#endif /* BSC_CELLXC */
#ifdef MPI
      integer  :: MPIerror
      real(dp) :: sbuffer(7), rbuffer(7)
#endif

C    Linres internal variables
      real(grid_p), pointer :: dvnoscf(:,:),dvnoscf_par(:),
     &                         dRhoscf(:,:),dRhoscf_par(:,:),
     &                         VLR(:,:),VLR_par(:,:), 
     &                         dummy_Vna(:), 
     &                         dVXCpot(:,:),
     &                         dRhoscf0_par(:,:,:), dRhoatm_par(:,:)  
      integer :: iaend, iainit,atder
C


#ifdef DEBUG
      call write_debug( '    PRE DHSCF' )
#endif

      if ( spin%H /= size(Dscf,dim=2) ) then
         call die('Spin components is not equal to options.')
      end if

      if ( debug_dhscf ) then
         write(*,debug_fmt) Node,'DM',
     &        (sqrt(sum(Dscf(:,ispin)**2)),ispin=1,spin%H)
         write(*,debug_fmt) Node,'H',
     &        (sqrt(sum(Hmat(:,ispin)**2)),ispin=1,spin%H)
      end if
      
!-------------------------------------------------------------------- BEGIN
C ----------------------------------------------------------------------
C Start of SCF iteration part
C ----------------------------------------------------------------------

C ----------------------------------------------------------------------
C     At the end of DHSCF_INIT, and also at the end of any previous
!     call to dhscf, we were in the UNIFORM distribution
!     Rhoatm, Rhopcc and Vna were in UNIFORM dist, sequential form
!     The index array endpht was in the QUADRATIC distribution
C ----------------------------------------------------------------------

#ifdef _TRACE_
      call MPI_Barrier( MPI_Comm_World, MPIerror )
      call MPItrace_restart( )
#endif
      call timer( 'DHSCF', 1 )
      call timer( 'DHSCF3', 1 )

      nullify( Vscf, Vscf_par, DRho, DRho_par,
     &         Vaux, Vaux_par, Chlocal, Totchar )
#ifdef BSC_CELLXC
      nullify( Vscf_gga, DRho_gga)
#endif /* BSC_CELLXC */

      volume = volcel(cell)

C-------------------------------------------------------------------------
      
      if (analyze_charge_density_only) then
         !! Use the functionality in the first block
         !! of the routine to get charge files and partial charges
         call setup_analysis_options()
      endif
      
      if (filesOut%vna .ne. ' ') then
        ! Uniform dist, sequential form
#ifdef NCDF_4
        if ( write_cdf ) then
           call cdf_save_grid(trim(slabel)//'.nc','Vna',1,ntml,Vna)
        else
           call write_rho( filesOut%vna, 
     &          cell, ntm, nsm, ntpl, 1, Vna)
           call write_grid_netcdf( cell, ntm, 1, ntpl, Vna, "Vna")
        end if
#else
        call write_rho( filesOut%vna,
     &       cell, ntm, nsm, ntpl, 1, Vna )
        call write_grid_netcdf( cell, ntm, 1, ntpl, Vna, "Vna")
#endif
      endif

C     Allocate memory for DRho using the UNIFORM data distribution
      call re_alloc( DRho, 1, ntpl, 1, nspin, 'DRho','dhscf' )

C Find number of diagonal spin values
      nsd  = min( nspin, 2 )
      if ( nsd == 1 ) then
         rnsd = 1._grid_p
      else
         rnsd = 1._grid_p / nsd
      end if

C ----------------------------------------------------------------------
C Find SCF electron density at mesh points. Store it in array DRho
C ----------------------------------------------------------------------
!
!     The reading routine works in the uniform distribution, in
!     sequential form
!
      if (present(use_rhog_in)) then
            use_rhog = use_rhog_in
         else
            use_rhog = .false.
      endif
      if (use_rhog) then
         ! fourier transform back into drho
         call rhofft(cell, ntml(1), ntml(2), ntml(3), ntm, nspin,
     $                  DRho,rhog_in,BACKWARD)

      else if (read_charge_cdf) then
         call read_grid_netcdf(ntm(1:3),nspin,ntpl,DRho,"Rho")
         read_charge_cdf = .false.
      else if (read_deformation_charge_cdf) then
         call read_grid_netcdf(ntm(1:3),nspin,ntpl,DRho,"DeltaRho")
         ! Add to diagonal components only
         do ispin = 1,nsd
            do ip= 1, ntpl
C             rhoatm and Drho are in sequential mode
              DRho(ip,ispin) = DRho(ip,ispin) + rhoatm(ip) * rnsd
            enddo
         enddo
         read_deformation_charge_cdf = .false.
      else
        ! Set the QUADRATIC distribution and allocate memory for DRho_par
        ! since the construction of the density from the DM and orbital
        ! data needs that distribution
        if (nodes.gt.1) then
           call setMeshDistr( QUADRATIC, nsm, nsp,
     &                        nml, nmpl, ntml, ntpl )
        endif
        call re_alloc( DRho_par, 1, ntpl, 1, nspin,
     &                 'DRho_par','dhscf' )       
        if (Spiral) then
          call rhoofdsp( norb, nmpl, maxnd, numd, listdptr, listd,
     &                   nspin, Dscf, DRho_par, nuo, nuotot, iaorb,
     &                   iphorb, isa, qspiral )
        else
          call rhoofd( norb, nmpl, maxnd, numd, listdptr, listd,
     &                 nspin, Dscf, DRho_par, 
     &                 nuo, nuotot, iaorb, iphorb, isa )
        endif
        ! DRHO_par is here in QUADRATIC, clustered form

C       Set the UNIFORM distribution again and copy DRho to it
        if (nodes.gt.1) then
           call setMeshDistr( UNIFORM, nsm, nsp, nml, nmpl, ntml, ntpl )
        endif

        do ispin = 1, nspin
          fsrc => DRho_par(:,ispin)
          fdst => DRho(:,ispin)
         ! Sequential to be able to write it out
         ! if nodes==1, this call will just reorder
          call distMeshData( QUADRATIC, fsrc,
     &                       UNIFORM, fdst, TO_SEQUENTIAL )
        enddo
        call de_alloc( DRho_par, 'DRho_par','dhscf' )

        if ( debug_dhscf ) then
           write(*,debug_fmt) Node,'DRho',
     &          (sqrt(sum(DRho(:,ispin)**2)),ispin=1,nspin)
        end if

        if (save_initial_charge_density) then
           ! This section is to be deprecated in favor
           ! of "analyze_charge_density_only"
           ! (except for the special name for the .nc file)
#ifdef NCDF_4
          if ( write_cdf ) then
             call cdf_save_grid(trim(slabel)//'.nc','RhoInit',nspin,
     &            ntml,DRho)
          else
             call write_rho( "RHO_INIT", cell, ntm, nsm, ntpl,
     $            nspin, DRho)
             call write_grid_netcdf( cell, ntm, nspin, ntpl, DRho,
     $            "RhoInit")
          end if
#else
          call write_rho( "RHO_INIT", cell, ntm, nsm, ntpl,
     $         nspin, DRho)
          call write_grid_netcdf( cell, ntm, nspin, ntpl, DRho,
     $         "RhoInit")
#endif
          call timer('DHSCF3',2)
          call timer('DHSCF',2)
          call bye("STOP after producing RHO_INIT from input DM")
        endif

      endif

      if (mix_charge) then
         ! Save fourier transform of charge density
         call rhofft(cell, ntml(1), ntml(2), ntml(3), ntm, nspin,
     $        DRho,rhog,FORWARD)
      endif
!
!     Proper place to integrate Hirshfeld and Voronoi code,
!     since we have just computed rhoatm and Rho.

      if (want_partial_charges) then
        ! The endpht array is in the quadratic distribution, so
        ! we need to use it for this...
        if (nodes.gt.1) then
           call setMeshDistr( QUADRATIC, nsm, nsp,
     &                        nml, nmpl, ntml, ntpl )
        endif
        call re_alloc( DRho_quad, 1, ntpl, 1, nspin,
     &                 'DRho_quad','dhscf' )       
        call re_alloc( rhoatm_quad, 1, ntpl,
     &                 'rhoatm_quad','dhscf' )
        ! Redistribute grid-density
        do ispin = 1, nspin
          fsrc => DRho(:,ispin)
          fdst => DRho_quad(:,ispin)
         ! if nodes==1, this call will just reorder
          call distMeshData( UNIFORM, fsrc,
     &                       QUADRATIC, fdst, TO_CLUSTER )
        enddo
        call distMeshData( UNIFORM, rhoatm,
     &                     QUADRATIC, rhoatm_quad, TO_CLUSTER )

        call compute_partial_charges(DRho_quad,rhoatm_quad,
     .                  nspin, iaorb, iphorb, 
     .                  isa, nmpl,dvol)

        call de_alloc(rhoatm_quad,'rhoatm_quad','dhscf')
        call de_alloc(Drho_quad,'DRho_quad','dhscf')
        if (nodes.gt.1) then
           call setMeshDistr( UNIFORM, nsm, nsp,
     &                        nml, nmpl, ntml, ntpl )
        endif
      endif

C ----------------------------------------------------------------------
C Save electron density
C ----------------------------------------------------------------------
      if (filesOut%rho .ne. ' ') then
        !  DRho is already using a uniform, sequential form
#ifdef NCDF_4
        if ( write_cdf ) then
           call cdf_save_grid(trim(slabel)//'.nc','Rho',nspin,ntml,
     &          DRho)
        else
           call write_rho( filesOut%rho, cell, ntm, nsm, ntpl, nspin,
     &          DRho )
           call write_grid_netcdf( cell, ntm, nspin, ntpl, DRho, "Rho")
        end if
#else
        call write_rho( filesOut%rho, cell, ntm, nsm, ntpl, nspin,
     &       DRho )
        call write_grid_netcdf( cell, ntm, nspin, ntpl, DRho, "Rho")
#endif
      endif
C-----------------------------------------------------------------------
C Save TD-electron density after every given number of steps- Rafi, Jan 2016
C-----------------------------------------------------------------------
      call write_tdrho(filesOut)
      if (filesOut%tdrho .ne. ' ') then
        !  DRho is already using a uniform, sequential form
        call write_rho( filesOut%tdrho, cell, ntm, nsm, ntpl, nspin,
     &                  DRho )
        call write_grid_netcdf( cell, ntm, nspin, ntpl, DRho, "TDRho")
      endif
C ----------------------------------------------------------------------
C Save the diffuse ionic charge and/or the total (ionic+electronic) charge
C ----------------------------------------------------------------------
      if (filesOut%psch .ne. ' ' .or. filesOut%toch .ne. ' ') then
C       Find diffuse ionic charge on mesh
        ! Note that the *OnMesh routines, except PhiOnMesh,
        ! work with any distribution, thanks to the fact that
        ! the ipa, idop, and indexp arrays are distro-specific
        call re_alloc( Chlocal, 1, ntpl, 'Chlocal', 'dhscf' )
        call LocalChargeOnMesh( na, isa, ntpl, Chlocal, indxua )
        ! Chlocal comes out in clustered form, so we convert it
        call reord( Chlocal, Chlocal, nml, nsm, TO_SEQUENTIAL )

        if ( debug_dhscf ) then
           write(*,debug_fmt) Node,'Chlocal',sqrt(sum(Chlocal**2))
        end if

C       Save diffuse ionic charge 
        if (filesOut%psch .ne. ' ') then
#ifdef NCDF_4
          if ( write_cdf ) then
             call cdf_save_grid(trim(slabel)//'.nc','Chlocal',1,ntml,
     &            Chlocal)
          else
             call write_rho( filesOut%psch, cell, ntm, nsm, ntpl, 1,
     &            Chlocal)
             call write_grid_netcdf( cell, ntm, 1, ntpl, Chlocal,
     &            'Chlocal' )
          end if
#else
          call write_rho( filesOut%psch, cell, ntm, nsm, ntpl, 1,
     &         Chlocal)
          call write_grid_netcdf( cell, ntm, 1, ntpl,
     &         Chlocal, 'Chlocal' )
#endif
        endif

C       Save total (ionic+electronic) charge 
        if ( filesOut%toch .ne. ' ') then
           ! *****************
           ! **  IMPORTANT  **
           ! The Chlocal array is re-used to minimize memory
           ! usage. In the this small snippet the Chlocal
           ! array will contain the total charge, and
           ! if the logic should change, (i.e. should Chlocal
           ! be retained) is the Totchar needed to be re-instantiated.
           ! *****************

!$OMP parallel default(shared), private(ispin,ip)
           do ispin = 1, nsd
!$OMP do 
           do ip = 1, ntpl
              Chlocal(ip) = Chlocal(ip) + DRho(ip,ispin)
           end do
!$OMP end do 
           end do
!$OMP end parallel

          ! See note above
#ifdef NCDF_4
           if ( write_cdf ) then
              call cdf_save_grid(trim(slabel)//'.nc','RhoTot',1,ntml,
     &             Chlocal)
           else
              call write_rho(filesOut%toch,cell,ntm,nsm,ntpl,1,Chlocal)
              call write_grid_netcdf( cell, ntm, 1, ntpl, Chlocal, 
     &             "TotalCharge")
           end if
#else
           call write_rho( filesOut%toch, cell, ntm, nsm, ntpl, 1,
     &          Chlocal )
           call write_grid_netcdf( cell, ntm, 1, ntpl,
     &          Chlocal, "TotalCharge")
#endif
        end if 
        call de_alloc( Chlocal, 'Chlocal', 'dhscf' )
      endif

C ----------------------------------------------------------------------
C Save the total charge (model core + valence) for Bader analysis
C ----------------------------------------------------------------------
      
      ! The test for toch guarantees that we are in "analysis mode"
      if (filesOut%toch .ne. ' ' .and. savebader) then
        call save_bader_charge()
      endif

C Find difference between selfconsistent and atomic densities

      !Both DRho and rhoatm are using a UNIFORM, sequential form
!$OMP parallel do default(shared), private(ispin,ip),
!$OMP&collapse(2)
      do ispin = 1,nsd
        do ip = 1,ntpl
          DRho(ip,ispin) = DRho(ip,ispin) - rhoatm(ip) * rnsd
        enddo
      enddo
!$OMP end parallel do

C ----------------------------------------------------------------------
C Save electron density difference
C ----------------------------------------------------------------------
      if (filesOut%drho .ne. ' ') then
#ifdef NCDF_4
        if ( write_cdf ) then
           call cdf_save_grid(trim(slabel)//'.nc','RhoDelta',nspin,ntml,
     &          DRho)
        else
           call write_rho( filesOut%drho, cell, ntm, nsm, ntpl, nspin,
     &          DRho )
           call write_grid_netcdf( cell, ntm, nspin, ntpl, DRho, 
     &          "DeltaRho")
        end if
#else
        call write_rho( filesOut%drho, cell, ntm, nsm, ntpl, nspin,
     &       DRho )
        call write_grid_netcdf( cell, ntm, nspin, ntpl,
     &       DRho, "DeltaRho")
#endif
      endif

      if (present(charge_density_only)) then 
         if (charge_density_only) then
            call timer('DHSCF3',2)
            call timer('DHSCF',2)
            call de_alloc( DRho, 'DRho', 'dhscf' )
            RETURN
         endif
      endif

! End of analysis section
! Can exit now, if requested
      
        if (analyze_charge_density_only) then 
            call timer('DHSCF3',2)
            call timer('DHSCF',2)
           call bye("STOP after analyzing charge from input DM")
        endif

!-------------------------------------------------------------
C     Transform spin density into sum and difference

      ! TODO Check for NC/SO:
      ! Should we diagonalize locally at every point first?
      if (nsd .eq. 2) then
!$OMP parallel do default(shared), private(rhotot,ip)
        do ip = 1,ntpl
          rhotot     = DRho(ip,1) + DRho(ip,2)
          DRho(ip,2) = DRho(ip,2) - DRho(ip,1)
          DRho(ip,1) = rhotot
        enddo
!$OMP end parallel do
      endif

C Add a background charge to neutralize the net charge, to
C model doped systems. It only adds the charge at points
C where there are atoms (i.e., not in vacuum).          
C First, call with 'task=0' to add background charge    
      if (doping_active) call doping_uniform(cell,ntpl,0,
     $     DRho(:,1),rhoatm)
      
      ! Add doping in cell (from ChargeGeometries/Geometry.Charge)
      ! Note that this routine will return immediately if no dopant is present
      call charge_add('+',cell, ntpl, DRho(:,1) )

C ----------------------------------------------------------------------
C Calculate the dipole moment
C ----------------------------------------------------------------------
      dipol(1:3) = 0.0_dp
      if (shape .ne. 'bulk') then

C Find center of system
        x0(1:3) = 0.0_dp
        do ia = 1,nua
          x0(1:3) = x0(1:3) + xa(1:3,ia) / nua
        enddo

C Find dipole
        ! This routine is distribution-blind
        ! and will reduce over all processors.
        call dipole( cell, ntm, ntml(1), ntml(2), ntml(3), nsm,
     &               DRho, x0, dipol )

        ! Orthogonalize dipole to bulk directions
        if (shape .eq. 'chain') then
          const = ddot(3,dipol,1,bcell,1) / ddot(3,bcell,1,bcell,1)
          dipol(1:3) = dipol(1:3) - const * bcell(1:3,1)
        else if (shape .eq. 'slab') then
          call cross( bcell(1,1), bcell(1,2), b1Xb2 )
          const = ddot(3,dipol,1,b1Xb2,1) / ddot(3,b1Xb2,1,b1Xb2,1)
          dipol(1:3) = const * b1Xb2(1:3)
       end if

       if ( TSmode ) then
         if ( N_elec > 1 ) then
          ! Orthogonalize dipole to electrode transport directions
          do ia = 1 , N_Elec
            x0 = Elecs(ia)%cell(:,Elecs(ia)%t_dir)
            const = ddot(3,dipol,1,x0,1) / ddot(3,x0,1,x0,1)
            dipol(1:3) = dipol(1:3) - const * x0
          end do
        else if ( (shape == 'molecule') .or. (shape == 'chain') ) then
          ! Only allow dipole correction for chains and molecules
          ! along the semi-infinite direciton.
          ! Note this is *only* for 1-electrode setups
          ! Note that since the above removes the periodic directions
          ! this should not do anything for 'chain' with the same semi-infinite
          ! direction
          x0 = Elecs(1)%cell(:,Elecs(1)%t_dir)
          const = ddot(3,dipol,1,x0,1) / ddot(3,x0,1,x0,1)
          dipol(1:3) = const * x0
        end if
       end if
          
      endif

C ----------------------------------------------------------------------
C     Find Hartree potential of DRho = rhoscf-rhoatm. Store it in Vaux
C ----------------------------------------------------------------------
!     Solve Poisson's equation
      call re_alloc( Vaux, 1, ntpl, 'Vaux', 'dhscf' )

      call poison( cell, ntml(1), ntml(2), ntml(3), ntm, DRho,
     &             DUscf, Vaux, DStres, nsm )

      if ( debug_dhscf ) then
         write(*,debug_fmt) Node,'Poisson',sqrt(sum(Vaux(:)**2))
      end if

      ! Vscf is in the UNIFORM, sequential form, and only using
      ! the first spin index

      ! We require that even the SIESTA potential is "fixed"
      ! NOTE, this will only do something if
      !   TS.Hartree.Fix is set
      call ts_hartree_fix( ntm, ntml, Vaux)
      
C Add contribution to stress from electrostatic energy of rhoscf-rhoatm
      if (istr .eq. 1) then
        stressl(1:3,1:3) = stressl(1:3,1:3) + DStres(1:3,1:3)
      endif

C ----------------------------------------------------------------------
C     Find electrostatic (Hartree) energy of full SCF electron density
C     using the original data distribution 
C ----------------------------------------------------------------------
      Uatm = Uharrs
      Uscf = 0._dp
!$OMP parallel do default(shared), private(ip), 
!$OMP&reduction(+:Uscf)
      do ip = 1, ntpl
        Uscf = Uscf + Vaux(ip) * rhoatm(ip)
      enddo
!$OMP end parallel do
      Uscf = Uscf * dVol + Uatm + DUscf

C Call doping with 'task=1' to remove background charge added previously
C The extra charge thus only affects the Hartree energy and potential,
C but not the contribution to Enascf ( = \Int_{Vna*\rho})
      if (doping_active) call doping_uniform(cell,ntpl,1,
     $     DRho(:,1),rhoatm)

      ! Remove doping in cell (from ChargeGeometries/Geometry.Charge)
      ! Note that this routine will return immediately if no dopant is present
      call charge_add('-',cell, ntpl, DRho(:,1) )

C ----------------------------------------------------------------------
C Add neutral-atom potential to Vaux
C ----------------------------------------------------------------------
      Enaatm = 0.0_dp
      Enascf = 0.0_dp
!$OMP parallel do default(shared), private(ip),
!$OMP&reduction(+:Enaatm,Enascf)
      do ip = 1, ntpl
        Enaatm   = Enaatm + Vna(ip) * rhoatm(ip)
        Enascf   = Enascf + Vna(ip) * DRho(ip,1)
        Vaux(ip) = Vaux(ip) + Vna(ip)
      enddo
!$OMP end parallel do
      Enaatm = Enaatm * dVol
      Enascf = Enaatm + Enascf * dVol

C ----------------------------------------------------------------------
C Add potential from external electric field (if present)
C ----------------------------------------------------------------------
      if ( acting_efield ) then
         if ( dipole_correction ) then
            field = get_field_from_dipole(dipol,cell)
            if (Node == 0) then
	       write(6,'(a,3f12.4,a)')
     $              'Dipole moment in unit cell   =', dipol/Debye, ' D'
	       write(6,'(a,3f12.6,a)')
     $              'Electric field for dipole correction =',
     $              field/eV*Ang, ' eV/Ang/e'
            end if
            ! The dipole correction energy has an extra factor
            ! of one half because the field involved is internal.
            ! See the paper by Bengtsson DOI:10.1103/PhysRevB.59.12301
            ! Hence we compute this part separately
            DUext = -0.5_dp * ddot(3,field,1,dipol,1)
         else
            field = 0._dp
            DUext = 0._dp
         end if
         ! Add the external electric field
         field = field + user_specified_field
         ! This routine expects a sequential array,
         ! but it is distribution-blind
         call add_potential_from_field( field, cell, nua, isa, xa,
     &                                  ntm, nsm, Vaux )
         ! Add energy of external electric field
         DUext = DUext - ddot(3,user_specified_field,1,dipol,1)
      endif

! ---------------------------------------------------------------------
!     Transiesta:
!     add the potential corresponding to the (possible) voltage-drop.
!     note that ts_voltage is not sharing the reord wih efield since
!     we should not encounter both at the same time.
! ---------------------------------------------------------------------
      if (TSmode.and.IsVolt.and.TSrun) then
         ! This routine expects a sequential array,
         ! in whatever distribution
#ifdef TRANSIESTA_VOLTAGE_DEBUG
!$OMP parallel workshare default(shared)
         Vaux(:) = 0._dp
!$OMP end parallel workshare
#endif
         call ts_voltage(cell, ntm, ntml, Vaux)
#ifdef TRANSIESTA_VOLTAGE_DEBUG
         call write_grid_netcdf( cell, ntm, 1, ntpl, Vaux, 
     &        "TransiestaHartreePotential")
         call timer('ts_volt', 3)
         call bye('transiesta debug for Hartree potential')
#endif
      endif

! ----------------------------------------------------------------------
! Add potential from user defined geometries (if present)
! ----------------------------------------------------------------------
      call hartree_add( cell, ntpl, Vaux )

C ----------------------------------------------------------------------
C     Save electrostatic potential
C ----------------------------------------------------------------------
      if (filesOut%vh .ne. ' ') then
        ! Note that only the first spin component is used
#ifdef NCDF_4
        if ( write_cdf ) then
           call cdf_save_grid(trim(slabel)//'.nc','Vh',1,ntml,
     &          Vaux)
        else
           call write_rho( filesOut%vh, cell, ntm, nsm, ntpl, 1, Vaux )
           call write_grid_netcdf( cell, ntm, 1, ntpl, Vaux, 
     &          "ElectrostaticPotential")
        end if
#else
        call write_rho( filesOut%vh, cell, ntm, nsm, ntpl, 1, Vaux )
        call write_grid_netcdf( cell, ntm, 1, ntpl,
     &       Vaux, "ElectrostaticPotential")
#endif
      endif

C     Get back spin density from sum and difference
      ! TODO Check for NC/SO:
      ! Should we diagonalize locally at every point first?
      if (nsd .eq. 2) then
!$OMP parallel do default(shared), private(ip,rhotot)
        do ip = 1, ntpl
          rhotot     = DRho(ip,1)
          DRho(ip,1) = 0.5_dp * (rhotot - DRho(ip,2))
          DRho(ip,2) = 0.5_dp * (rhotot + DRho(ip,2))
        enddo
!$OMP end parallel do
      endif

C ----------------------------------------------------------------------
#ifndef BSC_CELLXC
C Set uniform distribution of mesh points and find my processor mesh box
C This is the interface to JM Soler's own cellxc routine, which sets
C up the right distribution internally.
C ----------------------------------------------------------------------

      call jms_setMeshDistr( distrID=JDGdistr, nMesh=ntm, 
     .                   nNodesX=1, nNodesY=ProcessorY, nBlock=nsm )
      call myMeshBox( ntm, JDGdistr, myBox )

C ----------------------------------------------------------------------
#endif /* ! BSC_CELLXC */
C Exchange-correlation energy
C ----------------------------------------------------------------------
      call re_alloc( Vscf, 1, ntpl, 1, nspin, 'Vscf', 'dhscf' )

      if (npcc .eq. 1) then
         
!$OMP parallel do default(shared), private(ip,ispin), collapse(2)
       do ispin = 1,nsd
          do ip= 1, ntpl
             DRho(ip,ispin) = DRho(ip,ispin) +
     &            (rhopcc(ip)+rhoatm(ip)) * rnsd
          enddo
       enddo
!$OMP end parallel do

      else

!$OMP parallel do default(shared), private(ip,ispin), collapse(2)
       do ispin = 1,nsd
          do ip= 1, ntpl
             DRho(ip,ispin) = DRho(ip,ispin) +
     &            rhoatm(ip) * rnsd
          enddo
       enddo
!$OMP end parallel do

      end if

      ! Write the electron density used by cellxc
      if (filesOut%rhoxc .ne. ' ') then
#ifdef NCDF_4
        if ( write_cdf ) then
           call cdf_save_grid(trim(slabel)//'.nc','RhoXC',nspin,ntml,
     &          DRho )
        else
           call write_rho( filesOut%rhoxc, cell, ntm, nsm, ntpl, nspin,
     &          DRho )
           call write_grid_netcdf( cell, ntm, nspin, ntpl, DRho, 
     &          "RhoXC")
        end if
#else
        call write_rho( filesOut%rhoxc, cell, ntm, nsm, ntpl, nspin,
     &       DRho )
        call write_grid_netcdf( cell, ntm, nspin, ntpl, DRho, 
     &       "RhoXC")
#endif
      endif

!     Everything now is in UNIFORM, sequential form

      call timer("CellXC",1)
#ifdef BSC_CELLXC
      if (nodes.gt.1) then
         call setMeshDistr( LINEAR, nsm, nsp, nml, nmpl, ntml, ntpl )
      endif

      call re_alloc( Vscf_gga, 1, ntpl, 1, nspin, 'Vscf_gga', 'dhscf' )
      call re_alloc( DRho_gga, 1, ntpl, 1, nspin, 'DRho_gga', 'dhscf' )

      ! Redistribute all spin densities
      do ispin = 1, nspin
        fsrc => DRho(:,ispin)
        fdst => DRho_gga(:,ispin)
        call distMeshData( UNIFORM, fsrc, LINEAR, fdst, KEEP )
      enddo

      call bsc_cellxc( 0, 0, cell, ntml, ntml, ntpl, 0, aux3, nspin,
     &             DRho_gga, Ex, Ec, DEx, DEc, Vscf_gga,
     &             dummy_DVxcdn, stressl )
#endif /* BSC_CELLXC */

#ifndef BSC_CELLXC

        call cellXC( 0, cell, ntm, myBox(1,1), myBox(2,1),
     .                           myBox(1,2), myBox(2,2),
     .                           myBox(1,3), myBox(2,3), nspin,
     .             DRho, Ex, Ec, DEx, DEc, stressXC, Vscf )

#else /* BSC_CELLXC */

      if (nodes.gt.1) then
         call setMeshDistr( UNIFORM, nsm, nsp, nml, nmpl, ntml, ntpl )
      endif

      ! Redistribute to the Vxc array
      do ispin = 1,nspin
        fsrc => Vscf_gga(:,ispin)
        fdst => Vscf(:,ispin)
        call distMeshData( LINEAR, fsrc, UNIFORM, fdst, KEEP )
      enddo
#endif /* BSC_CELLXC */

      if ( debug_dhscf ) then
         write(*,debug_fmt) Node,'XC',
     &        (sqrt(sum(Vscf(:,ispin)**2)),ispin=1,nspin)
      end if
      
      
#ifndef BSC_CELLXC
!     Vscf is still sequential after the call to JMS's cellxc
#else /* BSC_CELLXC */
      call de_alloc( DRho_gga, 'DRho_gga', 'dhscf' )
      call de_alloc( Vscf_gga, 'Vscf_gga', 'dhscf' )
#endif /* BSC_CELLXC */

      Exc =  Ex +  Ec
      Dxc = DEx + DEc

      call timer("CellXC",2)

!     Vscf contains only Vxc, and is UNIFORM and sequential
!     Now we add up the other contributions to it, at 
!     the same time that we get DRho back to true DeltaRho form
!$OMP parallel default(shared), private(ip,ispin)

      ! Hartree potential only has diagonal components
      do ispin = 1,nsd
        if (npcc .eq. 1) then
!$OMP do
           do ip = 1,ntpl
              DRho(ip,ispin) = DRho(ip,ispin) - 
     &             (rhoatm(ip)+rhopcc(ip)) * rnsd
              Vscf(ip,ispin) = Vscf(ip,ispin) + Vaux(ip)
           enddo
!$OMP end do
        else
!$OMP do
           do ip = 1,ntpl
              DRho(ip,ispin) = DRho(ip,ispin) - rhoatm(ip) * rnsd
              Vscf(ip,ispin) = Vscf(ip,ispin) + Vaux(ip)
           enddo
!$OMP end do
        endif
      enddo
!$OMP end parallel

#ifndef BSC_CELLXC
      stress = stress + stressXC
#endif /* ! BSC_CELLXC */

C ----------------------------------------------------------------------
C     Save total potential
C ----------------------------------------------------------------------
      if (filesOut%vt .ne. ' ') then
#ifdef NCDF_4
        if ( write_cdf ) then
           call cdf_save_grid(trim(slabel)//'.nc','Vt',nspin,ntml,
     &          Vscf)
        else 
           call write_rho( filesOut%vt, cell, ntm, nsm, ntpl, nspin,
     &          Vscf )
           call write_grid_netcdf( cell, ntm, nspin, ntpl, Vscf, 
     &          "TotalPotential")
        end if
#else
        call write_rho( filesOut%vt, cell, ntm, nsm, ntpl, nspin, Vscf)
        call write_grid_netcdf( cell, ntm, nspin, ntpl,
     &       Vscf, "TotalPotential")
#endif
      endif

C ----------------------------------------------------------------------
C Print vacuum level
C ----------------------------------------------------------------------

      if (filesOut%vt/=' ' .or. filesOut%vh/=' ') then
        forall(ispin=1:nsd) 
     .    DRho(:,ispin) = DRho(:,ispin) + rhoatm(:) * rnsd
        call vacuum_level( ntpl, nspin, DRho, Vscf, 
     .                     np_vac, Vmax_vac, Vmean_vac )
        forall(ispin=1:nsd) 
     .    DRho(:,ispin) = DRho(:,ispin) - rhoatm(:) * rnsd
        if (np_vac>0 .and. Node==0) print'(/,a,2f12.6,a)', 
     .    'dhscf: Vacuum level (max, mean) =', 
     .    Vmax_vac/eV, Vmean_vac/eV, ' eV'
      endif

C ----------------------------------------------------------------------
C     Find SCF contribution to hamiltonian matrix elements
C ----------------------------------------------------------------------
      if (iHmat .eq. 1) then
         if (nodes.gt.1) then
            call setMeshDistr( QUADRATIC, nsm, nsp,
     &                         nml, nmpl, ntml, ntpl )
         endif

!        This is a work array, to which we copy Vscf
         call re_alloc( Vscf_par, 1, ntpl, 1, nspin,
     &                 'Vscf_par', 'dhscf' )

         do ispin = 1, nspin
           fsrc => Vscf(:,ispin)
           fdst => Vscf_par(:,ispin)
           call distMeshData( UNIFORM, fsrc,
     &                        QUADRATIC, fdst, TO_CLUSTER )
         enddo

         if (Spiral) then
            call vmatsp( norb, nmpl, dvol, nspin, Vscf_par, maxnd,
     &           numd, listdptr, listd, Hmat, nuo,
     &           nuotot, iaorb, iphorb, isa, qspiral )
         else
            call vmat( norb, nmpl, dvol, spin, Vscf_par, maxnd,
     &           numd, listdptr, listd, Hmat, nuo,
     &           nuotot, iaorb, iphorb, isa )
         endif

         call de_alloc( Vscf_par,  'Vscf_par', 'dhscf' )
         if (nodes.gt.1) then
!          Everything back to UNIFORM, sequential
           call setMeshDistr( UNIFORM, nsm, nsp,
     &                         nml, nmpl, ntml, ntpl )
         endif
      endif

C---------------------------------------------------------------------
C     LINRES: Find Pulay terms
C---------------------------------------------------------------------
      if (present(dHmat0)) then
        if (nodes.gt.1) then 
          call setMeshDistr( QUADRATIC, nsm, nsp,
     &                      nml, nmpl, ntml, ntpl )
        endif

!    This is a work array, to which we copy Vscf
        call re_alloc( Vscf_par, 1, ntpl, 1, nspin,
     &                 'Vscf_par', 'dhscf' )

        do ispin = 1, nspin
          fsrc => Vscf(:,ispin)
          fdst => Vscf_par(:,ispin)
          call distMeshData( UNIFORM, fsrc,
     &                        QUADRATIC, fdst, TO_CLUSTER )
        enddo

        call dvmat(norb, nmpl, dvol, nspin, Vscf_par, maxnd,
     &             numd, listdptr, listd, dHmat0,nuo,
     &             nuotot,iaorb, iphorb, isa, ialr)

        call de_alloc( Vscf_par,  'Vscf_par', 'dhscf' )
        if (nodes.gt.1) then
!          Everything back to UNIFORM, sequential
           call setMeshDistr( UNIFORM, nsm, nsp,
     &                         nml, nmpl, ntml, ntpl )
        endif
      endif
C---------------------------------------------------------------------

#ifdef MPI
C     Global reduction of Uscf/DUscf/Uatm/Enaatm/Enascf
#ifndef BSC_CELLXC
!     Note that Exc and Dxc are already reduced in the new cellxc
#endif /* ! BSC_CELLXC */
      sbuffer(1) = Uscf
      sbuffer(2) = DUscf
      sbuffer(3) = Uatm
      sbuffer(4) = Enaatm
      sbuffer(5) = Enascf
#ifdef BSC_CELLXC
      sbuffer(6) = Exc
      sbuffer(7) = Dxc
#else
      sbuffer(6:7) = 0._dp
#endif /* BSC_CELLXC */
      call MPI_AllReduce( sbuffer, rbuffer, 7, MPI_double_precision,
     &                     MPI_Sum, MPI_Comm_World, MPIerror )
      Uscf   = rbuffer(1)
      DUscf  = rbuffer(2)
      Uatm   = rbuffer(3)
      Enaatm = rbuffer(4)
      Enascf = rbuffer(5)
#ifdef BSC_CELLXC
      Exc    = rbuffer(6)
      Dxc    = rbuffer(7)
#endif /* BSC_CELLXC */
#endif /* MPI */

C     Add contribution to stress from the derivative of the Jacobian of ---
C     r->r' (strained r) in the integral of Vna*(rhoscf-rhoatm)
      if (istr .eq. 1) then
        do i = 1,3
          stress(i,i) = stress(i,i) + ( Enascf - Enaatm ) / volume
        enddo
      endif

C     Stop time counter for SCF iteration part
      call timer( 'DHSCF3', 2 )

C ----------------------------------------------------------------------
C     End of SCF iteration part
C ----------------------------------------------------------------------

      if (ifa.eq.1 .or. istr.eq.1) then
C ----------------------------------------------------------------------
C Forces and stress : SCF contribution
C ----------------------------------------------------------------------
C       Start time counter for force calculation part
        call timer( 'DHSCF4', 1 )

C       Find contribution of partial-core-correction
        if (npcc .eq. 1) then
          call reord( rhopcc, rhopcc, nml, nsm, TO_CLUSTER )
          call reord( Vaux, Vaux, nml, nsm, TO_CLUSTER )
          ! The partial core calculation only acts on
          ! the diagonal spin-components (no need to
          ! redistribute un-used elements)
          do ispin = 1, nsd
            call reord( Vscf(:,ispin), Vscf(:,ispin),
     &                  nml, nsm, TO_CLUSTER )
          enddo

          call PartialCoreOnMesh( na, isa, ntpl, rhopcc, indxua, nsd,
     &                            dvol, volume, Vscf, Vaux, Fal,
     &                            stressl, ifa.ne.0, istr.ne.0 )
   
          call reord( rhopcc, rhopcc, nml, nsm, TO_SEQUENTIAL )
          call reord( Vaux, Vaux, nml, nsm, TO_SEQUENTIAL )
          ! ** see above
          do ispin = 1, nsd
            call reord( Vscf(:,ispin), Vscf(:,ispin),
     &                  nml, nsm, TO_SEQUENTIAL )
          enddo

          if ( debug_dhscf ) then
             write(*,debug_fmt) Node,'PartialCore',
     &            (sqrt(sum(Vscf(:,ispin)**2)),ispin=1,nsd)
          end if
        endif

        if ( harrisfun) then
!         Forhar deals internally with its own needs
!         for distribution changes
#ifndef BSC_CELLXC
          call forhar( ntpl, nspin, nml, ntml, ntm, npcc, cell, 
#else /* BSC_CELLXC */
          call forhar( ntpl, nspin, nml, ntml, npcc, cell, 
#endif /* BSC_CELLXC */
     &                 rhoatm, rhopcc, Vna, DRho, Vscf, Vaux )
!         Upon return, everything is UNIFORM, sequential form
        endif

C     Transform spin density into sum and difference
        ! TODO NC/SO
        ! Should we perform local diagonalization?
        if (nsd .eq. 2) then
!$OMP parallel do default(shared), private(rhotot,ip)
          do ip = 1,ntpl
            rhotot     = DRho(ip,1) + DRho(ip,2)
            DRho(ip,2) = DRho(ip,2) - DRho(ip,1)
            DRho(ip,1) = rhotot
          enddo
!$OMP end parallel do
        endif

C       Find contribution of neutral-atom potential
        call reord( Vna, Vna, nml, nsm, TO_CLUSTER )
        call reord( DRho, DRho, nml, nsm, TO_CLUSTER )
        call NeutralAtomOnMesh( na, isa, ntpl, Vna, indxua, dvol, 
     &                          volume, DRho, Fal, stressl, 
     &                          ifa.ne.0, istr.ne.0 )
        call reord( DRho, DRho, nml, nsm, TO_SEQUENTIAL )
        call reord( Vna, Vna, nml, nsm, TO_SEQUENTIAL )

        if (nodes.gt.1) then
           call setMeshDistr( QUADRATIC, nsm, nsp,
     &                       nml, nmpl, ntml, ntpl )
        endif

        call re_alloc( Vscf_par, 1, ntpl, 1, nspin,
     &                 'Vscf_par', 'dhscf' )
        do ispin = 1, nspin
          fsrc => Vscf(:,ispin)
          fdst => Vscf_par(:,ispin)
          call distMeshData( UNIFORM, fsrc,
     &                       QUADRATIC, fdst, TO_CLUSTER )
        enddo

!       Remember that Vaux contains everything except Vxc
        call re_alloc( Vaux_par, 1, ntpl, 'Vaux_par', 'dhscf' )
        call distMeshData( UNIFORM, Vaux,
     &                     QUADRATIC, Vaux_par, TO_CLUSTER )

        call dfscf( ifa, istr, na, norb, nuo, nuotot, nmpl, nspin,
     &              indxua, isa, iaorb, iphorb, 
     &              maxnd, numd, listdptr, listd, Dscf, datm,
     &              Vscf_par, Vaux_par, dvol, volume, Fal, stressl )

        call de_alloc( Vaux_par, 'Vaux_par', 'dhscf' )
        call de_alloc( Vscf_par, 'Vscf_par', 'dhscf' )

        if (nodes.gt.1) then
          call setMeshDistr( UNIFORM, nsm, nsp, nml, nmpl, ntml, ntpl )
        endif


C       Stop time counter for force calculation part
        call timer( 'DHSCF4', 2 )
C ----------------------------------------------------------------------
C       End of force and stress calculation
C ----------------------------------------------------------------------
      endif

C-----------------------------------------------------------------------
C        Linres block: elements of the dynamical matrix
C----------------------------------------------------------------------
C variables initialization
      if (present(dDscf)) then

        call timer( 'DHSCF-LR', 1 )
        nullify (dvnoscf)
        call re_alloc( dvnoscf, 1, ntpl , 1, 3,'dvnoscf', 'dhscf' )
        dvnoscf=0.0_dp
C------------------------------------------------------------------------

C Compute gradients of density and atom -------------------------------
C drhoscf0= sum_mu,nu rho_mu,nu * grad phi_mu * phi_nu + cc
C drhoatm = 2 sum_mu rho_mu * grad phi_mu * phi_mu 
C drhopcc = similar as dVna for gradient/derivative
        atder=0 !LR flag to compute derivatives(0) or gradients(1)
        if ((present(LRfirst)).and. (LRfirst)) then

          if (present(dHmat0)) atder=1  !flag to compute gradients

C Compute gradient/derivative of the partial core correction density
          if (npcc .eq. 1) then
            call re_alloc(drhopcc,1,ntpl,1,3,'drhopcc',
     &                                  'linres_matrices')
            call reord( rhopcc, rhopcc, nml, nsm, TO_CLUSTER )
            call reord( Vaux, Vaux, nml, nsm, TO_CLUSTER )
            drhopcc(1:ntpl,1:3)=0.0_dp
            do ispin = 1, nsd
              call reord( Vscf(:,ispin), Vscf(:,ispin),
     &                nml, nsm, TO_CLUSTER )
            enddo
            do jx=1,3
              call reord( drhopcc(:,jx), drhopcc(:,jx),
     &                   nml, nsm, TO_CLUSTER )
            enddo
            call PartialCoreOnMesh( na, isa, ntpl, rhopcc, indxua,
     &                            nsd, dvol, volume, Vscf, Vaux, Fal,
     &                            stressl, .FALSE.,.FALSE.,
     &                            drhopcc, atder, ialr, nua)
            call reord( rhopcc, rhopcc, nml, nsm, TO_SEQUENTIAL )
            call reord( Vaux, Vaux, nml, nsm, TO_SEQUENTIAL )
            do jx=1,3
              call reord( drhopcc(:,jx), drhopcc(:,jx),
     &                   nml, nsm, TO_SEQUENTIAL )
            enddo
            do ispin = 1, nsd
              call reord( Vscf(:,ispin), Vscf(:,ispin),
     &                  nml, nsm, TO_SEQUENTIAL )
            enddo
          endif ! npcc calculated gradient/derivative of rhopcc (drhopcc)

! Allocate variables which are stored in linres_matrices 
! (UNIFORM distribution)
          call re_alloc(dRhoscf0,1,ntpl,1,nspin,1,3,
     &                                'dRhoscf0','linres_matrices')
          call re_alloc(dRhoatm,1,ntpl,1,3,'dRhoatm','linres_matrices')
          dRhoscf0 = 0.0_dp 
          dRhoatm = 0.0_dp 

          if (nodes.gt.1) then
            call setMeshDistr( QUADRATIC, nsm, nsp,
     &                        nml, nmpl, ntml, ntpl )
          endif

! Par dummy variables in QUADRATIC distribution
          nullify(dRhoscf0_par)
          call re_alloc(dRhoscf0_par,1,ntpl,1,nspin,1,3,
     &                                'dRhoscf0_par','dhscf')
          nullify(dRhoatm_par)
          call re_alloc(dRhoatm_par,1,ntpl,1,3,
     &                                 'dRhoatm_par','dhscf')

          call dfscf2(norb, nmpl, dvol, nspin, maxnd,
     &              numd, listdptr, listd,
     &              nuo, nuotot, iaorb, iphorb, isa,
     &              atder, ialr, datm, Dscf,
     &              dRhoatm_par, dRhoscf0_par)
! Densities are returned in a CLUSTERED QUADRATIC form 

! Set the UNIFORM distribution again and copy densities
          if (nodes.gt.1) then
             call setMeshDistr( UNIFORM, nsm, nsp, nml, nmpl,
     &                         ntml, ntpl )
          endif

          ! TODO, these loops will create temporary copies 
          ! due to the strided pointer. I.e. we might
          ! as well do a copy.

          do ispin = 1, nspin
            do jx =1,3
              fsrc =>dRhoscf0_par(:,ispin,jx)
              fdst =>dRhoscf0(:,ispin,jx)
              call distMeshData( QUADRATIC, fsrc,
     &                    UNIFORM, fdst, TO_SEQUENTIAL )
            enddo
          enddo 
          do jx = 1,3
            fsrc => dRhoatm_par(:,jx)
            fdst => dRhoatm(:,jx)
            call distMeshData( QUADRATIC, fsrc,
     &                      UNIFORM, fdst, TO_SEQUENTIAL )
          enddo

          call de_alloc(dRhoscf0_par,'dRhoscf0_par','dhscf')
          call de_alloc(dRhoatm_par,'dRhoatm_par','dhscf')

        endif !Present LRfirst
C----------------------------------------------------------------------

C Compute the perturbed Neutral Atom potential and added into dvnoscf---- 
        nullify(dummy_Vna)
        call re_alloc( dummy_Vna, 1, ntpl, 'dummy_Vna', 'dhscf' )
        call reord( DRho, DRho, nml, nsm, TO_CLUSTER )
        do jx = 1 , 3
           call reord( dvnoscf(:,jx), dvnoscf(:,jx), nml,
     &          nsm, TO_CLUSTER ) 
        end do
        call NeutralAtomOnMesh( na, isa, ntpl, dummy_Vna, indxua,
     &                          dVol, volume, DRho, Fal, stressl,
     &                          ifa.ne.0, istr.ne.0, dvnoscf, atder,
     &                          ialr, nua )
        call reord( DRho, DRho, nml, nsm, TO_SEQUENTIAL )
        call de_alloc( dummy_Vna, 'dummy_Vna','dhscf' )
        do jx = 1 , 3
           call reord( dvnoscf(:,jx), dvnoscf(:,jx), nml,
     &          nsm, TO_SEQUENTIAL)
        end do
C Neutral atom potential is stored into dvnoscf
C--------------------------------------------------------------------------

C-----------------------------------------------------------------------------
C Compute the perturbed Vhartree and Vxc terms
C Derivative Space loop for perturbed Hamiltonian
C-----------------------------------------------------------------------------

!$OMP parallel default(shared), private(ip,ispin)

! For the perturbed Vxc we need total density dRho
        if ( npcc == 1 ) then
!$OMP do collapse(2)
           do ispin = 1, nsd
              do ip = 1, ntpl
                 DRho(ip,ispin) = DRho(ip,ispin) + 
     &                (rhoatm(ip)+rhopcc(ip)) * rnsd
              end do
           end do
!$OMP end do
        else
!$OMP do collapse(2)
           do ispin = 1, nsd
              do ip = 1, ntpl
                 DRho(ip,ispin) = DRho(ip,ispin) + rhoatm(ip) * rnsd
              end do
           end do
!$OMP end do
        end if

!$OMP end parallel

        ! Loop on each Cartesian coordinate
        do jx = 1, 3

C Initialization of VL matrix. It will store perturbed potential
C Initialization of the perturbed density dRhoscf
C Both are in the UNIFORM distribution
          nullify( VLR )
          call re_alloc( VLR, 1, ntpl, 1, nspin,'VLR','dhscf' )
          nullify(dRhoscf)
          call re_alloc(dRhoscf,1,ntpl,1,nspin,'dRhoscf','dhscf')

C Compute the perturbed density--------------------------------------
C           RHOOFD-------------
            ! dn(r)/dr = Sum(mu,nu) 2*Rho(mu,nu)* dPhi(mu)/dr * 
            !Phi(nu)+ Sum(mu,nu) dRho(mu,nu)/dr * Phi(mu) * Phi(nu)
            ! first term computed in dfscf2, last term in rhoofd
          
          if (nodes.gt.1) then ! set QUADRATIC distribution
            call setMeshDistr( QUADRATIC, nsm, nsp,
     &                        nml, nmpl, ntml, ntpl )
          endif
          call re_alloc( dRho_par, 1, ntpl, 1, nspin,
     &                'dRho_par','dhscf' )

          call rhoofd( norb, nmpl, maxnd, numd, listdptr, listd,
     &              nspin, dDscf(:,:,jx), dRho_par, nuo,
     &              nuotot, iaorb, iphorb, isa )

C       Set the UNIFORM distribution again and copy DRho_par to dRhoscf
          if (nodes.gt.1) then
            call setMeshDistr( UNIFORM, nsm, nsp, nml, nmpl, ntml, ntpl)
          endif

          do ispin = 1, nspin
            call distMeshData( QUADRATIC, dRho_par(:,ispin),
     &                     UNIFORM, dRhoscf(:,ispin), TO_SEQUENTIAL )
          enddo
          call de_alloc( dRho_par, 'dRho_par','dhscf' )

          if (present(dHmat0)) then !gradient case
            dRhoscf(1:ntpl,1:nspin) = dRhoscf0(1:ntpl,1:nspin,jx)
          else !perturbed case: add rhoofd to dfscf2
            dRhoscf(1:ntpl,1:nspin) = dRhoscf(1:ntpl,1:nspin) +
     &                     dRhoscf0(1:ntpl,1:nspin,jx)
          endif
C-----------------------------------------------------------------------------

C-----------------------------------------------------------------------------
C Hartree potential drho-drhoatm
C-----------------------------------------------------------------------------

!$OMP parallel default(shared), private(ip,ispin,Rhotot)

C         Find the difference between selfconsistent and 
C         atomic densities
!$OMP do collapse(2)
          do ispin = 1 , nsd
            do ip = 1, ntpl
              dRhoscf(ip,ispin) = dRhoscf(ip,ispin) -
     &                       dRhoatm(ip,jx) * rnsd
            end do
          end do
!$OMP end do ! implicit barrier KEEP

C         Transform spin density into sum and difference ---
          if ( nsd == 2 ) then
!$OMP do
             do ip = 1 , ntpl
                Rhotot = dRhoscf(ip,1) + dRhoscf(ip,2)
                dRhoscf(ip,2) = dRhoscf(ip,2) - dRhoscf(ip,1)
                dRhoscf(ip,1) = Rhotot
             end do
!$OMP end do
          end if

!$OMP end parallel
          
C Hartree potential of DRHO - DRHOATM stored into VLR
          call poison( cell, ntml(1), ntml(2), ntml(3), ntm,
     &           dRhoscf, DUscf, VLR, DStres, nsm )
C Hartree potential stored in VLR
C--------------------------------------------------------------------

!$OMP parallel default(shared), private(ip,ispin,Rhotot)

C Adding Vna+VH
!$OMP do
          do ip = 1, ntpl
             VLR(ip,1) = VLR(ip,1) + dvnoscf(ip,jx)
             dvnoscf(ip,jx) = VLR(ip,1)
          end do
!$OMP end do

C Now, dvnoscf=VLR=Vna+VH. dvnoscf will be constant whereas in  
C VLR, the Vxc term will be added 
C     Get back spin density from sum and difference
          if ( nsd == 2 ) then
!$OMP do
             do ip = 1, ntpl
                Rhotot = dRhoscf(ip,1)
                dRhoscf(ip,1) = 0.5_dp *( Rhotot - dRhoscf(ip,2) )
                dRhoscf(ip,2) = 0.5_dp *( Rhotot + dRhoscf(ip,2) )
                VLR(ip,2) = VLR(ip,1)
             end do
!$OMP end do
          end if

C------------------------------------------------------------------------
C Adding final XC potential + add drhopcc to density
C------------------------------------------------------------------------
          if ( npcc == 1 ) then
!$OMP do collapse(2)
             do ispin = 1,nsd
                do ip = 1,ntpl
                   dRhoscf(ip,ispin) = dRhoscf(ip,ispin) +
     &                  (drhopcc(ip,jx) + dRhoatm(ip,jx)) * rnsd
                end do
             end do
!$OMP end do
          else
!$OMP do collapse(2)
             do ispin = 1,nsd
                do ip = 1,ntpl
                   dRhoscf(ip,ispin) = dRhoscf(ip,ispin) +
     &                  dRhoatm(ip,jx) * rnsd
                end do
             end do
!$OMP end do
          end if
!$OMP end parallel
          
!In order to keep the Vscf as out argument from cellXC, we create another
! variable
          call timer("CellXC-LR",1)
          nullify(dVXCpot)
          call re_alloc( dVXCpot, 1, ntpl, 1, nspin, 'dVXCpot','dhscf' )
          dVXCpot(:,:)=0.0_dp
          call cellXC( 0, cell, ntm, myBox(1,1), myBox(2,1),
     .                           myBox(1,2), myBox(2,2),
     .                           myBox(1,3), myBox(2,3), nspin,
     .             DRho, Ex, Ec, DEx, DEc, stressXC, dVXCpot, dRhoscf )

!$OMP parallel do default(shared), private(ip,ispin), collapse(2)
          do ispin = 1, nsd
             do ip = 1, ntpl
                VLR(ip,ispin) = VLR(ip,ispin) + dVXCpot(ip,ispin)
             end do
          end do
!$OMP end parallel do

          call de_alloc( dVXCpot, 'dVXCpot','dhscf' )
          call timer("CellXC-LR",2)

C Here, VLR contains Vxc+Vh+Vna
C------------------------------------------------------------------------

C Calculate the total perturbed Hamiltonian
       if ((.not.(present(dHmat0))).and.(present(LRfirst))) then

! Set the QUADRATIC DISTRIBUTION for VLR_par
         if (nodes.gt.1) then
            call setMeshDistr( QUADRATIC, nsm, nsp,
     &                        nml, nmpl, ntml, ntpl )
         endif
         nullify(VLR_par)
         call re_alloc( VLR_par, 1, ntpl, 1, nspin,
     &                'VLR_par','dhscf' )

         do ispin = 1, nspin
           fsrc => VLR(:,ispin)
           fdst => VLR_par(:,ispin)
           call distMeshData( UNIFORM, fsrc,
     &                       QUADRATIC, fdst, TO_CLUSTER )
         enddo
         call timer('vmat-LR',1)

         call vmat( norb, nmpl, dvol, spin, VLR_par,
     &         maxnd, numd, listdptr, listd, dHmat(:,:,jx),
     &         nuo, nuotot, iaorb, iphorb, isa)

         call timer('vmat-LR',2)

         call de_alloc( VLR_par,  'VLR_par', 'dhscf' )

! Back to the UNIFORM distribution
         if (nodes.gt.1) then
           call setMeshDistr( UNIFORM, nsm, nsp, nml, nmpl, ntml, ntpl )
         endif

       endif

C------------------------------------------------------------------------
C Compute the dynamical matrix
C------------------------------------------------------------------------
        if (((present(dHmat0)).and.(present(LRfirst)).and.(iai.eq.ialr))
     . .or.((.not.(present(dHmat0))).and.(.not.(present(LRfirst)))))then
            if (present(dHmat0)) then
              iainit = iai
              iaend = iaf
            else
              iainit = ialr
              iaend = ialr
            endif

!$OMP parallel default(shared), private(ip,ispin,Rhotot)

! To dynamat, we need the density dRhoscf and the unperturbed charge dRho=dRho-rhoatm

            if ( npcc == 1 ) then
!$OMP do collapse(2)
               do ispin = 1, nsd
                  do ip = 1, ntpl
                     dRhoscf(ip,ispin) = dRhoscf(ip,ispin) -
     &                    (drhopcc(ip,jx) + dRhoatm(ip,jx)) * rnsd
                  end do
               end do
!$OMP end do nowait
!$OMP do collapse(2)
               do ispin = 1, nsd
                  do ip = 1, ntpl
                     DRho(ip,ispin) = DRho(ip,ispin) -
     &                    (rhoatm(ip) + rhopcc(ip)) * rnsd
                  end do
               end do
!$OMP end do ! implicit barrier required !

            else
!$OMP do collapse(2)
               do ispin = 1, nsd
                  do ip = 1, ntpl
                     dRhoscf(ip,ispin) = dRhoscf(ip,ispin) -
     &                    dRhoatm(ip,jx) * rnsd
                  end do
               end do
!$OMP end do nowait
!$OMP do collapse(2)
               do ispin = 1, nsd
                  do ip = 1, ntpl
                     DRho(ip,ispin) = DRho(ip,ispin)
     &                    - rhoatm(ip) * rnsd
                  end do
               end do
!$OMP end do ! implicit barrier required !
            end if

            if ( nsd == 2 ) then
!$OMP do
               do ip = 1, ntpl
                  Rhotot = dRhoscf(ip,1) + dRhoscf(ip,2)
                  dRhoscf(ip,2) = dRhoscf(ip,2) - dRhoscf(ip,1)
                  dRhoscf(ip,1) = Rhotot

                  Rhotot = dRho(ip,1) + dRho(ip,2)
                  dRho(ip,2) = dRho(ip,2) - dRho(ip,1)
                  dRho(ip,1) = Rhotot
               end do
!$OMP end do
            endif
!$OMP end parallel

! SEQUENTIAL to CLUSTER DRho

!            call reord(DRho, DRho, nml, nsm, TO_CLUSTER)
            
! Set everything to QUADRATIC distribution
            if (nodes.gt.1) then
               call setMeshDistr( QUADRATIC, nsm, nsp,
     &              nml, nmpl, ntml, ntpl )
            end if

            nullify(VLR_par,dvnoscf_par,DRhoscf_par)

            call re_alloc( Vscf_par, 1, ntpl, 1, nspin,
     &                    'Vscf_par', 'dhscf' )
            call re_alloc( VLR_par, 1, ntpl, 1, nspin,
     &                    'VLR_par','dhscf' )
            call re_alloc( DRho_par, 1, ntpl, 1, nspin, 
     &                    'DRho_par','dhscf' ) !it will store the perturbed density
            call re_alloc( DRhoscf_par, 1, ntpl, 1, nspin,
     &                    'DRhoscf_par','dhscf' ) !it will store the perturbed density
            call re_alloc(dvnoscf_par,1,ntpl,'dvnoscf_par','dhscf')

            call distMeshData( UNIFORM, dvnoscf(:,jx),
     &              QUADRATIC, dvnoscf_par(:), TO_CLUSTER )

            do ispin = 1, nspin
               call distMeshData( UNIFORM, Vscf(:,ispin),
     &              QUADRATIC, Vscf_par(:,ispin), TO_CLUSTER )
               call distMeshData( UNIFORM, VLR(:,ispin),
     &              QUADRATIC, VLR_par(:,ispin), TO_CLUSTER )
               call distMeshData( UNIFORM, dRho(:,ispin),
     &              QUADRATIC, DRho_par(:,ispin), TO_CLUSTER )
               call distMeshData( UNIFORM, dRhoscf(:,ispin),
     &              QUADRATIC, DRhoscf_par(:,ispin), TO_CLUSTER )
              
               do ia = iainit, iaend
                  call dynamat(norb, nuo, na, nua, nuotot, nspin, ispin,
     &                 jx, nmpl, ntpl, iphorb, iaorb, ia, numd,
     &                 listd, listdptr, indxua, maxnd, isa, atder,
     &                 Dscf(:,ispin), dDscf(:,ispin,jx), dvol,
     &                 DRho_par(:,ispin), DRhoscf_par(:,ispin), Datm, 
     &                 Vscf_par(:,ispin),VLR_par(:,:), 
     &                 dvnoscf_par(:), dynmat(:,:,jx,ia)) 
               end do
            end do
            call de_alloc( Vscf_par, 'Vscf_par', 'dhscf' )
            call de_alloc( VLR_par,   'VLR_par', 'dhscf' )
            call de_alloc( DRho_par, 'DRho_par', 'dhscf' )
            call de_alloc( DRhoscf_par, 'DRhoscf_par', 'dhscf' )
            call de_alloc( dvnoscf_par, 'dvnoscf_par','dhscf')


C     Set the UNIFORM distribution again
            if (nodes.gt.1) then
               call setMeshDistr( UNIFORM, nsm, nsp, nml,
     &              nmpl, ntml, ntpl )
            end if

! Back from CLUSTER to SEQUENTIAL dvnoscf and DRho
            call reord(dvnoscf(:,jx), dvnoscf(:,jx), nml, nsm,
     &           TO_SEQUENTIAL)
!            call reord(DRho, DRho, nml, nsm, TO_SEQUENTIAL)

! Back to the total density (for the perturbed Vxc calculation)
            if (nspin .eq. 2) then
!$OMP do
              do ip = 1, ntpl
                Rhotot = dRho(ip,1)
                dRho(ip,1) = 0.5_dp *( Rhotot - dRho(ip,2) )
                dRho(ip,2) = 0.5_dp *( Rhotot + dRho(ip,2) )
              end do
!$OMP end do
            endif

            if ( npcc == 1 ) then
!$OMP parallel do default(shared), private(ip,ispin), collapse(2)
               do ispin = 1,nsd
                  do ip = 1,ntpl
                     DRho(ip,ispin) = DRho(ip,ispin) +
     &                    (rhoatm(ip)+rhopcc(ip)) * rnsd
                  end do
               end do
!$OMP end parallel do
            else
!$OMP parallel do default(shared), private(ip,ispin), collapse(2)
               do ispin = 1,nsd
                  do ip = 1,ntpl
                     DRho(ip,ispin) = DRho(ip,ispin)
     &                    + rhoatm(ip) * rnsd
                  end do
               end do
!$OMP end parallel do
            end if
          end if !Dynamat if
C----------------------------------------------------------------------------
          call de_alloc(VLR,'VLR','dhscf')
          call de_alloc(dRhoscf,'dRhoscf','dhscf')

       end do ! Cartesian loop
       
C-----------------------------------------------------------------------------
C End Cartesian loop
C-----------------------------------------------------------------------------
        if (present(dHmat0)) then
! When the gradients are computed, these variables are deallocated
! For non-first calculation, these variables are stored in linres_matrices module
           call de_alloc(dRhoatm, 'dRhoatm', 'linres_matrices')
           call de_alloc(dRhoscf0, 'dRhoscf0', 'linres_matrices')
           call de_alloc(drhopcc, 'drhopcc', 'linres_matrices')
        endif

        call de_alloc(dvnoscf,'dvnoscf','dhscf')
        call timer( 'DHSCF-LR', 2 )
      endif 
C---------------------------------------------------------------------
C      End of Linres block
C----------------------------------------------------------------------


!     We are in the UNIFORM distribution
!     Rhoatm, Rhopcc and Vna are in UNIFORM dist, sequential form
!     The index array endpht is in the QUADRATIC distribution

C     Stop time counter
      call timer( 'DHSCF', 2 )

C ----------------------------------------------------------------------
C     Free locally allocated memory
C ----------------------------------------------------------------------
      call de_alloc( Vaux, 'Vaux', 'dhscf' )
      call de_alloc( Vscf, 'Vscf', 'dhscf' )
      call de_alloc( DRho, 'DRho', 'dhscf' )

#ifdef DEBUG
      call write_debug( '    POS DHSCF' )
#endif
!------------------------------------------------------------------------ END
      CONTAINS

      subroutine save_bader_charge()
      use meshsubs, only: ModelCoreChargeOnMesh
#ifdef NCDF_4
      use siesta_options, only: write_cdf
      use m_ncdf_siesta, only: cdf_save_grid
#endif
      ! Auxiliary routine to output the Bader Charge
      !
      real(grid_p), pointer :: BaderCharge(:) => null()

      call re_alloc( BaderCharge, 1, ntpl, name='BaderCharge',
     &                 routine='dhscf' )

      ! Find a model core charge by re-scaling the local charge
      call ModelCoreChargeOnMesh( na, isa, ntpl, BaderCharge, indxua )
      ! It comes out in clustered form, so we convert it
      call reord( BaderCharge, BaderCharge, nml, nsm, TO_SEQUENTIAL)
      do ispin = 1,nsd
         BaderCharge(1:ntpl) = BaderCharge(1:ntpl) + DRho(1:ntpl,ispin)
      enddo
#ifdef NCDF_4
      if ( write_cdf ) then
         call cdf_save_grid(trim(slabel)//'.nc','RhoBader',1,ntml,
     &        BaderCharge)
      else
         call write_rho( trim(slabel)// ".BADER", cell,
     $        ntm, nsm, ntpl, 1, BaderCharge )
         call write_grid_netcdf( cell, ntm, 1, ntpl,
     $        BaderCharge, "BaderCharge")
      end if
#else
      call write_rho( trim(slabel)// ".BADER", cell,
     $     ntm, nsm, ntpl, 1, BaderCharge )
      call write_grid_netcdf( cell, ntm, 1, ntpl,
     $     BaderCharge, "BaderCharge")
#endif

      call de_alloc( BaderCharge, name='BaderCharge' )
      end subroutine save_bader_charge

      subroutine setup_analysis_options()
      !! For the analyze-charge-density-only case,
      !! avoiding any diagonalization
      
      use siesta_options, only: hirshpop, voropop
      use siesta_options, only: saverho, savedrho, saverhoxc
      use siesta_options, only: savevh, savevt, savevna
      use siesta_options, only: savepsch, savetoch
      
      want_partial_charges = (hirshpop .or. voropop) 

      if (saverho)   filesOut%rho   = trim(slabel)//'.RHO' 
      if (savedrho)  filesOut%drho  = trim(slabel)//'.DRHO'
      if (saverhoxc) filesOut%rhoxc = trim(slabel)//'.RHOXC'
      if (savevh)    filesOut%vh    = trim(slabel)//'.VH'  
      if (savevt)    filesOut%vt    = trim(slabel)//'.VT'  
      if (savevna)   filesOut%vna   = trim(slabel)//'.VNA' 
      if (savepsch)  filesOut%psch  = trim(slabel)//'.IOCH'
      if (savetoch)  filesOut%toch  = trim(slabel)//'.TOCH'
      
      end subroutine setup_analysis_options
      
      end subroutine dhscf

      subroutine delk_wrapper(isigneikr, norb, maxnd,
     &                        numd, listdptr, listd,
     &                        nuo,  nuotot, iaorb, iphorb, isa )

      use m_delk,  only  : delk  ! The real workhorse, similar to vmat

      use moreMeshSubs,   only : setMeshDistr
      use moreMeshSubs,   only : UNIFORM, QUADRATIC
      use parallel,       only : Nodes
      use mesh,           only : nsm, nsp

!
!     This is a wrapper to call delk, using some of the module
!     variables of m_dhscf, but from outside dhscf itself.
!
      integer                  :: isigneikr, 
     &                            norb, nuo, nuotot, maxnd,
     &                            iaorb(*), iphorb(*), isa(*), 
     &                            numd(nuo),     
     &                            listdptr(nuo), listd(maxnd)

! The dhscf module variables used are:
! nmpl
! dvol
! nml
! nmpl
! ntml
! ntpl
!
! Some of them might be put somewhere else (mesh?) to allow some
! of the kitchen-sink functionality of dhscf to be made more modular.
! For example, this wrapper might live independently if enough mesh
! information is made available to it.

C ----------------------------------------------------------------------
C Calculate matrix elements of exp(i \vec{k} \cdot \vec{r})
C ----------------------------------------------------------------------
        if (isigneikr .eq. 1 .or. isigneikr .eq. -1) then

          if (nodes.gt.1) then
            call setMeshDistr( QUADRATIC, nsm, nsp,
     &                         nml, nmpl, ntml, ntpl )
          endif

          call delk( isigneikr, norb, nmpl, dvol, maxnd,
     &               numd, listdptr, listd,
     &               nuo,  nuotot, iaorb, iphorb, isa )

          if (nodes.gt.1) then
!           Everything back to UNIFORM, sequential
            call setMeshDistr( UNIFORM, nsm, nsp,
     &                         nml, nmpl, ntml, ntpl )
          endif

        endif
      end subroutine delk_wrapper

      end module m_dhscf
