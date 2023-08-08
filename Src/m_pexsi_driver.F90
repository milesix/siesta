
! Tangled code
module m_pexsi_solver

 use precision, only  : dp

 implicit none

 real(dp), save :: prevDmax  ! For communication of max diff in DM in scf loop
                             ! used in the heuristics for N_el tolerance
 public :: prevDmax
#ifdef SIESTA__PEXSI
 public :: pexsi_solver

CONTAINS

! This version uses separate distributions for Siesta 
! (setup_H et al) and PEXSI.
!
subroutine pexsi_solver(iscf, no_u, no_l, nspin_in,  &
     maxnh, numh, listhptr, listh, H, S, qtot, DM, EDM, &
     ef, Entropy, temp, delta_Efermi)

    use fdf
    use parallel, only   : SIESTA_worker, BlockSize
    use parallel, only   : SIESTA_Group, SIESTA_Comm
    use m_mpi_utils, only: globalize_sum, globalize_max
    use m_mpi_utils, only: broadcast
    use units,       only: Kelvin, eV
    use m_redist_spmatrix, only: aux_matrix, redistribute_spmatrix
    use class_Distribution
    use alloc,             only: re_alloc, de_alloc
    use siesta_options,    only: dDtol
#ifdef MPI
    use mpi_siesta
#endif
use f_ppexsi_interface
use iso_c_binding
use m_pexsi, only: plan, pexsi_initialize_scfloop

#ifdef TRACING_SOLVEONLY
      use extrae_module
#endif

implicit          none

integer, intent(in)  :: iscf  ! scf step number
integer, intent(in)  :: maxnh, no_u, no_l, nspin_in
integer, intent(in), target  :: listh(maxnh), numh(no_l), listhptr(no_l)
real(dp), intent(in), target :: H(maxnh,nspin_in), S(maxnh)
real(dp), intent(in) :: qtot
real(dp), intent(out), target:: DM(maxnh,nspin_in), EDM(maxnh,nspin_in)
real(dp), intent(out)        :: ef  ! Fermi energy
real(dp), intent(out)        :: Entropy ! Entropy/k, dimensionless
real(dp), intent(in)         :: temp   ! Electronic temperature
real(dp), intent(in)         :: delta_Efermi  ! Estimated shift in E_fermi
integer        :: ih, i
integer        :: info
logical        :: write_ok
!------------
external         :: timer
integer          :: World_Comm, mpirank, ierr
!
real(dp)  :: temperature, numElectronExact
integer   :: norbs, scf_step
real(dp)  :: delta_Ef
!
integer   :: nspin
integer :: PEXSI_Pole_Group, PEXSI_Spatial_Group, World_Group
integer, allocatable :: pexsi_pole_ranks_in_world(:)
integer  :: nworkers_SIESTA
integer, allocatable :: siesta_ranks_in_world(:)
integer :: PEXSI_Pole_Group_in_World
integer, allocatable :: PEXSI_Pole_ranks_in_World_Spin(:,:)
integer :: PEXSI_Pole_Comm, PEXSI_Spatial_Comm, PEXSI_Spin_Comm
integer :: numNodesTotal
integer :: npPerPole
logical  :: PEXSI_worker
!
type(Distribution)   :: dist1
type(Distribution), allocatable, target   :: dist2_spin(:)
type(Distribution), pointer :: dist2
integer  :: pbs, color, spatial_rank, spin_rank
type(aux_matrix), allocatable, target :: m1_spin(:)
type(aux_matrix) :: m2
type(aux_matrix), pointer :: m1
integer :: nrows, nnz, nnzLocal, numColLocal
integer, pointer, dimension(:) ::  colptrLocal=> null(), rowindLocal=>null()
!
real(dp), pointer, dimension(:) :: &
        HnzvalLocal=>null(), SnzvalLocal=>null(),  &
        DMnzvalLocal => null() , EDMnzvalLocal => null(), &
        FDMnzvalLocal => null()
!
integer :: ispin, pexsi_spin
real(dp), save :: PEXSINumElectronToleranceMin, &
            PEXSINumElectronToleranceMax, &
            PEXSINumElectronTolerance
logical, save  :: first_call = .true.
real(dp), save :: muMin0, muMax0, mu
real(dp)       :: on_the_fly_tolerance
type(f_ppexsi_options) :: options
!
integer                :: isSIdentity
integer                :: verbosity
integer                :: inertiaMaxIter
!
real(dp), save         :: energyWidthInertiaTolerance
real(dp)               :: pexsi_temperature, two_kT
real(dp), allocatable :: numElectronSpin(:), numElectronDrvMuSpin(:)
integer :: numTotalPEXSIIter
integer :: numTotalInertiaIter
real(dp) :: numElectronDrvMuPEXSI, numElectronPEXSI
real(dp) :: numElectron_out, numElectronDrvMu_out
real(dp) :: deltaMu
real(dp)       :: bs_energy, eBandH, free_bs_energy
real(dp)       :: buffer1
real(dp), save :: previous_pexsi_temperature
!  --------  for serial compilation
#ifndef MPI
    call die("PEXSI needs MPI")
#else
!
! Our global communicator is a duplicate of the passed communicator
!
call MPI_Comm_Dup(MPI_Comm_DFT, World_Comm, ierr)
call mpi_comm_rank( World_Comm, mpirank, ierr )

! NOTE:  fdf calls will assign values to the whole processor set,
! but some other variables will have to be re-broadcast (see examples
! below)

call timer("pexsi", 1)  

if (SIESTA_worker) then

   ! rename some intent(in) variables, which are only
   ! defined for the Siesta subset

   norbs = no_u
   nspin = nspin_in
   scf_step = iscf
   delta_Ef = delta_Efermi
   numElectronExact = qtot 

   ! Note that the energy units for the PEXSI interface are arbitrary, but
   ! H, the interval limits, and the temperature have to be in the
   ! same units. Siesta uses Ry units.

   temperature      = temp

   if (mpirank==0) write(6,"(a,f10.2)") &
               "Electronic temperature (K): ", temperature/Kelvin
endif
!
call broadcast(norbs,comm=World_Comm)
call broadcast(scf_step,comm=World_Comm)
call broadcast(delta_Ef,comm=World_Comm)
call broadcast(numElectronExact,World_Comm)
call broadcast(temperature,World_Comm)
call broadcast(nspin,World_Comm)
! Imported from modules, but set only in Siesta side
call broadcast(prevDmax,comm=World_Comm)
call broadcast(dDtol,comm=World_Comm)
call MPI_Comm_Group(World_Comm,World_Group, ierr)
call MPI_Group_Size(SIESTA_Group, nworkers_SIESTA, ierr)
allocate(siesta_ranks_in_world(nworkers_SIESTA))
call MPI_Group_translate_ranks( SIESTA_Group, nworkers_SIESTA, &
     (/ (i,i=0,nworkers_SIESTA-1) /), &
     World_Group, siesta_ranks_in_world, ierr )
call newDistribution(dist1,World_Comm,siesta_ranks_in_world, &
                     TYPE_BLOCK_CYCLIC,BlockSize,"bc dist")
deallocate(siesta_ranks_in_world)
call MPI_Barrier(World_Comm,ierr)

call mpi_comm_size( World_Comm, numNodesTotal, ierr )

npPerPole  = fdf_get("PEXSI.np-per-pole",4)
if (nspin*npPerPole > numNodesTotal) &
          call die("PEXSI.np-per-pole is too big for MPI size")

! "Row" communicator for independent PEXSI operations on each spin
! The name refers to "spatial" degrees of freedom.
color = mod(mpirank,nspin)    ! {0,1} for nspin = 2, or {0} for nspin = 1
call MPI_Comm_Split(World_Comm, color, mpirank, PEXSI_Spatial_Comm, ierr)

! "Column" communicator for spin reductions
color = mpirank/nspin       
call MPI_Comm_Split(World_Comm, color, mpirank, PEXSI_Spin_Comm, ierr)

! Group and Communicator for first-pole team of PEXSI workers
!
call MPI_Comm_Group(PEXSI_Spatial_Comm, PEXSI_Spatial_Group, Ierr)
call MPI_Group_incl(PEXSI_Spatial_Group, npPerPole,   &
     (/ (i,i=0,npPerPole-1) /),&
     PEXSI_Pole_Group, Ierr)
call MPI_Comm_create(PEXSI_Spatial_Comm, PEXSI_Pole_Group,&
     PEXSI_Pole_Comm, Ierr)


call mpi_comm_rank( PEXSI_Spatial_Comm, spatial_rank, ierr )
call mpi_comm_rank( PEXSI_Spin_Comm, spin_rank, ierr )
PEXSI_worker = (spatial_rank < npPerPole)   ! Could be spin up or spin down

! PEXSI blocksize 
pbs = norbs/npPerPole


allocate(pexsi_pole_ranks_in_world(npPerPole))
call MPI_Comm_Group(World_Comm, World_Group, Ierr)

call MPI_Group_translate_ranks( PEXSI_Pole_Group, npPerPole, &
     (/ (i,i=0,npPerPole-1) /), &
     World_Group, pexsi_pole_ranks_in_world, ierr )

! Include the actual world ranks in the distribution object

allocate (PEXSI_Pole_ranks_in_World_Spin(npPerPole,nspin))
call MPI_AllGather(pexsi_pole_ranks_in_world,npPerPole,MPI_integer,&
     PEXSI_Pole_Ranks_in_World_Spin(1,1),npPerPole, &
     MPI_integer,PEXSI_Spin_Comm,ierr)

! Create distributions known to all nodes
allocate(dist2_spin(nspin))
do ispin = 1, nspin
   call newDistribution(dist2_spin(ispin), World_Comm, &
                        PEXSI_Pole_Ranks_in_World_Spin(:,ispin),  &
                        TYPE_PEXSI, pbs, "px dist")
enddo
deallocate(pexsi_pole_ranks_in_world,PEXSI_Pole_Ranks_in_World_Spin)
call MPI_Barrier(World_Comm,ierr)

pexsi_spin = spin_rank+1  ! {1,2}
! This is done serially on the Siesta side, each time
! filling in the structures in one PEXSI set

allocate(m1_spin(nspin))
do ispin = 1, nspin

   m1 => m1_spin(ispin)

   if (SIESTA_worker) then
      m1%norbs = norbs
      m1%no_l  = no_l
      m1%nnzl  = sum(numH(1:no_l))
      m1%numcols => numH
      m1%cols    => listH
      allocate(m1%vals(2))
      m1%vals(1)%data => S(:)
      m1%vals(2)%data => H(:,ispin)

   endif  ! SIESTA_worker

   call timer("redist_orbs_fwd", 1)

   ! Note that we cannot simply wrap this in a pexsi_spin test, as
   ! there are Siesta nodes in both spin sets.
   ! We must discriminate the PEXSI workers by the distribution info
   dist2 => dist2_spin(ispin)
   call redistribute_spmatrix(norbs,m1,dist1,m2,dist2,World_Comm)
   
   call timer("redist_orbs_fwd", 2)

   if (PEXSI_worker .and. (pexsi_spin == ispin) ) then

      nrows = m2%norbs          ! or simply 'norbs'
      numColLocal = m2%no_l
      nnzLocal    = m2%nnzl
      call MPI_AllReduce(nnzLocal,nnz,1,MPI_integer,MPI_sum,PEXSI_Pole_Comm,ierr)

      call re_alloc(colptrLocal,1,numColLocal+1,"colptrLocal","pexsi_solver")
      colptrLocal(1) = 1
      do ih = 1,numColLocal
         colptrLocal(ih+1) = colptrLocal(ih) + m2%numcols(ih)
      enddo

      rowindLocal => m2%cols
      SnzvalLocal => m2%vals(1)%data
      HnzvalLocal => m2%vals(2)%data

      call re_alloc(DMnzvalLocal,1,nnzLocal,"DMnzvalLocal","pexsi_solver")
      call re_alloc(EDMnzvalLocal,1,nnzLocal,"EDMnzvalLocal","pexsi_solver")
      call re_alloc(FDMnzvalLocal,1,nnzLocal,"FDMnzvalLocal","pexsi_solver")

      call memory_all("after setting up H+S for PEXSI (PEXSI_workers)",PEXSI_Pole_Comm)

   endif ! PEXSI worker
enddo

! Make these available to all
! (Note that the values are those on process 0, which is in the spin=1 set
! In fact, they are only needed for calls to the interface, so the broadcast
! could be over PEXSI_Spatial_Comm only.

call MPI_Bcast(nrows,1,MPI_integer,0,World_Comm,ierr)
call MPI_Bcast(nnz,1,MPI_integer,0,World_Comm,ierr)

call memory_all("after setting up H+S for PEXSI",World_comm)

if (first_call) then

! Initial guess of chemical potential and containing interval
! When using inertia counts, this interval can be wide.
! Note that mu, muMin0 and muMax0 are saved variables

   muMin0           = fdf_get("PEXSI.mu-min",-1.0_dp,"Ry")
   muMax0           = fdf_get("PEXSI.mu-max", 0.0_dp,"Ry")
   mu               = fdf_get("PEXSI.mu",-0.60_dp,"Ry")

   PEXSINumElectronToleranceMin =  &
         fdf_get("PEXSI.num-electron-tolerance-lower-bound",0.01_dp)
   PEXSINumElectronToleranceMax =  &
         fdf_get("PEXSI.num-electron-tolerance-upper-bound",0.5_dp)

   ! start with largest tolerance
   ! (except if overriden by user)
   PEXSINumElectronTolerance = fdf_get("PEXSI.num-electron-tolerance",&
                                       PEXSINumElectronToleranceMax)
   first_call = .false.
else
!
!  Here we could also check whether we are in the first scf iteration
!  of a multi-geometry run...
!
   ! Use a moving tolerance, based on how far DM_out was to DM_in
   ! in the previous iteration (except if overriden by user)

   call get_on_the_fly_tolerance(prevDmax,on_the_fly_tolerance)

   ! Override if tolerance is explicitly specified in the fdf file
   PEXSINumElectronTolerance =  fdf_get("PEXSI.num-electron-tolerance",&
                                        on_the_fly_tolerance)
endif
!
call f_ppexsi_set_default_options( options )

options%muPEXSISafeGuard = fdf_get("PEXSI.mu-pexsi-safeguard",0.05_dp,"Ry")
options%maxPEXSIIter = fdf_get("PEXSI.mu-max-iter",10)

isSIdentity = 0

options%numPole  = fdf_get("PEXSI.num-poles",40)
options%gap      = fdf_get("PEXSI.gap",0.0_dp,"Ry")

! deltaE is in theory the spectrum width, but in practice can be much smaller
! than | E_max - mu |.  It is found that deltaE that is slightly bigger
! than  | E_min - mu | is usually good enough.

options%deltaE     = fdf_get("PEXSI.delta-E",3.0_dp,"Ry") ! Lin: 10 Ry...

! Ordering flag:
!   1: Use METIS
!   0: Use PARMETIS/PTSCOTCH
options%ordering = fdf_get("PEXSI.ordering",1)

! Number of processors for symbolic factorization
! Only relevant for PARMETIS/PT_SCOTCH
options%npSymbFact = fdf_get("PEXSI.np-symbfact",1)

verbosity = fdf_get("PEXSI.verbosity",1)
options%verbosity = verbosity

call get_current_temperature(pexsi_temperature)
options%temperature = pexsi_temperature
!
!  Set guard smearing for later use
!
two_kT = 2.0_dp * pexsi_temperature

options%numElectronPEXSITolerance = PEXSINumElectronTolerance

! Stop inertia count if mu has not changed much from iteration to iteration.

options%muInertiaTolerance =  &
     fdf_get("PEXSI.inertia-mu-tolerance",0.05_dp,"Ry")

! One-sided expansion of interval if correct mu falls outside it
options%muInertiaExpansion =  &
     fdf_get("PEXSI.lateral-expansion-inertia",3.0_dp*eV,"Ry") 


! Other user options

! Maximum number of iterations for computing the inertia                                          
! in a given scf step (until a proper bracket is obtained)                                        
inertiaMaxIter   = fdf_get("PEXSI.inertia-max-iter",5)

! Energy-width termination tolerance for inertia-counting
! By default, it is the same as the mu tolerance, to match
! the criterion in the simple DFT driver
energyWidthInertiaTolerance =  &
     fdf_get("PEXSI.inertia-energy-width-tolerance", &
             options%muInertiaTolerance,"Ry")
if (scf_step == 1) then
   call pexsi_initialize_scfloop(PEXSI_Spatial_Comm,npPerPole,spatial_rank,info)
   call check_info(info,"initialize_plan")
endif
call f_ppexsi_load_real_hs_matrix(&
      plan,&
      options,&
      nrows,&
      nnz,&
      nnzLocal,&
      numColLocal,&
      colptrLocal,&
      rowindLocal,&
      HnzvalLocal,&
      isSIdentity,&
      SnzvalLocal,&
      info) 

call check_info(info,"load_real_sym_hs_matrix")


if (scf_step == 1) then
   call timer( "pexsi-symb", 1 )
   ! This is only needed for inertia-counting
   call f_ppexsi_symbolic_factorize_real_symmetric_matrix(&
        plan, &
        options,&
        info)
   call check_info(info,"symbolic_factorize_real_symmetric_matrix")

   call f_ppexsi_symbolic_factorize_complex_symmetric_matrix(&
        plan, &
        options,&
        info)
   call check_info(info,"symbolic_factorize_complex_symmetric_matrix")
   call timer( "pexsi-symb", 2 )
endif
options%isSymbolicFactorize = 0 ! We do not need it anymore
!
numTotalInertiaIter = 0

call timer("pexsi-solver", 1)

if (need_inertia_counting()) then
   call get_bracket_for_inertia_count( )  
   call do_inertia_count(plan,muMin0,muMax0,mu)
else
   !  Maybe there is no need for bracket, just for mu estimation
   call get_bracket_for_solver()
endif

numTotalPEXSIIter = 0
allocate(numElectronSpin(nspin),numElectronDrvMuSpin(nspin))

solver_loop: do
   if (numTotalPEXSIIter > options%maxPEXSIIter ) then
      ! Do not die immediately, and trust further DM normalization
      ! to fix the number of electrons for unstable cases
      ! call die("too many PEXSI iterations")
      if (mpirank == 0) then
         write(6,"(a)") " ** Maximum number of PEXSI-solver iterations reached without convergence"
         write(6,"(a)") " .... an attempt will be made to normalize the density-matrix"
         write(6,"(a)") " This will succeed or not depending on the normalization tolerance (see manual)"
      endif
   endif
   if(mpirank == 0) then
      write (6,"(a,f9.4,a,f9.5)") 'Computing DM for mu(eV): ', mu/eV, &
           ' Tol: ', PEXSINumElectronTolerance
      write (6,"(a,f9.4,f9.5)") 'Monitoring bracket: ', muMin0/eV, muMax0/eV
   endif
   
   call f_ppexsi_calculate_fermi_operator_real(&
        plan,&
        options,&
        mu,&
        numElectronExact,&
        numElectron_out,&
        numElectronDrvMu_out,&
        info)
   
   call check_info(info,"fermi_operator")
      
   ! Per spin
   numElectron_out = numElectron_out / nspin
   numElectronDrvMu_out =  numElectronDrvMu_out / nspin
   
   ! Gather the results for both spins on all processors
   
   call MPI_AllGather(numElectron_out,1,MPI_Double_precision,&
         numElectronSpin,1,MPI_Double_precision,PEXSI_Spin_Comm,ierr)
   call MPI_AllGather(numElectronDrvMu_out,1,MPI_Double_precision,&
         numElectronDrvMuSpin,1,MPI_Double_precision,PEXSI_Spin_Comm,ierr)
   
   numElectronPEXSI = sum(numElectronSpin(1:nspin))
   numElectronDrvMuPEXSI = sum(numElectronDrvMuSpin(1:nspin))
   
   if (mpirank == 0) then
      write(6,"(a,f10.4)") "Fermi Operator. mu: ", mu/eV
      if (nspin == 2) then
         write(6,"(a,2f10.4,a,f10.4)") "Fermi Operator. numElectron(Up,Down): ", &
                        numElectronSpin(1:nspin), " Total: ", numElectronPEXSI
         write(6,"(a,2f10.4,a,f10.4)") "Fermi Operator. dN_e/dmu(Up,Down): ", &
                        numElectronDrvMuSpin(1:nspin)*eV, " Total: ", numElectronDrvMuPEXSI*eV
      else
         write(6,"(a,f10.4)") "Fermi Operator. numElectron: ", numElectronPEXSI
         write(6,"(a,f10.4)") "Fermi Operator. dN_e/dmu: ", numElectronDrvMuPEXSI*eV
      endif
   endif

   numTotalPEXSIIter =  numTotalPEXSIIter + 1

   if (abs(numElectronPEXSI-numElectronExact) > PEXSINumElectronTolerance) then
   
      deltaMu = - (numElectronPEXSI - numElectronExact) / numElectronDrvMuPEXSI
      ! The simple DFT driver uses the size of the jump to flag problems:
      ! if (abs(deltaMu) > options%muPEXSISafeGuard) then
   
      if ( ((mu + deltaMu) < muMin0) .or. ((mu + deltaMu) > muMax0) ) then
         if (mpirank ==0) then
            write(6,"(a,f9.3)") "DeltaMu: ", deltaMu, " is too big. Falling back to IC"
         endif
   
         ! We choose for now to expand the bracket to include the jumped-to point
   
         muMin0 = min(muMin0,mu+deltaMu)
         muMax0 = max(muMax0,mu+deltaMu)
   
         call do_inertia_count(plan,muMin0,muMax0,mu)
   
         cycle solver_loop
   
      endif
      mu = mu + deltaMu
      cycle solver_loop
   else
      ! Converged
      if (mpirank == 0) then
         write(6,"(a,f10.4)") "PEXSI solver converged. mu: ", mu/eV
      endif
      exit solver_loop
   endif
end do solver_loop

deallocate(numElectronSpin,numElectronDrvMuSpin)
call timer("pexsi-solver", 2)


if( PEXSI_worker ) then
   call f_ppexsi_retrieve_real_dft_matrix(&
        plan,&
        DMnzvalLocal,&
        EDMnzvalLocal,&
        FDMnzvalLocal,&
        eBandH,&          ! Will not be available
        bs_energy,&
        free_bs_energy,&
        info)
   call check_info(info,"retrieve_real_symmetric_dft_matrix")

   if (nspin == 2) then
      ! The matrices have to be divided by two...
      DMnzvalLocal(:) = 0.5_dp * DMnzvalLocal(:)
      EDMnzvalLocal(:) = 0.5_dp * EDMnzvalLocal(:)
      !!! Watch out with this. Internals??
      FDMnzvalLocal(:) = 0.5_dp * FDMnzvalLocal(:)
   endif

endif

if ((mpirank == 0) .and. (verbosity >= 1)) then
   write(6,"(a,i3)") " #&s Number of solver iterations: ", numTotalPEXSIIter
   write(6,"(a,i3)") " #&s Number of inertia iterations: ", numTotalInertiaIter
   write(6,"(a,f12.5,f12.4,2x,a2)") "mu, N_e:", mu/eV, &
        numElectronPEXSI, "&s"
endif

if (PEXSI_worker) then

   free_bs_energy = 0.0_dp
   bs_energy = 0.0_dp
   eBandH = 0.0_dp
   do i = 1,nnzLocal
      free_bs_energy = free_bs_energy + SnzvalLocal(i) * &
           ( FDMnzvalLocal(i) ) 
      bs_energy = bs_energy + SnzvalLocal(i) * &
           ( EDMnzvalLocal(i) )
      eBandH = eBandH + HnzvalLocal(i) * &
           ( DMnzvalLocal(i) )
   enddo

   ! First, reduce over the Pole_comm

   call globalize_sum( free_bs_energy, buffer1, comm=PEXSI_Pole_Comm )
   free_bs_energy = buffer1 
   call globalize_sum( bs_energy, buffer1, comm=PEXSI_Pole_Comm )
   bs_energy = buffer1
   call globalize_sum( eBandH, buffer1, comm=PEXSI_Pole_Comm )
   eBandH = buffer1

   ! Now, reduce over both spins

   call globalize_sum( free_bs_energy, buffer1, comm=PEXSI_Spin_Comm )
   ! Note that we need an extra term: mu*N for the free energy
   free_bs_energy = buffer1 + mu*numElectronPEXSI
   call globalize_sum( bs_energy, buffer1, comm=PEXSI_Spin_Comm )
   bs_energy = buffer1
   call globalize_sum( eBandH, buffer1, comm=PEXSI_Spin_Comm )
   eBandH = buffer1

   ! This output block will be executed only if World's root node is
   ! in one of the leading pole groups. This might not be so

   if ((mpirank == 0) .and. (verbosity >= 2)) then
      write(6, "(a,f12.4)") "#&s Tr(S*EDM) (eV) = ", bs_energy/eV
      write(6,"(a,f12.4)") "#&s Tr(H*DM) (eV) = ", eBandH/eV
      write(6,"(a,f12.4)") "#&s Tr(S*FDM) + mu*N (eV) = ", (free_bs_energy)/eV
   endif

   ef = mu
   ! Note that we use the S*EDM version of the band-structure energy
   ! to estimate the entropy, by comparing it to S*FDM This looks
   ! consistent, but note that the EDM is not used in Siesta to
   ! estimate the total energy, only the DM (via the density) (that
   ! is, the XC and Hartree correction terms to Ebs going into Etot
   ! are estimated using the DM)

   Entropy = - (free_bs_energy - bs_energy) / temp

   ! ef and Entropy are now known to the leading-pole processes
endif ! PEXSI_worker


do ispin = 1, nspin

   m1 => m1_spin(ispin)

   if (PEXSI_worker .and. (pexsi_spin == ispin)) then
      ! Prepare m2 to transfer

      call de_alloc(FDMnzvalLocal,"FDMnzvalLocal","pexsi_solver")
      call de_alloc(colPtrLocal,"colPtrLocal","pexsi_solver")

      call de_alloc(m2%vals(1)%data,"m2%vals(1)%data","pexsi_solver")
      call de_alloc(m2%vals(2)%data,"m2%vals(2)%data","pexsi_solver")

      m2%vals(1)%data => DMnzvalLocal(1:nnzLocal)
      m2%vals(2)%data => EDMnzvalLocal(1:nnzLocal)

   endif

   ! Prepare m1 to receive the results
   if (SIESTA_worker) then
      nullify(m1%vals(1)%data)    ! formerly pointing to S
      nullify(m1%vals(2)%data)    ! formerly pointing to H
      deallocate(m1%vals)
      nullify(m1%numcols)         ! formerly pointing to numH
      nullify(m1%cols)            ! formerly pointing to listH
   endif

   call timer("redist_orbs_bck", 1)
   dist2 => dist2_spin(ispin)
   call redistribute_spmatrix(norbs,m2,dist2,m1,dist1,World_Comm)
   call timer("redist_orbs_bck", 2)

   if (PEXSI_worker .and. (pexsi_spin == ispin)) then
      call de_alloc(DMnzvalLocal, "DMnzvalLocal", "pexsi_solver")
      call de_alloc(EDMnzvalLocal,"EDMnzvalLocal","pexsi_solver")

      nullify(m2%vals(1)%data)    ! formerly pointing to DM
      nullify(m2%vals(2)%data)    ! formerly pointing to EDM
      deallocate(m2%vals)
      ! allocated in the direct transfer
      call de_alloc(m2%numcols,"m2%numcols","pexsi_solver")
      call de_alloc(m2%cols,   "m2%cols",   "pexsi_solver")
   endif

   ! We assume that Siesta's root node also belongs to one of the
   ! leading-pole PEXSI communicators.
   ! Note that by wrapping the broadcasts for SIESTA_workers we
   ! do not make ef and Entropy known to the non-leading PEXSI processes.
   
   if (SIESTA_worker) then
      call broadcast(ef,comm=SIESTA_Comm)
      call broadcast(Entropy,comm=SIESTA_Comm)
      ! In future, m1%vals(1,2) could be pointing to DM and EDM,
      ! and the 'redistribute' routine check whether the vals arrays are
      ! associated, to use them instead of allocating them.
      DM(:,ispin)  = m1%vals(1)%data(:)    
      EDM(:,ispin) = m1%vals(2)%data(:)    
      ! Check no_l
      if (no_l /= m1%no_l) then
         call die("Mismatch in no_l")
      endif
      ! Check listH
      if (any(listH(:) /= m1%cols(:))) then
         call die("Mismatch in listH")
      endif

      call de_alloc(m1%vals(1)%data,"m1%vals(1)%data","pexsi_solver")
      call de_alloc(m1%vals(2)%data,"m1%vals(2)%data","pexsi_solver")
      deallocate(m1%vals)
      ! allocated in the direct transfer
      call de_alloc(m1%numcols,"m1%numcols","pexsi_solver") 
      call de_alloc(m1%cols,   "m1%cols",   "pexsi_solver")

   endif
enddo
call timer("pexsi", 2)



call delete(dist1)
do ispin = 1, nspin
  call delete(dist2_spin(ispin))
enddo
deallocate(dist2_spin)
deallocate(m1_spin)

call MPI_Comm_Free(PEXSI_Spatial_Comm, ierr)
call MPI_Comm_Free(PEXSI_Spin_Comm, ierr)
call MPI_Comm_Free(World_Comm, ierr)

! This communicator was created from a subgroup.
! As such, it is MPI_COMM_NULL for those processes
! not in the subgroup (non PEXSI_workers). Only
! defined communicators can be freed
if (PEXSI_worker) then
   call MPI_Comm_Free(PEXSI_Pole_Comm, ierr)
endif

call MPI_Group_Free(PEXSI_Spatial_Group, ierr)
call MPI_Group_Free(PEXSI_Pole_Group, ierr)
call MPI_Group_Free(World_Group, ierr)
#endif

CONTAINS
    
subroutine do_inertia_count(plan,muMin0,muMax0,muInertia)
  use iso_c_binding, only : c_intptr_t
  use m_convergence

  integer(c_intptr_t)      :: plan
  real(dp), intent(inout)  :: muMin0, muMax0
  real(dp), intent(out)    :: muInertia

  real(dp)            ::   muMinInertia, muMaxInertia
  integer             ::   nInertiaRounds

  real(dp), parameter ::   eps_inertia = 0.1_dp
  type(converger_t)   ::   conv_mu
  logical             ::   bad_lower_bound, bad_upper_bound
  logical             ::   interval_problem, inertiaFlag, inertiaFlagOld
  real(dp)            ::   inertia_electron_width
  real(dp)            ::   inertia_original_electron_width
  real(dp)            ::   inertia_energy_width
  real(dp)            ::   muLower, muUpper
  integer             ::   numMinICountShifts, numShift
  integer             ::   numNodesSpatial

  real(dp), allocatable :: shiftVec(:), inertiaVec(:)
  real(dp), allocatable :: inertiaVec_out(:)

  integer :: imin, imax

  ! ADDITIONAL VARIABLES ADDED BY JAH
  real(dp) :: sigma
  integer :: tau
  real(dp) :: hsShift
  real(dp), allocatable :: muCandidate(:)
  real(dp), allocatable :: muCandidateTemp(:)
  integer :: inertiaElec
  integer :: muCandidateCount 
  real(dp), allocatable :: neLower(:)
  real(dp), allocatable :: neUpper(:)

  real(dp) :: muMin
  real(dp) :: muMax
  real(dp) :: rangeOld
  real(dp) :: muRange
  real(dp) :: rangeL
  real(dp) :: muMinTemp
  real(dp) :: muMaxTemp
  ! integer :: inertiaFlag ! ??
  real(dp) :: rangeNew
  real(dp) :: updateRange
  real(dp) :: muEstimate
  integer :: k, l

  real(dp) :: sig_mod ! sigma modifier
  real(dp) :: param_1 
  real(dp) :: param_2
  
  ! Minimum number of sampling points for inertia counts                                            
  numMinICountShifts = fdf_get("PEXSI.inertia-min-num-shifts", 10)
  call mpi_comm_size( PEXSI_Spatial_Comm, numNodesSpatial, ierr )
  numShift = numNodesSpatial/npPerPole
  do
     if (numShift < numMinICountShifts) then
        numShift = numShift + numNodesSpatial/npPerPole
     else
        exit
     endif
  enddo

  allocate(shiftVec(numShift), inertiaVec(numShift))
  allocate(inertiaVec_out(numShift))
  

  nInertiaRounds = 0
!----------------------START OF REFINE INTERVAL--------------------------------!
  refine_interval: do
      numTotalInertiaIter = numTotalInertiaIter + 1

      if (mpirank == 0) then
          write(6,*) "Beginning refine_interval..."
          write(6,*) "numTotalInertiaIter: ", numTotalInertiaIter
      endif

      
      ! What does this do? Is this necessary
      options%muMin0 = muMin0
      options%muMax0 = muMax0

      if (mpirank == 0) then
          write(6,*) "muMin0: ", muMin0
          write(6,*) "muMax0: ", muMax0
      endif
      
      if (mpirank == 0) then
         write (6,"(a,2f9.4,a,a,i4)") 'Calling inertiaCount: [', &
              muMin0/eV, muMax0/eV, "] (eV)", &
              " Nshifts: ", numShift
      endif
      
      call timer("pexsi-inertia-ct", 1)

      hsShift = (muMax0 - muMin0) / (numShift - 1)
      
      if (mpirank == 0) then
          ! NEEDS TO BE DECLARED ABOVE
          ! real(dp) :: hsShift
          write(6,*) "numShift: ", numShift
          write(6,*) "hsShift: ", hsShift
      endif 

      do i = 1, numShift
         shiftVec(i) = muMin0 + (i-1) * hsShift
      enddo

      ! write(6,*) "Printing shiftVec: ", shiftVec

      ! do k=1,size(shiftVec),12
      ! write(6, '(12F8.3)') (shiftVec(l), l=k, min(i+11,size(shiftVec))
      ! end do
      
      call f_ppexsi_inertia_count_real_matrix(&
           plan,&
           options,&
           numShift,&
           shiftVec,&
           inertiaVec_out,&
           info) 
      
      call check_info(info,"inertia-count")
      
      ! All-Reduce to add the (two) spin inertias
      ! so that all processors have the complete inertiaVec(:)
      call MPI_AllReduce(inertiaVec_out, inertiaVec, &
           numShift, MPI_Double_precision, &
           MPI_Sum, PEXSI_Spin_Comm, ierr)
      
      ! If nspin=1, each state is doubly occupied
      ! 
      inertiaVec(:) = 2 * inertiaVec(:) / nspin
      
      call timer("pexsi-inertia-ct", 2)
      
      interval_problem = .false.
      
      ! Why are these called bad?
      if(mpirank == 0) then
         bad_lower_bound = (inertiaVec(1) > (numElectronExact - 0.1)) 
         bad_upper_bound = (inertiaVec(numShift) < (numElectronExact + 0.1)) 
      endif
      
      call broadcast(bad_lower_bound,comm=World_Comm)
      call broadcast(bad_upper_bound,comm=World_Comm)
      
      if (bad_lower_bound) then
         interval_problem =  .true.
         muMin0 = muMin0 - options%muInertiaExpansion ! 0.5
         if (mpirank==0) then
            write (6,"(a,2f12.4,a,2f10.4)") 'Wrong inertia-count interval (lower end). Counts: ', &
                 inertiaVec(1), inertiaVec(numShift), &
                 ' New interval: ', muMin0/eV, muMax0/eV
         endif
      endif
      if (bad_upper_bound) then
         interval_problem =  .true.
         muMax0 = muMax0 + options%muInertiaExpansion ! 0.5
         if (mpirank==0) then
            write (6,"(a,2f12.4,a,2f10.4)") 'Wrong inertia-count interval (upper end). Counts: ', &
                 inertiaVec(1), inertiaVec(numShift), &
                 ' New interval: ', muMin0/eV, muMax0/eV
         endif
      endif
      
      ! If there is an interval problem keep going...
      if (interval_problem) then
         ! do nothing more, stay in loop
         cycle refine_interval
      endif

      ! Set inertiaFlag to be true indicating that inertia should continue
      ! until the flag is set to be false elsewhere in the code.
      inertiaFlag = .true.
      
      ! COMMENTED 7/25 BY J.A.H.
      ! nInertiaRounds = nInertiaRounds + 1
      ! 
      ! imin = 1; imax = numShift
      ! 
      ! do i = 1, numShift
      !    if (inertiaVec(i) < numElectronExact - eps_inertia) then
      !       imin = max(imin,i)
      !    endif
      !    if (inertiaVec(i) > numElectronExact + eps_inertia) then
      !       imax = min(imax,i)
      !    endif
      ! enddo
      ! muMaxInertia = shiftVec(imax)
      ! muMinInertia = shiftVec(imin)
      
      ! Get the band edges by interpolation
      ! muLower = interpolate(inertiaVec,shiftVec,numElectronExact-eps_inertia)
      ! muUpper = interpolate(inertiaVec,shiftVec,numElectronExact+eps_inertia)
      
      ! muInertia = 0.5_dp * (muUpper + muLower)

      ! NEEDS TO BE DECLARED ABOVE
      ! real(dp) :: sigma
      ! integer :: tau
      sig_mod = fdf_get("PEXSI.sig_mod", 1._dp)
      sigma = 3._dp * temperature 
      sigma = sigma * sig_mod
      tau = ceiling(sigma / hsShift)

      ! real(dp), allocatable :: muCandidate(:)
      ! real(dp), allocatable :: muCandidateTemp(:)
      ! integer :: inertiaElec
      ! integer :: muCandidateCount

      if (allocated(muCandidate)) deallocate(muCandidate)
      allocate(muCandidate(numShift))
      
      muCandidate = 0._dp
      muCandidateCount = 0._dp
      do i = 1, numShift
          if (inertiaVec(i) == numElectronExact) then
              muCandidateCount = muCandidateCount + 1
              muCandidate(muCandidateCount) = shiftVec(i) 
              inertiaElec = inertiaVec(i)
          end if
      end do

      !!!!!!!!!!! BEGIN TRIM muCandidate !!!!!!!!!!!!!
      ! Allocate temporary array if muCandidate is allocated
      if (allocated(muCandidate)) then
          allocate(muCandidateTemp(muCandidateCount))
          muCandidateTemp = muCandidate(1:muCandidateCount)
          deallocate(muCandidate)
          allocate(muCandidate(muCandidateCount))
          muCandidate = muCandidateTemp
          deallocate(muCandidateTemp)
      end if
      ! call move_alloc(muCandidateTemp, muCandidate)
      !!!!!!!!!!! END TRIM muCandidate !!!!!!!!!!!!!

      ! real(dp), allocatable ::  NeLower
      ! real(dp), allocatable ::  NeUpper

      if (allocated(neLower)) deallocate(neLower)
      allocate(neLower(numShift))
      if (allocated(neUpper)) deallocate(neUpper)
      allocate(neUpper(numShift))

      ! Let's modify this coefficient and see how this changes things
      ! 0.3 + 0.7 for first NeLower
      ! 0.7 + 0.3 for NeUpper
      ! Create parameters to go in front of these (inertieVec(..) + (inertiaVec(..))

      param_1 = fdf_get("PEXSI.param_1", 1._dp)
      param_2 = fdf_get("PEXSI.param_2", 1._dp)
      ! Check the indexing here
      do i = tau, numShift - 1
          ! param_1 + param_2 should sum to 2
          neLower(i + 1) = 0.5_dp * (param_1 * inertiaVec(i + 1 - tau) + param_2 * inertiaVec(i + 1))
          neUpper(i + 1 - tau) = 0.5_dp * (param_1 * inertiaVec(i + 1 - tau) + param_2 * inertiaVec(i + 1))
      end do

      imin = 1
      imax = numShift

      ! Did I implement this if/else if correctly?
      do i = 2, numShift - 1
          if ((NeUpper(i) < numElectronExact) .and. (NeUpper(i + 1) >= numElectronExact)) then
          imin = i
          else if ((NeLower(i) > numElectronExact) .and. (NeLower(i - 1) <= numElectronExact)) then
          imax = i
          end if
      end do


    ! real(dp) :: muMin
    ! real(dp) :: muMax
    muMin = muMin0
    muMax = muMax0

    ! real(dp) :: rangeOld
    rangeOld = abs(muMin - muMax)

    ! real(dp) :: muRange
    ! real(dp) :: rangeL
    ! real(dp) :: muMinTemp
    ! real(dp) :: muMaxTemp
    ! logical :: inertiaFlag ?? Is this already declared above somewhere??
    ! real(dp) :: rangeNew
    ! real(dp) :: updateRange
    ! real(dp) :: muEstimate

    if (size(muCandidate) > 0) then
        if (muCandidate(1) == muCandidate(size(muCandidate))) then
            muRange = max(sigma, hsShift)
            ! rangeL = abs(muMin - muMax)
            muMin = max(muMin, muCandidate(1) - sigma)
            muMax = min(muMax, muCandidate(1) + sigma)
        else
            muMinTemp = max(muMin, muCandidate(size(muCandidate)) - sigma)
            muMaxTemp = min(muMax, muCandidate(1) + sigma)
                if (muMinTemp - muMaxTemp > 0) then
                    inertiaFlag = .false.
                else if (muMin - muMax < 0) then
                    muMin = muMinTemp
                    muMax = muMaxTemp
                end if
        end if

        rangeNew = abs(muMin - muMax)
        updateRange = rangeNew - rangeOld
        muEstimate =  (muMin + muMax) / 2
        muInertia = muEstimate

        ! write (6, '(A, F6.2, A)') "The mu has been found with time ", timeInertloopEnd - timeInertiaStaA, " [s]."

        if (mpirank == 0) then
            write (6, '(A, F6.2, A)') "mu is estimated to be ", muEstimate, "."
            write (6, '(A, F6.2, A)') "muMax ", muMax, "."
            write (6, '(A, F6.2, A)') "muMin ", muMin, "."
            write (6, '(A, F6.2, A)') "computed electrons ", inertiaElec, "."
            write (6, '(A, F6.2, A)') "shift equals ", hsShift, "."
        end if


        ! Print statements:
        ! Time it took to find a mu estimate
        ! muMax
        ! muMin
        ! Computed number of electrons (inertiaElec)
        ! The shift (hsShift)

        ! where is muInertiaTolerance defined? (do I need to declare it?)
        ! Should be options%muInertiaTolerance

        ! write (6, *) "For debugging:"
        ! write (6, *) "options%muInertiaTolerance"
        ! write (6, *) options%muInertiaTolerance

       ! write (6, *) "muInertiaTolerance"
       ! write (6, *) muInertiaTolerance

        if ((muMax - muMin < 2._dp * options%muInertiaTolerance) &
        .or. (updateRange > -0.01_dp * options%muInertiaTolerance)) then
            inertiaFlag = .false.
        end if

        muMinInertia = muMin
        muMaxInertia = muMax
        muCandidateCount = 0
        muCandidate = 0._dp
    end if 

    !! End of large if statement

    if (size(muCandidate) == 0) then
        muMin = shiftVec(imin)
        muMax = shiftVec(imax)
        muMinInertia = shiftVec(imin)
        muMaxInertia = shiftVec(imax)
        muEstimate =  (muMin + muMax) / 2
        muInertia = muEstimate
        ! What is matrix_size?? 
        if ((imin == 1 .and. imax == numshift) &
        .or. (NeLower(imin) == 0 .and. NeUpper(imax) == nrows) &
        .or. (muMax - muMin < 2.0 * options%muInertiaTolerance)) then
           ! muMinInertia = shiftVec(imin)
           ! muMaxInertia = shiftVec(imax)
            inertiaFlag = .false.

            ! WRITE (*, '(A, E12.5, A)') "The mu has been found with time ", timeInertloopEnd - timeInertiaStaA, " [s]."

            if (mpirank == 0) then
               write(6, '(A, E12.5, A)') "The mu is estimated to be ", muEstimate, "."
               write(6, '(A, E12.5, A)') "The muMAX ", muMaxInertia, "."
               write(6, '(A, E12.5, A)') "The muMin ", muMinInertia, "."
               write(6, '(A, I8, A)') "The computed electrons ", inertiaElec, "."
               write(6, '(A, E12.5, A)') "The shift equals ", hsShift, "."
            end if

            exit
        end if 

        if (mpirank == 0) then
            write (6, *) 
            write (6, *) "Inertia Counting"
            write (6, '(A, F6.2, A, F6.2, A)') "(muMin, muMax)   = (", muMin, ", ", muMax, ")"
            write (6, '(A, I6)') "numShift           = ", numShift
            write (6, *)
        end if

    end if

nInertiaRounds = nInertiaRounds + 1
      
      if (mpirank == 0) then
         write (6,"(a,i3,f10.4,i3,f10.4)") 'imin, muMinInertia, imax, muMaxInertia: ',&
                imin, muMinInertia/eV, imax, muMaxInertia/eV
         write (6,"(a,2f10.4,a,f10.4)") 'muLower, muUpper: ', muLower/eV, muUpper/eV, &
              ' mu estimated: ', muInertia/eV
      endif
        
        if (mpirank==0) then
      
           inertia_energy_width = (muMaxInertia - muMinInertia)
           ! Note that this is the width of the starting interval...
           inertia_original_electron_width = (inertiaVec(numShift) - inertiaVec(1))
           inertia_electron_width = (inertiaVec(imax) - inertiaVec(imin))
      
           write (6,"(a,2f9.4,a,f9.4,3(a,f10.3))") ' -- new bracket (eV): [', &
                muMinInertia/eV, muMaxInertia/eV,  &
                "] estimated mu: ", muInertia/eV, &
                " Nel width: ", inertia_electron_width, &
                " (Base: ", inertia_original_electron_width, &
                " ) E width: ", inertia_energy_width/eV
      
           if (nInertiaRounds == 1) then
              call reset(conv_mu)
              call set_tolerance(conv_mu,options%muInertiaTolerance)
           endif
           call add_value(conv_mu, muInertia)
      
      
           ! COMMENTED 7/27 BY JAH
          !logical :: inertiaFlagOld
           inertiaFlagOld = .true.
           !inertiaNumElectronTolerance is replace with options%numElectronPEXSITolerance 
           if (inertia_original_electron_width < options%numElectronPEXSITolerance) then
              write (6,"(a)") 'Leaving inertia loop: electron tolerance'
              inertiaFlagOld = .false.
           endif
           ! if (inertia_original_electron_width < inertiaNumElectronTolerance) then
           !    write (6,"(a)") 'Leaving inertia loop: electron tolerance'
           !    inertiaFlagOld = .false.
           ! endif
      !!$     if (inertia_electron_width < inertiaMinNumElectronTolerance) then
      !!$        write (6,"(a)") 'Leaving inertia loop: minimum workable electron tolerance'
      !!$        inertiaFlagOld = .false.
      !!$     endif
      
           ! This is the first clause of Lin's criterion
           ! in the simple DFT driver. The second clause is the same as the next one
           ! when the energy-width tolerance is the same as the mu tolerance (my default)
           ! We renamed inertiaFlag (from old code) to inertiaFlagOld
           ! I am not sure about the basis for this

           ! COMMENTED 7/27 BY JAH
           ! if (abs(muMaxInertia -numElectronExact) < eps_inertia ) then
           !    write (6,"(a,f12.6)") "Leaving inertia loop: |muMaxInertia-N_e|: ", &
           !         abs(muMaxInertia -numElectronExact)
           !    inertiaFlagOld = .false.
           ! endif
           ! if (inertia_energy_width < energyWidthInertiaTolerance) then
           !    write (6,"(a,f12.6)") 'Leaving inertia loop: energy width tolerance: ', &
           !     energyWidthInertiaTolerance/eV
           !    inertiaFlagOld = .false.
           ! endif
           ! if (is_converged(conv_mu)) then
           !    write (6,"(a,f12.6)") 'Leaving inertia loop: mu tolerance: ', options%muInertiaTolerance/eV
           !    inertiaFlagOld = .false.
           ! endif
           ! if (nInertiaRounds == inertiaMaxIter) then
           !    write (6,"(a)") 'Leaving inertia loop: too many rounds'
           !    inertiaFlagOld = .false.
           ! endif
        endif


        ! COMMENTED 7/27 BY JAH
        ! call broadcast(inertiaFlagOld,comm=World_Comm)
      
        ! if (inertiaFlagOld) then
        !    ! stay in loop
        !    ! These values should be guarded, in case the refined interval
        !    ! is too tight. Use 2*kT
        !    ! 
        !    muMin0 = muMinInertia - two_kT
        !    muMax0 = muMaxInertia + two_kT
        ! else
        !    exit refine_interval
        ! endif

        ! call broadcast(inertiaFlagOld,comm=World_Comm)

        call broadcast(inertiaFlag,comm=World_Comm)
      
        if (inertiaFlag .and. inertiaFlagOld) then
           ! stay in loop
           ! These values should be guarded, in case the refined interval
           ! is too tight. Use 2*kT
           ! 
           muMin0 = muMinInertia
           muMax0 = muMaxInertia
        else
           exit refine_interval
        endif
      
  enddo refine_interval

  ! Figure out what else needs to be deallocated
  if (allocated(shiftVec)) deallocate(shiftVec)
  if (allocated(inertiaVec)) deallocate(inertiaVec)
  if (allocated(muCandidate)) deallocate(muCandidate)
  if (allocated(muCandidateTemp)) deallocate(muCandidateTemp)
  if (allocated(neLower)) deallocate(neLower)
  if (allocated(neUpper)) deallocate(neUpper)


  ! deallocate(inertiaVec_out) hmm where is this used?
   
 end subroutine do_inertia_count
!----------------------END OF REFINE INTERVAL----------------------------------!

!
! This routine encodes the heuristics to compute the
! tolerance dynamically.
!
subroutine get_on_the_fly_tolerance(dDmax,tolerance)
real(dp), intent(in)  :: dDmax
real(dp), intent(out) :: tolerance

real(dp) :: tolerance_preconditioner
real(dp) :: tolerance_target_factor, tolerance_exp
real(dp), save :: previous_tolerance
logical :: new_algorithm

new_algorithm = fdf_get("PEXSI.dynamical-tolerance",.false.)
!
!
if (new_algorithm) then

!   By default, the tolerance goes to the (minimum) target 
!   at a level 5 times dDtol

   tolerance_target_factor = fdf_get("PEXSI.tolerance-target-factor",5.0_dp)

!
!  This can range in a (0.5,2.0) interval, approximately

   tolerance_preconditioner = fdf_get("PEXSI.tolerance-preconditioner",1.0_dp)

   if (scf_step > 1 ) then

      tolerance_exp = log10(dDmax/(tolerance_target_factor*dDtol))
      ! 
  !   range = log10(PEXSINumElectronToleranceMax/PEXSINumElectronToleranceMin)
      tolerance_exp = max(tolerance_exp,0.0_dp)*tolerance_preconditioner
      tolerance = PEXSINumElectronToleranceMin * 10.0_dp**tolerance_exp
      tolerance = min(tolerance,PEXSINumElectronToleranceMax)

      if (tolerance > previous_tolerance) then
         if (mpirank==0) write(6,"(a,f10.2)") &
              "Will not raise PEXSI solver tolerance to: ", &
              tolerance
         tolerance = previous_tolerance
      endif
      previous_tolerance = tolerance
   else
      ! No heuristics for now for first step
      ! Note that this should really change in MD or geometry optimization
      previous_tolerance = huge(1.0_dp)
      tolerance = PEXSINumElectronToleranceMax

   endif
else
   tolerance = Max(PEXSINumElectronToleranceMin, &
                              Min(dDmax*1.0, PEXSINumElectronToleranceMax))
endif

if (mpirank==0) write(6,"(a,f10.2)") &
     "Current PEXSI solver tolerance: ", tolerance

end subroutine get_on_the_fly_tolerance

!------------------------------------------------------------------
! This function will determine whether an initial inertia-counting
! stage is needed, based on user input and the level of convergence
!
! Variables used through host association for now:
!
!      scf_step
!      prevDmax, safe_dDmax_NoInertia
!
! Some logging output is done, so this function is not pure.

function need_inertia_counting() result(do_inertia_count)
logical :: do_inertia_count

real(dp) :: safe_dDmax_NoInertia
integer  :: isInertiaCount, numInertiaCounts

! Use inertia counts?
! The use of this input variable is deprecated. Warn the user
! only if there is a disagreement.

isInertiaCount = fdf_get("PEXSI.inertia-count",-1)
! For how many scf steps?
numInertiaCounts = fdf_get("PEXSI.inertia-counts",3)

if ((isInertiaCount == 0) .and. (numInertiaCounts > 0)) then 
   if (mpirank == 0) write(6,"(a,i4)")  &
        "Warning: Inertia-counts turned off by legacy parameter" // &
        " PEXSI.inertia-count"
   numInertiaCounts = 0
endif

safe_dDmax_NoInertia = fdf_get("PEXSI.safe-dDmax-no-inertia",0.05)

do_inertia_count = .false.

write_ok = ((mpirank == 0) .and. (verbosity >= 1))

if (numInertiaCounts > 0) then
  if (scf_step .le. numInertiaCounts) then
     if (write_ok) write(6,"(a,i4)")  &
      "&o Inertia-count step scf_step<numIC", scf_step
     do_inertia_count = .true.
  endif
else  if (numInertiaCounts < 0) then
   if (scf_step <= -numInertiaCounts) then
      if (write_ok) write(6,"(a,i4)") &
           "&o Inertia-count step scf_step<-numIC ", scf_step
      do_inertia_count = .true.
   else if (prevDmax > safe_dDmax_NoInertia) then
      if (write_ok) write(6,"(a,i4)") &
           "&o Inertia-count step as prevDmax > safe_Dmax ", scf_step
      do_inertia_count = .true.
   endif
endif

end function need_inertia_counting

!---------------------------------------------------------------
!  Chooses the proper interval for the call to the driver
!  in case we need a stage of inertia counting  
!
subroutine get_bracket_for_inertia_count()

 real(dp)       :: safe_width_ic
 real(dp)       :: safe_dDmax_Ef_inertia

 safe_width_ic = fdf_get("PEXSI.safe-width-ic-bracket",4.0_dp*eV,"Ry")
 safe_dDmax_Ef_Inertia = fdf_get("PEXSI.safe-dDmax-ef-inertia",0.1)

write_ok = ((mpirank == 0) .and. (verbosity >= 1))

 ! Proper bracketing                                                           
 if (scf_step > 1) then
   if (prevDmax < safe_dDmax_Ef_inertia) then
      ! Shift brackets using estimate of Ef change from previous iteration 
      !                                                                    
      if (write_ok) write(6,"(a)") &
         "&o Inertia-count bracket shifted by Delta_Ef"
      ! This might be risky, if the final interval of the previous iteration   
      ! is too narrow. We should broaden it by o(kT)                           
      ! The usefulness of delta_Ef is thus debatable...                        

      muMin0 = muMin0 + delta_Ef - two_kT
      muMax0 = muMax0 + delta_Ef + two_kT
   else
      ! Use a large enough interval around the previous estimation of   
      ! mu (the gap edges are not available...)  
      if (write_ok) write(6,"(a)") "&o Inertia-count safe bracket"
!      muMin0 = min(muLowerEdge - 0.5*safe_width_ic, muMinInertia)
      muMin0 = min(mu - 0.5*safe_width_ic, muMin0)
!      muMax0 = max(muUpperEdge + 0.5*safe_width_ic, muMaxInertia)
      muMax0 = max(mu + 0.5*safe_width_ic, muMax0)
   endif
 else
    if (write_ok) write(6,"(a)") &
       "&o Inertia-count called with iscf=1 parameters"
 endif
end subroutine get_bracket_for_inertia_count

subroutine get_bracket_for_solver()

    real(dp)       :: safe_width_solver
    real(dp)       :: safe_dDmax_Ef_solver

safe_width_solver = fdf_get("PEXSI.safe-width-solver-bracket",4.0_dp*eV,"Ry")
safe_dDmax_Ef_solver = fdf_get("PEXSI.safe-dDmax-ef-solver",0.05)

write_ok = ((mpirank == 0) .and. (verbosity >= 1))

! Do nothing for now
! No setting of  muMin0 and muMax0 yet, pending clarification of flow

  if (scf_step > 1) then
     if (prevDmax < safe_dDmax_Ef_solver) then
        if (write_ok) write(6,"(a)") "&o Solver mu shifted by delta_Ef"
        mu = mu + delta_Ef
     endif
     ! Always provide a safe bracket around mu, in case we need to fallback
     ! to executing a cycle of inertia-counting
     if (write_ok) write(6,"(a)") "&o Safe solver bracket around mu"
     muMin0 = mu - 0.5*safe_width_solver
     muMax0 = mu + 0.5*safe_width_solver
  else
     if (write_ok) write(6,"(a)") "&o Solver called with iscf=1 parameters"
     ! do nothing. Keep mu, muMin0 and muMax0 as they are inherited
  endif
end subroutine get_bracket_for_solver

!------------------------------------------------------
! If using the "annealing" feature, this routine computes
! the current temperature to use in the PEXSI solver
!
subroutine get_current_temperature(pexsi_temperature)
  real(dp), intent(out) :: pexsi_temperature

 logical  :: use_annealing
 real(dp) :: annealing_preconditioner, temp_factor
 real(dp) :: annealing_target_factor

 use_annealing = fdf_get("PEXSI.use-annealing",.false.)
 if (use_annealing) then
   annealing_preconditioner = fdf_get("PEXSI.annealing-preconditioner",1.0_dp)
!   By default, the temperature goes to the target at a level 10 times dDtol
   annealing_target_factor = fdf_get("PEXSI.annealing-target-factor",10.0_dp)

   if (scf_step > 1 ) then

      ! Examples for target_factor = 10, dDtol=0.0001:
      ! prevDmax=0.1, preconditioner=1, factor=3
      ! prevDmax=0.1, preconditioner=2, factor=5
      ! prevDmax=0.1, preconditioner=3, factor=7
      ! prevDmax<=0.001, factor = 1
      ! prevDmax<0.001, factor = 1

      temp_factor = (log10(prevDmax/(annealing_target_factor*dDtol)))
      temp_factor = 1 + annealing_preconditioner * max(0.0_dp, temp_factor)

      pexsi_temperature = temp_factor * temperature
      if (pexsi_temperature > previous_pexsi_temperature) then
         if (mpirank==0) write(6,"(a,f10.2)") &
              "Will not raise PEXSI temperature to: ", &
              pexsi_temperature/Kelvin
         pexsi_temperature = previous_pexsi_temperature
      endif
      previous_pexsi_temperature = pexsi_temperature
   else
      ! No heuristics for now for first step
      previous_pexsi_temperature = huge(1.0_dp)
      pexsi_temperature = temperature
      !   Keep in mind for the future if modifying T at the 1st step
      !      previous_pexsi_temperature = pexsi_temperature
   endif
else
      pexsi_temperature = temperature
endif
if (mpirank==0) write(6,"(a,f10.2)") &
     "Current PEXSI temperature (K): ", pexsi_temperature/Kelvin
end subroutine get_current_temperature

function interpolate(xx,yy,x) result(val)
!
! Interpolate linearly in the (monotonically increasing!) arrays xx and yy
!
integer, parameter :: dp = selected_real_kind(10,100)

real(dp), intent(in) :: xx(:), yy(:)
real(dp), intent(in) :: x
real(dp)             :: val

integer :: i, n

n = size(xx)
if (size(yy) /= n) call die("Mismatch in array sizes in interpolate")

if ( (x < xx(1)) .or. (x > xx(n))) then
   call die("Interpolate: x not in range")
endif

do i = 2, n
   if (x <= xx(i)) then
      val = yy(i-1) + (x-xx(i-1)) * (yy(i)-yy(i-1))/(xx(i)-xx(i-1))
      exit
   endif
enddo

end function interpolate

subroutine check_info(info,str)
integer, intent(in) :: info
character(len=*), intent(in) :: str

    if(mpirank == 0) then
       if (info /= 0) then
          write(6,*) trim(str) // " info : ", info
          call die("Error exit from " // trim(str) // " routine")
       endif
      call pxfflush(6)
    endif       
end subroutine check_info 

end subroutine pexsi_solver
#endif
end module m_pexsi_solver
! End of tangled code
