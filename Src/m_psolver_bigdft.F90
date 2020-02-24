! 
! Copyyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!  Interface module to call BigDFT solver instead of standar Siesta 
!  poison subroutine. It requires compilation of Siesta with additional
!  libraries.
!  More info at:
!  http://bigdft.org/Wiki/index.php?title=BigDFT_website 
! 
!  Created by Pablo Lopez-Tarifa and Daniel Sanchez Portal @CFM(2018)
!  Optimizations (mainly memory) done by Nick Papior (2020)
module m_psolver_bigdft

  use precision,      only : dp, grid_p
  use parallel,       only : Node, Nodes, ionode, ProcessorY
  use fdf 
  use units, only : Ang, pi
  implicit none

  logical, save         :: firsttime = .true.

  public :: poisson_bigdft

contains

  !> \brief General purpose of the subroutine poisson_bigdft
  !!
  !! This subroutine handles the full hartree potential calculation
  !! given by Psolver package of BigDFT. 
  !!
  !! The potential is calculated in two scenarios:
  !!
  !!    In the case of a serial run:
  !!
  !!    S.1 onetothreedarray, transforms array Rho from 1D to 3D
  !!         (compulsory by BigDFT functions).
  !!    S.2 pkernel_init and pkernel_set, set the coulomb kernel
  !!         object.
  !!    S.3 electrostatic_Solver, call the Hartree potential calculation.
  !!    S.4 threetoonedarray, transforms back V from 3D to 1D.
  !!
  !!    In the case of a parallel run:
  !!
  !!    P.1 call gather_rho, gather the partial Rho arrays carried by
  !!        each of the processors into a global Rho array.
  !!    P.2 local_grids, the information about the local grid extensions
  !!        is collected in this array and futher used in point P.4.
  !!    P.3 onetothreedarray, transforms array total Rho from 1D to 3D.
  !!    P.4 pkernel_init and pkernel_set, set the coulomb kernel
  !!         object.
  !!    P.5 electrostatic_Solver, call the Hartree potential calculation.
  !!    P.6 spread_potential, subroutine that spreads back the total
  !!        potential just calculated into the processors , using the
  !!        grid division stored in local_grids.
  !!
  !! 3. Memory deallocation.
  !!
  !> \brief Further improvements
  !! We are currently working on:
  !!
  !! 1. Cleaning the overal implementation, this is a draft version.
  !!
  !! 2. Parsing the stress tensor to Siesta dhscf.F.
  !!
  !! 3. Design of a fdf block to declare pkernel implicit solvent parameters.
  !!    As for this version options are handled by 'simple' input fdf keywords.
  !!
  ! 
  ! #################################################################### !
  subroutine poisson_bigdft(cell, lgrid, grid, ntm, RHO, V,  eh)
    ! #################################################################### !
    !                                                                      !
    !   This subroutine sets up the kernel object for a bigDFT Poisson     !
    !   solver calculation of the Hartree potential.                       !
    !                                                                      !
    ! #################################################################### !

    use units,          only : eV
    use siesta_geom,    only : shape, na_u, na_s, xa
    use alloc,          only : re_alloc, de_alloc
    use Poisson_Solver, only : coulomb_operator,  Electrostatic_Solver, &
        pkernel_set, pkernel_set_epsilon
    use siesta_options, only : bigdft_cavity
    use mesh,           only : nsm
    use PStypes,        only : pkernel_init, pkernel_free, pkernel_get_radius
    use yaml_output
    use f_utils
    use dictionaries, dict_set => set
    use mpi_siesta

    ! PSolver uses atlab to handle geometry etc.
    use at_domain

    real( dp), intent(in)      :: cell(3,3)        ! Unit-cell vectors.
    integer, intent(in)        :: lgrid(3)          ! Local mesh divisions.
    integer, intent(in)        :: grid(3)           ! Total mesh divisions.
    real( grid_p), intent(in)  :: RHO(:)            ! Input density.
    integer, intent(in)         :: ntm
    real( grid_p), intent(out)  :: V(:)              ! Output potential.
    real( dp), intent(out)      :: eh                 ! Electrostatic energy.
    real( dp), dimension(1,3), save :: hgrid        ! Uniform mesh spacings in the three directions.
    integer                     :: i, j, k, iatoms
    real(dp)                   :: einit
    type(dictionary), pointer  :: dict               ! Input parameters.
    integer                     :: isf_order, fd_order
    type(coulomb_operator)     :: pkernel
    logical                     :: pkernel_verbose 

    real(grid_p), pointer       :: Paux_1D(:)        ! 1D auxyliary array.
    real(grid_p), pointer       :: Paux_3D(:,:,:)  ! 3D auxyliary array.
    integer :: ngrid
    integer                     :: World_Comm, mpirank, ierr

    real(dp)                    :: radii(na_u) 
    character(len=2)          :: atm_label(na_u) 
    real(dp)                    :: delta, factor 
    integer, allocatable        :: local_grids(:,:)
    type(domain) :: dom

    call timer( 'BigDFT_solv', 1)

    ! f_util and dictionary initialisations:
    call f_lib_initialize()
    call dict_init(dict)

    nullify(Paux_3D)
    nullify(Paux_1D)

    ! Total number of mesh-points
    ngrid = product(grid)

    allocate(local_grids(3,Nodes))
    call re_alloc(Paux_1D, 1, ngrid, &
        'Paux_1D', 'poison_bigdft' )

    ! Calculation of cell angles and grid separation for psolver:
    if ( firsttime ) then
      do i = 1, 3
        hgrid(1,i)= sqrt(cell(1, i)**2 + cell(2,i)**2 + cell(3,i)**2) / grid(i)
      end do
    end if

    ! Setting of dictionary variables, done at every calling:  
    call set_bigdft_variables(dom, cell, dict, pkernel_verbose)

#ifdef MPI
!!! Collection of RHO: 

    ! RHO is collected to a global array.
    call gather_rho(RHO, grid, 2, lgrid(1) * lgrid(2) * lgrid(3), Paux_1D)

    ! Distribute (now all have the same grid)
    call MPI_Bcast(Paux_1D, ngrid, &
        MPI_double_precision, 0, MPI_COMM_WORLD, ierr) 

    call MPI_AllGather(lgrid, 3, MPI_Integer, &
        local_grids(1,1), 3, MPI_Integer, MPI_COMM_WORLD, ierr)

    ! PSolver interface requires a 3D array
    ! We will use pointer tricks
    call pointer_1D_3D(Paux_1D, grid(1), grid(2), grid(3), Paux_3D)
#else
    call pointer_1D_3D(Rho, grid(1), grid(2), grid(3), Paux_3D)
#endif

    pkernel = pkernel_init(Node, Nodes, dict, dom, grid, hgrid)

    call pkernel_set( pkernel, verbose=pkernel_verbose)

    ! Implicit solvent?
    if ( bigdft_cavity ) then

      ! Set a radius for each atom in Angstroem (UFF, Bondi, Pauling, etc ...)
      call set_label(atm_label)
      do iatoms = 1 , na_u
        radii(iatoms) = pkernel_get_radius(pkernel, atname=atm_label(iatoms))
      end do

      ! Multiply for a constant prefactor (see the Soft-sphere paper JCTC 2017)
      factor = pkernel%cavity%fact_rigid
      do iatoms = 1 , na_u
        radii(iatoms) = factor * radii(iatoms) * Ang
      end do
      call pkernel_set_epsilon(pkernel, nat=na_u, rxyz=xa, radii=radii)
    end if

    ! Free calling dictionary
    call dict_free(dict)

    call Electrostatic_Solver(pkernel, Paux_3D, ehartree=eh)

    ! Energy must be devided by the total number of nodes because the way 
    ! Siesta computes the Hartree energy. Indeed, BigDFT provides the value
    ! alredy distributed. 
    eh = 2.0_dp * eh / Nodes
#ifdef MPI
    ! Now we spread back the potential on different nodes. 
    call timer('spread_pot', 1)
    Paux_1D = 2.0_dp * Paux_1D - sum(2.0_dp * Paux_1D) / ngrid
    call spread_potential(Paux_1D, grid, 2, ngrid, local_grids, V) 
    call timer('spread_pot', 2)
#else
    V(:) = 2.0_dp * Paux_1D - sum(2.0_dp * Paux_1D) / ngrid
#endif

    ! Ending calculation
    call pkernel_free( pkernel)
    call f_lib_finalize_noreport() 
    call de_alloc( Paux_1D, 'Paux', 'poisson_bigdft' )
    deallocate(local_grids)

    firsttime = .false.

    call timer( 'BigDFT_solv', 2)
    
  end subroutine poisson_bigdft

  subroutine set_bigdft_variables(dom, cell, dict, pkernel_verbose)

    use dictionaries, dict_set => set
    use at_domain
    use futile
    use siesta_geom,    only : shape
    use siesta_options, only : bigdft_isf_order, bigdft_fd_order,    &
        bigdft_verbose, bigdft_delta,         &
        bigdft_fact_rigid, bigdft_cavity_type,&
        bigdft_cavity, bigdft_radii_type,     &
        bigdft_gps_algorithm, bigdft_epsilon, &
        bigdft_gammaS, bigdft_alphaS,         &
        bigdft_betaV, bigdft_atomic_radii

    ! #################################################################### !
    ! Subroutine to set-up BigDFT Poisson solver variables and boundary    !
    ! condition.                                                           !
    ! #################################################################### !
    type(domain), intent(inout) :: dom
    real(dp), intent(in) :: cell(3,3)
    type(dictionary), pointer, intent(out) :: dict         ! Input parameters.
    logical, intent(out)  :: pkernel_verbose 
    ! #################################################################### !

    ! Initialize domain
    select case ( shape )
    case ( 'molecule' )
      dom = domain_new(ATOMIC_UNITS, [FREE_BC, FREE_BC, FREE_BC], abc=cell)
      if ( firsttime .and. IONode ) then
        write(6,'(/,a,/,a/,a,/,a,/,a,/,a)')               &
            'bigdft_solver: initialising kernel for a molecular system',&
            '#================= WARNING bigdft_solver ===============#',&
            '                                                         ',&
            ' Molecular system must be place at the center of the box,',&
            ' otherwise box-edge effects might lead to wrong Hartree  ',&
            ' energy values.                                          ',&
            '#=======================================================#'
      end if
    case ( 'slab' )
      dom = domain_new(ATOMIC_UNITS, [PERIODIC_BC, FREE_BC, PERIODIC_BC], abc=cell)
      if ( firsttime .and. IONode ) then
        write(6,'(/,a,/,a/,a,/,a,/,a,/,a)')                 &
            'bigdft_solver: initialising kernel for a slab system',     &
            '#================= WARNING bigdft_solver ===============#',&
            ' Empty direction of the periodic system has to be set    ',&
            ' along the Y-axis.                                       ',&
            ' Additionally, simulation box must be large enough to    ',&
            ' to avoid box edges effects in the Hartree potential.    ',&
            '#=======================================================#'
      end if
    case ( 'bulk' )
      dom = domain_new(ATOMIC_UNITS, [PERIODIC_BC, PERIODIC_BC, PERIODIC_BC], abc=cell)
      if ( firsttime .and. IONode ) then
        write(6,"(a)")                                     &
            'bigdft_solver: initialising kernel for a periodic system'
      end if
    case ( 'wire' )
      dom = domain_new(ATOMIC_UNITS, [FREE_BC, FREE_BC, PERIODIC_BC], abc=cell)
      if ( firsttime .and. IONode ) then
        write(6,'(/,a,/,a/,a,/,a,/,a,/,a)')                 &
            'bigdft_solver: initialising kernel for a 1D-wire system  ',&
            '#================= WARNING bigdft_solver ===============#',&
            ' Empty directions of the wire system has to be set       ',&
            ' along the X and Y-axis.                                 ',&
            ' Additionally, simulation box must be large enough to    ',&
            ' to avoid box edges effects in the Hartree potential.    ',&
            '#=======================================================#'
      end if
    case default
      ! This should be enough to break PSolver ;)
      dom = domain_null()
    end select

    ! Data distribution will be always global:
    call set( dict//'setup'//'global_data', .true.) ! Hardwired, it cannot be otherwise.

    ! Set kernel and cavity variables: 
    call set( dict//'kernel'//'isf_order', bigdft_isf_order)
    call set( dict//'environment'//'fd_order', bigdft_fd_order)  
    call set( dict//'setup'//'verbose', bigdft_verbose)
    if( bigdft_cavity) then
      call set( dict//'environment'//'delta', bigdft_delta)
      call set( dict//'environment'//'fact_rigid', bigdft_fact_rigid)
      call set( dict//'environment'//'cavity', bigdft_cavity_type) 
      call set( dict//'environment'//'radii_set', bigdft_radii_type)
      call set( dict//'environment'//'gps_algorithm', bigdft_gps_algorithm)
      call set( dict//'environment'//'epsilon', bigdft_epsilon)
      call set( dict//'environment'//'gammaS', bigdft_gammaS)
      call set( dict//'environment'//'alphaS', bigdft_alphaS)
      call set( dict//'environment'//'betaV', bigdft_betaV) 
      call set( dict//'environment'//'atomic_radii', bigdft_atomic_radii)
    endif
    pkernel_verbose = bigdft_verbose   ! Verbose from setup is different from pkernel_init
    !$      call set(dict//'environment'//'pi_eta','0.6') ??

    ! Printing some info and warnings:
    if ( firsttime ) then
      if ( IONode ) then
        write(6,"(a)") 'bigdft_solver: Poisson solver variables'
        write(6,"(2a)")        'Solver boundary condition', shape
        write(6,"(a, L2)")    'verbosity:', bigdft_verbose
        write(6,"(a)")        'cavity type:', bigdft_cavity_type 
        write(6,"(a, I3)")    'isf_order:', bigdft_isf_order 
        write(6,"(a, I3)")    'fd_order:', bigdft_fd_order 
        if( bigdft_cavity ) then
          write(6,"(a, F10.5)") 'delta:', bigdft_delta 
          write(6,"(a, F10.5)") 'fact_rigid:', bigdft_fact_rigid 
          write(6,"(a, F10.5)") 'epsilon:', bigdft_epsilon 
          write(6,"(a, F10.5)") 'gammaS:', bigdft_gammaS
          write(6,"(a, F10.5)") 'alphaS:', bigdft_alphaS
          write(6,"(a, F10.5)") 'betaV:', bigdft_betaV
          write(6,"(a)")        'radii_set:', bigdft_radii_type
          write(6,"(a)")        'gps_algorithm:', bigdft_gps_algorithm
        endif
      end if
    end if

  end subroutine set_bigdft_variables

  subroutine gather_rho(rho, mesh, nsm, maxp, rho_total)

    ! Based on write_rho subroutine writen by J.Soler July 1997.           !
    ! Modifications are meant to send information to rho_total array       !
    ! istead to a *.RHO file.                                              !
    !                                                                      !
    ! #################################################################### !

    use alloc,          only : re_alloc, de_alloc
    use mpi_siesta

    integer, intent( in)          ::     nsm         
    integer, intent( in)          ::     mesh( 3)
    integer, intent( in)          ::     maxp
    real( grid_p), intent( in)    ::     rho( maxp)
    real( grid_p), pointer , intent( out) ::     rho_total( :) 
    integer :: i, ip, iu, is, np, BlockSizeY, BlockSizeZ, ProcessorZ
    integer :: meshnsm( 3), NRemY, NRemZ, iy, iz, izm, Ind, Ind_rho, ir
    integer :: Ind2, MPIerror, Status( MPI_Status_Size), BNode, NBlock
    real(grid_p), pointer :: bdens( :) => null()

#ifdef DEBUG
    call write_debug( 'ERROR: write_rho not ready yet' )
    call write_debug( fname )
#endif

    ! Work out density block dimensions
    if (mod(Nodes,ProcessorY).gt.0)                               &
        call die('ERROR: ProcessorY must be a factor of the' //  &
        ' number of processors!')
    ProcessorZ = Nodes/ ProcessorY
    BlockSizeY =((((mesh( 2)/ nsm)- 1)/ ProcessorY) + 1) *nsm

    call re_alloc(bdens, 1, BlockSizeY*mesh(1), 'bdens', 'write_rho')

    np = mesh(1) * mesh(2) * mesh(3)

    meshnsm(1) = mesh(1) / nsm
    meshnsm(2) = mesh(2) / nsm
    meshnsm(3) = mesh(3) / nsm

    Ind = 0
    Ind_rho = 1

    ! Loop over Z dimension of processor grid
    do iz = 1, ProcessorZ

      BlockSizeZ = meshnsm(3) / ProcessorZ
      NRemZ = meshnsm(3) - BlockSizeZ * ProcessorZ
      if ( iz - 1 < NRemZ ) BlockSizeZ = BlockSizeZ + 1
      BlockSizeZ = BlockSizeZ * nsm

      ! Loop over local Z mesh points
      do izm = 1, BlockSizeZ

        ! Loop over blocks in Y mesh direction
        do iy = 1, ProcessorY

          ! Work out size of density sub-matrix to be transfered
          BlockSizeY = meshnsm(2) / ProcessorY
          NRemY = meshnsm(2) - BlockSizeY * ProcessorY
          if ( iy-1 < NRemY ) BlockSizeY = BlockSizeY + 1
          BlockSizeY = BlockSizeY * nsm
          NBlock = BlockSizeY * mesh(1)

          ! Work out which node block is stored on
          BNode = (iy -1) * ProcessorZ + iz - 1
          if ( BNode == 0 .and. Node == BNode ) then
            do ir = 1, BlockSizeY
              do i =  1, mesh(1)
                rho_total(Ind_rho) = rho(Ind+i)
                Ind_rho = Ind_rho + 1
              end do
              Ind = Ind + mesh(1)
            end do
          else if ( Node == 0 ) then
            call MPI_Recv( bdens, NBlock, MPI_grid_real, BNode,  &
                1, MPI_Comm_World, Status,  MPIerror)
          else if ( Node == BNode ) then
            call MPI_Send(rho(Ind + 1), NBlock, MPI_grid_real,&
                0, 1, MPI_Comm_World, MPIerror)
            Ind = Ind + NBlock
          end if

          if ( BNode /= 0 ) then
            if ( Node == 0 ) then
              Ind2 = 0
              do ir = 1, BlockSizeY
                do i =  1 ,  mesh( 1) 
                  rho_total(Ind_rho) = bdens(Ind2+i)
                  Ind_rho = Ind_rho + 1
                end do
                Ind2 = Ind2 + mesh( 1)
              end do
            end if
          end if
        end do

      end do
    end do
    
    ! Deallocate density buffer memory
    call de_alloc( bdens, 'bdens', 'write_rho' )

  end subroutine gather_rho

  subroutine spread_potential( V_total, mesh, nsm, maxp,           &
      local_grids, V_local)

    ! #################################################################### !
    !                                                                      !
    !     Based on read_rho writen by J.Soler July 1997.                   !
    !     This subroutine distributes a potential from a total array       !
    !     into all nodes according to local_grids meshes.                  !
    !                                                                      !
    ! #################################################################### !
    ! *************************** INPUT **********************************
    ! real    V_total(maxp)   : Total potential. 
    ! integer nsm             : Number of sub-mesh points per mesh point
    !                           (not used in this version)
    ! integer maxp            : First dimension of array rho
    ! ************************** OUTPUT **********************************
    ! real*8  cell(3,3)       : Lattice vectors
    ! integer mesh(3)         : Number of mesh divisions of each
    !                           lattice vector
    ! real    rho(maxp) : Electron density
    ! #################################################################### !

    use alloc,          only : re_alloc, de_alloc
    use parallelsubs, only: HowManyMeshPerNode
    use mpi_siesta

    real(grid_p), pointer,intent(in) :: V_total(:) 
    real(grid_p), intent(out) :: V_local(:) 
    integer, intent(in)      ::     nsm, maxp
    integer, intent(in)      ::     mesh(3)
    integer, intent(in) ::  local_grids(:,:)
    real(grid_p)             ::     V_aux(maxp)
    external          io_assign, io_close, memory
    integer    ip, iu, is, np, ns,                   & 
        BlockSizeY, BlockSizeZ, ProcessorZ,   &
        meshnsm(3), npl, NRemY, NRemZ,        &
        iy, iz, izm, Ind, Ind_Vaux,  ir

    integer    Ind2, MPIerror, Request, meshl(3),     &
        Status(MPI_Status_Size), BNode, NBlock
    logical    ltmp
    real(grid_p), pointer :: bdens(:) => null()
    logical    baddim, found

    ! Work out density block dimensions
    if (mod(Nodes,ProcessorY).gt.0)                               &
        call die('ERROR: ProcessorY must be a factor of the' //  &
        ' number of processors!')
    ProcessorZ = Nodes / ProcessorY
    BlockSizeY = (((mesh(2) / nsm) - 1) / ProcessorY + 1) * nsm

    call re_alloc(bdens, 1, BlockSizeY*mesh(1), 'bdens', &
        'spread_potential' )

    np = mesh(1) * mesh(2) * mesh(3)

    ! Get local dimensions
    meshnsm(:) = mesh(:) / nsm

    call HowManyMeshPerNode(meshnsm,Node,Nodes,npl,meshl)

    ! Check dimensions
    baddim = .false.

    if ( npl .gt. maxp ) baddim = .true.

    ! Globalise dimension check
    call MPI_AllReduce(baddim,ltmp,1,MPI_logical,MPI_Lor,  &
        MPI_Comm_World,MPIerror)
    baddim = ltmp
    if ( baddim ) then
      call die("Dim of array V_total too small in spread_potential")
    endif

    Ind = 0
    Ind_Vaux = 1

    ! Loop over Z mesh direction
    do iz = 1, ProcessorZ

      ! Work out number of mesh points in Z direction
      BlockSizeZ = meshnsm(3) / ProcessorZ
      NRemZ = meshnsm(3) - BlockSizeZ * ProcessorZ
      if ( iz - 1 < NRemZ ) BlockSizeZ = BlockSizeZ + 1
      BlockSizeZ = BlockSizeZ * nsm

      ! Loop over local Z mesh points
      do izm = 1, BlockSizeZ

        ! Loop over blocks in Y mesh direction
        do iy = 1, ProcessorY

          ! Work out size of density sub-matrix to be transfered
          BlockSizeY = meshnsm(2) / ProcessorY
          NRemY = meshnsm(2) - BlockSizeY * ProcessorY
          if ( iy -1 < NRemY ) BlockSizeY = BlockSizeY + 1
          BlockSizeY = BlockSizeY * nsm
          NBlock = BlockSizeY * mesh(1)

          ! Work out which node block is stored on
          BNode = (iy - 1) * ProcessorZ + iz - 1
          if (BNode == 0 .and. Node == BNode ) then
            ! If density sub-matrix is local Node 0 then just read it in
            do ir = 1, BlockSizeY
              do ip = 1, mesh( 1) 
                V_aux(Ind + ip) = V_total(Ind_Vaux)
                Ind_Vaux = Ind_Vaux + 1
              enddo
              Ind = Ind + mesh( 1)
            enddo

          elseif ( Node == 0 ) then

            ! If this is Node 0 then read and send density sub-matrix
            Ind2 = 0
            do ir = 1, BlockSizeY
              do ip = 1, mesh(1)
                bdens(Ind2 + ip) = V_total(Ind_Vaux)
                Ind_Vaux = Ind_Vaux + 1
              end do
              Ind2 = Ind2 + mesh( 1)
            end do
            call MPI_Send(bdens, NBlock, MPI_grid_real, BNode, &
                1, MPI_Comm_World, MPIerror)

          else if ( Node == BNode ) then
            
            ! If this is the Node where the density sub-matrix is, then receive
            call MPI_Recv(V_aux(Ind + 1), NBlock, &
                MPI_grid_real, 0, 1, MPI_Comm_World, Status, MPIerror)
            Ind = Ind + NBlock

          end if

        end do

      end do

    end do

    V_local(:) = V_aux(1:product(local_grids(:,Node+1)))

    ! Deallocate density buffer memory
    call de_alloc( bdens, 'bdens', 'spread_potential' )

  end subroutine spread_potential

  
  subroutine set_label(atm_label)

    ! Assignation of atom labels. 
    use siesta_geom,    only : na_u, isa
    use chemical,       only : atomic_number

    character(len=2), intent(out), dimension(na_u) :: atm_label
    integer :: aux_Zatom(na_u)             ! Auxiliary array containing Z_atom for each atom. 
    integer :: i

    do i = 1, na_u
      aux_Zatom(i) = atomic_number(isa(i))
    end do

    do i = 1, na_u
      select case ( aux_Zatom(i) )
      case ( 1 )
        atm_label(i) = "H "
      case ( 2 )
        atm_label(i) = "He"
      case ( 3 )
        atm_label(i) = "Li"
      case ( 4 )
        atm_label(i) = "Be"
      case ( 5 )
        atm_label(i) = "B "
      case ( 6 )
        atm_label(i) = "C "
      case ( 7 )
        atm_label(i) = "N "
      case ( 8 )
        atm_label(i) = "O "
      case ( 9 )
        atm_label(i) = "F "
      case ( 10 )
        atm_label(i) = "Ne"
      case ( 11 )
        atm_label(i) = "Na"
      case ( 12 )
        atm_label(i) = "Mg"
      case ( 13 )
        atm_label(i) = "Al"
      case ( 14 )
        atm_label(i) = "Si"
      case ( 15 )
        atm_label(i) = "P "
      case ( 16 )
        atm_label(i) = "S "
      case ( 17 )
        atm_label(i) = "Cl"
      case ( 18 )
        atm_label(i) = "Ar"
      case ( 19 )
        atm_label(i) = "K "
      case ( 20 )
        atm_label(i) = "Ca"
      case ( 21 )
        atm_label(i) = "Sc"
      case ( 22 )
        atm_label(i) = "Ti"
      case ( 23 )
        atm_label(i) = "V "
      case ( 24 )
        atm_label(i) = "Cr"
      case ( 25 )
        atm_label(i) = "Mn"
      case ( 26 )
        atm_label(i) = "Fe"
      case ( 27 )
        atm_label(i) = "Co"
      case ( 28 )
        atm_label(i) = "Ni"
      case ( 29 )
        atm_label(i) = "Cu"
      case ( 30 )
        atm_label(i) = "Zn"
      case ( 31 )
        atm_label(i) = "Ga"
      case ( 32 )
        atm_label(i) = "Ge"
      case ( 33 )
        atm_label(i) = "As"
      case ( 34 )
        atm_label(i) = "Se"
      case ( 35 )
        atm_label(i) = "Br"
      case ( 36 )
        atm_label(i) = "Kr"
      case ( 37 )
        atm_label(i) = "Rb"
      case ( 38 )
        atm_label(i) = "Sr"
      case ( 39 )
        atm_label(i) = "Y "
      case ( 40 )
        atm_label(i) = "Zr"
      case ( 41 )
        atm_label(i) = "Nb"
      case ( 42 )
        atm_label(i) = "Mo"
      case ( 43 )
        atm_label(i) = "Tc"
      case ( 44 )
        atm_label(i) = "Ru"
      case ( 45 )
        atm_label(i) = "Rh"
      case ( 46 )
        atm_label(i) = "Pd"
      case ( 47 )
        atm_label(i) = "Ag"
      case ( 48 )
        atm_label(i) = "Cd"
      case ( 49 )
        atm_label(i) = "In"
      case ( 50 )
        atm_label(i) = "Sn"
      case ( 51 )
        atm_label(i) = "Sb"
      case ( 52 )
        atm_label(i) = "Te"
      case ( 53 )
        atm_label(i) = "I "
      case ( 54 )
        atm_label(i) = "Xe"
      case ( 55 )
        atm_label(i) = "Cs"
      case ( 56 )
        atm_label(i) = "Ba"
      case ( 57 )
        atm_label(i) = "La"
      case ( 58 )
        atm_label(i) = "Ce"
      case ( 59 )
        atm_label(i) = "Pr"
      case ( 60 )
        atm_label(i) = "Nd"
      case ( 61 )
        atm_label(i) = "Pm"
      case ( 62 )
        atm_label(i) = "Sm"
      case ( 63 )
        atm_label(i) = "Eu"
      case ( 64 )
        atm_label(i) = "Gd"
      case ( 65 )
        atm_label(i) = "Tb"
      case ( 66 )
        atm_label(i) = "Dy"
      case ( 67 )
        atm_label(i) = "Ho"
      case ( 68 )
        atm_label(i) = "Er"
      case ( 69 )
        atm_label(i) = "Tm"
      case ( 70 )
        atm_label(i) = "Yb"
      case ( 71 )
        atm_label(i) = "Lu"
      case ( 72 )
        atm_label(i) = "Hf"
      case ( 73 )
        atm_label(i) = "Ta"
      case ( 74 )
        atm_label(i) = "W "
      case ( 75 )
        atm_label(i) = "Re"
      case ( 76 )
        atm_label(i) = "Os"
      case ( 77 )
        atm_label(i) = "Ir"
      case ( 78 )
        atm_label(i) = "Pt"
      case ( 79 )
        atm_label(i) = "Au"
      case ( 80 )
        atm_label(i) = "Hg"
      case ( 81 )
        atm_label(i) = "Tl"
      case ( 82 )
        atm_label(i) = "Pb"
      case ( 83 )
        atm_label(i) = "Bi"
      case ( 84 )
        atm_label(i) = "Po"
      case ( 85 )
        atm_label(i) = "At"
      case ( 86 )
        atm_label(i) = "Rn"
      case ( 87 )
        atm_label(i) = "Fr"
      case ( 88 )
        atm_label(i) = "Ra"
      case ( 89 )
        atm_label(i) = "Ac"
      case ( 90 )
        atm_label(i) = "Th"
      case ( 91 )
        atm_label(i) = "Pa"
      case ( 92 )
        atm_label(i) = "U "
      case ( 93 )
        atm_label(i) = "Np"
      case ( 94 )
        atm_label(i) = "Pu"
      case ( 95 )
        atm_label(i) = "Am"
      case ( 96 )
        atm_label(i) = "Cm"
      case ( 97 )
        atm_label(i) = "Bk"
      case ( 98 )
        atm_label(i) = "Cf"
      case ( 99 )
        atm_label(i) = "Es"
      case ( 100 )
        atm_label(i) = "Fm"
      case ( 101 )
        atm_label(i) = "Md"
      case ( 102 )
        atm_label(i) = "No"
      case ( 103 )
        atm_label(i) = "Lr"
      case ( 104 )
        atm_label(i) = "Rf"
      case ( 105 )
        atm_label(i) = "Db"
      case ( 106 )
        atm_label(i) = "Sg"
      case ( 107 )
        atm_label(i) = "Bh"
      case ( 108 )
        atm_label(i) = "Hs"
      case ( 109 )
        atm_label(i) = "Mt"
      case ( 110 )
        atm_label(i) = "Ds"
      end select

    end do

  end subroutine set_label

  subroutine pointer_1D_3D(A1D, n1, n2, n3, A3D)
    integer, intent(in) :: n1, n2, n3
    real(grid_p), target, intent(in) :: A1D(n1, n2, n3)
    real(grid_p), pointer :: A3D(:,:,:)
    A3D => A1D(:,:,:)
  end subroutine pointer_1D_3D

end module m_psolver_bigdft
