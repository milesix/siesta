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
!
      module m_psolver_bigdft

!
      use precision,      only : dp, grid_p
      use parallel,       only : Node, Nodes, ionode, ProcessorY
      use fdf 
      use units, only : Ang, pi
      implicit none
!
      logical, save         :: firsttime = .true.
!
      public :: poisson_bigdft
!
      CONTAINS
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
!
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
!
! #################################################################### !
!
      real( dp), intent( in)      :: cell( 3, 3)        ! Unit-cell vectors.
      integer, intent( in)        :: lgrid( 3)          ! Local mesh divisions.
      integer, intent( in)        :: grid( 3)           ! Total mesh divisions.
      real( grid_p), intent( in)  :: RHO( :)            ! Input density.
      integer, intent(in)         :: ntm
      real( grid_p), intent(out)  :: V( :)              ! Output potential.
      real( dp), intent(out)      :: eh                 ! Electrostatic energy.
      real( dp), dimension( 1, 3), save :: hgrid        ! Uniform mesh spacings in the three directions.
      integer                     :: i, j, k, iatoms
      real( dp)                   :: einit
      type( dictionary), pointer  :: dict               ! Input parameters.
      integer                     :: isf_order, fd_order
      type( coulomb_operator)     :: pkernel
      logical                     :: pkernel_verbose 
!
      real(grid_p), pointer       :: Paux_1D( :)        ! 1D auxyliary array.
      real(grid_p), pointer       :: Paux_3D( :, :, :)  ! 3D auxyliary array.
      integer                     :: World_Comm, mpirank, ierr
!
      real(dp), save              :: alpha, beta, gamma ! Box axis angles.  
      real(dp)                    :: angdeg( 3)
      real(dp)                    :: radii( na_u) 
      character( len= 1)          :: boundary_case 
      character( len= 2)          :: atm_label( na_u) 
      real(dp)                    :: delta, factor 
      integer, allocatable        :: local_grids( :)
!
! #################################################################### !
!
      call timer( 'BigDFT_solv', 1)
! 
! f_util and dictionary initialisations:
! 
      call f_lib_initialize()
      call dict_init(dict)
!
! Allocation:
!
      nullify( Paux_3D)
      nullify( Paux_1D)
!
      allocate( local_grids( Nodes* 3))
      call re_alloc( Paux_3D, 1, grid( 1), 1, grid( 2), 1, grid( 3),   &
                     'Paux_3D', 'poison_bigdft' )
      call re_alloc( Paux_1D, 1, grid( 1)* grid( 2)* grid( 3),         &
                     'Paux_1D', 'poison_bigdft' )
!
! Calculation of cell angles and grid separation for psolver:
!
      if( firsttime) then
!
        call set_cell_angles( cell, angdeg) 
! 
        alpha= angdeg( 1)/ 180.0_dp* pi
        beta = angdeg( 2)/ 180.0_dp* pi
        gamma= angdeg( 3)/ 180.0_dp* pi
!
        do i= 1, 3
          hgrid( 1, i)=( sqrt( cell( 1, i)**2 +cell( 2, i)**2 +        &
                               cell( 3, i)**2)/ grid(i))
        enddo
      endif
!
! Setting of dictionary variables, done at every calling:  
!
      call set_bigdft_variables( dict, boundary_case, pkernel_verbose)
!
#ifdef MPI
!
!!! Collection of RHO: 
!
! RHO is collected to a common array.
!
!        call timer('gather_rho',1)
        call gather_rho( RHO, grid, 2, lgrid(1) * lgrid(2) * lgrid(3), Paux_1D)
!
        call MPI_BCAST( Paux_1D, grid(1) * grid(2) * grid(3),         &
                        MPI_double_precision, 0, MPI_COMM_WORLD, ierr) 
!        call timer('gather_rho',2)
!
!!! Local grids: 
!
! ... and grids on each of the nodes is store. 
! 
        local_grids(1 + node * 3: 3 * (node + 1)) = lgrid( :)
!
        call MPI_Gather( local_grids(1+ node * 3: 3 * (node+1)), 3,    &
             MPI_integer, local_grids(1 + node * 3: 3 *(node + 1)), 3,&
             MPI_integer, 0, MPI_COMM_WORLD, ierr)
!
        call MPI_BCAST( local_grids, Nodes * 3, MPI_integer, 0,        &
                        MPI_COMM_WORLD, ierr)
#endif
!
#ifdef MPI
!
!!! Total density passes from a 1D to a 3D-aaray: (compulsory for BigDFT)
!
! shame on me! I should find more elegant way to make it. 
!
        call onetothreedarray( Paux_1D, Paux_3D, grid)
#else
        call onetothreedarray( RHO, Paux_3D, grid)
#endif
!
!!! BigDFT solver starts now
!
!
        pkernel=pkernel_init( Node, Nodes, dict, boundary_case,        &
                         (/ grid( 1), grid( 2), grid( 3)/), hgrid,     &
                        alpha_bc= alpha, beta_ac= beta, gamma_ab= gamma)
!
        call pkernel_set( pkernel, verbose=pkernel_verbose)
!
! Implicit solvent?
!
        if( bigdft_cavity) then
!
! Set a radius for each atom in Angstroem (UFF, Bondi, Pauling, etc ...)
!
          call set_label( atm_label)
          do iatoms= 1, na_u
            radii( iatoms) = pkernel_get_radius( pkernel, atname= atm_label( iatoms))
          end do
!
! Multiply for a constant prefactor (see the Soft-sphere paper JCTC 2017)
!
          factor= pkernel%cavity%fact_rigid
          do iatoms= 1, na_u
            radii( iatoms)= factor* radii( iatoms)* Ang
          enddo
          call pkernel_set_epsilon( pkernel, nat= na_u, rxyz= xa, radii= radii)
        endif
!        
        call dict_free( dict)
!
!        call timer( 'PSolver', 1)
        call Electrostatic_Solver( pkernel, Paux_3D, ehartree= eh)
!        call timer('PSolver',2)
!
! Energy must be devided by the total number of nodes because the way 
! Siesta computes the Hartree energy. Indeed, BigDFT provides the value
! alredy distributed. 
!
        eh= 2.0_dp* eh/ Nodes
        Paux_3D= 2.0_dp* Paux_3D-( sum( 2.0_dp* Paux_3D)/              &
                 (grid( 1) * grid( 2) * grid( 3)))
!
#ifdef MPI
! 
! we go back from a 3D to a 1D array: 
!
        call threetoonedarray( Paux_1D, Paux_3D, grid) 
!
! Now we spread back the potential on different nodes. 
!
        call timer('spread_pot',1)
        call spread_potential( Paux_1D, grid, 2, grid( 1) * grid( 2) * &
                   grid( 3), local_grids, V) 
        call timer('spread_pot',2)
# else 
        call threetoonedarray(V, Paux_3D, grid) 
#endif
!
! Ending calculation:
!
        call pkernel_free( pkernel)
        call f_lib_finalize_noreport() 
        call de_alloc( Paux_1D, 'Paux', 'poisson_bigdft' )
        call de_alloc( Paux_3D, 'Paux', 'poisson_bigdft' )
        deallocate( local_grids)
 
      call timer( 'BigDFT_solv', 2)
!
      firsttime=.false.
!
! #################################################################### !
      end subroutine poisson_bigdft 
! #################################################################### !
!
! #################################################################### !
      subroutine set_bigdft_variables(dict, boundary_case, pkernel_verbose)
! #################################################################### !
      use dictionaries, dict_set => set
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
      type( dictionary), pointer, intent(out) :: dict         ! Input parameters.
      character( len= 1), intent(out)         :: boundary_case 
      logical, intent(out)  :: pkernel_verbose 
! #################################################################### !
!
! Data distribution will be always global:
!
      call set( dict//'setup'//'global_data', .true.) ! Hardwired, it cannot be otherwise.
!
! Select the boundary condition:
!
      select case(shape)
       case( 'molecule') 
        boundary_case= "F"
       case( 'slab')
        boundary_case= "S"
       case( 'bulk')
        boundary_case= "P"
!      1D WIRE SYSTEM ===NOT IMPLEMENTED YET=== 
!       case( 'wire')  
!        boundary_case= "W"
      end select 
! 
! Set kernel and cavity variables: 
!
!
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
      pkernel_verbose= bigdft_verbose   ! Verbose from setup is different from pkernel_init
!$      call set(dict//'environment'//'pi_eta','0.6') ??
! 
! Printing some info and warnings:
!
!#ifdef DEBUG
      if (firsttime) then 
        if( boundary_case.eq. 'F') then  
          if (ionode) write(6,'(/,a,/,a/,a,/,a,/,a,/,a)')               &
            'bigdft_solver: initialising kernel for a molecular system',&
            '#================= WARNING bigdft_solver ===============#',&
            '                                                         ',&
            ' Molecular system must be place at the center of the box,',&
            ' otherwise box-edge effects might lead to wrong Hartree  ',&
            ' energy values.                                          ',&
            '#=======================================================#'
      elseif (boundary_case.eq. 'S') then  
        if (ionode) write(6,'(/,a,/,a/,a,/,a,/,a,/,a)')                 &
            'bigdft_solver: initialising kernel for a slab system',     &
            '#================= WARNING bigdft_solver ===============#',&
            ' Empty direction of the periodic system has to be set    ',&
            ' along the Y-axis.                                       ',&
            ' Additionally, simulation box must be large enough to    ',&
            ' to avoid box edges effects in the Hartree potential.    ',&
            '#=======================================================#'
      elseif( boundary_case.eq. 'P') then  
        if (ionode) write(6,"(a)")                                     &
          'bigdft_solver: initialising kernel for a periodic system'
!      1D WIRE SYSTEM ===NOT IMPLEMENTED YET=== 
!        elseif (boundary_case .eq. 'W') then                          & 
!          'bigdft_solver: initialising kernel for a 1D-wire system' 
      endif 
!
        if (ionode) then 
          write(6,"(a)") 'bigdft_solver: Poisson solver variables'
          write(6,* )    'Solver boundary condition: ', boundary_case 
          write(6,"(a, L2)")    ' verbosity:', bigdft_verbose
          write(6,"(a)")        ' cavity type:', bigdft_cavity_type 
          write(6,"(a, I3)")    ' isf_order:', bigdft_isf_order 
          write(6,"(a, I3)")    ' fd_order:', bigdft_fd_order 
          if( bigdft_cavity) then
            write(6,"(a, F10.5)") ' delta:', bigdft_delta 
            write(6,"(a, F10.5)") ' fact_rigid:', bigdft_fact_rigid 
            write(6,"(a, F10.5)") ' epsilon:', bigdft_epsilon 
            write(6,"(a, F10.5)") ' gammaS:', bigdft_gammaS
            write(6,"(a, F10.5)") ' alphaS:', bigdft_alphaS
            write(6,"(a, F10.5)") ' betaV:', bigdft_betaV
            write(6,"(a)")        ' radii_set:', bigdft_radii_type
            write(6,"(a)")        ' gps_algorithm:', bigdft_gps_algorithm
          endif
        endif
      endif 
!#endif
! #################################################################### !
      end subroutine set_bigdft_variables
! #################################################################### !
      subroutine set_cell_angles( cell, celang)
! #################################################################### !
!     Subroutine to calculate the cell angles. Roughly copied from     !
!     outcell.f of E. Artacho, December 1997.
! #################################################################### !
! #################################################################### !
      real( dp), intent( in)      :: cell( 3, 3)        ! Unit-cell vectors.
      real( dp), intent( out)     :: celang( 3)
      real( dp)                   :: cellm( 3)
      integer                     :: i
! #################################################################### !
      do i = 1,3
        cellm( i) = dot_product( cell(:,i),cell(:,i))
        cellm( i) = sqrt( cellm( i))
      enddo
!
      celang( 1)= dot_product( cell(:,1),cell(:,2))
      celang( 1)= acos( celang( 1)/(cellm( 1)*cellm( 2)))* 180._dp/ pi
      celang( 2)= dot_product( cell( :,1), cell( :,3))
      celang( 2)= acos( celang( 2)/(cellm(1)*cellm( 3)))* 180._dp/ pi
      celang( 3)= dot_product( cell( :,2),cell(:,3))
      celang( 3)= acos( celang( 3)/( cellm( 2)* cellm( 3)))* 180._dp/ pi
! #################################################################### !
      end subroutine set_cell_angles
! #################################################################### !
! #################################################################### !
      subroutine onetothreedarray(rho, Paux, grid)
! #################################################################### !
! Utility to copy 1D array to output 3D array.                         !
! Arguments:                                                           !
! rho(in)   :: 1D array.                                               !
! Paux(out) :: 3D array.                                               !
! grid      :: 3D array dimiensions.                                   !
! #################################################################### !
      integer, intent(in)       :: grid(3) 
      real(grid_p), intent(in)  :: rho(:)
      real(grid_p), pointer,intent(out) :: Paux(:, :, :)
      integer                   :: i, j, k, w
! #################################################################### !
!
      w = 1
      do i = 1, grid( 1)
!
       do j = 1, grid( 2)
!
        do k = 1, grid(3)

         Paux( i, j, k) = rho( w)
!
         w = w + 1
!
        enddo
!
       enddo
!
      enddo
! #################################################################### !
      end subroutine onetothreedarray 
! #################################################################### !
!
! #################################################################### !
      subroutine threetoonedarray(rho, Paux, grid)
! #################################################################### !
! Utility to copy 3D array to output 1D array.                         !
! Arguments:                                                           !
! rho(out) :: 1D array.                                                !
! Paux(in) :: 3D array.                                                !
! grid     :: 3D array dimiensions.                                    !
! #################################################################### !
      integer, intent(in)               :: grid( 3) 
      real(grid_p), intent(out)         :: rho( :)
      real(grid_p), pointer,intent(in)  :: Paux(:, :, :)
      integer                           :: i, j, k, w
! #################################################################### !
!
      w = 1
      do i = 1, grid( 1)
!
       do j = 1, grid( 2)
!
        do k = 1, grid( 3)
!
         rho( w)=Paux( i, j, k)
!
         w = w + 1
        enddo
!
       enddo
!
      enddo
!
! #################################################################### !
      end subroutine threetoonedarray
! #################################################################### !
!
! #################################################################### !
      subroutine gather_rho(rho, mesh, nsm, maxp, rho_total)
!                                                                      !
! Based on write_rho subroutine writen by J.Soler July 1997.           !
! Modifications are meant to send information to rho_total array       !
! istead to a *.RHO file.                                              !
!                                                                      !
! #################################################################### !
!
      use alloc,          only : re_alloc, de_alloc
      use mpi_siesta
!
      integer, intent( in)          ::     nsm         
      integer, intent( in)          ::     mesh( 3)
      integer, intent( in)          ::     maxp
      real( grid_p), intent( in)    ::     rho( maxp)
      real( grid_p), pointer , intent( out) ::     rho_total( :) 
      integer  i, ip, iu, is, np, BlockSizeY, BlockSizeZ, ProcessorZ
      integer  meshnsm( 3), NRemY, NRemZ, iy, iz, izm, Ind, Ind_rho, ir
      integer  Ind2, MPIerror, Request, Status( MPI_Status_Size),      &
               BNode, NBlock
      real(grid_p), pointer         :: bdens( :) => null()
      real(grid_p), pointer         :: temp( :) => null()
! #################################################################### !
!
#ifdef DEBUG
        call write_debug( 'ERROR: write_rho not ready yet' )
        call write_debug( fname )
#endif
!
! Work out density block dimensions
!
      if (mod(Nodes,ProcessorY).gt.0)                               &
           call die('ERROR: ProcessorY must be a factor of the' //  &
           ' number of processors!')
!
      ProcessorZ= Nodes/ ProcessorY
      BlockSizeY=((((mesh( 2)/ nsm)- 1)/ ProcessorY) + 1) *nsm
!
      call re_alloc(bdens, 1, BlockSizeY*mesh(1), 'bdens', 'write_rho')
      call re_alloc( temp, 1, mesh(1), 'temp', 'write_rho' )
!
      np = mesh( 1) * mesh( 2) * mesh( 3)

      meshnsm( 1) = mesh( 1) / nsm
      meshnsm( 2) = mesh( 2) / nsm
      meshnsm( 3) = mesh( 3) / nsm
!
      Ind = 0
      Ind_rho= 1
!
! Loop over Z dimension of processor grid
!
      do iz = 1, ProcessorZ
!
        BlockSizeZ = (meshnsm( 3) / ProcessorZ)
        NRemZ = meshnsm( 3) - BlockSizeZ * ProcessorZ
        if (iz - 1 .lt. NRemZ) BlockSizeZ = BlockSizeZ + 1
        BlockSizeZ = BlockSizeZ * nsm
!
! Loop over local Z mesh points
!
          do izm = 1, BlockSizeZ
!
! Loop over blocks in Y mesh direction
!
            do iy = 1, ProcessorY
!
! Work out size of density sub-matrix to be transfered
!
              BlockSizeY = (meshnsm(2)/ProcessorY)
              NRemY = meshnsm(2) - BlockSizeY*ProcessorY
              if (iy-1 .lt. NRemY) BlockSizeY = BlockSizeY + 1
              BlockSizeY = BlockSizeY * nsm
              NBlock = BlockSizeY * mesh( 1)
!
! Work out which node block is stored on
!
              BNode = ( iy -1) * ProcessorZ + iz - 1
              if (BNode .eq. 0 .and. Node .eq. BNode) then
                do ir = 1, BlockSizeY
                  temp( 1: mesh(1)) =                                  &
                        real( rho( Ind + 1: Ind + mesh( 1)), kind = dp)
!
                  do i =  1, mesh( 1) 
                    rho_total( Ind_rho) = temp( i)
                    Ind_rho = Ind_rho + 1
                  enddo 
!
                  Ind = Ind + mesh( 1)
                enddo
!
               elseif (Node.eq.0) then
! If this is Node 0 then recv and write density sub-matrix
!
                 call MPI_IRecv( bdens, NBlock, MPI_grid_real, BNode,  &
                      1, MPI_Comm_World, Request, MPIerror)
!
                 call MPI_Wait( Request, Status, MPIerror)
!
                elseif ( Node .eq. BNode) then
!
! If this is the Node where the density sub-matrix is, then send
                  call MPI_ISend( rho( Ind + 1), NBlock, MPI_grid_real,&
                       0, 1, MPI_Comm_World, Request, MPIerror)
                  call MPI_Wait( Request, Status, MPIerror)
                  Ind = Ind + NBlock
                endif

                if (BNode .ne. 0) then
!
                  call MPI_Barrier(MPI_Comm_World,MPIerror)
                  if (Node .eq. 0) then
                    Ind2 = 0
                  do ir = 1, BlockSizeY
                     temp( 1 :mesh( 1)) =                               & 
                       real(bdens( Ind2 + 1: Ind2 +mesh( 1)), kind = dp)
                    do i =  1, mesh( 1) 
                      rho_total( Ind_rho) = temp( i)
                      Ind_rho = Ind_rho + 1
                    enddo 
!
                  Ind2 = Ind2 + mesh( 1)
!
                enddo
!
              endif
            endif
!
          enddo

        enddo
!
      enddo
!
! Deallocate density buffer memory
!
      call de_alloc( bdens, 'bdens', 'write_rho' )
      call de_alloc( temp, 'temp', 'write_rho' )
! #################################################################### !
      end subroutine gather_rho
! #################################################################### !
! #################################################################### !
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
!
      use alloc,          only : re_alloc, de_alloc
      use parallelsubs, only: HowManyMeshPerNode
      use mpi_siesta
!
! #################################################################### !
!
      real(grid_p), pointer,intent(in) :: V_total(:) 
      real(grid_p), intent(out) :: V_local(:) 
      integer, intent(in)      ::     nsm, maxp
      integer, intent(in)      ::     mesh(3)
      integer, intent(in) ::  local_grids(:)
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
! #################################################################### !
!
! Work out density block dimensions
!
      if (mod(Nodes,ProcessorY).gt.0)                               &
           call die('ERROR: ProcessorY must be a factor of the' //  &
           ' number of processors!')
      ProcessorZ = Nodes / ProcessorY
      BlockSizeY = ((((mesh( 2) / nsm) - 1) / ProcessorY) + 1) * nsm
!
      call re_alloc( bdens, 1, BlockSizeY*mesh(1), 'bdens',         &
                    'spread_potential' )
!
      np = mesh( 1) * mesh( 2) * mesh( 3)
!
! Get local dimensions
!
      meshnsm( 1) = mesh( 1) / nsm
      meshnsm( 2) = mesh( 2) / nsm
      meshnsm( 3) = mesh( 3) / nsm

      call HowManyMeshPerNode(meshnsm,Node,Nodes,npl,meshl)

! Check dimensions
      baddim = .false.
!
      if ( npl .gt. maxp) baddim = .true.
!
! Globalise dimension check
!
      call MPI_AllReduce(baddim,ltmp,1,MPI_logical,MPI_Lor,  &
           MPI_Comm_World,MPIerror)
      baddim = ltmp
!
      if ( baddim) then
        call die("Dim of array V_total too small in spread_potential")
      endif
!
      Ind = 0
      Ind_Vaux = 1
!
! Loop over Z mesh direction
!
      do iz = 1, ProcessorZ
!
! Work out number of mesh points in Z direction
!
        BlockSizeZ = (meshnsm( 3)/ProcessorZ)
        NRemZ = meshnsm( 3) - BlockSizeZ * ProcessorZ
        if (iz - 1 .lt. NRemZ) BlockSizeZ = BlockSizeZ + 1
        BlockSizeZ = BlockSizeZ * nsm
! Loop over local Z mesh points
        do izm = 1, BlockSizeZ
! Loop over blocks in Y mesh direction
          do iy = 1, ProcessorY
! Work out size of density sub-matrix to be transfered
             BlockSizeY = (meshnsm( 2) / ProcessorY)
             NRemY = meshnsm( 2) - BlockSizeY * ProcessorY
             if (iy -1 .lt. NRemY) BlockSizeY = BlockSizeY + 1
             BlockSizeY = BlockSizeY * nsm
             NBlock = BlockSizeY * mesh( 1)
! Work out which node block is stored on
             BNode = (iy - 1) * ProcessorZ + iz - 1
             if (BNode .eq. 0 .and. Node .eq. BNode) then
! If density sub-matrix is local Node 0 then just read it in
               do ir = 1, BlockSizeY
                 do ip = 1, mesh( 1) 
                   V_aux( Ind + ip) = V_total( Ind_Vaux)
                   Ind_Vaux = Ind_Vaux + 1
                 enddo
                 Ind = Ind + mesh( 1)
               enddo
!
               elseif (Node .eq. 0) then
!
! If this is Node 0 then read and send density sub-matrix
                 Ind2 = 0
                 do ir = 1, BlockSizeY
                   do ip = 1, mesh(1)
                     bdens( Ind2 + ip) = V_total( Ind_Vaux)
                     Ind_Vaux = Ind_Vaux + 1
                   enddo
                   Ind2 = Ind2 + mesh( 1)
                 enddo
                 call MPI_ISend(bdens, NBlock, MPI_grid_real, BNode,   &
                       1, MPI_Comm_World, Request, MPIerror)
                 call MPI_Wait(Request, Status, MPIerror)

                elseif (Node .eq. BNode) then
! If this is the Node where the density sub-matrix is, then receive
                  call MPI_IRecv( V_aux( Ind + 1), NBlock,             &
                        MPI_grid_real, 0, 1, MPI_Comm_World,           &
                        Request, MPIerror)
                  call MPI_Wait( Request, Status, MPIerror)
                  Ind = Ind + NBlock
!
                endif
!
                if( BNode .ne. 0) then
                  call MPI_Barrier( MPI_Comm_World, MPIerror)
                endif
!
              enddo

           enddo

         enddo

      V_local( :) = V_aux( 1: local_grids( 1 + 3 * node) *             &
               local_grids( 2 + 3 * node) * local_grids( 3 + 3 * node))
!
! Deallocate density buffer memory
!
      call de_alloc( bdens, 'bdens', 'spread_potential' )
!
      end subroutine spread_potential 
!
!-----------------------------------------------------------------------
! 
      subroutine set_label( atm_label)
!
! Assignation of atom labels. 
!----------------------------------------------------------------------!
      use siesta_geom,    only : na_u, isa
      use chemical,       only : atomic_number
!
      character( len= 2), intent( out), dimension( na_u) :: atm_label
      integer :: aux_Zatom( na_u)             ! Auxiliary array containing Z_atom for each atom. 
      integer :: i
!**********************************************************************!
!
       do i=1, na_u
         aux_Zatom( i) = atomic_number( isa(i))
       enddo
!
       do i= 1, na_u
         select case(aux_Zatom( i))
           case( 1) 
            atm_label( i)="H "
           case( 2)
            atm_label( i)="He"
           case( 3)
            atm_label( i)="Li"
           case( 4)
            atm_label( i)="Be"
           case( 5)
            atm_label( i)="B "
           case( 6)     
            atm_label( i)="C "
           case( 7)     
            atm_label( i)="N "
           case( 8)     
            atm_label( i)="O "
           case( 9)   
            atm_label( i)="F "
           case( 10)   
            atm_label( i)="Ne"
           case( 11) 
            atm_label( i)="Na"
           case( 12)   
            atm_label( i)="Mg"
           case( 13) 
            atm_label( i)="Al"
           case( 14)   
            atm_label( i)="Si"
           case( 15)    
            atm_label( i)="P "
           case( 16)    
            atm_label( i)="S "
           case( 17)   
            atm_label( i)="Cl"
           case( 18)  
            atm_label( i)="Ar"
           case( 19)   
            atm_label( i)="K "
           case( 20)  
            atm_label( i)="Ca"
           case( 21)
            atm_label( i)="Sc"
           case( 22)        
            atm_label( i)="Ti"
           case( 23)        
            atm_label( i)="V "
           case( 24)        
            atm_label( i)="Cr"
           case( 25)        
            atm_label( i)="Mn"
           case( 26)        
            atm_label( i)="Fe"
           case( 27)       
            atm_label( i)="Co"
           case( 28)       
            atm_label( i)="Ni"
           case( 29)       
            atm_label( i)="Cu"
           case( 30)       
            atm_label( i)="Zn"
           case( 31)       
            atm_label( i)="Ga"
           case( 32)       
            atm_label( i)="Ge"
           case( 33)       
            atm_label( i)="As"
           case( 34)       
            atm_label( i)="Se"
           case( 35)       
            atm_label( i)="Br"
           case( 36)       
            atm_label( i)="Kr"
           case( 37)       
            atm_label( i)="Rb"
           case( 38)       
            atm_label( i)="Sr"
           case( 39)        
            atm_label( i)="Y "
           case( 40)        
            atm_label( i)="Zr"
           case( 41)        
            atm_label( i)="Nb"
           case( 42)        
            atm_label( i)="Mo"
           case( 43)        
            atm_label( i)="Tc"
           case( 44)        
            atm_label( i)="Ru"
           case( 45)        
            atm_label( i)="Rh"
           case( 46)        
            atm_label( i)="Pd"
           case( 47)        
            atm_label( i)="Ag"
           case( 48)        
            atm_label( i)="Cd"
           case( 49)        
            atm_label( i)="In"
           case( 50)        
            atm_label( i)="Sn"
           case( 51)        
            atm_label( i)="Sb"
           case( 52)        
            atm_label( i)="Te"
           case( 53)        
            atm_label( i)="I "
           case( 54)        
            atm_label( i)="Xe"
           case( 55)        
            atm_label( i)="Cs"
           case( 56)        
            atm_label( i)="Ba"
           case( 57)        
            atm_label( i)="La"
           case( 58)        
            atm_label( i)="Ce"
           case( 59)        
            atm_label( i)="Pr"
           case( 60)        
            atm_label( i)="Nd"
           case( 61)        
            atm_label( i)="Pm"
           case( 62)        
            atm_label( i)="Sm"
           case( 63)        
            atm_label( i)="Eu"
           case( 64)        
            atm_label( i)="Gd"
           case( 65)        
            atm_label( i)="Tb"
           case( 66)        
            atm_label( i)="Dy"
           case( 67)        
            atm_label( i)="Ho"
           case( 68)        
            atm_label( i)="Er"
           case( 69)        
            atm_label( i)="Tm"
           case( 70)        
            atm_label( i)="Yb"
           case( 71)        
            atm_label( i)="Lu"
           case( 72)        
            atm_label( i)="Hf"
           case( 73)        
            atm_label( i)="Ta"
           case( 74)        
            atm_label( i)="W "
           case( 75)        
            atm_label( i)="Re"
           case( 76)        
            atm_label( i)="Os"
           case( 77)        
            atm_label( i)="Ir"
           case( 78)        
            atm_label( i)="Pt"
           case( 79)        
            atm_label( i)="Au"
           case( 80)        
            atm_label( i)="Hg"
           case( 81)        
            atm_label( i)="Tl"
           case( 82)        
            atm_label( i)="Pb"
           case( 83)        
            atm_label( i)="Bi"
           case( 84)        
            atm_label( i)="Po"
           case( 85)        
            atm_label( i)="At"
           case( 86)        
            atm_label( i)="Rn"
           case( 87)        
            atm_label( i)="Fr"
           case( 88)        
            atm_label( i)="Ra"
           case( 89)        
            atm_label( i)="Ac"
           case( 90)        
            atm_label( i)="Th"
           case( 91)        
            atm_label( i)="Pa"
           case( 92)        
            atm_label( i)="U "
           case( 93)        
            atm_label( i)="Np"
           case( 94)        
            atm_label( i)="Pu"
           case( 95)        
            atm_label( i)="Am"
           case( 96)        
            atm_label( i)="Cm"
           case( 97)        
            atm_label( i)="Bk"
           case( 98)        
            atm_label( i)="Cf"
           case( 99)        
            atm_label( i)="Es"
           case( 100)       
            atm_label( i)="Fm"
           case( 101)       
            atm_label( i)="Md"
           case( 102)       
            atm_label( i)="No"
           case( 103)       
            atm_label( i)="Lr"
           case( 104)       
            atm_label( i)="Rf"
           case( 105)       
            atm_label( i)="Db"
           case( 106)       
            atm_label( i)="Sg"
           case( 107)       
            atm_label( i)="Bh"
           case( 108) 
            atm_label( i)="Hs"  
           case( 109) 
            atm_label( i)="Mt"  
           case( 110) 
            atm_label( i)="Ds"  
         end select 

       enddo
!**********************************************************************!
      end subroutine set_label 
!**********************************************************************!
!
!-----------------------------------------------------------------------
! 
      end module m_psolver_bigdft
