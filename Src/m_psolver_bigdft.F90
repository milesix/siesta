! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!> \brief General purpose of the poisson_bigdft module 
!! 
!!  Interface module to call BigDFT solver instead of standar Siesta 
!!  poison. It requires compilation of Siesta with additional
!!  BigDFT libraries.
!!  More info at:
!!  http://bigdft.org/Wiki/index.php?title=BigDFT_website 
!!
!!  Created by Pablo Lopez-Tarifa and Daniel Sanchez Portal @CFM(2018)
!
      module m_psolver_bigdft

!
      use precision,      only : dp, grid_p
      use parallel,       only : Node, Nodes, ionode, ProcessorY
#ifdef MPI
      use mpi_siesta
#endif
!
      implicit none
      public :: poisson_bigdft
!
      CONTAINS
!
!> \brief General purpose of the subroutine poisson_bigdft
!!
!! This subroutine handles the full hartree potential calculation 
!! given by Psolver package of BigDFT. The main points are:
!!
!! 1. The system is classified accoring its symmetry (shape). 
!! 
!! 2. The potential is calculated in two scenarios:
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
!! 2. Inclusion of other simmetries. This has been already worked-out 
!!    in serial so hard work is not expected here. 
!! 
!! 3. Design of an additional subroutine (poisson_bigdft_enviroment ?)
!!    to insert an implicit solvent effect. New arguments as atom kinds
!!    and possitions must to be passed through. 
!!
!! 4. Hardwired solution in the case of shape=molecule (i.e. gas-phase), 
!!    to move the molecule to the centre of the box. The same issue 
!!    has to be solved to 2D periodic systems. 
!!    
!!  
! #################################################################### !
      subroutine poisson_bigdft(cell, lgrid, grid, ntm, RHO, V,  eh)
! #################################################################### !
!                                                                      !
!   This subroutine sets up the kernel object for a bigDFT Poisson     !
!   solver calculation of the Hartree potential.                       !
!                                                                      !
! #################################################################### !
!
      use units,          only: eV
      use siesta_geom,    only: shape
      use alloc,          only : re_alloc, de_alloc
      use Poisson_Solver, only: coulomb_operator,  electrostatic_solver, &
                                pkernel_set
      use mesh,           only : nsm
      use PStypes,        only: pkernel_init, pkernel_free
      use yaml_output
      use futile
      use f_utils
      use dictionaries, dict_set => set
!
! #################################################################### !
!
      real(dp), intent(in)      :: cell( 3, 3)
      integer, intent(in)       :: lgrid( 3)       ! local mesh divisions 
      integer, intent(in)       :: grid( 3)        ! total mesh divisions 
      real(grid_p), intent(in)  :: RHO( :)         ! input density 
      integer, intent(in)       :: ntm
      real(grid_p), intent(out) :: V( :)           ! output potential 
      real(dp), intent(out)     :: eh              ! electrostatic energy 
!
      real(dp), dimension(1,3)  :: hgrid !uniform mesh spacings in the three directions
      integer                   :: i, j, k
      real(dp)                  :: pi
      real(dp)                  :: einit
      type(dictionary), pointer :: dict            ! input parameters
      integer                   :: isf_order, fd_order
      type(coulomb_operator)    :: pkernel
!
      real(grid_p), pointer :: Paux_3D( :, :, :)  ! 3D auxyliary array 
      real(grid_p), pointer :: Paux_1D( :)        ! !D auxyliary array 
      real(dp) :: alpha,beta,gamma !possible angles of the box axis
      real(dp) :: angdeg(3)
      integer               :: World_Comm, mpirank, ierr
      integer, allocatable  :: local_grids( :)
      parameter( pi = 4.D0 * DATAN( 1.D0))
!
! #################################################################### !
! 
! f_util initialisation:
! 
      call f_lib_initialize()
!
! Grid separation for psolver:
!
!
      do i=1,3
        hgrid(1,i)=(sqrt(cell(1,i)**2+cell(2,i)**2+cell(3,i)**2)       &
                /grid(i))
      enddo
!
! 
! Initialise futil dictionary  
!
!#ifdef
!      call MPI_Init( error )
!#endif
      call dict_init(dict)
      dict=>dict_new()
      call set(dict//'setup'//'global_data',.true.)
!      call set(dict//'setup'//'taskgroup_size', 2)
!      call pkernel_lowpez()
!
      nullify(Paux_3D)
      nullify(Paux_1D)
!
      allocate( local_grids( Nodes*3))
      call re_alloc( Paux_3D, 1, grid(1), 1, grid(2), 1, grid(3),     &
                     'Paux_3D', 'poison_bigdft' )
      call re_alloc( Paux_1D, 1, grid(1)*grid(2)*grid(3),             &
                     'Paux_1D', 'poison_bigdft' )
!
#ifdef MPI
!
!!! Collection of RHO: 
!
! RHO is collectedi to a common array.
        call gather_rho( RHO, grid, 2, lgrid(1) * lgrid(2) * lgrid(3), Paux_1D)
!
        call MPI_BCAST( Paux_1D, grid(1) * grid(2) * grid(3),         &
                        MPI_double_precision, 0, MPI_COMM_WORLD, ierr) 
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
        call onetothreedarray( Paux_1D, Paux_3D, grid)
#else
        call onetothreedarray( RHO, Paux_3D, grid)
#endif
!
!     MOLECULAR ISOLATED SYSTEM
!
      if (shape .eq. 'molecule') then  
!
!!! BigDFT solver starts now:
!
        if (ionode) write(6,"(a)")                                     &
            'bigdft_solver: initialising kernel for a molecular system'
!
        pkernel=pkernel_init( Node, Nodes, dict, 'F', (/grid(1),       &
                              grid(2), grid(3)/), hgrid)
!
        call pkernel_set( pkernel)
!
        call electrostatic_solver( pkernel, Paux_3D, ehartree = eh)
!
! Enery must be devided by the total number of nodes because the way 
! Siesta computes the Hartree energy. Indeed, BigDFT provides the value
! alredy distributed. 
!
        eh = eh / Nodes
!
#ifdef MPI
! 
! we go back from a 3D to a 1D array: 
!
        call threetoonedarray( Paux_1D, Paux_3D, grid) 
!
! Now we spread back the potential on different nodes. 
!
        call spread_potential( Paux_1D, grid, 2, grid(1) * grid(2) * &
                   grid(3), local_grids, V) 
# else 
        call threetoonedarray(V, Paux_3D, grid) 
#endif
!
      endif
!
!     2D PERIODIC SYSTEM 
!
      if (shape .eq. 'slab') then  
       if (ionode) write(6,"(a)")                                   &
          'bigdft_solver: initialising kernel for a slab system'
      endif
!
!     3D PERIODIC SYSTEM 
!
      if (shape .eq. 'bulk') then  
       if (ionode) write(6,"(a)")                                   &
          'bigdft_solver: initialising kernel for a periodic system'
      endif
!
!      1D WIRE SYSTEM ===NOT IMPLEMENTED YET=== 
!
!       if (shape .eq. 'bulk') then  
!        if (Node.eq.0) write(6,"(a)") 
!     $    'bigdft_solver: initialising kernel for a periodic system'
!        endif

!
! Ending calculation:
!
        call dict_free(dict)
        call pkernel_free(pkernel)
        call f_lib_finalize()
        call de_alloc( Paux_1D, 'Paux', 'dhscf_init' )
        call de_alloc( Paux_3D, 'Paux', 'dhscf_init' )
!
      end subroutine poisson_bigdft 
!
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
!
      integer, intent(in)          ::     nsm         
      integer, intent(in)          ::     mesh(3)
      integer, intent(in)          ::     maxp
      real(grid_p), intent(in)     ::     rho(maxp)
      real(grid_p), pointer , intent(out) ::     rho_total(:) 
      integer  i, ip, iu, is, np, BlockSizeY, BlockSizeZ, ProcessorZ
      integer  meshnsm(3), NRemY, NRemZ, iy, iz, izm, Ind, Ind_rho, ir
#ifdef MPI
      integer    Ind2, MPIerror, Request, Status(MPI_Status_Size), BNode, NBlock
#endif
      real(grid_p), pointer     :: bdens(:) => null()
      real(grid_p),     pointer :: temp(:) => null()
! #################################################################### !
#ifdef MPI
!
! Work out density block dimensions
!
      if (mod(Nodes,ProcessorY).gt.0)                               &
           call die('ERROR: ProcessorY must be a factor of the' //  &
           ' number of processors!')
!
      ProcessorZ = Nodes / ProcessorY
      BlockSizeY = ((((mesh( 2) / nsm) - 1)/ProcessorY) + 1) *nsm
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
#else
!
      call die('Error gather_rho does not work in serial runs') 
!
#endif
! #################################################################### !
      end subroutine gather_rho
! #################################################################### !
! #################################################################### !
      subroutine spread_potential( V_total, mesh, nsm, maxp,           &
                     local_grids, V_local)
! #################################################################### !
!                                                                      !
!     Based on read_rho writen by J.Soler July 1997,                   !
!     this subroutine takes a potential from a total array and         !
!     distributes it into all nodes according to local_grids meshes.   !
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
#ifdef MPI
      integer    Ind2, MPIerror, Request, meshl(3),     &
                 Status(MPI_Status_Size), BNode, NBlock
#endif
      logical    ltmp
      real(grid_p), pointer :: bdens(:) => null()
      logical    baddim, found
! #################################################################### !
#ifdef MPI
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
                if (BNode .ne. 0) then
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
#else
!
      call die('Error spread_potential does not work in serial runs') 
!
#endif
      end subroutine spread_potential 
!
!-----------------------------------------------------------------------
! 
      end module m_psolver_bigdft
