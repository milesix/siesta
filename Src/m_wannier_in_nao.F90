!
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!> \brief General purpose of the m_wannier_in_nao module:
!! Here we write the expansion of a Wannier function, 
!! output of the WANNIER-90 code, in the numerical atomic orbital (NAO)
!! basis set of SIESTA.
!! 
!! The analytical expressions can be found in
!!
!! <https://personales.unican.es/junqueraj/JavierJunquera_files/Notes/Wannier/wannier_in_nao.pdf>
!!

module m_wannier_in_nao

  use precision,      only: dp                ! Real double precision type

  implicit none

!
! Variables related with the coefficients of the wavefunctions and
! eigenvalues at the Wannier90 k-point mesh
!
  complex(dp), pointer :: coeffs_wan_nao(:,:) => null()
                                             ! Coefficients of the 
                                             !   Wannier functions in a basis
                                             !   of NAO
                                             !   First  index: Wannier function
                                             !   Second index: NAO in the 
                                             !       supercell
  complex(dp), pointer :: coeffs_opt(:,:,:) => null()

  CONTAINS

  subroutine wannier_in_nao
    use units,          only: Ang            ! Conversion factor from Bohr to Ang
    use parallel,       only: Nodes          ! Total number of Nodes
    use parallel,       only: Node           ! Local Node
    use parallel,       only: IONode         ! Input/output node
    use parallel,       only: BlockSize      ! Input/output node
    use parallelsubs,   only: GetNodeOrbs    ! Calculates the number of orbitals stored 
                                             !    on the local Node.
    use parallelsubs,   only: LocalToGlobalOrb
                                             ! Converts an orbital index in the 
                                             !    local frame to the global frame
    use parallelsubs,   only: GlobalToLocalOrb
                                             ! Converts an orbital index in the
                                             !    global frame to the local frame
                                             !    if the orbital is local to this node. 
                                             !    Otherwise the pointer is
                                             !    return as zero and can therefore 
                                             !    be used to test whether the
                                             !  orbital is local or not.
    use siesta_geom,    only: na_s           ! Number of atoms in the supercell
    use siesta_geom,    only: xa             ! Atomic positions for all the
                                             !   atoms in the supercell
                                             !   In Bohrs
    use siesta_geom,    only: isa            ! Species index of each atom
    use siesta_geom,    only: scell          ! Supercell lattice vectors
    use atomlist,       only: no_s           ! Number of atomic orbitals in the
                                             !   supercell
    use atomlist,       only: no_u           ! Number of atomic orbitals in the
                                             !   unit cell
    use atomlist,       only: indxuo         ! Index of equivalent orbital in 
                                             !   the unit cell
    use atomlist,       only: iaorb          ! Atomic index of each orbital
    use atomlist,       only: iphorb         ! Orbital index of each orbital 
                                             !   in its atom
    use atomlist,       only: lasto          ! Position of last orbital
                                             !   of each atom
    use atmfuncs,       only: labelfis       ! Atomic label
    use atmfuncs,       only: cnfigfio       ! Principal quantum number of the
                                             !   shell to what the orbital 
                                             !   belongs (for polarization 
                                             !   orbitals the quantum number 
                                             !   corresponds to the shell which
                                             !   is polarized by the orbital io)
    use atmfuncs,       only: symfio         ! Symmetry of the orbital
    use w90_parameters, only: num_kpts       ! Total number of k-points that
                                             !   will be used within
                                             !   Wannier90
    use w90_parameters, only: kpt_latt       ! Coordinates of the k-points
                                             !   that will be used within
                                             !   Wannier90
                                             !   In fractional units, i. e.
                                             !   in units of the reciprocal
                                             !   lattice vectors
                                             !   First  index: component
                                             !   Second index: k-vector
    use w90_parameters, only: recip_lattice  ! Reciprocal lattice vectors
                                             !   as used inside
                                             !   Wannier90
                                             !   The factor 2.0 * pi
                                             !   is included
                                             !   BE CAREFUL: In WANNIER90,the
                                             !   order of the indices is
                                             !   inverted with respect to
                                             !   SIESTA
                                             !   First  index: vector
                                             !   Second index: component
                                             !   In Angstroms^-1
    use w90_parameters, only: num_proj       ! Number of projections
    use w90_parameters, only: num_bands      ! Number of bands Wannierized
    use w90_parameters, only: u_matrix       ! Unitary rotations from the 
                                             !   optimal subspace to the
                                             !   optimally smooth states.
    use w90_parameters, only: u_matrix_opt   ! Unitary rotations from the 
    use w90_parameters, only: gamma_only     ! Only the gamma point will be
                                             !   used within Wannier90?
    use w90_parameters, only : disentanglement ! logical value
                                             !   .true.:
                                             !   disentanglement active
                                             !   .false.:
                                             !   disentanglement inactive
    use w90_constants,  only: bohr_angstrom_internal  
                                             ! Conversion factor from 
                                             !   Bohr to Ang
    use w90_parameters, only: wannier_centres
                                             ! Centres of the Wannier functions
    use w90_in_siesta_types,       only: num_proj_local
                                             ! Number of projectors that will be 
                                             !   handled in the local node
    use m_switch_local_projection, only: coeffs
                                             ! Coefficients of the wavefunctions
                                             !   First  index: orbital
                                             !   Second index: band
                                             !   Third  index: k-point
    use m_switch_local_projection, only: seedname
                                             ! Name of the file where the 
                                             !   WANNIER90 code
                                             !   reads or dumps the
                                             !   information.
    use writewave,      only: writew         ! Subroutine to dump the 
                                             !   coefficients of the 
                                             !   wavefunction
                                             !   in a .WFSX file
    use writewave,      only: setup_wfs_list ! Subroutine to setup the
                                             !   list of wave functions 
                                             !   (in this case, Wannier func.)
                                             !   to be plotted
    use writewave,      only: wfs_filename   ! Name of the file where the 
                                             !   coefficients of the Wanniers
                                             !   will be dumped
    use writewave,      only: wwf            ! Write wave functions?

!
!   Variables related with the search for neighbour orbitals
!   around a given Wannier center
!
    use neighbour,     only: maxnna          ! Maximum number of neighbours
    use neighbour,     only: iana=>jan       ! Atom-index of neighbours
    use neighbour,     only: r2ki=>r2ij      ! Squared distances to neighbors
    use neighbour,     only: xki=>xij        ! Vector from a given atom
                                             !   to its neighbours
    use neighbour,     only: x0              ! Position of the point around
                                             !   which we are going to compute
                                             !   the neighbours
    use neighbour,     only: mneighb         ! Subroutine to compute the
                                             !   neighbour list


!
! Allocation/Deallocation routines
!
    use m_switch_local_projection, only: coeffs
                                             ! Coefficients of the wavefunctions
    use alloc,                     only: re_alloc  
                                             ! Reallocation routines
    use alloc,                     only: de_alloc  
                                             ! Deallocation routines

!
! Termination routines
!
    
    use sys,                       only: die

#ifdef MPI
    use mpi_siesta
    use alloc,                     only: de_alloc  
                                             ! Deallocation routines
    use m_orderbands,              only: which_band_in_node  
                                             ! Given a node and a local index,
                                             !   this array gives the
                                             !   global index of the band
                                             !   stored there
    use m_orderbands,              only: sequential_index_included_bands
                                             ! Sequential number of the
                                             !   bands included for
                                             !   wannierization
                                             !   (the bands are listed
                                             !   in order of incremental
                                             !   energy)
    use m_switch_local_projection, only: nincbands_loc
                                             ! Number of bands for
                                             !   wannierization
                                             !   after excluding bands
                                             !   in the local Node
    use parallelsubs,              only: set_blocksizedefault
#endif

! 
! Internal variables 
! 
    integer  :: iu                           ! Logical unit
    integer  :: ik                           ! Counter for loop on k-points
    integer  :: iorb                         ! Counter for loop on atomic 
                                             !   orbitals
    integer  :: iorb_global                  ! Counter for loop on atomic 
                                             !   orbitals (global index)
    integer  :: iorb_local                   ! Counter for loop on atomic 
                                             !   orbitals (local index)
    integer  :: iband                        ! Counter for loop on bands
    integer  :: iproj                        ! Counter for loop on projectors
    integer  :: iproj_global                 ! Counter for loop on projectors
    integer  :: iprojn                       ! Counter for loop on projectors
    integer  :: iprojm                       ! Counter for loop on projectors
    integer  :: index_orbital_unit           ! Equivalent index in the unit cell
                                             !   of an orbital in the supercell
    integer  :: index_atom                   ! Index of the atom to which an
                                             !   atomic orbital belongs to
    integer  :: nk                           ! Number of k-points written in
                                             !   the .WFSX file
                                             !   Since the number of Wanniers 
                                             !   functions to be plotted do not 
                                             !   depend in k-points,
                                             !   it is set to 1
    integer  :: maxspn                       ! Maximum number of spin components
    logical  :: gamma                        ! Whether only the Gamma-point is 
                                             !    sampled.
                                             !    Since the dependency of the 
                                             !    phase is already included in 
                                             !    the coefficients, we set it
                                             !    up to .true.
    real(dp) :: kpoint(3)                    ! Coordinates of the k-point
    real(dp) :: kxmu                         ! Dot product of the k-point and
                                             !   the position of the atom in the
                                             !   supercell
    real(dp) :: ckxmu                        ! Cosine of kxmu
    real(dp) :: skxmu                        ! Sine of kxmu
    real(dp) :: kdummy(3)                    ! Dummy variable for the k-points
    real(dp), pointer :: psi(:,:,:)          ! Dummy variable to store the 
                                             !   coefficients of the Wanniers
                                             !   in a basis of atomic orbitals
                                             !   to call writew
    real(dp), target :: aux(2,no_s*5)        ! Dummy variable that will play
                                             !   the role of the eigenvalues
                                             !   in the call to writew
    real(dp) :: rangemax                     ! Maximum module of the supercell
                                             !   lattice vectors
    real(dp) :: vecmod                       ! Module of the supercell
                                             !   lattice vectors
    integer  :: ivector                      ! Counter for supercell lattice 
                                             !   vectors
    integer  :: nuo                          ! Number of orbitals stored in the 
                                             !   local node
    integer  :: nneig                        ! Maximum number of neighbours
    integer  :: ina                          ! Counter for loop on neighbours
    integer  :: ia                           ! Index of neighbour atom
    integer  :: is                           ! Species of neighbour atom
    real(dp) :: rki                          ! Distance to neighbour
    integer  :: no_local                     ! Number of orbitals in the 
                                             !   supercell
                                             !   stored in the local node
    integer  :: blocksizeincbands_wannier    !  Maximum number of bands
                                             !   considered for per node
    integer  :: blocksize_save               !  Maximum number of bands
    logical, dimension(:), pointer :: listed
    logical, dimension(:), pointer :: listedall 
                                             ! List of orbitals of the supercell
                                             !   stored in the local node
    complex(dp), dimension(:,:), pointer :: psiloc => null() 
                                             ! Coefficients of the wave
                                             !   function (in complex format)
#ifdef MPI
    integer     :: i
    integer     :: iband_global              ! Global index for a band
    integer     :: iband_sequential          ! Sequential index of the band
    integer     :: MPIerror
    complex(dp), dimension(:,:), pointer :: auxloc => null()
                                             ! Temporal array for the
                                             !   the global reduction of Amnmat
    integer, external :: numroc
#endif

!!   For debugging
!#ifdef MPI
!    do iproj = 1, num_proj
!      write(6,'(a,3i5,3f12.5)')                                     &
! &      'Nodes, Node, iproj, wannier_centres = ',                   &
! &      Nodes, Node, iproj, wannier_centres(:,iproj) 
!    enddo 
!    call MPI_barrier(MPI_Comm_world,i)
!    call die()
!#endif
!!   End debugging

! Find maximum range and maximum number of KB projectors
    rangemax = 0.0_dp
    do ivector = 1, 3
      vecmod = sqrt(dot_product(scell(:,ivector),scell(:,ivector)))
      rangemax = max(rangemax,vecmod)
    enddo
    rangemax = rangemax / 2.0_dp

!!   For debugging
!    write(6,'(a,2i5,f12.5)')                        &
! &    'Nodes, Node, rangemax = ',                   &
! &     Nodes, Node, rangemax
!!    call MPI_barrier(MPI_Comm_world,i)
!!    call die()
!!   End debugging

!   Calculate the number of orbitals in the unit cell stored in the local node
    call GetNodeOrbs(no_u,Node,Nodes,nuo)

!   Calculate the number of orbitals in the supercell stored in the local node
    call GetNodeOrbs(no_s,Node,Nodes,no_local)

!   Calculate the number of projections stored in the local node
    call set_blocksizedefault( Nodes, num_proj,     &
 &                             blocksizeincbands_wannier )
    num_proj_local = numroc( num_proj, blocksizeincbands_wannier, &
 &                           Node, 0, Nodes )

!!   For debugging
!    write(6,'(a,6i5)')                              &
! &    'Nodes, Node, no_u, nuo, no_s, no_local = ',  &
! &     Nodes, Node, no_u, nuo, no_s, no_local
!    write(6,'(a,5i5)')                              &
! &    'Nodes, Node, blocksizeincbands_wannier, num_proj, num_proj_local = ',  &
! &     Nodes, Node, blocksizeincbands_wannier, num_proj, num_proj_local 
!    blocksize_save = BlockSize
!    BlockSize = blocksizeincbands_wannier
!    do iproj = 1, num_proj_local
!      call LocalToGlobalOrb(iproj, Node, Nodes, iproj_global)
!      write(6,'(a,4i5)')' Nodes, Node, iproj, iproj_global = ', &
! &                        Nodes, Node, iproj, iproj_global 
!    enddo
!    BlockSize = blocksize_save
!    call MPI_barrier(MPI_Comm_world,i)
!    call die()
!!   End debugging

!   Initialize the list of visited orbitals
    nullify( listed )
    call re_alloc( listed, 1, no_s, 'listed', 'wannier_in_nao' )
    listed(:) = .false.
    nullify( listedall )
    call re_alloc( listedall, 1, no_s, 'listedall', 'wannier_in_nao' )
    listedall(:) = .false.

!   Make list of all orbitals needed for this node
    do iorb = 1, no_local
      ! we need this process's orbitals...
      call LocalToGlobalOrb(iorb,Node,Nodes,iorb_global)
      listedall(iorb_global) = .true.
!!      For debugging
!      write(6,'(a,4i5,l5)')                              &
! &      'Nodes, Node, iorb, iorb_global, listedall = ',  &
! &       Nodes, Node, iorb, iorb_global, listedall(iorb_global)
!!      End debugging
    enddo

!   Allocate the array where the coefficients of the Wannier functions
!   in a basis of Numerical Atomic Orbitals will be stored
    nullify( coeffs_wan_nao )
    call re_alloc( coeffs_wan_nao,                                  &
 &                 1, num_proj_local,                               &
 &                 1, no_s,                                         &
 &                 name='coeffs_wan_nao', routine='wannier_in_nao')
    coeffs_wan_nao = cmplx(0.0_dp,0.0_dp,kind=dp)

!   Allocate the array where the coefficients of the 
!   bands that will be wannierized
!   If there is no-disentanglement required, the coefficients are
!   the same that come out from the diagonalization.
!   If a disentanglement procedure is required 
!   (more bands that Wannier functions),
!   then we follow the recipe described in Sec. III A of the paper by
!   I. Souza et al. Phys. Rev. B 65 035109 (2001)
!   To find the N-dimensional subspace if the number of bands for
!   a k-point, N_k, is larger than N, we have to multiply the
!   coefficients by a matrix that comes from Wannier90
    nullify( coeffs_opt )
    call re_alloc( coeffs_opt,                                      &
 &                 1, no_u,                                         &
 &                 1, num_proj,                                     &
 &                 1, num_kpts,                                     &
 &                 name='coeffs_opt', routine='wannier_in_nao')
    coeffs_opt = cmplx(0.0_dp,0.0_dp,kind=dp)

!   Allocate memory related with a local variable where the coefficients
!   of the eigenvector at the k-point will be stored
!   Only num_bands are retained for wannierization, that is why the
!   second argument is made equal to num_bands
    call re_alloc( psiloc, 1, no_u, 1, num_bands,     &
 &                'psiloc', 'wannier_in_nao' )

#ifdef MPI
!   Store the local bands in this node on a complex variable
    do ik = 1, num_kpts

!     Initialize the local coefficient matrix for every k-point
      psiloc(:,:) = cmplx(0.0_dp, 0.0_dp, kind=dp)

      do iband = 1, nincbands_loc
        iband_global     = which_band_in_node(Node,iband)
        iband_sequential = sequential_index_included_bands(iband_global)

!!       For debugging
!        write(6,'(a,7i5)')                        &
! &       'Nodes, Node, ik, nbands_loc, iband, iband_global, iband_sequential=',&
! &        Nodes, Node, ik, nincbands_loc, iband, iband_global, iband_sequential
!!       End debugging

        do iorb = 1, no_u
          psiloc(iorb,iband_sequential) = coeffs(iorb,iband,ik)
        enddo 

      enddo 
!     Allocate workspace array for global reduction
      call re_alloc( auxloc, 1, no_u, 1, num_bands,                  &
 &                   name='auxloc', routine='wannier_in_nao' )
!     Global reduction of auxloc matrix
      auxloc(:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
      call MPI_AllReduce( psiloc(1,1), auxloc(1,1),        &
 &                        no_u*num_bands,                  &
 &                        MPI_double_complex,MPI_sum,MPI_Comm_World,MPIerror )
!     After this reduction, all the nodes know the coefficients of the
!     wave function for the point ik, for all the bands and for all atomic
!     orbitals
      psiloc(:,:) = auxloc(:,:)

      if( disentanglement ) then
        do iorb = 1, no_u
          do iprojn = 1, num_proj
            do iband = 1, num_bands
              coeffs_opt(iorb,iprojn,ik) = coeffs_opt(iorb,iprojn,ik) +  &
 &              u_matrix_opt(iband,iprojn,ik) * psiloc(iorb,iband)
            enddo 
          enddo
        enddo
      else
        coeffs_opt(:,:,ik) = psiloc(:,:) 
      endif

!!     For debugging
!      write(6,'(a,i5,3f12.5)')      &
! &      ' ik, kpt_latt(:,ik) = ',   &
! &        ik, kpt_latt(:,ik) 
!      if( disentanglement ) then
!        do iprojn = 1, num_proj
!          do iband = 1, num_bands
!            write(6,'(a,3i5,2f12.5)')                & 
! &           'Node, iband, iprojn, u_matrix_opt = ', &
! &            Node, iband, iprojn, u_matrix_opt(iband,iprojn,ik)
!          enddo 
!        enddo 
!      endif
!      do iband = 1, num_bands
!        do iorb = 1, no_u
!          write(6,'(3i5,2f12.5)') Node, iband, iorb, psiloc(iorb,iband)
!        enddo 
!      enddo 
!      do iprojn = 1, num_proj
!        do iorb = 1, no_u
!          write(6,'(a,3i5,2f12.5)')                                     &
! &          'Node, iprojn, iorb, coeffs_opt(iorb,iprojn,ik) = ',   &
! &          Node, iprojn, iorb, coeffs_opt(iorb,iprojn,ik)
!        enddo 
!      enddo 
!      call MPI_barrier(MPI_Comm_world,i)
!      call die()
!!     End debugging

    enddo    ! End loop on k-points
#else
    do ik = 1, num_kpts
      do iband = 1, num_bands
        do iorb = 1, no_u
          psiloc(iorb,iband) = coeffs(iorb,iband,ik)
        enddo 
      enddo
      if( disentanglement ) then
        do iorb = 1, no_u
          do iprojn = 1, num_proj
            do iband = 1, num_bands
              coeffs_opt(iorb,iprojn,ik) = coeffs_opt(iorb,iprojn,ik) +  &
 &              u_matrix_opt(iband,iprojn,ik) * psiloc(iorb,iband)
            enddo 
          enddo
        enddo
      else
        coeffs_opt(:,:,ik) = psiloc(:,:) 
      endif
    enddo
#endif

!   Initialize neighb subroutine
    nneig = 0
    x0(:) = 0.0_dp
    call mneighb( scell, rangemax, na_s, xa, 0, 0, nneig )
!!   For debugging
!    write(6,'(a,4i5)')                        &
! &    'Nodes, Node, nneig = ',                  &
! &     Nodes, Node, nneig
!    call MPI_barrier(MPI_Comm_world,i)
!    call die()
!!   End debugging

!   Loop over all the local Wannier functions
!   All processes will be doing this loop over Wanniers
    do iprojn = 1, num_proj_local
      blocksize_save = BlockSize
      BlockSize = blocksizeincbands_wannier
      call LocalToGlobalOrb(iprojn, Node, Nodes, iproj_global)
      BlockSize = blocksize_save

!     Search for all the atoms neighbour to that particular Wannier function.
!     The relative vector xki between the Wannier and the atom is centered
!     on the Wannier
!
!     The Wannier centers are assumed to be given in Angstroms.
!     we transform them to Bohrs
      x0(:) = wannier_centres(:,iproj_global) * Ang
      call mneighb( scell, rangemax, na_s, xa, 0, 0, nneig )

      if ( nneig .gt. maxnna )                                     &
 &       call die('wannier_in_nao: insufficient array shapes; see neighb(..)')

!     From the position of the center of the Wannier
!     with respect the origin, and
!     the position of the atom with respect the Wannier,
!     compute the position of the atom
      do ina = 1, nneig
        ia  = iana(ina)
        is  = isa(ia)
        xki(:,ina) = xki(:,ina) + x0(:)

!!       For debugging
!        write(6,'(a,5i5)')                                   &
! &        'Nodes, Node, iprojn, iproj_global, nneig = ' ,    &
! &         Nodes, Node, iprojn, iproj_global, nneig 
!        rki = sqrt(r2ki(ina))
!        ia  = iana(ina)
!        is  = isa(ia)
!        write(6,'(a,6i5,4f12.5)')                            &
! &        'Nodes, Node, iprojn, ina, r2ki, rki, ia, is = ',   &
! &         Nodes, Node, iprojn, ina, ia, is, xki(:,ina), rki
!!       End debugging

!       Loop on the atomic orbitals of the neighbour atom
        do iorb = lasto(ia-1)+1, lasto(ia) 

!         Identify the index of the atom in the unit cell
          index_orbital_unit = indxuo(iorb)

!          call GlobalToLocalOrb(iorb,Node,Nodes,iorb_local)

!!         Only calculate if needed locally in our MPI process
!          if (.not. listedall(iorb)) CYCLE

!         Up to here, we know:
!         - the projector that will be expressed in a basis of NAO
!         - the atomic orbital for which the coefficient will be found
!         - the position of the atom where the atomic orbital is centered
!         Now, we perform the sum on k-point in the Equation written
!         in the header of the subroutine

          do ik = 1, num_kpts
!           Compute the coordinates of the k-point (in Ang^-1)
            kpoint(:) = kpt_latt(1,ik) * recip_lattice(1,:) +      &
 &                      kpt_latt(2,ik) * recip_lattice(2,:) +      &
 &                      kpt_latt(3,ik) * recip_lattice(3,:) 
!           Transform the coordinates of the k-point to Bohr^-1
            kpoint(:) = kpoint(:) * bohr_angstrom_internal

!           Compute the dot product between the k-point and the
!           atomic position
            kxmu = kpoint(1) * xki(1,ina) +                        &
 &                 kpoint(2) * xki(2,ina) +                        &
 &                 kpoint(3) * xki(3,ina) 
            ckxmu = dcos(kxmu)
            skxmu = dsin(kxmu)

            do iprojm = 1, num_proj
              coeffs_wan_nao(iprojn,iorb) =                        &
 &              coeffs_wan_nao(iprojn,iorb) +                      &
 &              u_matrix(iprojm,iproj_global,ik)            *      &
 &              coeffs_opt(index_orbital_unit,iprojm,ik)    *      &
 &              cmplx(ckxmu,skxmu,kind=dp) 
            enddo 

          enddo ! End loop on k-points
        enddo   ! End loop on atomic orbitals
      enddo     ! End loop on neighbour atoms 
    enddo       ! End loop on projections

!   Divide by the number of k-points
    coeffs_wan_nao = coeffs_wan_nao / num_kpts

!!   For debugging
!    do iproj = 1, num_proj_local
!      do iorb = 1, no_s
!        write(6,'(a,5i5,5f12.5)')   &
! &        'Node, Nodes, iproj, iorb, indxuo, coeffs_wan_nao(iproj,iorb) = ',  &
! &         Node, Nodes, iproj, iorb, indxuo(iorb), coeffs_wan_nao(iproj,iorb),&
! &         xa(:,iaorb(iorb)) + wannier_centres(:,iproj) * Ang
!      enddo 
!    enddo 
!    call MPI_barrier(MPI_Comm_world,i)
!    call die()
!   End debugging

!   Allocate the array where the coefficients of the Wannier functions
!   in a basis of Numerical Atomic Orbitals will be stored
    nullify( psi )
    call re_alloc( psi,                                             &
 &                 1, 2,                                            &
 &                 1, no_s,                                         &
 &                 1, num_proj_local,                               &
 &                 name='psi', routine='wannier_in_nao')
    psi = 0.0_dp

!!    For debugging
!    do ik = 1, num_kpts
!      do iproj = 1, num_proj
!        do iorb = 1, no_u
!          write(6,'(a,3i5,4f12.5)')' ik, iproj, iorb, coeff = ', &
! &          ik, iproj, iorb, coeffs(iorb,iproj,ik),              &
! &          u_matrix(iproj,iproj,ik)
!        enddo 
!      enddo 
!    enddo 
!!    End debugging


!!   For debugging
!    do iorb = 1, no_s
!      write(6,'(a,3i5,3f12.5,3f12.5)')'iorb, indxuo, iaorb = ',   &
! &      iorb, indxuo(iorb), iaorb(iorb), xa(:,iaorb(iorb)),       &
! &      xafold(:,iaorb(iorb)) 
!! &      iatfold(:,iaorb(iorb))
!    enddo
!     do iproj = 1, num_proj_local
!       do iorb = 1, no_s
!!!         if( indxuo(iorb) .eq. 15) then
!!         if( aimag(coeffs_wan_nao(iproj,iorb)) .gt. 1.d-5 ) then
!         if( Node .eq. 1) then
!         write(6,'(a,4i5,5f12.5)') ' Node, Nodes, iproj, iorb, coeffs_wan_nao = ', &
! &         Node, Nodes, iproj, iorb, coeffs_wan_nao(iproj,iorb), xa(:,iaorb(iorb))
!         endif
!!         endif
!       enddo 
!     enddo 
!     call MPI_barrier(MPI_Comm_world,i)
!     call die()
!!   End debugging

!   Set up the variables to call writew
    do iproj = 1, num_proj_local
      do iorb = 1, no_s
        psi(1,iorb,iproj) = real(coeffs_wan_nao(iproj,iorb))
        psi(2,iorb,iproj) = aimag(coeffs_wan_nao(iproj,iorb))
!!       For debugging
!        if( abs(psi(1,iorb,iproj)) .gt. 1.d-5) then
!        if (iproj .eq. 4 ) then
!        write(6,'(a,2i5,8f12.5)')' iproj, iorb, psi = ', &
! &        iproj, iorb, psi(1,iorb,iproj), psi(2,iorb,iproj), xa(:,iaorb(iorb)), xafold(:,iaorb(iorb)) - xa(:,1)
!        endif 
!        if (iproj .eq. 5 ) then
!        write(6,'(a,2i5,8f12.5)')' iproj, iorb, psi = ', &
! &        iproj, iorb, psi(1,iorb,iproj), psi(2,iorb,iproj), xa(:,iaorb(iorb)), xafold(:,iaorb(iorb)) - xa(:,2)
!        endif 
!        endif
!!       End debugging
      enddo 
    enddo

    aux          = 0.0_dp
    kdummy       = 0.0_dp
    wfs_filename = trim(seedname)//".WFSX"
    wwf          = .false.
    nk           = 1
    gamma        = .false.
    maxspn       = 1

! 
!   Open the WFSX file and print the header of the file
!   with information of the atomic orbitals
!   This was done in the subroutine wwave when the coefficients of the 
!   wave functions at particular k-points are required
!
    if( IONode ) then
      call io_assign( iu )
      open(iu, file=wfs_filename,form="unformatted",status='unknown')

      write(iu) nk, gamma
      write(iu) maxspn
      write(iu) no_s
      write(iu) (iaorb(iorb),labelfis(isa(iaorb(iorb))),                    &
 &              iphorb(iorb), cnfigfio(isa(iaorb(iorb)),iphorb(iorb)),      &
 &              symfio(isa(iaorb(iorb)),iphorb(iorb)), iorb=1,no_s)
      call io_close( iu )
    endif

!!     For debugging
!      write(6,'(a,i5,l5)')' nk, gamma = ', nk, gamma
!      write(6,'(a,i5)')   ' maxspn = ', maxspn
!      write(6,'(a,i5)')   ' no_s   = ', no_s
!      do iorb = 1, no_s
!        write(6,*) iorb, iaorb(iorb),labelfis(isa(iaorb(iorb))),            &
! &                 iphorb(iorb), cnfigfio(isa(iaorb(iorb)),iphorb(iorb)),   &
! &                 symfio(isa(iaorb(iorb)),iphorb(iorb))
!      enddo 
!!     End debugging

    call setup_wfs_list(nk,no_s,1,num_proj,.false.,.false.)

    call writew( no_s, num_proj, 1, kdummy, 1, &
 &               aux, psi, gamma )


  end subroutine wannier_in_nao

end module m_wannier_in_nao
